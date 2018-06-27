package di_blois.onePassClustering;

import di_blois.clustering.RandPerm;
import di_blois.comparator.Comparator;
import di_blois.comparator.Minkowski;
import di_blois.dataset.Data;
import org.apache.commons.math3.distribution.TDistribution;
import di_blois.reader.IReader;
import di_blois.reader.SimpleDataReader;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.atomic.AtomicInteger;

/*
 * Implements a one pass di_blois.clustering algorithm based on various statistical bounds
 * The objective is to evaluate these bounds in the context of di_blois.clustering
 * in terms of max. data set size that is handled, computation speed, rand index, confusion errors
 * ... 
 */
public class BoundedClustering<D> {

	// Constants
	public static final int HOEFFDING_BOUND = 0;
	public static final int BERNSTEIN_BOUND = 1;
	public static final int STUDENT_BOUND = 2;
	public final double error_probability = 0.05;
    //Spawning thread pool
    private final ThreadPoolExecutor pool;
	
	
	// Variables
	private RandPerm r; // Random number generator
	private Random rand;
	private List<List<Integer>> clusters;
	private List<D> data;
	private long comparisons;
	private long comptemplate;
	private long ambiguities;
	private int bound_type = BERNSTEIN_BOUND;
	private double range_value;
	private double mindist = Double.MAX_VALUE;
	
	
	
	/** Default class constructor */
	public BoundedClustering(int bound_type){
		this.bound_type = bound_type;
		this.r = new RandPerm();
		this.rand = new Random();
		this.clusters = new ArrayList<>();
		this.comparisons = 0;
		this.comptemplate = 0;
		this.ambiguities = 0;
		pool = (ThreadPoolExecutor) Executors.newFixedThreadPool(Math.max(2, Runtime.getRuntime().availableProcessors() - 2));
	}
	
	/** Main di_blois.clustering function: scan the di_blois.dataset and assigns data objects to clusters */
	public void cluster(List<D> dataset, Comparator<D> comp, double learning_sample, double size){

		// init variables
		this.data = dataset;
		int count = dataset.size();
		int[] order = r.randperm(dataset.size());

		
		// template computation
		double template = 0;
		int nbrenc = (int) (learning_sample * count);
		this.comptemplate = nbrenc;
		
		double maxdistance = Double.MIN_VALUE;
		double mindistance = Double.MAX_VALUE;
		
		double d;
		for (int k = 0; k < nbrenc; k++){
			//get two different elements at random
			int i = r.nextInt(count);
			int j = i;
			while (i == j) j = r.nextInt(count);
			//compare them and add their distance to the threshold
			try {
				d = comp.compare(this.data.get(i), this.data.get(j));
				template += d;
				if (d > maxdistance) maxdistance = d;
				else if (d < mindistance) mindistance = d;
			} catch (Exception e) {
				System.err.println("Error while comparing" + i + " and " + j);
				e.printStackTrace();
			}
		}
		
		// normalize threshold to mean of the seen exemples
		template = template / nbrenc;

		// range value evaluation: as overestimation is not problematic, 
		// estimated range = 2 times previously evaluated range
		// avoid computing the exact range --> quadratic complexity
		
		range_value = 2* (maxdistance - mindistance);
		
		// cluster building
		
		for (int i = 0; i < dataset.size(); i++){
			// return nearest cluster and set global variable mindist
			if (i%1000==0) {
				System.out.println(i + " points done. With " + comparisons + " comparisons.");
				comparisons = 0;
			}
			
			int cluster = race(order[i], clusters, comp, 0.5);
			
			// at this point, mindist is set to the distance to 'cluster'
			
			// assignment to nearest cluster
			if (cluster == -1 || mindist > template){
				// build a new clustersynchronized
				List<Integer> c = new ArrayList<Integer>();
				c.add(order[i]);
				this.clusters.add(c);
			} else {
				// join nearest cluster
				clusters.get(cluster).add(order[i]);
			}
		}
		
		// finalize cluster
		if (size > 0){
			ArrayList<Integer> others = new ArrayList<Integer>();

			int minSize = (int)(size * count);

			// Suppression effective des clusters
			int c = 0;
			while (c < clusters.size()){
				if (clusters.get(c).size() < minSize){
					others.addAll(clusters.get(c));
					clusters.remove(c);
				} else c++;
			}
			
			if (others.size() > 0) clusters.add(others);
		}

        //Shutting down the pool
        pool.shutdown();
	}
	
	/** Race the clusters */
	private int race(int obj, List<List<Integer>> clusters, Comparator<D> cp, double epsilon_multiplier){
		
		// initially, all clusters are in the race
		boolean[] inrace = new boolean[clusters.size()];
		Arrays.fill (inrace,true);
				
		int still_in_race = clusters.size();
				
		double[] lower_bounds = new double[clusters.size()];
		double[] upper_bounds = new double[clusters.size()];
		Arrays.fill(upper_bounds, Double.MAX_VALUE);
		double[] mean_values = new double[clusters.size()];
		
		// Optimization on distance tables
		
		double[][] distances = new double[clusters.size()][];
		int ind = 0;
		int maxsize = 0;
		for (List<Integer> l : clusters){
			distances[ind] = new double[l.size()];
			if (l.size() > maxsize) maxsize = l.size();
			ind ++; 
		}



		int current_winner = 0;

		for (int i = 0 ; i < maxsize ; i++) {

		    // We setup a latch allowing us to signal the main thread we have finished
            CountDownLatch latch = new CountDownLatch(clusters.size());

			for (int c = 0; c < clusters.size(); c++){

                int finalC = c, finalI = i;
                pool.execute(() -> {
                    List<Integer> clust = clusters.get(finalC);
                    double epsilon;

                    if (inrace[finalC] && (finalI < clust.size())) {

                        //sample another data
                        int d = clust.get(rand.nextInt(clust.size()));

                        try {
                            distances[finalC][finalI] = cp.compare(this.data.get(obj), this.data.get(d));
                            mean_values[finalC] = empiricalMean(distances[finalC], finalI +1);

                            if (bound_type == HOEFFDING_BOUND) {
                                epsilon = epsilon_multiplier * hoeffdingsEpsilon(mean_values[finalC], finalI +1, error_probability, range_value);

                            } else if (bound_type == BERNSTEIN_BOUND){
                                //use bernstein bound
                                double std = empiricalStd(distances[finalC], finalI +1, mean_values[finalC]);
                                epsilon = epsilon_multiplier * bernsteinEpsilon(mean_values[finalC], std, finalI +1, error_probability, range_value);
                            } else {
                                double std = unbiasedEmpiricalStd(distances[finalC], finalI +1, mean_values[finalC]);
                                epsilon = epsilon_multiplier * studentsEpsilon(std, finalI +1, 1 - error_probability, range_value);
                            }


                            lower_bounds[finalC] = mean_values[finalC] - epsilon;
                            upper_bounds[finalC] = mean_values[finalC] + epsilon;


                        } catch (Exception e) {
                            System.err.println("Error while comparing" + obj + " and " + d);
                            e.printStackTrace();
                        }
                    }
                    latch.countDown();
                });

			} // end for each cluster
            comparisons += clusters.size();
            //Wait for all comparisons to finish
            try {latch.await();}
            catch (InterruptedException e) {
			    System.err.printf("Fatal error compute thread interrupted, aborting commutation for %s. %n", data.get(obj));
			    return 0;
            }

            for (int c = 0; c < clusters.size(); c++) {
                if (upper_bounds[c] < upper_bounds[current_winner])
                    current_winner = c;
            }


            // eliminate all loosers
			for (int c = 0; c < clusters.size(); c++) {
				if (inrace[c] && (lower_bounds[c] > upper_bounds[current_winner])) {
					inrace[c] = false;
					still_in_race--;
				}
			}

			if (still_in_race < 2){
				//if there is no cluster yet, stop
				//if only one is left, we have a winner
				break;
			}
		} // end for each datum


		// float mindist = Float.MAX_VALUE;
		mindist = Double.MAX_VALUE;
		
		switch (still_in_race) {
			case 0:
				// we are starting with the first data; create the first cluster
				// MOD : if # clusters = 0 then current_winner is set to -1 to avoid first cluster assignment to be considered an error
				if (clusters.size() == 0) current_winner = -1;
				break;
			case 1:
				// we have a clear winner ; use its distance; not used
				mindist = mean_values[current_winner];
				break;
			default:
				this.ambiguities ++;
				// if the race ended without a clear winner, use mean_values to decide
				for (int c = 0 ; c < clusters.size(); c++) {
					if (inrace[c] && (mean_values[c] < mindist)) {
						mindist = mean_values[c];
						current_winner = c;
					}
				}
				break;
		}
		
		return current_winner;
	}



	public List<List<Integer>> getClusters() {
		return clusters;
	}

	

	
	public int getNbClusters(){
		return clusters.size();
	}
	
	/** rand index computation */
	public double rand(int [] partition, List<String> theLabels){
        double res = 0d;
        double cpt = 0d;
        for (int i=0; i < this.data.size(); i++){
            for (int j=i+1; j < this.data.size(); j++){
                boolean b1 = (partition[i] == partition[j]);
                boolean b2 = theLabels.get(i).equalsIgnoreCase(theLabels.get(j));
                if ((b1 && b2) || (!b1 && !b2)) res ++;
                cpt ++;
            }
        }
        return (res/cpt);
    }
	
	/** di_blois.clustering building */
    public int[] createPartition(){
    	int[] part = new int[this.data.size()];
    	for (int c = 0; c < clusters.size(); c++){
    		for (int k = 0; k < clusters.get(c).size(); k++) part[clusters.get(c).get(k)] = c;
    	}
    	return part;
    }
	
    /** Returns all comparisons: template + racing comparisons */
	public long getTotalNbComparisons(){
		return this.comparisons + this.comptemplate;
	}
	
	/** Returns only racing comparisons */
	public long getRacingNbComparisons(){
		return this.comparisons;
	}
    
	/** Returns the number of ambiguous assignment */
	public long getAmbiguities(){
		return this.ambiguities;
	}
	
	/**
	 * Calcul de la matrice de contingence a partir de listes
	 * @param partition : partition decouverte
	 * @param partition2 : partition theorique
	 * @return
	 */
	private int[][] contingence(List<String> partition, int[] partition2){
		
		// Construction de dictionnaires pour les deux partitions :
		// au cas ou certains clusters seraient vides et pour traduire les string en Integer
		Map<String, Integer> dico1 = new HashMap<String, Integer>();
		Map<Integer, Integer> dico2 = new HashMap<Integer, Integer>();
		
		int count = 0;
		for (String value : partition){
			if (!dico1.containsKey(value)){
				dico1.put(value, count);
				count ++;
			}
		}
		count = 0;
		for (int value : partition2){
			if (!dico2.containsKey(value)){
				dico2.put(value, count);
				count ++;
			}
		}
		int nbclasses = dico1.size();
		int nbclusters = dico2.size();
		
		int[][] contingence = new int [nbclusters][nbclasses];
		
		for (int i = 0; i < partition.size(); i++){
			contingence[dico2.get(partition2[i])][dico1.get(partition.get(i))]++;
		}
		
		return contingence;
	}
	
	private static int max_row(int[] row){
		int max = row[0];
		for (int val : row)
			if (val > max) max = val;
		return max;
	}

	public List<List<D>> getPartition(){
		List<List<D>> res = new ArrayList<>(clusters.size());
		clusters.forEach(cluster -> {
			List<D> tmp = new ArrayList<>(cluster.size());
			res.add(tmp);
			cluster.forEach(index -> tmp.add(data.get(index)));
		});
		return res;
	}
	
	// Return 1/N * sum_{k in clusters} max_{c in classes} |k inter c|
	public double purity(List<String> partition, int[] partition2){
		double sum = 0;
		int[][] cont = contingence(partition, partition2);
		for (int[] tab : cont){
			sum += max_row(tab);
		}
		return sum / partition.size();
	}

}

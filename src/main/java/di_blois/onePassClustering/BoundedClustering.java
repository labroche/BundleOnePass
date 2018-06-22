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
		r = new RandPerm();
	}
	
	/** Main di_blois.clustering function: scan the di_blois.dataset and assigns data objects to clusters */
	public void cluster(List<D> dataset, Comparator<D> comp, double learning_sample, double size){
		
		
		// init variables
		this.data = dataset;
		int count = dataset.size();
		int[] order = r.randperm(dataset.size());
		clusters = new ArrayList<List<Integer>>();
		rand = new Random();
		this.comparisons = 0;
		this.comptemplate = 0;
		this.ambiguities = 0;
		
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
		
		// unused - get exact value range
		/*
		for (int i = 0 ; i < count ; i++) {
			for (int j = i+1 ; j < count ; j++) {
				d = 0;
				try {
					d = cp.compare(this.dat.get(i), this.dat.get(j));
					if (d > maxdistance) maxdistance = d;
					//eventually check for minimal distance between non-equal points also
					if (d < mindistance) mindistance = d;
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		*/
		
		// range value evaluation: as overestimation is not problematic, 
		// estimated range = 2 times previously evaluated range
		// avoid computing the exact range --> quadratic complexity
		
		range_value = 2* (maxdistance - mindistance);
		
		// cluster building
		
		for (int i = 0; i < dataset.size(); i++){
			// return nearest cluster and set global variable mindist
			
			int cluster = race(order[i], clusters, comp, 1.0);
			
			// at this point, mindist is set to the distance to 'cluster'
			
			// assignment to nearest cluster
			if (cluster == -1 || mindist > template){
				// build a new cluster
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
		double epsilon;

		int i;
		for (i = 0 ; i < maxsize ; i++) {
			for (int c = 0; c < clusters.size(); c++){
				
				List<Integer> clust = clusters.get(c);
				
				if (inrace[c] && (i < clust.size())) {
					
					//sample another data
					int d = clust.get(rand.nextInt(clust.size()));
					
					try {
						distances[c][i] = cp.compare(this.data.get(obj), this.data.get(d));
						
						mean_values[c] = empiricalMean(distances[c], i+1);
								
						if (bound_type == HOEFFDING_BOUND) {
							epsilon = epsilon_multiplier * hoeffdingsEpsilon(mean_values[c], i+1, error_probability, range_value);
							
						} else if (bound_type == BERNSTEIN_BOUND){
							//use bernstein bound
							double std = empiricalStd(distances[c], i+1, mean_values[c]);
							epsilon = epsilon_multiplier * bernsteinEpsilon(mean_values[c], std, i+1, error_probability, range_value);							
						} else {
							double std = unbiasedEmpiricalStd(distances[c], i+1, mean_values[c]);
							epsilon = epsilon_multiplier * studentsEpsilon(std, i+1, 1 - error_probability, range_value);
						}
								
								
						lower_bounds[c] = mean_values[c] - epsilon;
						upper_bounds[c] = mean_values[c] + epsilon;
						if (upper_bounds[c] < upper_bounds[current_winner]) current_winner = c;
					} catch (Exception e) {
						System.err.println("Error while comparing" + obj + " and " + d);
						e.printStackTrace();
					}
				}
			} // end for each cluster

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
	
	/**
	 * Computes the confidence interval for the given sample using Hoeffding's bound
	 * Main idea : P(|E_true - E_est| > epsilon) < 2e^{-2n*epsilon�/B�}
	 *    
	 * @param errorprobability - the probability that the real mean is outside 
	 * 								the bounds for this sample 
	 * @param value_range - the largest possible range of the observations
	 * @return a two-element array containing min and max 
	 */
	
	private double hoeffdingsEpsilon (double mean, int nelements, double errorprobability, double value_range) {
		return value_range*Math.sqrt(Math.log(2/errorprobability)/(2*nelements));
	}
	
	
	/**
	 * Computes the confidence interval for the given sample using Student's bound
	 * Main idea: the sample distribution function is normal with unknown variance
	 * Compared to a normal distribution bound, Student introduces a correction factor
	 *    
	 * @param errorprobability - the probability that the real mean is outside 
	 * 								the bounds for this sample 
	 * @param value_range - the largest possible range of the observations
	 * @return a two-element array containing min and max 
	 */
	private double studentsEpsilon (double std, int nelements, double significance, double value_range){
		double val = 0;
		
		if (nelements == 1) { 
			val = value_range; // there is only one point in the cluster: epsilon is set to the maximum possible (range)
		} else {
			TDistribution tDist = new TDistribution(nelements - 1);
			double conf = tDist.inverseCumulativeProbability((1.0 + significance) / 2);
		
			// the t-value has to be scaled with sigma = estimated_std_dev / sqrt(nelements) - see wikipedia Student t test
			val = conf * (std / Math.sqrt(nelements));
		}
		
		// System.out.println("Student (" + nelements + "): " + val);
		
		return val;
	}
	
	
	/**
	 * empirical Bernstein bound
	 * if the possible range is large w.r.t. the sample's standard deviation, the resulting
	 * bound is tighter than Hoeffding's
	 * @param errorprobability - the probability that the real mean is outside 
	 * 								the bounds for this sample 
	 * @param value_range - the largest possible range of the observations
	 * @return  
	 */
	private double bernsteinEpsilon (double mean, double std, int nelements, double errorprobability, double value_range) {
		double log3errorprob = Math.log(3/errorprobability);
		double epsilon = std*Math.sqrt(
				(2*log3errorprob)/nelements) + 3*value_range*log3errorprob/nelements;
		return epsilon;
	}
	
	/**
	 * computes the empirical mean of the current distance distribution
	 * @return
	 */
	private static double empiricalMean (double[] d, int nelements) {
		double mean = 0f;
		// incremental version possible ;)
		for (int i = 0 ; i < nelements ; i++) {
			mean += d[i];
		}
		return mean/nelements;
	}

	/**
	 * computes the empirical standard deviation of the current distance distribution
	 * @param mean
	 * @return
	 * warning: when there is 0 elements, std = 0 and when there is one element std = 0?
	 */
	private static double empiricalStd (double[] d, int nelements, double mean) {
		double tmp;
		double std = 0d;
		for (int i = 0 ; i < nelements ; i++) {
			tmp = d[i] - mean;
			std += tmp*tmp; 
		}
		return Math.sqrt(std/nelements);
	}
	
	/**
	 * computes the unbiased empirical standard deviation of the current distance distribution
	 * needed for Student's bound computation
	 * @param mean
	 * @return
	 * warning: when there is 1 or less element, std = 0
	 */
	private static double unbiasedEmpiricalStd (double[] d, int nelements, double mean) {
		double tmp;
		double std = 0d;
		for (int i = 0 ; i < nelements ; i++) {
			tmp = d[i] - mean;
			std += tmp*tmp; 
		}
		return Math.sqrt(std/(nelements - 1));
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
	
	// Return 1/N * sum_{k in clusters} max_{c in classes} |k inter c|
	public double purity(List<String> partition, int[] partition2){
		double sum = 0;
		int[][] cont = contingence(partition, partition2);
		for (int[] tab : cont){
			sum += max_row(tab);
		}
		return sum / partition.size();
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int nbrun = 100;
		int bound_type = BoundedClustering.BERNSTEIN_BOUND;
		
		// parsing arguments
		switch (args.length) {
		case 0:
			break;
		case 1:
			if (args[0].equals("B")) bound_type = BoundedClustering.BERNSTEIN_BOUND;
			else if (args[0].equals("H")) bound_type = BoundedClustering.HOEFFDING_BOUND;
			else bound_type = BoundedClustering.STUDENT_BOUND;
			break;
		case 2:
			if (args[0].equals("B")) bound_type = BoundedClustering.BERNSTEIN_BOUND;
			else if (args[0].equals("H")) bound_type = BoundedClustering.HOEFFDING_BOUND;
			else bound_type = BoundedClustering.STUDENT_BOUND;
			try{
				Integer iter = new Integer(args[1]);
				nbrun = iter.intValue();
			} catch (Exception e) {
				e.printStackTrace();
			}
			break;
		default: break;
		}
		
		String path = "../dat/data/";
		String outputpath = "../res/";
		String filename = "liste_alea2.txt";

		try{
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(path + filename)));
			String line;

			Minkowski cp = new Minkowski(2);
			
			
			BoundedClustering<float[]> algo;
			IReader<float[]> sdr = new SimpleDataReader();
			Data<float[]> dat = null;
			
			DateFormat dateFormat = new SimpleDateFormat("ddMMyyyy_HHmmss");
			Date date = new Date();
			String outputfilename = dateFormat.format(date);
			
			// BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputpath + outputfilename + ".txt")));
			// FileWriter out = new FileWriter(outputfilename);
			// out.close();
			
			while ((line = br.readLine()) != null){
				if (line.startsWith("#")) continue;
				String[] tmp = line.split(" ");
				dat = sdr.read(path + tmp[0]);
				double[] rand = new double[nbrun];
				double[] pur = new double[nbrun];
				double[] time = new double[nbrun];
				double[] clusts = new  double[nbrun];
				double[] comp = new double[nbrun];
				double[] amb = new double[nbrun];
				
				long debut, fin;
				
				for (int test = 0; test < nbrun; test ++){
					debut = System.currentTimeMillis();
					algo = new BoundedClustering<float[]>(bound_type);
					// params: List<D> di_blois.dataset, Comparator<D> comp, double learning_sample, double size
					algo.cluster(dat.getData(), cp, 0.1, 0d);
					fin = System.currentTimeMillis();
					rand[test] = algo.rand(algo.createPartition(), dat.getLabels());
					pur[test] = algo.purity(dat.getLabels(), algo.createPartition());
					time[test] = fin - debut;
					clusts[test] = algo.getNbClusters();
					comp[test] = algo.getTotalNbComparisons();
					amb[test] = algo.getAmbiguities();
				}
				
				double mean = empiricalMean(rand, nbrun);
				String s = tmp[0] + ";rand;" + mean + ";" + empiricalStd(rand, nbrun, mean) + ";time;";
				mean = empiricalMean(time, nbrun);
				s += mean + ";" + empiricalStd(time, nbrun, mean) + ";clusters;";
				mean = empiricalMean(clusts, nbrun);
				s += mean + ";" + empiricalStd(clusts, nbrun, mean) + "\n";
				System.out.print(s);
				FileWriter fw = new FileWriter(outputpath + outputfilename, true);
				fw.write(s);
				fw.close();
				// bw.write(s);
			}
			// bw.close();
		} catch (Exception e){e.printStackTrace();}

	}

}

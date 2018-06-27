package di_blois.onePassClustering;

import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import static di_blois.onePassClustering.BoundedClustering.BERNSTEIN_BOUND;
import static di_blois.onePassClustering.BoundedClustering.HOEFFDING_BOUND;
import static di_blois.onePassClustering.Bounds.*;

public class BopWorker<D> extends Worker<CompareTask<D>>{

    volatile boolean run = true;
    Random rd;

    @Override
    public void run() {
        rd = ThreadLocalRandom.current();
        while (run){
            CompareTask currentTask = null;
            try {
                currentTask = provider.take();
            } catch (InterruptedException e) {
                continue;
            }

            double epsilon;

            if (currentTask.inrace[currentTask.c] && (currentTask.i < currentTask.cluster.size())) {

                //sample another data
                int d = (int) currentTask.cluster.get(rd.nextInt(currentTask.cluster.size()));

                try {
                    currentTask.distances[currentTask.c][currentTask.i] = currentTask.cp.compare(currentTask.dataset.get(currentTask.obj), currentTask.dataset.get(d));
                    currentTask.mean_values[currentTask.c] = empiricalMean(currentTask.distances[currentTask.c], currentTask.i +1);

                    if (currentTask.bound == HOEFFDING_BOUND) {
                        epsilon = currentTask.epsilon_multiplier * hoeffdingsEpsilon(currentTask.mean_values[currentTask.c], currentTask.i +1, currentTask.error_probability, currentTask.range_value);

                    } else if (currentTask.bound == BERNSTEIN_BOUND){
                        //use bernstein bound
                        double std = empiricalStd(currentTask.distances[currentTask.c], currentTask.i +1, currentTask.mean_values[currentTask.c]);
                        epsilon = currentTask.epsilon_multiplier * bernsteinEpsilon(currentTask.mean_values[currentTask.c], std, currentTask.i +1, currentTask.error_probability, currentTask.range_value);
                    } else {
                        double std = unbiasedEmpiricalStd(currentTask.distances[currentTask.c], currentTask.i +1, currentTask.mean_values[currentTask.c]);
                        epsilon = currentTask.epsilon_multiplier * studentsEpsilon(std, currentTask.i +1, 1 - currentTask.error_probability, currentTask.range_value);
                    }


                    currentTask.lower_bounds[currentTask.c] = currentTask.mean_values[currentTask.c] - epsilon;
                    currentTask.upper_bounds[currentTask.c] = currentTask.mean_values[currentTask.c] + epsilon;


                } catch (Exception e) {
                    System.err.println("Error while comparing" + currentTask.obj + " and " + d);
                    e.printStackTrace();
                }
            }
        }

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
}

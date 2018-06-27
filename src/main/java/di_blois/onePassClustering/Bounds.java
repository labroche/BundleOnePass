package di_blois.onePassClustering;

import org.apache.commons.math3.distribution.TDistribution;

public class Bounds {
    /**
     * Computes the confidence interval for the given sample using Hoeffding's bound
     * Main idea : P(|E_true - E_est| > epsilon) < 2e^{-2n*epsilon�/B�}
     *
     * @param errorprobability - the probability that the real mean is outside
     * 								the bounds for this sample
     * @param value_range - the largest possible range of the observations
     * @return a two-element array containing min and max
     */

    public static double hoeffdingsEpsilon (double mean, int nelements, double errorprobability, double value_range) {
        return value_range*Math.sqrt(Math.log(2/errorprobability)/(2*nelements));
        // 2/errorprobability)/(2*nelements == 1/(errorprobability*nelement)
    }


    /**
     * Computes the confidence interval for the given sample using Student's bound
     * Main idea: the sample distribution function is normal with unknown variance
     * Compared to a normal distribution bound, Student introduces a correction factor
     *
     * @param value_range - the largest possible range of the observations
     * @return a two-element array containing min and max
     */
    public static double studentsEpsilon (double std, int nelements, double significance, double value_range){
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
    public static double bernsteinEpsilon (double mean, double std, int nelements, double errorprobability, double value_range) {
        double log3errorprob = Math.log(3/errorprobability);
        double epsilon = std*Math.sqrt(
                (2*log3errorprob)/nelements) + 3*value_range*log3errorprob/nelements;
        return epsilon;
    }

}

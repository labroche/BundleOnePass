package di_blois.onePassClustering;

import di_blois.comparator.Comparator;

import java.util.List;

public class CompareTask<D> {
    public final int bound, i, c, obj;
    public final List<Integer> cluster;
    public final boolean[] inrace;
    public final double[][] distances;
    public final double[] mean_values, upper_bounds, lower_bounds;
    public final Comparator<D> cp;
    public final List<D> dataset;
    public final double error_probability, epsilon_multiplier, range_value;

    public CompareTask(int bound, int i, int c, int obj, List<Integer> cluster, boolean[] inrace, double[][] distances, double[] mean_values, double[] upper_bounds, double[] lower_bounds, Comparator<D> cp, List<D> dataset, double error_probability, double epsilon_multiplier, double range_value) {
        this.bound = bound;
        this.i = i;
        this.c = c;
        this.obj = obj;
        this.cluster = cluster;
        this.inrace = inrace;
        this.distances = distances;
        this.mean_values = mean_values;
        this.upper_bounds = upper_bounds;
        this.lower_bounds = lower_bounds;
        this.cp = cp;
        this.dataset = dataset;
        this.error_probability = error_probability;
        this.epsilon_multiplier = epsilon_multiplier;
        this.range_value = range_value;
    }
}

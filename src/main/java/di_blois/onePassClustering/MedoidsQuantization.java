package di_blois.onePassClustering;

import di_blois.comparator.Comparator;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

public class MedoidsQuantization<D> {

    Comparator<D> metric;
    List<D> medoids;
    List<List<D>> partition;
    double minSep;

    double[] avgSep = new double[100]; // For each medoid average distance between itself and all other medoids
    double[] distance = new double[100];

    public MedoidsQuantization(Comparator<D> metric, double minSep) {
        this.metric = metric;
        this.minSep = minSep;
        medoids = new ArrayList<>(10);
        partition = new ArrayList<>();
    }

    public void cluster(Stream<D> data){
        Iterator<D> it = data.iterator();

        newCulster(it.next());

        it.forEachRemaining(d -> {
            for (int i = 0; i < medoids.size(); i++) {
                distance[i] = metric.compare(medoids.get(i), d);
            }

            int chosen = minI(distance, medoids.size());
            double chosenValue = distance[chosen];
            if (chosenValue < minSep){
                partition.get(chosen).add(d);
                updateMedoid(chosen, d);
            }else {
                newCulster(d);
            }


        });
    }

    private void updateMedoid(int chosen, D point) {
        if (medoids.size() == 1)
            return;
        double sepForP = avg(distance, medoids.size(), chosen);
        if (sepForP > avgSep[chosen]) {
            setMedoid(chosen, point);
        }
    }

    private void setMedoid(int cluster, D point) {
        medoids.set(cluster, point);
        updateAvgSep();
    }

    private void updateAvgSep() {
        if (medoids.size() > 1){
            for (int i = 0; i < medoids.size(); i++) {
                double sum = 0;
                for (int j = 0; j < medoids.size(); j++) {
                    if (j != i)
                        sum += metric.compare(medoids.get(i), medoids.get(j));
                }

                avgSep[i] = sum/(medoids.size());
            }
        }
        else {
            return;
            //We don't update medoids until we have at least 2 clusters
        }

    }

    public List<List<D>> getPartition() {
        return partition;
    }

    private void newCulster(D medoid){
        partition.add(new ArrayList<>(500));
        partition.get(partition.size() - 1).add(medoid);
        medoids.add(medoid);
        if (medoids.size()>= distance.length){
            double[] tmp = new double[distance.length + 100];
            System.arraycopy(distance, 0, tmp, 0, distance.length);
            distance = tmp;

            tmp = new double[avgSep.length + 100];
            System.arraycopy(avgSep, 0, tmp, 0, avgSep.length);
            avgSep = tmp;
        }
        updateAvgSep();
    }

    private int minI(double[] a, int bound){
        int j = 0;
        double best = a[0];
        for (int i = 0; i < bound; i++) {
            if (a[i] < best){
                best = a[i];
                j = i;
            }
        }
        return j;
    }

    private double avg(double[] a, int bound, int except){
        double sum = 0;
        for (int i = 0; i < bound; i++) {
            if (i != except)
                sum += a[i];
        }
        return sum/(bound - 1);
    }
}

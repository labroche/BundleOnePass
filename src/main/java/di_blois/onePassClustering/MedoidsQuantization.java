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
    }

    public void cluster(Stream<D> data){
        Iterator<D> it = data.iterator();

        newCulster(it.next());
        partition.add(new ArrayList<>(500));

        it.forEachRemaining(d -> {
            for (int i = 0; i < medoids.size(); i++) {
                try {
                    distance[i] = metric.compare(medoids.get(i), d);
                } catch (Exception e) {
                    distance[i] = Double.MAX_VALUE;
                    e.printStackTrace(); //TODO Create garbage stack for those points
                }
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
                    if (j != i) {
                        try {
                            sum += metric.compare(medoids.get(i), medoids.get(j));
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }

                avgSep[i] = sum/(medoids.size());
            }
        }
        else {
            //TODO what do we do if only one cluster !!!
        }

    }

    public List<List<D>> getPartition() {
        return partition;
    }

    private void newCulster(D medoid){
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

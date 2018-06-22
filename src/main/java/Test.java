
import di_blois.clustering.FCM;
import di_blois.comparator.L1;
import di_blois.comparator.Minkowski;
import di_blois.onePassClustering.BOPBundles;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import utilities.misc.Nd4jUtils;
import utilities.misc.Pair;

import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class Test {
    public static void main(String[] args) throws Exception {

        List<float[]> dataset = Nd4jUtils.loadCSVAsListNative(args[0], ";", 0, false); //60 sec - 85 sec for ND4J ...
        List<String> ids = Nd4jUtils.loadIds(args[0], ";", 0, "spotify:track:");

        List<Pair<String, float[]>> items = new ArrayList<>(dataset.size());
        for (int i = 0; i < dataset.size(); i++) {
            items.add(new Pair<>(ids.get(i), dataset.get(i)));
        }
        ids.clear();
        dataset.clear();
        System.out.println("Data loaded");

        FCM fcm = new FCM(7, 1.2);
        INDArray tmp = Nd4j.create(10000, items.get(0).right.length);
        for (int i = 0; i < 10000; i++) {
            tmp.putRow(i,Nd4j.create(items.get(i).right));
        }
        System.out.println("Begin bootstrap");
        fcm.fit(tmp);

        List<List<Pair<String, float[]>>> partition = new ArrayList<>();
        int[] labels = fcm.getLabels();
        IntSummaryStatistics stat = Arrays.stream(labels).summaryStatistics();
        for (int i = 0; i < stat.getMax() + 1; i++) {
            partition.add(new ArrayList<>());
        }

        for (int i = 0; i < labels.length; i++) {
            partition.get(labels[i]).add(items.get(i));
        }
        System.out.println("Bootstrap Done");

        L1 dis = new L1();
        double max_radius = 0;
        for (List<Pair<String, float[]>> cluster : partition){
            for (int i = 0; i < cluster.size(); i++) {
                for (int j = i; j < cluster.size(); j++) {
                    double d = dis.compare(cluster.get(i), cluster.get(j));
                    if (max_radius < d)
                        max_radius = d;
                }
            }
        }


        BOPBundles<Pair<String, float[]>> bop = new BOPBundles<>(BOPBundles.HOEFFDING_BOUND);
        bop.bootstrap(partition);
        System.out.println("BOP begins");
        long tstart = System.currentTimeMillis();
        bop.cluster(items.stream(), new L1(), 0.05, 10000);
        long tend = System.currentTimeMillis();
        System.out.println("BOP done");
        System.out.println("Exec Time was " + (tend - tstart)/1000 + " seconds");
        System.out.println("Found " + bop.getNbClusters() + " clusters.");
        System.out.println(bop.getTotalNbComparisons());

        Nd4jUtils.writePartition(bop.getPartition(), args[1]);
    }


}

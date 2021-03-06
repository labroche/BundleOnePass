import di_blois.clustering.FCM;
import di_blois.comparator.L1;
import di_blois.onePassClustering.BOPBundles;
import di_blois.onePassClustering.BoundedClustering;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import utilities.misc.Nd4jUtils;
import utilities.misc.Pair;

import java.util.*;
import java.util.stream.Stream;

public class TestOriginal {
    public static void main(String[] args) throws Exception{


        List<float[]> dataset = Nd4jUtils.loadCSVAsListNative(args[0], ";", 0, false);
        List<String> ids = Nd4jUtils.loadIds(args[0], ";", 0, "spotify:track:");

        List<Pair<String, float[]>> items = new ArrayList<>(dataset.size());
        for (int i = 0; i < dataset.size(); i++) {
            items.add(new Pair<>(ids.get(i), dataset.get(i)));
        }
        ids.clear();
        dataset.clear();
        System.out.println("Data loaded");

        BoundedClustering<Pair<String, float[]>> bop = new BoundedClustering<>(0);
        System.out.println("BOP begins");
        long tstart = System.currentTimeMillis();
        bop.cluster(items.subList(0, 15000), new L1(), 0.10, 0.01);
        long tend = System.currentTimeMillis();
        System.out.println("BOP done");
        System.out.println("Exec Time was " + (tend - tstart)/1000 + " seconds");
        System.out.println("Found " + bop.getNbClusters() + " clusters.");
        System.out.println(bop.getTotalNbComparisons());
        System.out.println("Cluster size:");
        bop.getPartition().forEach(cluster -> System.out.println(cluster.size()));
        Nd4jUtils.writePartition(bop.getPartition(), args[1]);


    }
}

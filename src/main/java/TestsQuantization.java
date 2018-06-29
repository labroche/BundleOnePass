import di_blois.comparator.L1;
import di_blois.onePassClustering.BoundedClustering;
import di_blois.onePassClustering.MedoidsQuantization;
import utilities.misc.Nd4jUtils;
import utilities.misc.Pair;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TestsQuantization {




    public static void main(String[] args) throws FileNotFoundException {
        List<float[]> dataset = Nd4jUtils.loadCSVAsListNative("/home/alex/stage/lda_tracks.csv", " ", 0, false);
        List<String> ids = Nd4jUtils.loadIds("/home/alex/stage/lda_tracks.csv", " ", 0, "");

        List<Pair<String, float[]>> items = new ArrayList<>(dataset.size());
        for (int i = 0; i < dataset.size(); i++) {
            //float[] tmp = new float[dataset.get(i).length - 1];
            //System.arraycopy(dataset.get(i), 1, tmp, 0, dataset.get(i).length - 1);
            items.add(new Pair<>(ids.get(i), dataset.get(i)));
            //System.err.println(ids.get(i)+"|"+ Arrays.toString(tmp));
        }
        ids.clear();
        dataset.clear();
        System.out.println("Data loaded");


        MedoidsQuantization<Pair<String, float[]>> mq = new MedoidsQuantization<>(new L1(), 0.17);
        mq.cluster(items.stream());

        System.out.printf("Found %d clusters.%n", mq.getPartition().size());

        Nd4jUtils.writePartition(mq.getPartition(), "/home/alex/stage/full_partition.txt");


    }
}

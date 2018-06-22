package di_blois.clustering;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.api.ops.impl.accum.AMax;
import org.nd4j.linalg.api.ops.impl.indexaccum.IMax;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.BiFunction;

import static java.lang.Double.NaN;
import static java.lang.Math.pow;

/**
 * Implementation of FCM (Fuzzy C-Means) Using the ND4J library
 */
public class FCM<D>{
    public static final int THREADS = 4;
    protected int k, max_iter = 100;
    protected double epsilon = 1e-3, m;
    protected INDArray data, membership, centers, distance, membership_old;
    protected int[] shape;
    protected Random rd;
    protected final double p;

    public FCM(int k, double m) {
        this.k = k;
        this.m = m;
        rd = new Random();
        p = 1 / (m - 1);
    }



    public INDArray fit(INDArray input){
        data = input;
        shape = new int[]{data.rows(), data.columns()};
        membership = new NDArray(shape[0], k);
        centers = new NDArray(k, shape[1]);

        //Randomly picking starting centers near data points
        randomInitCenters();
        updateMembership();
        membership_old = membership;

        int iter;
        for (iter = 0; iter < max_iter; ++iter){
            updateCenters();
            updateMembership();
            if (checkConvergence())
                break;
            membership_old = membership;
        }
        return membership;
    }

    private void updateCenters() {
        centers = Nd4j.zeros(k, shape[1]);
        ExecutorService executor = Executors.newCachedThreadPool();
        for(int j = 0; j < k; ++j){
            int finalJ = j;
            executor.execute(() -> {
                double sum = 0;
                for (int i = 0; i < shape[0]; i++) {
                    double p = pow(membership.getDouble(i, finalJ), m);
                    centers.getRow(finalJ).addi(data.getRow(i).mul(p));//Using addi will add 'In Place'
                    sum += p;
                }
                centers.getRow(finalJ).divi(sum);//same here we divide in place
            });
        }
        executor.shutdown();
        try {executor.awaitTermination(10, TimeUnit.MINUTES);}
        catch (InterruptedException e) {e.printStackTrace();}
    }

    protected void updateMembership() {

        /*
        distance = Nd4j.create(shape[0], k);
        for (int i = 0; i < shape[0]; ++i){
            INDArray point = data.getRow(i);
            for (int j = 0; j < k; ++j)
                distance.putScalar(i, j, point.squaredDistance(centers.getRow(j)) );
        }
        */
        distance = Nd4j.create(shape[0], k);
        ExecutorService executor = Executors.newCachedThreadPool();
        for (int j = 0; j < k; j++) {
            int finalJ = j;
            executor.execute(() -> {
                for (int i = 0; i < shape[0]; i++) {
                    distance.putScalar(i, finalJ, data.getRow(i).squaredDistance(centers.getRow(finalJ)) );
                }
            });
        }
        executor.shutdown();
        try {executor.awaitTermination(10, TimeUnit.MINUTES);}
        catch (InterruptedException e) {e.printStackTrace();}

        membership = Nd4j.zeros(shape[0], k);


            for (int i = 0; i < shape[0]; i++) {
                innerLoop(i);
            }


    }

    private void innerLoop(int i) {
        int warning = -1;
        for (int j = 0; j < k; j++) {
            double inner = innerTerm(i, j);
            if (inner == NaN) {
                warning = j;
                break;
            }
            membership.putScalar(i, j, inner);
        }
        if (warning != -1){
            for (int j = 0; j < k; j++) {
                if (j == warning){
                    membership.putScalar(i, j, 1.0);
                }else {
                    membership.putScalar(i, j, 0.0);
                }
            }
        }
    }

    private double innerTerm(int i, int j){
        double sum = 0, num = 0;
        for (int l = 0; l < k; l++) {
            double dis = distance.getDouble(i, l);
            if (dis == 0)
                return NaN;
            double t = pow(1/dis, p);
            sum += t;
            if (l == j)
                num = t;
        }
        return num/sum;
    }

    /**
     * We check if the bigest change in the membership matrix is below our epsilon threshold
     * @return Should we stop the algorithm
     */
    protected boolean checkConvergence(){
        INDArray diff = membership.sub(membership_old);
        double maxDiff = Nd4j.getExecutioner().execAndReturn(new AMax(diff)).getFinalResult().doubleValue();
        return  maxDiff < epsilon;
    }

    protected void randomInitCenters() {
        //We choose different starting points
        Set<Integer> indexes = new HashSet<>();
        while (indexes.size() < k)
            indexes.add(rd.nextInt(shape[0] - 1));
        int i = 0;
        for (Integer index : indexes){
            //We should NOT start on the point itself
            centers.putRow(i, data.getRow(index).add(data.getRow(index).mul(rd.nextDouble()-0.5)));
            ++i;
        }
    }

    public void setMax_iter(int max_iter) {
        this.max_iter = max_iter;
    }

    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }

    /**
     * Returns the labels found by the algorithm for each point in the data set (ie which cluster they are in the most)
     * @return label for each data point
     */
    public int[] getLabels(){
        int[] tmp = new int[shape[0]];
        for (int i = 0; i < shape[0]; i++) {
            tmp[i] = Nd4j.getExecutioner().execAndReturn(new IMax(membership.getRow(i))).getFinalResult();
        }
        return tmp;
    }

    public INDArray getCenters() {
        return centers;
    }


    class Range{
        final int left, right;

        public Range(int left, int right) {
            this.left = left;
            this.right = right;
        }
    }
}

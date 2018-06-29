package di_blois.comparator;

import utilities.misc.Pair;

/**
 * @author Alexandre
 * L1 norm, no check to go faster
 */
public class L1 implements Comparator<Pair<String,float[]>>{
    @Override
    public double compare(Pair<String,float[]> dat1, Pair<String,float[]> dat2) {
        float sum = 0;
        for (int i = 0; i < dat1.right.length; i++) {
            sum += Math.abs(dat1.right[i] - dat2.right[i]);
        }
        return sum/dat1.right.length;
    }
}

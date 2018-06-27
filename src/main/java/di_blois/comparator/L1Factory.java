package di_blois.comparator;

import utilities.misc.Pair;

public class L1Factory implements ComparatorFactory<Pair<String,float[]>> {
    @Override
    public Comparator<Pair<String, float[]>> build() {
        return new L1();
    }
}

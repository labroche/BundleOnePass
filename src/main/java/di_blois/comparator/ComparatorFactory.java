package di_blois.comparator;

public interface ComparatorFactory<D> {
    Comparator<D> build();
}

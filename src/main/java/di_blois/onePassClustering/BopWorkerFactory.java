package di_blois.onePassClustering;

public class BopWorkerFactory<D> implements WorkerFactory<Worker<CompareTask<D>>> {
    @Override
    public Worker<CompareTask<D>> build() {
        return new BopWorker();
    }
}

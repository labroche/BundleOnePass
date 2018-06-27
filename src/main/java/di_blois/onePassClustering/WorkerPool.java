package di_blois.onePassClustering;

import java.util.ArrayList;
import java.util.concurrent.LinkedBlockingQueue;

public class WorkerPool<W extends Worker<T>, T> {
    int workerNumber;
    LinkedBlockingQueue<T> queue;
    ArrayList<W> workers;


    public WorkerPool(int workerNumber, WorkerFactory<W> factory) {
        this.workerNumber = workerNumber;
        queue = new LinkedBlockingQueue<>();
        workers = new ArrayList<>(workerNumber);
        for (int i = 0; i < workerNumber; i++) {
            W w = factory.build();
            w.registerJobQueue(queue);
            w.setName("Worker_" + i);
            workers.add(w);
        }

        workers.forEach(Thread::start);

    }

    public void submit(T task){
        queue.offer(task);
    }

    public boolean isDone(){
        return queue.isEmpty();
    }
}

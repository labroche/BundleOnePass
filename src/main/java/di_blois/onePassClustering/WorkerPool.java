package di_blois.onePassClustering;

import java.util.ArrayList;
import java.util.concurrent.LinkedBlockingQueue;

public class WorkerPool<W extends Thread, T> {
    int workerNumber;
    LinkedBlockingQueue<T> queue;
    ArrayList<W> workers;



}

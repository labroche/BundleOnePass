package di_blois.onePassClustering;


import java.util.concurrent.LinkedBlockingQueue;

public abstract class Worker<T> extends Thread{
    public LinkedBlockingQueue<T> provider;

    public void registerJobQueue(LinkedBlockingQueue<T> queue){
        provider = queue;
    }
}

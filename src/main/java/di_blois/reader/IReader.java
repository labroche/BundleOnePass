package di_blois.reader;

import di_blois.dataset.Data;

public interface IReader<D> {
	public Data<D> read(String filename);
	public void open(String filename);
	public void close();
	public Data<D> next(int n);
}

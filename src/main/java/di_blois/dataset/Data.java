package di_blois.dataset;

import java.util.ArrayList;
import java.util.List;

// TODO : il manque la gestion des labels vides ou inexistants
public class Data<D> {
	
	List<D> dat;
	List<String> labels;
	int count, natt;
	
	public Data(List<D> d, List<String> l){
		this.dat = d;
		this.labels = l;
		this.count = d.size();
		// natt ne peut pas etre fixe
	}
	
	
	public int size(){return count;}	// # donnees
    public int dim(){return natt;}		// # attributs
    
    public List<String> getLabels(){return this.labels;}	// Recuperation des etiquettes ou null
    public List<D> getData(){return this.dat;}				// Recuperation des donnees

    public D get(int index){
    	if (index < dat.size())
    		return dat.get(index); 
    	else return null;
    }
    
    public void set(int index, D value){
    	if (index < dat.size())
    		this.dat.set(index, value);
    	else this.dat.add(index, value);
    }
    
    public void setLabel(int index, String value){
    	if (index < labels.size())
    		this.labels.set(index, value);
    	else this.labels.add(index, value);
    }
    
    public void addAllData(Data<D> b){
    	this.dat.addAll(b.getData());
    	this.labels.addAll(b.getLabels());
    	this.count = this.count + b.count;
    }
    
    public void addSomeData(Data<D> b, int[] index){
    	for (int ind : index){
    		this.dat.add(b.get(ind));
        	this.labels.add(b.getLabels().get(ind));
    	}
    	this.count = this.count + index.length;
    }
    
    public void addSomeData(Data<D> b, int start){
    	this.dat.addAll(b.getData().subList(start, b.getData().size() - 1));
    	this.labels.addAll(b.getLabels().subList(start, b.getData().size() - 1));
    	this.count = this.count + b.getData().size() - start;	
    }
    
    public Data<D> getSubDataSet(int[] index){
    	List<D> tmpdat = new ArrayList<D>();
    	List<String> tmplab = new ArrayList<String>();
    	for (int ind : index){
    		tmpdat.add(this.dat.get(ind));
        	try{
        		tmplab.add(this.getLabels().get(ind));
        	} catch (Exception e)
        	{
        		// Rien on renverra une liste vide dans ce cas
        		// TODO : verifier que cela fonctionne
        	}
    	}
    	return new Data<D>(tmpdat, tmplab);
    }
    
    public Data<D> getSubDataSet(int start, int length){
    	List<D> tmpdat = this.dat.subList(start, start + length);
    	try{
    		List<String> tmplab = this.labels.subList(start, start + length);
    		return new Data<D>(tmpdat, tmplab);
    	} catch (Exception e){
    		return new Data<D>(tmpdat, new ArrayList<String>());
    	}
    	
    }
}

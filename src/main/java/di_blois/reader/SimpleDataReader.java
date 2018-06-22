package di_blois.reader;

import di_blois.dataset.Data;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class SimpleDataReader implements IReader<float[]> {
	
	final String sep = " ";
	
	List<float[]> data;			// Liste des donnees codees sous forme de listes
	List<String> labels;		// Liste des labels
	
	BufferedReader br;			// Acces au flux du fichier
	int count;					// Compteur du nombre de lignes
	String line;				// Derniere ligne parcourue
	
	/**
	 * Fermeture du flux de donnees
	 */
	public void close() {
		try{
			if (br != null) br.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public Data<float[]> next(int n) {
		this.data.clear();
		this.labels.clear();
		try {
			int nb = 0; 						// Nombre de lignes lues
			boolean fini = false;
			while (nb < n && !fini){
				line = this.br.readLine();
				if (line !=null){
					traitementLigne(line, sep);
					count ++;
					nb++;
				} else {
					fini = true;
				}
			}
			return new Data<float[]>(this.data, this.labels);
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * Ouverture du flux et comptage du nombre de donnees
	 */
	public void open(String filename) {
		try{
			this.count = 0;
			this.data = new ArrayList<float[]>();
			this.labels = new ArrayList<String>();
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
			while ((line = br.readLine())!=null){
	            count ++;
	        }
			br.close();
			// Remise du flux au debut
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
		} catch (Exception e) {e.printStackTrace();}
	}

	public Data<float[]> read(String filename) {
		try{
			this.count = 0;
			this.data = new ArrayList<float[]>();
			this.labels = new ArrayList<String>();
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
			
			while ((line = br.readLine())!=null){
				traitementLigne(line, sep);
	            count ++;
	        }
			return new Data<float[]>(this.data, this.labels);
		} catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Lecture d'une ligne d'un fichier "data"
	 * Hypotheses :
	 * - Les valeurs sont separees par le caractere sep
	 * la derniere valeur est le label de la classe
	 */
	private void traitementLigne(String line, String sep){
		String[] tmp = line.split(sep);
		int i;
		float[] vect = new float[tmp.length - 1];
		for (i = 0; i < tmp.length - 1; i++){
			vect[i] = new Float(tmp[i]).floatValue();
		}
		
		this.data.add(vect);
		this.labels.add(tmp[tmp.length - 1]);
	}
}

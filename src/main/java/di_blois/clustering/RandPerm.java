package di_blois.clustering;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class RandPerm extends Random {
	private long seed;
	
	/**
	 * generated serialVersionUID
	 */
	private static final long serialVersionUID = 5804666636792594861L;

	public RandPerm () {
		seed = System.nanoTime();
	}
	
	public RandPerm (long seed) {
		this.seed = seed;
	}
	
	protected int next(int nbits){
		// Not thread-safe!
		long x = this.seed;
	    x ^= (x << 21);
	    x ^= (x >>> 35);
	    x ^= (x << 4);
	    this.seed = x;
	    x &= ((1L << nbits) -1);
	    return (int) x;
	}
	
/*	public long nextLong(){
		// Not thread-safe!
		long x = this.seed;
	    x ^= (x << 21);
	    x ^= (x >>> 35);
	    x ^= (x << 4);
	    this.seed = x;
	    return x;
	}
*/	
	/**
	 * returns a random permutation using a uniform distribution
	 * of length integers from 0 to length-1
	 * @param length
	 */
	public int[] randperm (int length) {
		int[] arr = new int[length];
		for (int i=1;i<length;i++){
			arr[i]=i;
		}
		shuffle (arr);
		return arr;
	}
	
	
	public List<Integer> randpermL (int length) {
		List<Integer> arr = new ArrayList<Integer>(length);
		for (int i=0;i<length;i++){
			arr.add(i);
		}
		Collections.shuffle(arr);
		return arr;
	}
	
	/**
	 * shuffles a boolean array in random order
	 * all permutations are equally probable
	 * @param arr
	 */
	public void shuffle (boolean[] arr) {
		boolean tmp;
		int j;
		for (int i=1;i<arr.length;i++){
			j = nextInt (i+1);
			tmp = arr[i];
			arr[i]=arr[j];
			arr[j]=tmp;
		}		
	}

	/**
	 * shuffles an integer array in random order
	 * all permutations are equally probable
	 * @param arr
	 */
	public void shuffle (int[] arr) {
		int tmp,j;
		for (int i=1;i<arr.length;i++){
			j = nextInt (i+1);
			tmp = arr[i];
			arr[i]=arr[j];
			arr[j]=tmp;
		}		
	}
	
	/**
	 * shuffles a block inside an array in random order
	 * @param arr 			the array
	 * @param fromindex		the index to start with (included)
	 * @param toindex		the index to end with (included)
	 */
	public void shuffle (int[] arr, int fromindex, int toindex) {
		int tmp,j;
		for (int i=1;i<toindex-fromindex+1;i++) {
			j = nextInt (i+1);
			tmp = arr[fromindex+i];
			arr[fromindex+i] = arr[fromindex+j];
			arr[fromindex+j] = tmp;
		}
	}
	
	public void uniformScramble (int[] arr, int param) {
		int[] indizes = new int[arr.length];
		int curpos = 0;
		for (int i=0;i<arr.length;i++) {
			if (nextInt(1000)<param)
				indizes[curpos++] = i;
		}
		int tmp,j;
		for (int i=1;i<curpos;i++) {
			j = nextInt (i+1);
			tmp = arr[indizes[i]];
			arr[indizes[i]] = arr[indizes[j]];
			arr[indizes[j]] = tmp;
		}
	}
	
	/**
	 * shuffles a double array in random order
	 * all permutations are equally probable
	 * @param arr
	 */
	public void shuffle (double[] arr) {
		double tmp;
		int j;
		for (int i=1;i<arr.length;i++){
			j = nextInt (i+1);
			tmp = arr[i];
			arr[i]=arr[j];
			arr[j]=tmp;
		}		
	}

	/**
	 * shuffles a double array in random order
	 * all permutations are equally probable
	 * @param arr
	 */
	public void shuffle (float[] arr) {
		float tmp;
		int j;
		for (int i=1;i<arr.length;i++){
			j = nextInt (i+1);
			tmp = arr[i];
			arr[i]=arr[j];
			arr[j]=tmp;
		}		
	}

	/**
	 * shuffles an array of type <T> in random order
	 * all permutations are equally probable
	 * @param <T>
	 * @param arr
	 */
	public <T> void shuffle (T[] arr) {
		T tmp;
		int j;
		for (int i=1;i<arr.length;i++){
			j = nextInt (i+1);
			tmp = arr[i];
			arr[i]=arr[j];
			arr[j]=tmp;
		}		
	}
}

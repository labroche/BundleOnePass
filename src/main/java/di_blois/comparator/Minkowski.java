package di_blois.comparator;


public class Minkowski implements Comparator<float[]> {
	/**
	 *puissance : si order = 1 alors distance de Manhattan, si order = 2 alors distance Euclidienne
	 */
	int order;
	
	public Minkowski(int order){
		this.order = order;
	}
	
	/**
	 * Fonction de comparaison de vecteurs numeriques
	 * @param v1 : premier vecteur
	 * @param v2 : second vecteur
	 * @return la distance de Minkowski entre les 2 vecteurs  l'ordre specifie
	 */
	public double compare(float[] v1, float[] v2){
		if (v1.length != v2.length){
			System.err.println("Error in comparisons vectors must be the same length");
			return -1;
		} else {
			double somme = 0;
			for (int i = 0; i < v1.length; i++){
				if (order == 1){
					somme += Math.abs(v1[i] - v2[i]);
				} else {
					somme += Math.pow(v1[i] - v2[i], order);
				}
			}
			return Math.pow(somme, 1.0/order);
		}
	}
}

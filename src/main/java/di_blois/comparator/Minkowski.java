package di_blois.comparator;


public class Minkowski implements Comparator<float[]> {
	
	int order;
	
	public Minkowski(int order){
		this.order = order;
	}
	
	/**
	 * Fonction de comparaison de vecteurs numeriques
	 * @param vector1 : premier vecteur
	 * @param vector2 : second vecteur
	 * @param order : puissance : si order = 1 alors distance de Manhattan, si order = 2 alors distance Euclidienne
	 * @return la distance de Minkowski entre les 2 vecteurs  l'ordre specifie
	 */
	public double compare(float[] v1, float[] v2) throws Exception{
		if (v1.length != v2.length){
			throw new Exception("Vectors must contain the same number of numerical values");
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

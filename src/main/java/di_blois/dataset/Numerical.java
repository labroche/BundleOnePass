package di_blois.dataset;

import java.util.List;

public class Numerical extends Data<float[]> {

	public Numerical(List<float[]> d, List<String> l) {
		super(d, l);
		if (!d.isEmpty()) natt = d.get(0).length;
	}

}

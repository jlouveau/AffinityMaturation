
package tip.maam3;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import tip.math.TipVector;

/**
 * Represents a vaccination strategy in an affinity maturation simulation.
 *
 * <p>
 * A vaccination strategy is characterized by its number of vaccines and a list
 * of vaccines.
 * 
 */

public class VaccinationStrategy {
	private final int vaccineCount;
	private final List<Vaccine> vaccines;

	private VaccinationStrategy(int vaccineCount, List<Vaccine> vaccines) {
		this.vaccineCount = vaccineCount;
		this.vaccines = Collections.unmodifiableList(vaccines);
	}

	/**
	 * Creates a new vaccination strategy and transforms the vaccineConc into
	 * equalAntigenConc.
	 *
	 * @param vsProp
	 *            the governing vaccination strategy properties.
	 *
	 * @return the new founder antigen.
	 */
	public static List<Vaccine> createVaccines(VaccinationStrategyProp vsProp, double vaccineConcentration) {

		List<Vaccine> vaccines = new ArrayList<Vaccine>();
		if (vsProp.getVaccineCount() > 0) {
			for (int indexVaccine = 0; indexVaccine < vsProp.getVaccineCount(); indexVaccine++) {
				int antigenNb = (int) vsProp.getAntigenCounts().get(indexVaccine);
				vaccines.add(Vaccine.createVaccine(antigenNb, (int) vsProp.getEndTimes().get(indexVaccine),
						(int) vsProp.getTimeSteps().get(indexVaccine),
						(int) vsProp.getMutationalDistances().get(indexVaccine), vaccineConcentration/antigenNb));
			}
		} else {
			throw new IllegalStateException("no vaccine was given.");
		}
		return vaccines;
	}
}

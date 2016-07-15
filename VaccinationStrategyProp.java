
package tip.maam3;

import tip.app.TipProperties;
import tip.math.TipVector;

/**
 * Defines global vaccination strategy properties:
 *
 * <ul>
 * <li>{@code VaccinationStrategy.vaccineCount}: Number of vaccines.</li>
 * <li>{@code VaccinationStrategy.antigenCount}: Number of antigens per vaccine
 * .</li>
 * <li>{@code VaccinationStrategy.endTime}: end time of each vaccine.</li>
 * <li>{@code VaccinationStrategy.timeStep}: time step > 0 if vaccine is
 * sequential.</li>
 * <li>{@code VaccinationStrategy.mutationalDistance}: mutational distance
 * between antigens -1 if antigens are randomly created.</li>
 * </ul>
 *
 * <p>
 * The properties must be assigned before accessing the global instance.
 */

public class VaccinationStrategyProp {
	private final int vaccineCount;
	private final TipVector antigenCounts;
	private final TipVector endTimes;
	private final TipVector timeSteps;
	private final TipVector mutationalDistances;

	private static VaccinationStrategyProp instance = null;

	private VaccinationStrategyProp() {
		this.vaccineCount = TipProperties.getRequiredInt("VaccinationStrategy.vaccineCount");

		this.antigenCounts = TipVector.parseCSV(TipProperties.getRequired("VaccinationStrategy.antigenCounts")).lock();

		this.endTimes = TipVector.parseCSV(TipProperties.getRequired("VaccinationStrategy.endTimes")).lock();

		this.timeSteps = TipVector.parseCSV(TipProperties.getRequired("VaccinationStrategy.timeSteps")).lock();

		this.mutationalDistances = TipVector
				.parseCSV(TipProperties.getRequired("VaccinationStrategy.mutationalDistances")).lock();

	}

	/**
	 * Returns the global instance.
	 *
	 * <p>
	 * The properties must be assigned (by loading the governing properties
	 * file) prior to accessing the global instance.
	 *
	 * @return the global instance.
	 *
	 * @throws RuntimeException
	 *             unless all required properties have been assigned in the
	 *             virtual machine.
	 */
	public static VaccinationStrategyProp instance() {
		if (instance == null)
			instance = new VaccinationStrategyProp();

		return instance;
	}

	/**
	 * Returns the number of vaccines to be injected.
	 *
	 * @return the number of vaccines to be injected.
	 */
	public int getVaccineCount() {
		return vaccineCount;
	}

	/**
	 * Returns the number of antigens for each vaccine.
	 *
	 * @return the number of antigens for each vaccine.
	 */
	public TipVector getAntigenCounts() {
		return antigenCounts;
	}

	/**
	 * Returns the end time of each vaccine.
	 *
	 * @return the end time of each vaccine.
	 */
	public TipVector getEndTimes() {
		return endTimes;
	}

	/**
	 * Returns the time step between antigens for each vaccine.
	 *
	 * @return the time step between antigens for each vaccine.
	 */
	public TipVector getTimeSteps() {
		return timeSteps;
	}

	/**
	 * Returns the mutational distance between antigens for each vaccine.
	 *
	 * @return the mutational distance between antigens for each vaccine.
	 */
	public TipVector getMutationalDistances() {
		return mutationalDistances;
	}

}

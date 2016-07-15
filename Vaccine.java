
package tip.maam3;

import java.util.ArrayList;
import java.util.List;

import tip.math.TipVector;

/**
 * Represents a vaccine in an affinity maturation simulation.
 *
 * <p>
 * A vaccine is characterized by the number of antigens, the time where the
 * vaccine stops having any effect, the time step between injections for a given
 * vaccine (time step = 0 if cocktail) and whether the antigens are randomly
 * chosen (mutational distance = -1) or if they are closely related with a
 * specified mutational distance. the concentration is at the time of injection
 * (different from the antigen concentration over time).
 * 
 */
public final class Vaccine {
	private final int antigenCount;
	private final int endTime;
	private final int timeStep;
	/** timeStep is 0 if the vaccine is a cocktail */
	private final double mutationalDistance;
	private double equalAntigenConcentration;

	private Vaccine(int antigenCount, int endTime, int timeStep, double mutationalDistance, double equalAntigenConcentration) {
		this.antigenCount = antigenCount;
		this.endTime = endTime;
		this.timeStep = timeStep;
		this.mutationalDistance = mutationalDistance;
		this.equalAntigenConcentration = equalAntigenConcentration;
	}

	/**
	 * Creates the vaccine.
	 *
	 * @param antigenCount
	 *            the number of antigens.
	 *
	 * @param timeStep
	 *            the time step between injection of 2 antigens.
	 * 
	 * @param mutationalDistance
	 *            the mutational distance between 2 antigens.
	 * 
	 * @param concentration.
	 *
	 * @return the vaccine.
	 */

	public static Vaccine createVaccine(int antigenCount, int endTime, int timeStep, double mutationalDistance,
			double equalAntigenConcentration) {
		Vaccine vaccine = new Vaccine(antigenCount, endTime, timeStep, mutationalDistance, equalAntigenConcentration);
		if (antigenCount < 0) {
			throw new IllegalStateException("no antigen was given.");
		}
		return vaccine;
	}

	/**
	 * Test if the vaccine is sequential
	 * 
	 * @param antigenCount
	 *            the number of antigens.
	 * 
	 * @param timeStep
	 *            the time step between injection of 2 antigens.
	 * 
	 * @return boolean
	 */

	public static boolean IsSequential(Vaccine vaccine) {
		return (vaccine.getAntigenCount() > 1 && vaccine.getTimeStep() > 0);
	}

	/**
	 * Creates the antigens for a vaccine.
	 *
	 * @param antigenCount
	 *            the number of antigens.
	 *
	 * @param timeStep
	 *            the time step between injection of 2 antigens.
	 * 
	 * @param mutationalDistance
	 *            the mutational distance between 2 antigens.
	 *
	 * @return a list containing the antigens composing the vaccine.
	 */

	public static List<Antigen> createAntigens(int antigenCount, double mutationalDistance, double equalAntigenConc,
			AntigenProp agProp) {
		// Check that mutationalDistance isn't too big
		if (mutationalDistance > ((int) Math.pow(2, agProp.getVariableLength())) / antigenCount) {
			throw new IllegalStateException(
					"the mutational distance is too big given the length and number of antigens");
		}

		List<Antigen> antigens = new ArrayList<Antigen>();
		antigens.addAll(createCocktail(antigenCount, mutationalDistance, equalAntigenConc, agProp));

		// CHECK SIZE of List<Antigen>
		if (antigens.size() != antigenCount) {
			throw new IllegalStateException("the number of antigens doesn't fit the properties of the vaccine.");
		}
		return antigens;
	}

	/**
	 * Creates the antigens for a cocktail. A founder antigen is created first.
	 * If antigenCount > 1, a variant is created from the founder by randomly
	 * adding mutations. The total number of mutations between the two strains
	 * is given by the mutational distance. The following antigens are creating
	 * from the last antigen in the list and before adding them, we check that
	 * there are no duplicates.
	 *
	 * @param antigenCount
	 *            the number of antigens.
	 *
	 * @param mutationalDistance
	 *            the mutational distance between 2 antigens.
	 *
	 * @return a list containing the antigens composing the cocktail.
	 */

	public static List<Antigen> createCocktail(int antigenCount, double mutationalDistance, double equalAntigenConc,
			AntigenProp agProp) {
		List<Antigen> antigens = new ArrayList<Antigen>();

		Antigen founder = Antigen.createRandom(agProp, equalAntigenConc);
		antigens.add(founder);
		if (antigenCount > 1) {
			antigens.addAll(
					Antigen.createVariants(antigenCount - 1, agProp, equalAntigenConc, mutationalDistance, founder));
		}
		return antigens;
	}

	/**
	 * Returns the number of antigens for each vaccine.
	 *
	 * @return the number of antigens for each vaccine.
	 */
	public int getAntigenCount() {
		return antigenCount;
	}

	/**
	 * Returns the end time of each vaccine.
	 *
	 * @return the end time of each vaccine.
	 */
	public int getEndTime() {
		return endTime;
	}

	/**
	 * Returns the time step between antigens for each vaccine.
	 *
	 * @return the time step between antigens for each vaccine.
	 */
	public int getTimeStep() {
		return timeStep;
	}

	/**
	 * Returns the mutational distance between antigens for each vaccine.
	 *
	 * @return the mutational distance between antigens for each vaccine.
	 */
	public double getMutationalDistance() {
		return mutationalDistance;
	}

	/**
	 * Returns the concentration of the antigens for all vaccines.
	 *
	 * @return the concentration of the antigens for all vaccines.
	 */
	public double getEqualAntigenConcentration() {
		return equalAntigenConcentration;
	}
}

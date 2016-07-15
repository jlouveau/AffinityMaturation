
package tip.maam3;

import java.util.ArrayList;
import java.util.List;

import tip.math.TipVector;

/**
 * Represents an antigen in an affinity maturation simulation. An antigen is the
 * list of epitopes introduced in the GC between the end of the previous
 * injection ( + Vaccine.getTimeStep() * index/Vaccine.getAntigenCount() ) and
 * the following injection
 * <p>
 * The epitopes are immutable but the concentration may change during the course
 * of affinity maturation.
 */
public final class Antigen {
	private final Epitope epitope;
	private double concentration;

	private Antigen(Epitope epitope, double concentration) {
		this.epitope = epitope;
		this.concentration = concentration;
	}

	/**
	 * Creates a new conserved antigen.
	 *
	 * @param agProp
	 *            the governing antigen properties.
	 *
	 * @param concentration
	 *            the antigen concentration.
	 *
	 * @return the new conserved antigen.
	 */
	public static Antigen createConserved(AntigenProp agProp, double concentration) {
		Epitope epitope = Epitope.createConserved(agProp);

		return new Antigen(epitope, concentration);
	}

	/**
	 * Creates a new variant antigen.
	 *
	 * @param variantCount
	 *            the number of variants to be created.
	 *
	 * @param agProp
	 *            the governing antigen properties.
	 * 
	 * @param concentration
	 *            the antigen concentration.
	 *
	 * @param mutationalDistance
	 *            the number of mutations between consecutive antigens
	 * 
	 * @param founder
	 *            the "wild type" antigen from which the variants are created.
	 *
	 * @return the new variant antigen.
	 */
	public static List<Antigen> createVariants(int variantCount, AntigenProp agProp, double concentration,
			double mutationalDistance, Antigen founder) {
		List<Antigen> variants = new ArrayList<Antigen>();

		Epitope epitope = Epitope.createVariant(mutationalDistance, founder);
		variants.add(new Antigen(epitope, concentration));

		if (variantCount > 1) {
			int index = 0;
			do {
				epitope = Epitope.createVariant(mutationalDistance, variants.get(index));
				if (!Epitope.hasDuplicate(epitope, founder) && !Epitope.hasDuplicateInList(epitope, variants)) {
					variants.add(new Antigen(epitope, concentration));
					index++;
				}
			} while (index < variantCount - 1);
		}
		return variants;
	}

	/**
	 * Creates a new random antigen.
	 *
	 * @param agProp
	 *            the governing antigen properties.
	 *
	 * @param concentration
	 *            the antigen concentration.
	 *
	 * @return the new random antigen.
	 */
	public static Antigen createRandom(AntigenProp agProp, double concentration) {
		Epitope epitope = Epitope.createRandom(agProp);

		return new Antigen(epitope, concentration);
	}

	/**
	 * Returns all epitopes to which a B cell receptor may bind (in an
	 * unmodifiable list).
	 *
	 * @return all epitopes to which a B cell receptor may bind.
	 */
	public Epitope getEpitope() {
		return epitope;
	}

	/**
	 * Returns the concentration of this antigen in the germinal center.
	 *
	 * @return the concentration of this antigen in the germinal center.
	 */
	public final double getConcentration() {
		return concentration;
	}

	/**
	 * Assigns the concentration of this antigen in the germinal center.
	 *
	 * @param concentration
	 *            the concentration of this antigen in the germinal center.
	 */
	public final void setConcentration(double concentration) {
		this.concentration = concentration;
	}

	@Override
	public String toString() {
		return String.format("Antigen(%s, %f)", getEpitope().toString(), getConcentration());
	}
}

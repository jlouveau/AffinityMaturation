
package tip.maam3;

import java.util.ArrayList;
import java.util.Collection;

import org.apache.commons.collections.primitives.ArrayDoubleList;

import tip.math.BitVector;
import tip.math.GeneralizedExtremeValueDistribution;

/**
 * Represents a B cell receptor as a vector coord of fixed length
 * receptorLength. Receptors are immutable: once created, they are fixed
 * forever.
 */

public final class Receptor {
	private final double[] coord;

	private Receptor(double[] coord) {
		this.coord = coord;
	}

	/**
	 * Creates a receptor for founder B cells.
	 *
	 * <p>
	 * Randomly generate a BitVector of length equal to the receptorLength.
	 * Check that the affinity of the receptor for the conserved antigen is
	 * bigger than the activationAffinity threshold.
	 *
	 * @param bcprop
	 *            BCell prop
	 *
	 * @param activationAffinity
	 *            the binding affinity low-threshold (in units of kT) between
	 *            the founder cell receptor and a conserved antigen.
	 *
	 * @return a randomly generated receptor with affinity above the
	 *         low-threshold.
	 */
	public static Receptor createFounder(BCellProp bcProp, Epitope founderEpitope) {
		double activationEnergy = bcProp.getActivationThreshold();
		double max = bcProp.getMax();
		double min = bcProp.getMin();
		Receptor founderReceptor = generateRandomReceptor(bcProp);
		int randomIndex = 0;
		do {
			randomIndex = (int) (Math.random() * bcProp.getReceptorLength());
			founderReceptor.getCoord()[randomIndex] = (max - min) * Math.random() + min;
// CHANGE IT so that it doesn't start from scratch?
		} while (founderReceptor.calculateEnergy(founderEpitope) < activationEnergy);

		return founderReceptor;

	}

	private static Receptor generateRandomReceptor(BCellProp bcProp) {
		double[] coord = new double[bcProp.getReceptorLength()];
		double max = bcProp.getMax();
		double min = bcProp.getMin();

		for (int i = 0; i < bcProp.getReceptorLength(); i++) {
			double RANDOM = Math.random();
			coord[i] = (max - min) * RANDOM + min;
		}
		return new Receptor(coord);
	}

	private boolean validatePair(Epitope epitope) {
		//System.out.println(paratope.getReceptorLength());
		//System.out.println(epitope.getLength());
		if (epitope == null){
			System.out.println("epitope is null");
		}
		return (this.getReceptorLength() == epitope.getLength());
			//throw new IllegalStateException("epitope and paratope lengths not equal.");
	}

	/**
	 * Calculates the binding energy between the receptor and the epitope as the sum
	 * of the product between the values of the receptor and the epitope at
	 * every bit.
	 * 
	 * @param epitope
	 * 
	 * @return a sum(V(k)*A(k))
	 */
	public double calculateEnergy(Epitope epitope) {
		BCellProp bcProp = BCellProp.instance();
		if (this.validatePair(epitope)){
			int length = bcProp.getReceptorLength();

			double energy = 0.0;
			int[] epitopeToIntegers = epitope.getCoord().formatToInteger();
			for (int i = 0; i < length; i++) {
				energy = energy + this.getCoord()[i] * (double) epitopeToIntegers[i];
			}
			return energy;
		}
		else{
			System.out.println(this.getReceptorLength());
			System.out.println(epitope.getLength());
			throw new IllegalStateException("epitope and paratope lengths not equal.");
			//return Double.NaN;
		}
	}

	/**
	 * Creates a new receptor with mutations around this receptor leaving this
	 * receptor unchanged.
	 *
	 * <p>
	 * The coordinates of the new receptor are a copy of the coordinates of the
	 * receptor with one mutated bit.
	 * 
	 * @param parent
	 *            is the parent B cell
	 *
	 * @return a new receptor with a random mutation away from the receptor of
	 *         the parent BCell
	 */
	public static Receptor mutate(BCell parent, Epitope target, BCellProp bcProp, AntigenProp agProp) {
		double[] newCoord = parent.getReceptor().getCoord();
		double sigma = bcProp.getSigma();
		double mu = bcProp.getMu();
		double kappa = bcProp.getKappa();
		double alpha = bcProp.getAlpha();
		double mutationMax = bcProp.getMutationMax();
		double mutationMin = bcProp.getMutationMin();
		int consLength = agProp.getConservedLength();

		int randomNum = (int) (Math.random() * (parent.getReceptor().getReceptorLength()));
		double delta = mutateSite(target, randomNum, kappa, sigma, mu);
		newCoord[randomNum] = newCoord[randomNum] + delta;

		if (randomNum >= consLength){
			//change variable site
			int randomCons = (int) (Math.random() * consLength);
			newCoord[randomCons] = newCoord[randomCons] - alpha * delta;
		}
		
		if (newCoord[randomNum] < mutationMin){
			newCoord[randomNum] = mutationMin;
		}
		if (newCoord[randomNum] > mutationMax){
			newCoord[randomNum] = mutationMax;
		}
		
		return new Receptor(newCoord);
	}

	public static double mutateSite(Epitope target, int siteIndex, double kappa, double sigma, double mu) {
		if (!(target.getCoord().get(siteIndex))) { //if V(k) == 1
			return -GeneralizedExtremeValueDistribution.gevdSampling(kappa, sigma, mu);
		} else {
			return GeneralizedExtremeValueDistribution.gevdSampling(kappa, sigma, mu);
		}
	}

	/**
	 * Finds the neutralization target from a collection of epitopes: the
	 * epitope to which this receptor binds most strongly.
	 *
	 * @param epitopes
	 *            the epitopes to search.
	 *
	 * @return the epitope to which the receptor binds most strongly.
	 */
	public static Epitope findTarget(Collection<Epitope> epitopes, Receptor paratope) {
		double targetEnergy = Double.NEGATIVE_INFINITY;
		Epitope targetEpitope = null;

		for (Epitope epitope : epitopes) {
			double energy = paratope.calculateEnergy(epitope);

			if (energy > targetEnergy) {
				targetEnergy = energy;
				targetEpitope = epitope;
			}
		}

		return targetEpitope;
	}

	@Override
	public String toString() {
		return String.format("Receptor(%s)", getCoord().toString());
	}

	/**
	 * Returns the number of bits of the receptor.
	 *
	 * @return the number of bits of the receptor.
	 */
	public int getReceptorLength() {
		return coord.length;
	}

	/**
	 * Returns the coordinates of this receptor.
	 *
	 * @return the coordinates of this receptor.
	 */
	public double[] getCoord() {
		return coord;
	}

}

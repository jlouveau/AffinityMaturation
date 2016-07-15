
package tip.maam3;

import tip.app.TipProperties;

/**
 * Defines global B cell properties:
 *
 * <ul>
 * <li>{@code BCell.activationThreshold}: low-threshold affinity for activation
 * of B cell.</li>
 * <li>{@code BCell.mutationProbability}: Probability that daughter cells mutate
 * after division.</li>
 * <li>{@code BCell.mutationAffinity}: Average change in binding free energy
 * upon mutation (in units of kT).</li>
 * </ul>
 *
 * <p>
 * The properties must be assigned before accessing the global instance.
 */
public final class BCellProp {
	private final int terminationMultiple;
	private final double activationThreshold;
	private final double mutationProbability;
	private final double max;
	private final double min;
	private final double mutationMax;
	private final double mutationMin;
	private final double alpha;
	private final double sigma;
	private final double mu;
	private final double kappa;

	private static BCellProp instance = null;

	private BCellProp() {
		this.terminationMultiple = TipProperties.getRequiredInt("BCell.terminationMultiple");
		this.activationThreshold = TipProperties.getRequiredDouble("BCell.activationThreshold");
		this.mutationProbability = TipProperties.getRequiredDouble("BCell.mutationProbability");
		this.max = TipProperties.getRequiredDouble("BCell.max");
		this.min = TipProperties.getRequiredDouble("BCell.min");
		this.mutationMax = TipProperties.getRequiredDouble("BCell.mutationMax");
		this.mutationMin = TipProperties.getRequiredDouble("BCell.mutationMin");
		this.alpha = TipProperties.getRequiredDouble("BCell.alpha");
		this.sigma = TipProperties.getRequiredDouble("BCell.sigma");
		this.mu = TipProperties.getRequiredDouble("BCell.mu");
		this.kappa = TipProperties.getRequiredDouble("BCell.kappa");
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
	public static BCellProp instance() {
		if (instance == null)
			instance = new BCellProp();

		return instance;
	}

	/**
	 * Returns the low-threshold affinity for activation of B cell
	 *
	 * @return the low-threshold affinity for activation of B cell.
	 */
	public double getActivationThreshold() {
		return activationThreshold;
	}

	/**
	 * Returns the probability that daughter cells mutate after division.
	 *
	 * @return the probability that daughter cells mutate after division.
	 */
	public double getMutationProbability() {
		return mutationProbability;
	}

	public double getAlpha() {
		return alpha;
	}

	public double getSigma() {
		return sigma;
	}

	public double getMu() {
		return mu;
	}

	public double getKappa() {
		return kappa;
	}

	public double getMin() {
		return min;
	}

	public double getMax() {
		return max;
	}

	public double getMutationMax() {
		return mutationMax;
	}

	public double getMutationMin() {
		return mutationMin;
	}

	/**
	 * Returns the length of the B cell receptor.
	 *
	 * @return the length of the B cell receptor.
	 */
	public int getReceptorLength() {
		AntigenProp agProp = AntigenProp.instance();
		return agProp.getLength();
	}

	public int getTerminationMultiple() {
		return terminationMultiple;
	}
}

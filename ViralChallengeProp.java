
package tip.maam3;

import tip.app.TipProperties;
import tip.math.TipVector;

/**
 * Defines global viral challenge properties:
 *
 * <ul>
 * <li>{@code ViralChallenge.challengerCount}: Number of challengers to present.
 * </li>
 * <li>{@code ViralChallenge.affinityThreshold}: Affinity threshold for
 * challenger neutralization.</li>
 * </ul>
 *
 * <p>
 * The properties must be assigned before accessing the global instance.
 */
public final class ViralChallengeProp {
	private final int challengerCount;
	private final double affinityThreshold;

	private static ViralChallengeProp instance = null;

	private ViralChallengeProp() {
		this.challengerCount = TipProperties.getRequiredInt("ViralChallenge.challengerCount");

		this.affinityThreshold = TipProperties.getRequiredDouble("ViralChallenge.affinityThreshold");
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
	public static ViralChallengeProp instance() {
		if (instance == null)
			instance = new ViralChallengeProp();

		return instance;
	}

	/**
	 * Returns the number of challengers to present.
	 *
	 * @return the number of challengers to present.
	 */
	public int getChallengerCount() {
		return challengerCount;
	}

	/**
	 * Returns the affinity threshold for challenger neutralization.
	 *
	 * @return the affinity threshold for challenger neutralization.
	 */
	public double getAffinityThreshold() {
		return affinityThreshold;
	}
}

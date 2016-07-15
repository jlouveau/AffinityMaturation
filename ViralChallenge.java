
package tip.maam3;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

import org.apache.commons.collections.primitives.ArrayDoubleList;

import tip.app.TipLog;
import tip.math.DoubleUtil;
import tip.math.MeanEstimator;
import tip.math.MedianEstimator;
import tip.math.TipVector;

/**
 * Computes antibody breadth for a collection of challenging antigens as a
 * function of neutralization threshold. The breadth is defined as the fraction
 * of non-conserved epitopes to which the receptor binds more strongly than a
 * specified threshold.
 */
public final class ViralChallenge {
	private final double breadth;
	private final MeanEstimator meanAffinity;
	private final MedianEstimator medianAffinity;

	private static Collection<Antigen> challengers = null;

	private ViralChallenge(double breadth, MeanEstimator meanAffinity, MedianEstimator medianAffinity) {
		this.breadth = breadth;
		this.meanAffinity = meanAffinity;
		this.medianAffinity = medianAffinity;
	}

	/**
	 * Challenge a B cell with a collection of antigens.
	 *
	 * <p>
	 * The challenger properties are read from the global properties object.
	 *
	 * @param bcell is the B cell to challenge. Its affinities for each
	 * challenger is determined and compared to the threshold to determine 
	 * the fraction of neutralized challengers which is the breadth.
	 *
	 * @return the results of the challenge.
	 */
	public static ViralChallenge challenge(BCell bcell) {
		double[] affinities = computeAffinities(bcell);
		double breadth = computeBreadth(affinities, getThreshold());

		return new ViralChallenge(breadth, MeanEstimator.estimate(affinities), MedianEstimator.estimate(affinities));
	}

	private static double[] computeAffinities(BCell bcell) {
		ArrayDoubleList affinities = new ArrayDoubleList();

		for (Antigen challenger : getChallengers()) {
			Epitope epitope = challenger.getEpitope();
			affinities.add(bcell.getReceptor().calculateEnergy(epitope));
		}

		return affinities.toArray();
	}

	/**
	 * Returns the fraction of viral challengers 
	 * which bind to a given bcell with affinity stronger than the threshold. 
	 * 
	 * @param affinities is a vector that stored the affinity between the 
	 *  Bcell and all challenger antigens.
	 * 
	 * @return
	 */

	private static double computeBreadth(double[] affinities, double threshold) {
		int neutralized = 0;

		for (double affinity : affinities){
			//System.out.println(affinity);
			if (affinity >= threshold){
				++neutralized;
			}
		}
		double breadthBcell = DoubleUtil.ratio(neutralized, affinities.length);
		if (breadthBcell > 1 || breadthBcell < 0){
			System.out.println(breadthBcell);
		}

		return breadthBcell;
	}

	/**
	 * Returns the antigen challengers.
	 *
	 * @return the antigen challengers (in an unmodifiable collection).
	 */
	public static Collection<Antigen> getChallengers() {
		if (challengers == null)
			challengers = createChallengers();

		return challengers;
	}

	private static Collection<Antigen> createChallengers() {
		TipLog.info("Creating viral challengers...");

		int challengerCount = ViralChallengeProp.instance().getChallengerCount();
		ArrayList<Antigen> challengerList = new ArrayList<Antigen>(challengerCount);

		while (challengerList.size() < challengerCount) {
			challengerList.add(Antigen.createRandom(AntigenProp.instance(), 0.0));
		}
		return Collections.unmodifiableList(challengerList);
	}

	/**
	 * Returns the affinity threshold.
	 *
	 * @return the affinity threshold.
	 */
	public static double getThreshold() {
		return ViralChallengeProp.instance().getAffinityThreshold();
	}

	/**
	 * Returns the computed breadth .
	 *
	 * @return the breadth depends on the 
	 *         affinity threshold.
	 */
	public double getBreadth() {
		return breadth;
	}

	/**
	 * Returns the mean affinity statistics.
	 *
	 * @return the mean affinity statistics.
	 */
	public MeanEstimator getMeanAffinity() {
		return meanAffinity;
	}

	/**
	 * Returns the median affinity statistics.
	 *
	 * @return the median affinity statistics.
	 */
	public MedianEstimator getMedianAffinity() {
		return medianAffinity;
	}
}

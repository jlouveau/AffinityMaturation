
package tip.maam3;

import java.io.File;

import tip.app.TipProperties;
import tip.dist.RealDistribution;
import tip.dist.RealDistributionType;
import tip.maam.FDCPresentation;
import tip.math.TipVector;

/**
 * Defines global germinal center properties:
 *
 * <ul>
 * <li>{@code GerminalCenter.germlineCount}: Number of unique founder B cells
 * (germlines) before replication.</li>
 * <li>{@code GerminalCenter.replicationFactor}: Number of B cells produced from
 * each germline by replication.</li>
 * <li>{@code GerminalCenter.cycleLimit}: Maximum number of mutation/selection
 * cycles.</li>
 * <li>{@code GerminalCenter.concentration}: concentration of vaccine.</li>
 * <li>{@code GerminalCenter.antigenPresentation}: Type of antigen presentation
 * by follicular dendritic cells.</li>
 * <li>{@code GerminalCenter.tcellSelectionRate}: Fraction of B cells receiving
 * T cell help at each maturation cycle.</li>
 * <li>{@code GerminalCenter.memorySelectionRate}: Fraction of B surviving B
 * cells that are selected to exit the germinal center and become memory cells.
 * </li>
 * <li>{@code GerminalCenter.plasmaSelectionRate}: Fraction of B surviving B
 * cells that are selected to exit the germinal center and become plasma cells,
 * but only if they exceeed an affinity threshold for epitope binding.</li>
 * <li>{@code GerminalCenter.plasmaAffinityThreshold}: Developing B cells must
 * bind at least one epitope with this affinity or greater to become plasma
 * cells.</li>
 * <li>
 * </ul>
 *
 * <p>
 * The properties must be assigned before accessing the global instance.
 */
public final class GerminalCenterProp {
	private final int germlineCount;
	private final int replicationFactor;
	private final int cycleLimit;
	private final int numberOfAntigensPresented;

	private final double memorySelectionRate;
	private final double plasmaSelectionRate;
	private final double plasmaAffinityThreshold;
	private final double tcellSelectionRate;
	private final double concentration;

	private static GerminalCenterProp instance = null;

	private GerminalCenterProp() {
		this.germlineCount = TipProperties.getRequiredInt("GerminalCenter.germlineCount");
		this.replicationFactor = TipProperties.getRequiredInt("GerminalCenter.replicationFactor");
		this.cycleLimit = TipProperties.getRequiredInt("GerminalCenter.cycleLimit");
		this.numberOfAntigensPresented = TipProperties.getRequiredInt("GerminalCenter.numberOfAntigensPresented");

		this.concentration = TipProperties.getRequiredDouble("GerminalCenter.concentration");		
		this.tcellSelectionRate = TipProperties.getRequiredDouble("GerminalCenter.tcellSelectionRate");
		this.memorySelectionRate = TipProperties.getRequiredDouble("GerminalCenter.memorySelectionRate");
		this.plasmaSelectionRate = TipProperties.getRequiredDouble("GerminalCenter.plasmaSelectionRate");
		this.plasmaAffinityThreshold = TipProperties.getRequiredDouble("GerminalCenter.plasmaAffinityThreshold");
	}

	/**
	 * Returns the global instance.
	 *
	 * <p>
	 * The properties must be assigned (e.g., by calling {@code load}) prior to
	 * accessing the global instance.
	 *
	 * @return the global instance.
	 *
	 * @throws RuntimeException
	 *             unless all required properties have been assigned in the
	 *             virtual machine.
	 */
	public static GerminalCenterProp instance() {
		if (instance == null)
			instance = new GerminalCenterProp();

		return instance;
	}

	/**
	 * Loads global antigen properties from a file.
	 *
	 * @param propFile
	 *            name of the governing properties file.
	 *
	 * @throws RuntimeException
	 *             unless the properties file contains valid values for all
	 *             required properties.
	 */
	public static void load(String propFile) {
		load(new File(propFile));
	}

	/**
	 * Loads global antigen properties from a file.
	 *
	 * @param propFile
	 *            the governing properties file.
	 *
	 * @throws RuntimeException
	 *             unless the properties file contains valid values for all
	 *             required properties.
	 */
	public static void load(File propFile) {
		TipProperties.loadFile(propFile);
	}

	/**
	 * Returns the number of unique founder B cells (germlines) before
	 * replication.
	 *
	 * @return the number of unique founder B cells (germlines) before
	 *         replication.
	 */
	public int getGermlineCount() {
		return germlineCount;
	}

	/**
	 * Returns the number of B cells produced from each germline by replication.
	 *
	 * @return the number of B cells produced from each germline by replication.
	 */
	public int getReplicationFactor() {
		return replicationFactor;
	}

	/**
	 * Returns the maximum number of mutation/selection cycles.
	 *
	 * @return the maximum number of mutation/selection cycles.
	 */
	public int getCycleLimit() {
		return cycleLimit;
	}

	/**
	 * Returns the fraction of B surviving B cells that are selected to exit the
	 * germinal center and become memory cells.
	 *
	 * @return the memory cell selection rate.
	 */
	public double getMemorySelectionRate() {
		return memorySelectionRate;
	}

	/**
	 * Returns the fraction of B surviving B cells that are selected to exit the
	 * germinal center and become plasma cells, but only if they exceeed an
	 * energy threshold for variable or conserved epitopes.
	 *
	 * @return the plasma cell selection rate.
	 */
	public double getPlasmaSelectionRate() {
		return plasmaSelectionRate;
	}

	/**
	 * Returns the affinity threshold for plasma cell selection: Developing B
	 * cells must bind at least one epitope with this affinity or better to
	 * become plasma cells.
	 *
	 * @return the affinity threshold for plasma cell selection.
	 */
	public double getPlasmaAffinityThreshold() {
		return plasmaAffinityThreshold;
	}
	
	public double getTcellSelectionRate(){
		return tcellSelectionRate;
	}

	/**
	 * Returns the vaccine concentration.
	 *
	 * @return the vaccine concentration.
	 */
	public double getConcentration() {
		return concentration;
	}

	/**
	 * Returns the type of antigen presentation by follicular dendritic cells.
	 *
	 * @return the type of antigen presentation by follicular dendritic cells.
	 */
	public int getNumberOfAntigensPresented() {
		return numberOfAntigensPresented;
	}
}

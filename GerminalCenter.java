
package tip.maam3;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections.primitives.ArrayDoubleList;

import tip.app.TipLog;
import tip.app.TipProperties;
import tip.math.DoubleUtil;
import tip.math.MeanEstimator;
import tip.math.MedianEstimator;
import tip.math.TipRandom;
import tip.math.TipVector;
import tip.util.ListUtil;

public final class GerminalCenter {
	private final double vaccineConcentration;
	private boolean isExtinguished;

	// List of vaccines from Vaccination Strategy
	List<Vaccine> vaccines = new ArrayList<Vaccine>();

	private final List<Antigen> antigens = new ArrayList<Antigen>();

	private final Collection<BCell> founderCells = new HashSet<BCell>(); 
	// Only the few germline B cells

	private final Collection<BCell> memoryCells = new HashSet<BCell>();
	private final Collection<BCell> plasmaCells = new HashSet<BCell>();

	// Unique epitopes from all antigens: must be a set to eliminate
	// potential duplicate conserved epitopes...
	private final List<Epitope> antigenEpitopes = new ArrayList<Epitope>();

	// List of B-cell lists: activeCells.get(k) contains the actively
	// maturing B cells at the kth cycle...
	private final List<List<BCell>> activeCells = new ArrayList<List<BCell>>();

	// Nearest epitope for each plasma cell...
	private final Map<BCell, Epitope> targetEpitopes = new HashMap<BCell, Epitope>();

	// Viral challenge to each plasma cell...
	private final Map<BCell, ViralChallenge> viralChallenges = new HashMap<BCell, ViralChallenge>();

	private final TipRandom random = TipRandom.instance();
	private final VaccinationStrategyProp vsProp = VaccinationStrategyProp.instance();
	private final AntigenProp agProp = AntigenProp.instance();
	private final GerminalCenterProp gcProp = GerminalCenterProp.instance();

	private int cycleIndex;

	// Creation must happen through public static methods...
	private GerminalCenter() {
		this.vaccineConcentration = gcProp.getConcentration();
		this.isExtinguished = false;
	}

	/**
	 * Creates a new germinal center and simulates the affinity maturation
	 * process.
	 *
	 * @return the germinal center, after affinity maturation has terminated.
	 *
	 * @throws RuntimeException
	 *             unless the antigen, B cell, and germinal center properties
	 *             have been assigned in the VM.
	 */
	public static GerminalCenter run(PrintWriter trialWriter) {
		GerminalCenter gc = new GerminalCenter();

		return gc.runAM(trialWriter);
	}

	/**
	 * Returns the number of cycles executed during affinity maturation.
	 *
	 * @return the number of cycles executed during affinity maturation.
	 */
	public int getCycleCount() {
		return cycleIndex;
	}

	/**
	 * Returns the concentration of each vaccine.
	 *
	 * @return the concentration of each vaccine.
	 */
	public static double getVaccineConc() {
		return GerminalCenterProp.instance().getConcentration();
	}

	/**
	 * Returns the antigens present in this germinal center.
	 *
	 * @return the antigens present in this germinal center (in an unmodifiable
	 *         list).
	 */
	public Collection<Antigen> getAntigens() {
		return Collections.unmodifiableCollection(antigens);
	}

//	/**
//	 * Returns the unique epitopes present in this germinal center.
//	 *
//	 * @return the unique epitopes present in this germinal center.
//	 */
//	public Collection<Epitope> getEpitopes() {
//		return Collections.unmodifiableCollection(antigenEpitopes);
//	}

	/**
	 * Returns the founder (germline) B cells that seeded this germinal center.
	 *
	 * @return the founder (germline) B cells that seeded this germinal center
	 *         (in an unmodifiable collection).
	 */
	public Collection<BCell> getFounderCells() {
		return Collections.unmodifiableCollection(founderCells);
	}

	/**
	 * Returns the unique memory cells produced by affinity maturation in this
	 * germinal center.
	 *
	 * @return the unique memory cells produced by affinity maturation in this
	 *         germinal center (in an unmodifiable collection).
	 */
	public Collection<BCell> getMemoryCells() {
		return Collections.unmodifiableCollection(memoryCells);
	}

	/**
	 * Returns the unique plasma cells produced by affinity maturation in this
	 * germinal center.
	 *
	 * @return the unique plasma cells produced by affinity maturation in this
	 *         germinal center (in an unmodifiable collection).
	 */
	public Collection<BCell> getPlasmaCells() {
		return Collections.unmodifiableCollection(plasmaCells);
	}

	/**
	 * Returns the final generation of active cells (from the last step of
	 * affinity maturation); will be empty if the germinal center extinguished.
	 *
	 * @return the final generation of active cells (unmodifiable list).
	 */
	public List<BCell> getActiveCells() {
		return Collections.unmodifiableList(activeCells.get(activeCells.size() - 1));
	}

	/**
	 * Computes the plasma production rate of the affinity maturation: the
	 * number of plasma cells produced for each replicated founder.
	 *
	 * @return the plasma production rate of the affinity maturation.
	 */
	public double computePlasmaProductionRate() {
		return DoubleUtil.ratio(plasmaCells.size(), activeCells.get(0).size());
	}

	/**
	 * Computes the active production rate of the affinity maturation: the
	 * number of active cells produced for each replicated founder.
	 *
	 * @return the active production rate of the affinity maturation.
	 */
	public double computeActiveProductionRate() {
		return DoubleUtil.ratio(getActiveCells().size(), activeCells.get(0).size());
	}

	/**
	 * Computes the mean generation of plasma cells produced by affinity
	 * maturation.
	 *
	 * @return the mean generation of plasma cells.
	 */
	public MeanEstimator computeMeanPlasmaGeneration() {
		return MeanEstimator.estimate(collectPlasmaGeneration());
	}

	/**
	 * Computes the median generation of plasma cells produced by affinity
	 * maturation.
	 *
	 * @return the median generation of plasma cells.
	 */
	public MedianEstimator computeMedianPlasmaGeneration() {
		return MedianEstimator.estimate(collectPlasmaGeneration());
	}

	private double[] collectPlasmaGeneration() {
		ArrayDoubleList values = new ArrayDoubleList();

		for (BCell plasmaCell : plasmaCells)
			values.add(plasmaCell.getGeneration());

		return values.toArray();
	}

	/**
	 * Computes the mean number of mutations in plasma cells produced by
	 * affinity maturation.
	 *
	 * @return the mean mutation count of plasma cells.
	 */
	public MeanEstimator computeMeanPlasmaMutationCount() {
		return MeanEstimator.estimate(collectPlasmaMutationCount());
	}

	/**
	 * Computes the median number of mutations in plasma cells produced by
	 * affinity maturation.
	 *
	 * @return the median mutation count of plasma cells.
	 */
	public MedianEstimator computeMedianPlasmaMutationCount() {
		return MedianEstimator.estimate(collectPlasmaMutationCount());
	}

	private double[] collectPlasmaMutationCount() {
		ArrayDoubleList values = new ArrayDoubleList();

		for (BCell plasmaCell : plasmaCells)
			values.add(plasmaCell.getMutationCount());

		return values.toArray();
	}

	/**
	 * Computes the mean affinity between plasma cells and their nearest
	 * epitope.
	 *
	 * @return the estimated mean affinity between plasma cells and their
	 *         nearest epitope.
	 */
	public MeanEstimator computeMeanTargetAffinity() {
		return MeanEstimator.estimate(collectTargetAffinity());
	}

	/**
	 * Computes the median affinity between plasma cells and their nearest
	 * epitope.
	 *
	 * @return the estimated median affinity between plasma cells and their
	 *         nearest epitope.
	 */
	public MedianEstimator computeMedianTargetAffinity() {
		return MedianEstimator.estimate(collectTargetAffinity());
	}

	private double[] collectTargetAffinity() {
		ArrayDoubleList values = new ArrayDoubleList();

		for (BCell plasmaCell : plasmaCells)
			values.add(plasmaCell.getReceptor().calculateEnergy(targetEpitopes.get(plasmaCell)));

		return values.toArray();
	}

	// /**
	// * Computes the mean specificity of plasma cells and their nearest
	// * epitopes.
	// *
	// * @return the estimated mean specificity.
	// */
	// public MeanEstimator computeMeanTargetSpecificity() {
	// return MeanEstimator.estimate(collectTargetSpecificity());
	// }
	//
	// /**
	// * Computes the median specificity of plasma cells and their nearest
	// * epitopes.
	// *
	// * @return the estimated median specificity.
	// */
	// public MedianEstimator computeMedianTargetSpecificity() {
	// return MedianEstimator.estimate(collectTargetSpecificity());
	// }

	// private double[] collectTargetSpecificity() {
	// ArrayDoubleList values = new ArrayDoubleList();
	//
	// for (BCell plasmaCell : plasmaCells)
	// values.add(Receptor.computeSpecificity(plasmaCell.getReceptor(),
	// targetEpitopes.get(plasmaCell)));
	//
	// return values.toArray();
	// }

	/**
	 * Computes the mean plasma cell neutralization breadth.
	 *
	 * @param index
	 *            the index of the viral challenge neutralization threshold.
	 *
	 * @return the estimated mean plasma cell neutralization breadth for the
	 *         given threshold.
	 */
	public MeanEstimator computeMeanBreadth(int index) {
		return MeanEstimator.estimate(collectBreadth(index));
	}

	/**
	 * Computes the median plasma cell neutralization breadth.
	 *
	 * @param index
	 *            the index of the viral challenge neutralization threshold.
	 *
	 * @return the estimated median plasma cell neutralization breadth for the
	 *         given threshold.
	 */
	public MedianEstimator computeMedianBreadth(int index) {
		return MedianEstimator.estimate(collectBreadth(index));
	}

	private double[] collectBreadth(int index) {
		ArrayDoubleList values = new ArrayDoubleList();

		for (BCell plasmaCell : plasmaCells)
			values.add(viralChallenges.get(plasmaCell).getBreadth());

		return values.toArray();
	}


	/**
	 * Maps each plasma cell to its nearest epitope.
	 *
	 * @return a map containing the nearest epitope for each plasma cell.
	 */
	public Map<BCell, Epitope> mapTargets() {
		return Collections.unmodifiableMap(targetEpitopes);
	}

	/**
	 * Maps each plasma cell to its viral challenge.
	 *
	 * @return a map containing the viral challenge for each plasma cell.
	 */
	public Map<BCell, ViralChallenge> mapChallenges() {
		return Collections.unmodifiableMap(viralChallenges);
	}

	private GerminalCenter runAM(PrintWriter trialWriter) {

		createVaccines();

		initialization(); // suppose that the first vaccine is a cocktail

		collectEpitopes();//stores the epitopes of the antigens of the first vaccine

		generateFounders(vaccines.get(0));

		// Special first cycle: replication only...
		cycleIndex = 0;
		replicateFounders();

		if (!vaccines.isEmpty()) {
			// ExecuteCycle() for the founder vaccine which is not sequential
			// -> cycle index is 1 or more
			do {
				executeCycle(trialWriter);
			} while (continueMaturation(vaccines.get(0)));
			antigenEpitopes.clear();
			
			if (vsProp.getVaccineCount() > 1) {
				// ExecuteCycle for the other vaccines which can be cocktails or
				// sequentials
				for (int vaccineIndex = 1; vaccineIndex < vsProp.getVaccineCount(); vaccineIndex++) {

					if (!Vaccine.IsSequential(vaccines.get(vaccineIndex))) {
						// if cocktail, then all antigens injected together
						generateCocktail(vaccines.get(vaccineIndex));
						collectEpitopes();
						do {
							executeCycle(trialWriter);
						} while (continueMaturation(vaccines.get(vaccineIndex)));
					}
					
					else {
						// if sequential
						int stepSize = vaccines.get(vaccineIndex).getTimeStep();
						int antigenNb = vaccines.get(vaccineIndex).getAntigenCount();
						int startCycleIndex = vaccines.get(vaccineIndex - 1).getEndTime();
						antigens.clear();
						antigenEpitopes.clear();

						List<Antigen> antigensForVaccine = new ArrayList<Antigen>();
						antigensForVaccine.addAll(
								Vaccine.createAntigens(antigenNb, vaccines.get(vaccineIndex).getMutationalDistance(),
										vaccines.get(vaccineIndex).getEqualAntigenConcentration(), agProp));

						for (int i = 0; i < antigenNb - 1; i++) {
							antigens.add(antigensForVaccine.get(i));
							collectEpitopes();
							do {
								executeCycle(trialWriter);
							} while (continueMaturation(vaccines.get(vaccineIndex))
									&& cycleIndex < startCycleIndex + stepSize * (i + 1));
							antigens.clear();
							antigenEpitopes.clear();
						}
						// The last antigen must stay in the GC until we reach
						// the gc cycle limit of 240 or end of vaccine.
						antigens.add(antigensForVaccine.get(antigenNb - 1));
						do {
							executeCycle(trialWriter);
						} while (continueMaturation(vaccines.get(vaccineIndex)));
					}
				}
			}
		} else {
			throw new IllegalStateException("No Vaccine was loaded");
		}
		return this;
	}

	private void createVaccines() {
		vaccines.addAll(VaccinationStrategy.createVaccines(vsProp, getVaccineConc()));
		// Check size
		if (vsProp.getVaccineCount() != vaccines.size()) {
			throw new IllegalStateException("the number of vaccines loaded is different from vaccineCount");
		}
		if (vaccines.isEmpty()) {
			throw new IllegalStateException("no vaccine was given");
		}
	}

	private void initialization() {
		if (!vaccines.isEmpty()) {
			if (!Vaccine.IsSequential(vaccines.get(0))) {
				generateCocktail(vaccines.get(0));
			} else {
				throw new IllegalStateException("the first vaccine should be a cocktail");
			}
		}
	}

	/**
	 * Creates a list of antigens with length equal to the antigenCount for the
	 * vaccine and mutationalDistance =< 2^(length of the epitope) /
	 * antigenCount. The antigens are created as a cocktail. What counts is
	 * their number and how far away from each other they are. If antigens is
	 * filled and the function is called, antigens is cleared before new
	 * antigens can be added.
	 * 
	 * @param vaccine
	 */
	private void generateCocktail(Vaccine vaccine) {
		if (!antigens.isEmpty()) {
			antigens.clear();
		}
		antigens.addAll(Vaccine.createAntigens(vaccine.getAntigenCount(), vaccine.getMutationalDistance(),
					vaccine.getEqualAntigenConcentration(), agProp));
	}

	private void collectEpitopes() {
		for (Antigen antigen : antigens)
			antigenEpitopes.add(antigen.getEpitope());
	}

	private void generateFounders(Vaccine vaccine) {
		int founderCount = gcProp.getGermlineCount();
		
		while (founderCells.size() < founderCount) {
			BCell founderCell = BCell.createFounder(selectAntigen());
			founderCells.add(founderCell);
		}
		for (BCell founderCell : founderCells){
			founderCell.findMaxEpitopeinList(antigenEpitopes);
		}
		//System.out.println("founders generated");
		//System.exit(1);
	}

	private void replicateFounders() {
		//
		// Must be called on the first cycle...
		//
		if (!activeCells.isEmpty())
			throw new IllegalStateException("Active cells are already present.");

		activeCells.add(newActiveList());

		for (BCell founder : founderCells)
			replicateFounder(founder);
	}

	private static List<BCell> newActiveList() {
		//
		// We use a linked list because we want efficient removal by
		// the iterator during the selection phases...
		//
		return new LinkedList<BCell>();
	}

	private void replicateFounder(BCell founder) {
		int replicationFactor = gcProp.getReplicationFactor();

		for (int replicationIndex = 0; replicationIndex < replicationFactor; replicationIndex++)
			activeCells.get(0).add(BCell.replicate(founder));
	}

	private void executeCycle(PrintWriter trialWriter) {
		cycleIndex++;

		divideActiveCells();
		TipLog.debug("After cell division:   %10d", activeCount());
		
		firstSurvivalSignal();
		TipLog.debug("After antigen binding: %10d", activeCount());

		competeTcellHelp();
		TipLog.debug("After T cell help:     %10d", activeCount());

		recycle();
		TipLog.debug("After recyling:        %10d", activeCount());
		
		
		writeTrialDetailLine(trialWriter);
	}
	
	private void writeTrialDetailLine(PrintWriter trialWriter) {
		DecimalFormat formatter = new DecimalFormat("##0.0#####");

		trialWriter.print(cycleIndex);
		trialWriter.print(",");
		trialWriter.print(activeCount());
		//System.out.println("active B cells " + activeCount());
		trialWriter.print(",");
		if(activeCount() < 1){
			trialWriter.print("all die");
			trialWriter.print(",");
			trialWriter.print("all die");
			trialWriter.print(",");
			trialWriter.print("all die");
			trialWriter.print(",");
			trialWriter.print(gcProp.getConcentration());
			trialWriter.print(",");
			trialWriter.print("all die");
			trialWriter.print(",");
			trialWriter.print("all die");
			trialWriter.print(",");
			trialWriter.print("all die");
			}
		else{
			trialWriter.print(formatter.format(computeMean(resolveDeleteriousMutationsCount(nowActive()))));
			trialWriter.print(",");
			trialWriter.print(formatter.format(computeMean(resolveBeneficialMutationsCount(nowActive())))); 
			trialWriter.print(",");
			trialWriter.print(formatter.format(computeMean(resolveMutationsCount((nowActive()))))); 

			trialWriter.print(",");
			trialWriter.print(gcProp.getConcentration());
			trialWriter.print(",");
			trialWriter.print(formatter.format(computeMean(resolveEnergyChange(nowActive()))));
			trialWriter.print(",");
			trialWriter.print(formatter.format(computeVariance(resolveEnergyChange(nowActive()))));
			trialWriter.print(",");
			trialWriter.print(formatter.format(resolveTotalAffinityChange(nowActive())));
		}
		trialWriter.println();
	}
	
	private ArrayDoubleList resolveDeleteriousMutationsCount(List<BCell> activeBcells){
		ArrayDoubleList array = new ArrayDoubleList();
		array.clear();
		for (BCell bcell : activeBcells){
			array.add(bcell.getDeleteriousMutCount());
		}
		return array;
	}
	
	private ArrayDoubleList resolveBeneficialMutationsCount(List<BCell> activeBcells){
		ArrayDoubleList array = new ArrayDoubleList();
		array.clear();
		for (BCell bcell : activeBcells){
			array.add(bcell.getBeneficialMutCount());
		}
		return array;
	}
	
	private ArrayDoubleList resolveMutationsCount(List<BCell> activeBcells){
		ArrayDoubleList array = new ArrayDoubleList();
		array.clear();
		for (BCell bcell : activeBcells){
			array.add(bcell.getMutationCount());
		}
		return array;
	}
	
	private ArrayDoubleList resolveEnergyChange(List<BCell> activeBcells){
		ArrayDoubleList array = new ArrayDoubleList();
		array.clear();
		for (BCell bcell : activeBcells){
			array.add(bcell.getEnergyChange());
			//System.out.println(bcell.getEnergyChange());
		}
		return array;
	}
	
	private double resolveTotalAffinityChange(List<BCell> activeBcells){
		double totalAffinity = BCell.resolveTotalAffinity(activeBcells);
		
		List<BCell> founderCellsList = new ArrayList<BCell>(founderCells);
		for (BCell found : founderCellsList){
			//System.out.println(java.util.Arrays.toString(found.getReceptor().getCoord()));
			//System.out.println(found.getMaxEnergy());
		}
		//System.out.println(founderCellsList.size());
		double totalAffinityInit = BCell.resolveTotalAffinity(founderCellsList);
		
		//System.out.println("initial Affinity " + totalAffinityInit);
		//System.out.println(totalAffinity);
		
		return totalAffinity/totalAffinityInit;
	}
	

	/**
	 * Computes the mean of a list of doubles
	 *
	 * @return the mean of a list of doubles.
	 */
	public double computeMean(ArrayDoubleList values) {
		double[] out = new double[values.size()];
		values.toArray(out);
		return MeanEstimator.estimate(out).getMean();
	}
	
	/**
	 * Computes the variance of a list of doubles
	 *
	 * @return the variance of a list of doubles.
	 */
	public double computeVariance(ArrayDoubleList values) {
		double[] out = new double[values.size()];
		values.toArray(out);
		return Math.pow(MeanEstimator.estimate(out).getSD(), 2);
	}


	private int activeCount() {
		return nowActive().size();
	}

	private List<BCell> nowActive() {
		return activeCells.get(cycleIndex);
	}

	private void divideActiveCells() {
		if (cycleIndex == 0)
			throw new IllegalStateException("Cannot divide on first cycle.");

		List<BCell> parents = activeCells.get(cycleIndex - 1);

		List<BCell> daughters = newActiveList();
		List<BCell> intermediates = new ArrayList<BCell>();

		// 1st division
		for (BCell parent : parents){
			intermediates.addAll(BCell.divide(parent));
		}

		// 2nd division
		for (BCell interm : intermediates){
			daughters.addAll(BCell.divide(interm));
		}

		activeCells.add(cycleIndex, daughters);
	}

	private void firstSurvivalSignal() {
		Iterator<BCell> iterator = nowActive().iterator();

		while (iterator.hasNext()) {
			firstSurvivalSignal(iterator);
		}
	}

	private void bindOneAntigen(BCell activeCell) {
		Antigen seenAntigen = selectAntigen();
		activeCell.bindAntigen(seenAntigen);
		activeCell.setEpitopeKey(seenAntigen.getEpitope());
	}

	private void bindAllAntigens(BCell activeCell) {
		for (Antigen antigen : antigens){
			activeCell.bindAntigen(antigen);
		}
		activeCell.setEpitopeKey(activeCell.getMaxEpitope());//epitopeKey is the maxEpitope in this case
		//activeCell.setEpitopeKey(Receptor.findTarget(antigenEpitopes, activeCell.getReceptor())); 
	}

	private void firstSurvivalSignal(Iterator<BCell> iterator) {
		BCell activeCell = iterator.next();
		activeCell.findMaxEpitopeinList(antigenEpitopes);//updates the maxEpitope and maxEnergy for the activeCell
		int numberOfPresentedAntigens = gcProp.getNumberOfAntigensPresented();
		qtyAntigens(activeCell, numberOfPresentedAntigens);

		double ratioOfOccupiedSites = activeCell.getAntigenQty() / (activeCell.getAntigenQty() +1); 

		if (ratioOfOccupiedSites > 1 || ratioOfOccupiedSites < 0) {
			throw new IllegalStateException("ratioOfOccupiedSites not within  0 and 1");
		}
		//System.out.println(activeCell.getAntigenQty());
		
		if (!random.accept(ratioOfOccupiedSites))
            iterator.remove();
	}

	private Antigen selectAntigen() {
		return ListUtil.select(antigens, random);
	}

	// Updates antigenQty for the given B cell
	private void qtyAntigens(BCell activeCell, int numberOfPresentedAntigens) {
		if (numberOfPresentedAntigens != 1) {
			bindAllAntigens(activeCell);
		} 
		else {
			bindOneAntigen(activeCell);
		}
	}

	private TipVector computeBoundAntigenQtyForAllBcells() {
		TipVector sum = null;

		for (BCell bcell : nowActive())
			sum = sum.add(bcell.getAntigenQty());

		return sum;
	}

	private void competeTcellHelp() {
		//
		//Use a hard threshold. The T cells select the B cells
		// with highest probability of binding antigens
        //
        // Keep the fraction R with the highest amount, where R is the
        // T cell selection rate...
        //
        double selectionRate = gcProp.getTcellSelectionRate();
        double deletionRate  = 1.0 - selectionRate;
        int    deletionCount = (int) (activeCount() * deletionRate);

        if (deletionCount > 0) {
            //
            // Sort B cells by the amount of internalized antigen
            // (with the lowest amount first) and then remove the
            // required number from the list...
            //
            Collections.sort(nowActive(), BCell.ANTIGEN_QTY_COMPARATOR);
            nowActive().subList(0, deletionCount).clear();
        }
        
//        
//		Iterator<BCell> iterator = nowActive().iterator();
//
//		while (iterator.hasNext())
//			competeTcellHelp(iterator);
//

	}

//	private void competeTcellHelp(Iterator<BCell> iterator) {
//		BCell activeCell = iterator.next();
//		double antigenQty = activeCell.getAntigenQty();
//		double averageQtyOtherBcells = (computeBoundAntigenQtyForAllBcells() - antigenQty) / (activeCount() - 1);
//
//		double probaTcellHelp = antigenQty / (antigenQty + averageQtyOtherBcells / gcProp.getConcentration());
//
//		// Check it's a probability
//		if (probaTcellHelp > 1 || probaTcellHelp < 0) {
//			throw new IllegalStateException("probaTcellHelp not within  0 and 1");
//		}
//
//		if (!random.accept(probaTcellHelp))
//			iterator.remove();
//	}

	private void recycle() {
		Iterator<BCell> iterator = nowActive().iterator();

		while (iterator.hasNext()) {
			BCell activeCell = iterator.next();

			if (selectPlasma(activeCell)) {
				iterator.remove();
				addPlasmaCell(activeCell);
			} else if (selectMemory(activeCell)) {
				iterator.remove();
				addMemoryCell(activeCell);
			}
		}
	}

	private boolean selectPlasma(BCell activeCell) {
		double selectionRate = gcProp.getPlasmaSelectionRate();
		double affinityThreshold = gcProp.getPlasmaAffinityThreshold();

		return (activeCell.getMaxEnergy() >= affinityThreshold && random.accept(selectionRate));
	}

	private boolean selectMemory(BCell activeCell) {
		//
		// In our model, the selection probability is independent of
		// the active cell characteristics...
		//
		double selectionRate = gcProp.getMemorySelectionRate();
		return (random.accept(selectionRate));
	}

	private void addPlasmaCell(BCell plasmaCell) {
		plasmaCells.add(plasmaCell);
		targetEpitopes.put(plasmaCell, plasmaCell.getMaxEpitope());
		//targetEpitopes.put(plasmaCell, Receptor.findTarget(antigenEpitopes, plasmaCell.getReceptor()));
		viralChallenges.put(plasmaCell, ViralChallenge.challenge(plasmaCell));
	}

	private void addMemoryCell(BCell memoryCell) {
		//
		// Not much to do with memory cells yet, just keep track...
		//
		memoryCells.add(memoryCell);
	}

	public boolean IsExtinguished() {
		return this.isExtinguished;
	}

	private boolean continueMaturation(Vaccine vaccine) {
		int activeSize = activeCount();
		BCellProp bcProp = BCellProp.instance();
		int initialSize = bcProp.getTerminationMultiple()*activeCells.get(0).size();

		if (cycleIndex >= gcProp.getCycleLimit()) {
			TipLog.info("Reached cycle limit.");
			logState();
			return false;
		}

		if (activeSize == 0) {
			TipLog.info("GC extinguished.");
			logState();
			this.isExtinguished = true;
			return false;
		}

		if (activeSize > initialSize) {
			TipLog.info("Maturation complete (initial size = %d).", initialSize);
			logState();
			return false;
		}

		if (cycleIndex >= vaccine.getEndTime()) {
			TipLog.info("End of vaccine");
			logState();
			return false;
		}

		return true;
	}

	private void logState() {
		TipLog.info("CYCLE:  %8d", cycleIndex);
		TipLog.info("ACTIVE: %8d", activeCount());
		TipLog.info("PLASMA: %8d", plasmaCells.size());
		TipLog.info("MEMORY: %8d", memoryCells.size());
	}

	public static void main(String[] args) throws IOException {
		String propName = args[0];
		File propFile = new File(propName);

		TipProperties.loadFile(propFile);
		TipRandom.initialize(TipProperties.getRequiredInt("AMDriver.randomSeed"));
		run(null);
	}
}

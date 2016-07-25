
package tip.maam3;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

import tip.maam.MutationType;
import tip.math.TipRandom;
import tip.math.TipVector;
import tip.util.ObjectUtil;

/**
 * Represents a B cell in an affinity maturation simulation.
 *
 * <p>
 * The receptor shape, parent cell, and founder cell are immutable, but the
 * bound antigen amount may change during affinity maturation.
 */
public final class BCell {
	private final BCell parent;
	private final BCell founder;
	private final Receptor receptor;
	private final int generation;
	private final int mutationCount;
	private final int deleteriousMutationCount;
	private final int beneficialMutationCount;

	private Epitope epitopeKey; // the antigen that we use for the mutation
								// process

	// Probability to bind and internalize one or several antigens 
	private double antigenQty = 0.0;

	// Strongest binder among the antigens encountered...
	private Epitope maxEpitope = null;

	// Energy for the strongest binder...
	private double maxEnergy = Double.NEGATIVE_INFINITY;

	// Computed on demand and cached...
	private int hashCode = 0;

	/**
	 * Places B cells in INCREASING ORDER (least first) of bound antigen.
	 */
	public static final Comparator<BCell> ANTIGEN_QTY_COMPARATOR = new AntigenQtyComparator();

	private static class AntigenQtyComparator implements Comparator<BCell> {
		@Override
		public int compare(BCell cell1, BCell cell2) {
			if  (cell1.getAntigenQty() == cell2.getAntigenQty())
				return 1;
			else
				return 0;
		}
	}

	private BCell(Receptor receptor, Epitope epitopeKey) {
		//
		// Creates a founder B cell...
		//
		this.parent = null;
		this.founder = this;
		this.receptor = receptor;
		this.generation = 0;
		this.mutationCount = 0;
		this.beneficialMutationCount = 0;
		this.deleteriousMutationCount = 0;
		this.epitopeKey = epitopeKey;
	}


	private BCell(BCell parent, Receptor receptor) {
		//
		// These attributes must be assigned first...
		//
		this.parent = parent;
		this.receptor = receptor;

		// Then the remainder are derived from the parent and receptor...
		this.founder = resolveFounder();
		this.generation = resolveGeneration();
		this.mutationCount = resolveMutationCount();
		this.beneficialMutationCount = resolveBeneficialMutationCount();
		this.deleteriousMutationCount = resolveDeleteriousMutationCount();
		
		//check mutation >= beneficial + deleterious
		checkMutations();
		
		this.epitopeKey = resolveEpitopeKey();
	}

	private void checkMutations(){
		if(mutationCount < beneficialMutationCount + deleteriousMutationCount){
			throw new IllegalStateException("total number of mutations smaller than the sum of beneficical and deleterious mutations");
		}
	}
	
	private Epitope resolveEpitopeKey() {
		if (parent == null)
			return null;
		else
			return parent.getEpitopeKey();
	}


	private BCell resolveFounder() {
		if (parent == null)
			return this;
		else
			return parent.getFounder();
	}

	private int resolveGeneration() {
		if (parent == null)
			return 0;
		else
			return parent.getGeneration() + 1;
	}

	private int resolveMutationCount() {
		if (parent == null)
			return 0;
		else if (receptor.equals(parent.getReceptor()))
			return parent.getMutationCount();
		else
			return parent.getMutationCount() + 1;
	}
	
	private int resolveDeleteriousMutationCount(){
		if (parent == null)
			return 0;
		else if (this.getReceptor().calculateEnergy(parent.getEpitopeKey()) 
				< parent.getReceptor().calculateEnergy(parent.getEpitopeKey())){
			return parent.getDeleteriousMutCount() + 1;
		}
		else 
			return parent.getDeleteriousMutCount();
	}

	private int resolveBeneficialMutationCount(){
	if (parent == null)
		return 0;
	else if (this.getReceptor().calculateEnergy(parent.getEpitopeKey()) 
			> parent.getReceptor().calculateEnergy(parent.getEpitopeKey())){
		return parent.getBeneficialMutCount() + 1;
	}
	else 
		return parent.getBeneficialMutCount();
}
	/**
	 * Creates a founder (germline) B cell according to the global antigen and B
	 * cell properties.
	 * 
	 * @param selected antigen from the initial vaccine.
	 *
	 * @return a founder (germline) B cell.
	 */
	public static BCell createFounder(Antigen antigen) {
		BCellProp bcProp = BCellProp.instance();
		Epitope epitope = antigen.getEpitope();
		return new BCell(Receptor.createFounder(bcProp, epitope), epitope);
	}

	/**
	 * Sorts a collection of B cells by generation.
	 *
	 * @param cells
	 *            a collection of B cells.
	 *
	 * @return a list {@code G} for which {@code G.get(k)} returns another list
	 *         containing all B cells from generation {@code k};
	 *         {@code G.size() - 1} is the index of the highest generation.
	 */
	public static List<List<BCell>> sortGenerations(Collection<BCell> cells) {
		List<List<BCell>> generations = new ArrayList<List<BCell>>();

		for (BCell cell : cells) {
			//
			// Add new lists as new generations are encountered for
			// the first time...
			//
			int genIndex = cell.getGeneration();

			while (generations.size() <= genIndex) {
				generations.add(new ArrayList<BCell>());
			}
			generations.get(genIndex).add(cell);
		}

		return generations;
	}

	/**
	 * Creates an identical daughter cell (without mutation).
	 *
	 * @return an identical daughter cell.
	 */
	public static BCell replicate(BCell parent) {
		return new BCell(parent, parent.getReceptor());
	}

	/**
	 * Creates a mutated daughter cell.
	 *
	 * @param mutationAffinity
	 *            the average change in binding affinity to accompany the
	 *            mutation.
	 *
	 * @return a mutated daughter cell.
	 *
	 * @see Receptor#mutate(double)
	 */
	public static BCell mutate(BCell parent, Epitope target, BCellProp bcProp, AntigenProp agProp) {
		return new BCell(parent, Receptor.mutate(parent, target, bcProp, agProp));
	}

	/**
	 * Simulates cell division with stochastic mutation.
	 *
	 * <p>
	 * This method creates two initially identical daughter cells and then
	 * subjects the daughter cells to stochastic mutation. Mutations occur with
	 * a probability specified by the {@code
	 * BCellProp} class. If a mutation occurs, its type (silent, lethal, or
	 * affinity-affecting) is selected at random and the outcome is applied
	 * accordingly.
	 *
	 * @return a list containing the daughter cells; the list may have length 0,
	 *         1, or 2 since lethal mutations may occur.
	 */
	public static List<BCell> divide(BCell parent) {
		List<BCell> daughters = new ArrayList<BCell>();

		BCell daughter1 = createDaughter(parent);
		BCell daughter2 = createDaughter(parent);

		if (daughter1 != null)
			daughters.add(daughter1);

		if (daughter2 != null)
			daughters.add(daughter2);

		return daughters;
	}

	private static BCell createDaughter(BCell parent) {
		if (acceptMutation())
			return mutateCell(parent);
		else
			return replicate(parent);
	}

	private static boolean acceptMutation() {
		return TipRandom.instance().accept(BCellProp.instance().getMutationProbability());
	}

	private static BCell mutateCell(BCell parent) {
		MutationType mutationType = MutationType.select();

		switch (mutationType) {
		case AFFINITY_AFFECTING:
			return mutate(parent, parent.getEpitopeKey(), BCellProp.instance(), AntigenProp.instance());

		case LETHAL:
			return null;

		case SILENT:
			return replicate(parent);
		}

		throw new IllegalStateException("Unknown mutation type.");
	}

	/**
	 * Simulates the binding of this B cell to an antigen; the amount of bound
	 * antigen is increased by the amount bound during this encounter.
	 *
	 * @param antigen
	 *            the presented antigen.
	 */
	public void bindAntigen(Antigen antigen) {
		Receptor Bcellreceptor = this.getReceptor(); 
		double energy = Bcellreceptor.calculateEnergy(antigen.getEpitope());
		double bindRate = computeBindingRate(energy, antigen.getConcentration());
		antigenQty = antigenQty + bindRate;
	}
	
	public void findMaxEpitopeinList(Collection<Epitope> epitopes){
		double energy = Double.NEGATIVE_INFINITY;
		for (Epitope epitope : epitopes){
			energy = receptor.calculateEnergy(epitope);
			if (energy > maxEnergy){
				maxEpitope = epitope;
				maxEnergy = energy;
			}
		}
	}

	private double computeBindingRate(double energy, double concentration) {
		BCellProp bcProp = BCellProp.instance();
		return concentration * computeBindingConstant(energy); 
		//return concentration * computeBindingConstant(energy - bcProp.getActivationThreshold());
	}

	private double computeBindingConstant(double energy) {
		return Math.exp(energy);
	}
	
	public double getEnergyChange(){
		double energy = this.getMaxEnergy();
		double energyInit = this.getFounder().getMaxEnergy();//getReceptor().calculateEnergy(this.getFounder().getMaxEpitope());
		if (energyInit !=0){
			return energy/energyInit;
		}
		else{
			return Double.NaN;//energy;
		}
	}
	
	public static double resolveTotalAffinity(List<BCell> bcells){
		double affinity = 0;
		for (BCell bcell : bcells){
			affinity = affinity + bcell.computeBindingConstant(bcell.getMaxEnergy());
		}
		return affinity;
	}

	/**
	 * Returns the founder cell from which this B cell is derived.
	 *
	 * @return the founder cell from which this B cell is derived.
	 */
	public BCell getFounder() {
		return founder;
	}

	/**
	 * Identifies founder cells.
	 *
	 * @return {@code true} iff this cell is a founder cell.
	 */
	public boolean isFounder() {
		return founder == this;
	}

	/**
	 * Returns the number of cell divisions that have occurred since the founder
	 * cell.
	 *
	 * @return the number of cell divisions that have occurred since the founder
	 *         cell.
	 */
	public int getGeneration() {
		return generation;
	}

	/**
	 * Returns the number of mutations this B cell has accumulated from its
	 * originating founder cell.
	 *
	 * @return the number of mutations this B cell has accumulated from its
	 *         originating founder cell.
	 */
	public int getMutationCount() {
		return mutationCount;
	}

	/**
	 * Returns the parent of this B cell ({@code null} for founder cells).
	 *
	 * @return the parent of this B cell ({@code null} for founder cells).
	 */
	public BCell getParent() {
		return parent;
	}

	/**
	 * Returns the receptor shape for this B cell.
	 *
	 * @return the receptor shape for this B cell.
	 */
	public Receptor getReceptor() {
		return receptor;
	}

	/**
	 * Returns the (dimensionless) quantity of antigen internalized by this B
	 * cell.
	 *
	 * @return the (dimensionless) quantity of antigen internalized by this B
	 *         cell.
	 */
	public double getAntigenQty() {
		return antigenQty;
	}

	/**
	 * Returns the energy for the strongest binding epitope.
	 *
	 * @return the energy for the strongest binding epitope.
	 */
	public double getMaxEnergy() {
		return maxEnergy;
	}

	/**
	 * Returns the strongest binding epitope among those encountered by this B
	 * cell.
	 *
	 * @return the strongest binding epitope among those encountered by this B
	 *         cell.
	 */
	public Epitope getMaxEpitope() {
		return maxEpitope;
	}

	public Epitope getEpitopeKey() {
		return epitopeKey;
	}

	/**
	 * Extracts the receptors from a collection of B cells.
	 *
	 * @param cells
	 *            the B cells to process.
	 *
	 * @return a list where element {@code k} is the receptor for cell {@code k}
	 *         in the input collection.
	 */
	public static List<Receptor> getReceptors(Collection<BCell> cells) {
		List<Receptor> receptors = new ArrayList<Receptor>(cells.size());

		for (BCell cell : cells)
			receptors.add(cell.getReceptor());

		return receptors;
	}

	@Override
	public boolean equals(Object that) {
		return (that instanceof BCell) && equalsBCell((BCell) that);
	}

	private boolean equalsBCell(BCell that) {
		return ObjectUtil.equals(this.parent, that.parent) && ObjectUtil.equals(this.receptor, that.receptor)
				&& this.generation == that.generation && this.mutationCount == that.mutationCount;
	}

	@Override
	public int hashCode() {
		if (hashCode == 0)
			hashCode = computeHashCode();

		return hashCode;
	}

	private int computeHashCode() {
		return Arrays.asList(parent, receptor, generation, mutationCount).hashCode();
	}

	@Override
	public String toString() {
		return String.format("BCell(%d, %s)", generation, receptor.toString());
	}

	public void setEpitopeKey(Epitope seenEpitope) {
		this.epitopeKey = seenEpitope;
	}

	public void setEpitopeMax(Epitope maxEnergyEpitope) {
		this.maxEpitope = maxEnergyEpitope;
	}
	
	public int getDeleteriousMutCount() {
		return deleteriousMutationCount;
	}

	public int getBeneficialMutCount() {
		return beneficialMutationCount;
	}
}

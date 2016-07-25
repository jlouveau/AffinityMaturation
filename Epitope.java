
package tip.maam3;

import java.util.List;

import tip.math.BitVector;
import tip.math.TipRandom;

/**
 * Represents an antigen epitope as a vector of 0 and 1;
 * 
 * 1 for unmutated bits and -1 for mutated bits; the first conservedLength bits
 * are always 1. the following variableLength bits can be mutated.
 *
 * <p>
 * Epitopes are immutable: once created, their properties are fixed forever.
 */
public final class Epitope {
	private final BitVector coord;
	private final int length;
	private final int conservedLength;

	private Epitope(BitVector coord, int length, int conservedLength) {
		this.coord = coord;
		this.length = length;
		this.conservedLength = conservedLength;
	}

	/**
	 * Creates a conserved epitope for the founder antigen. All bits are "false"
	 * = unmutated
	 *
	 * @param agProp.
	 *
	 * @return the founder epitope.
	 */
	public static Epitope createConserved(AntigenProp agProp) {
		BitVector coord = new BitVector(agProp.getLength());

		return new Epitope(coord, agProp.getLength(), agProp.getConservedLength());
	}

	/**
	 * Creates a random epitope. Starts with all bits unmutated or false. then
	 *
	 * @param agProp.
	 *
	 * @return the founder epitope.
	 */
	public static Epitope createRandom(AntigenProp agProp) {
		BitVector conservedCoord = new BitVector(agProp.getConservedLength());
//		System.out.println("conserved part " + conservedCoord.toString());
		BitVector variableCoord = BitVector.random(agProp.getVariableLength(), TipRandom.instance());
//		System.out.println("variable part " + variableCoord.toString());
		BitVector coord = BitVector.cat(conservedCoord, variableCoord);

		// Check size
		if (coord.length() != agProp.getLength()) {
			throw new IllegalStateException("problem when concatenating conserved and variable bits");
		}

		return new Epitope(coord, agProp.getLength(), agProp.getConservedLength());
	}

	/**
	 * Creates a variant epitope from a copy of a founder epitope. Then add
	 * mutationalDistance mutations.
	 *
	 * @param epitopeDim
	 *            the dimensionality of the shape space.
	 *
	 * @return the conserved epitope.
	 */
	public static Epitope createVariant(double mutationalDistance, Antigen founder) {
		int nbMutations = 0;
		BitVector newCoord = founder.getEpitope().getCoord().copy();

		do {
			int randomNum = founder.getEpitope().getConservedLength() + TipRandom.instance().nextInt(founder.getEpitope().getVariableLength());
					//+ (int) (Math.random() * (founder.getEpitope().getVariableLength()));
			newCoord.flip(randomNum);
			nbMutations++;
		} while (nbMutations < mutationalDistance);

		return new Epitope(newCoord, founder.getEpitope().getLength(), founder.getEpitope().getConservedLength());
	}

	/**
	 * Check if the variant epitope is equal to the epitope of a given antigen.
	 * 
	 * @param variant
	 *            epitope
	 * 
	 * @param antigen
	 * 
	 * @return boolean for epitope equal to the epitope of the antigen
	 */
	public static boolean hasDuplicate(Epitope variant, Antigen antigen) {
		BitVector coordinates = antigen.getEpitope().getCoord();
		if (coordinates.equals(variant.getCoord())) {
			return true;
		}
		return false;
	}

	public static boolean hasDuplicateInList(Epitope variant, List<Antigen> antigens) {
		for (int i = 0; i < antigens.size(); i++) {
			if (hasDuplicate(variant, antigens.get(i))) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Returns the total number of bits in this epitope.
	 *
	 * @return the total number of bits in this epitope.
	 */
	public int getLength() {
		return length;
	}

	/**
	 * Returns the number of conserved bits for this epitope.
	 *
	 * @return the number of conserved bits for this epitope.
	 */
	public int getConservedLength() {
		return conservedLength;
	}

	/**
	 * Returns the indexes of the conserved bits in this epitope.
	 *
	 * @return the indexes of the conserved bits in this epitope.
	 */
	public int[] getConservedIndexes() {
		int[] indexes = new int[getConservedLength()];

		for (int k = 0; k < indexes.length; k++)
			indexes[k] = k;

		return indexes;
	}

	/**
	 * Returns the number of variable bits for this epitope.
	 *
	 * @return the number of variable bits for this epitope.
	 */
	public int getVariableLength() {
		return getLength() - getConservedLength();
	}

	/**
	 * Returns the indexes of the variable bits in this epitope.
	 *
	 * @return the indexes of the variable bits in this epitope.
	 */
	public int[] getVariableIndexes() {
		int[] indexes = new int[getVariableLength()];

		for (int k = 0; k < indexes.length; k++)
			indexes[k] = k + getConservedLength();

		return indexes;
	}

	/**
	 * Returns the coordinates of this epitope.
	 *
	 * @return the coordinates of this epitope.
	 */
	public BitVector getCoord() {
		return coord;
	}

	@Override
	public String toString() {
		return String.format("Epitope(%s)", getCoord().toString());
	}
}

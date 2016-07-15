
package tip.maam3;

import tip.app.TipProperties;

/**
 * Defines global antigen properties:
 *
 * <ul>
 * <li>{@code Antigen.conservedLength}: Number of unmutated bits on the antigen.
 * </li>
 * <li>{@code Antigen.variableLength}: Number of variable bits.</li>
 * </ul>
 *
 * <p>
 * The properties must be assigned before accessing the global instance.
 */
public final class AntigenProp {
	private final int conservedLength;
	private final int variableLength;

	private static AntigenProp instance = null;

	private AntigenProp() {
		this.conservedLength = TipProperties.getRequiredInt("Antigen.conservedLength");
		this.variableLength = TipProperties.getRequiredInt("Antigen.variableLength");

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
	public static AntigenProp instance() {
		if (instance == null)
			instance = new AntigenProp();

		return instance;
	}

	/**
	 * Returns the total number of bits.
	 *
	 * @return the total number of bits.
	 */
	public int getLength() {
		return conservedLength + variableLength;
	}

	/**
	 * Returns the number of conserved bits.
	 *
	 * @return the number of conserved bits.
	 */
	public int getConservedLength() {
		return conservedLength;
	}

	/**
	 * Returns the number of variable bits.
	 *
	 * @return the number of variable bits.
	 */
	public int getVariableLength() {
		return variableLength;
	}
}

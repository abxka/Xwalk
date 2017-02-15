/*
 * (C) 2010 Abdullah Kahraman
 *
 * This software is part of the open-source project "Xwalk". You can use this
 * software under the terms of the
 * Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 * (http://creativecommons.org/licenses/by-nc-sa/3.0/).
 * This means that you
 * 1.) can copy, modify, distribute the software
 * 2.) must give credit to the author
 * 3.) must not use this work for commercial purposes
 * 4.) must license derivative works under the same or a similar license.
 *
 */

package structure.matter.parameter;

import structure.constants.Constants.ElementTypes;

/**
 * Element types as found in the periodic table.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public enum Element {
    // ORGANIC PROTEIN ATOMS
    CARBON("C", ElementTypes.ORGANIC), NITROGEN("N", ElementTypes.ORGANIC),
    OXYGEN("O", ElementTypes.ORGANIC), SULPHUR("S", ElementTypes.ORGANIC),
    HYDROGEN("H", ElementTypes.ORGANIC),

    // ORGANIC SMALL MOLECULE ATOMS
    PHOSPHOR("P", ElementTypes.ORGANIC),

    // METAL ATOMS
    POTASSIUM("K", ElementTypes.METAL), MAGNESIUM("MG", ElementTypes.METAL),
    MANGANESE("MN", ElementTypes.METAL), CALCIUM("CA", ElementTypes.METAL),
    IRON("FE", ElementTypes.METAL), SODIUM("NA", ElementTypes.METAL),
    ZINC("ZN", ElementTypes.METAL), CADMIUM("CD", ElementTypes.METAL),
    URANIUM("U", ElementTypes.METAL), LUTETIUM("LU", ElementTypes.METAL),
    LITHIUM("LI", ElementTypes.METAL), SILICON("SI", ElementTypes.METALLOID),
    COPPER("CU", ElementTypes.METAL), ALUMINUM("AL", ElementTypes.METAL),
    MERCURY("HG", ElementTypes.METAL), BERYLLIUM("BE", ElementTypes.METAL),

    // NON-METAL ATOMS
    IODINE("I", ElementTypes.NON_METAL),FLUORINE("F", ElementTypes.NON_METAL),
    CHLORINE("CL", ElementTypes.NON_METAL),
    BROMINE("BR", ElementTypes.NON_METAL),
    SELENIUM("SE", ElementTypes.NON_METAL);

/*    // SPECIAL ELEMENT TYPES FOR PARSE
    CARBON_AROMATIC("CB", ElementTypes.ORGANIC),
    HYDROGEN_ALIPHATIC("HL", ElementTypes.ORGANIC),

    // SPECIAL ELEMENT TYPES FOR CHARMM
    HYDROGEN_BOUND_TO_CARBON("HC", ElementTypes.ORGANIC),
    HYDROGEN_BOUND_TO_OXYGEN_OR_NITROGEN("HO", ElementTypes.ORGANIC),
    HYDROGEN_BOUND_TO_SULPHUR("HS", ElementTypes.ORGANIC);
*/
    //--------------------------------------------------------------------------

    /**
     * One or two letter code of the element.
     */
    private String elementSymbol;
    /**
     * Element type of the element, i.e. whether organic, metallic etc.
     */
    private ElementTypes elementType;
    //--------------------------------------------------------------------------

    /**
     * Constructor.
     * @param symbol
     *    - One or two letter code of the element.
     * @param type
     *    - Element type of the element
     */
    Element(final String symbol, final ElementTypes type) {
        this.elementSymbol = symbol;
        this.elementType = type;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the symbol of the element as found in a periodic table.
     * @return String object holding the symbol of the element.
     */
    public String getSymbol() {
        return this.elementSymbol;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the type of the element.
     * @return ElementTypes object.
     */
    public ElementTypes getType() {
        return this.elementType;
    }
}

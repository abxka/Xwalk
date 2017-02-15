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

/**
 * Supported atom enum types for cross-linking.
 * @author Abdullah
 * @version 0.1
 * @since 0.1
 *
 */
public enum AtomType {
    CARBON(Element.CARBON, "C"),
    CARBON_ALPHA(Element.CARBON, "CA"),
    CARBON_BETA(Element.CARBON, "CB"),
    CARBON_GAMMA(Element.CARBON, "CG"),
    CARBON_DELTA(Element.CARBON, "CD"),
    CARBON_DELTA1(Element.CARBON, "CD1"),
    CARBON_DELTA2(Element.CARBON, "CD2"),
    CARBON_GAMMA1(Element.CARBON, "CG1"),
    CARBON_GAMMA2(Element.CARBON, "CG2"),
    CARBON_EPSILON(Element.CARBON, "CE"),
    CARBON_EPSILON1(Element.CARBON, "CE1"),
    CARBON_EPSILON2(Element.CARBON, "CE2"),
    CARBON_EPSILON3(Element.CARBON, "CE3"),
    CARBON_Z(Element.CARBON, "CZ"),
    CARBON_Z2(Element.CARBON, "CZ2"),
    CARBON_Z3(Element.CARBON, "CZ3"),
    CARBON_H2(Element.CARBON, "CH2"),

    NITROGEN(Element.NITROGEN, "N"),
    NITROGEN_DELTA1(Element.NITROGEN, "ND1"),
    NITROGEN_DELTA2(Element.NITROGEN, "ND2"),
    NITROGEN_EPSILON(Element.NITROGEN, "NE"),
    NITROGEN_EPSILON1(Element.NITROGEN, "NE1"),
    NITROGEN_EPSILON2(Element.NITROGEN, "NE2"),
    NITROGEN_H1(Element.NITROGEN, "NH1"),
    NITROGEN_H2(Element.NITROGEN, "NH2"),
    NITROGEN_Z(Element.NITROGEN, "NZ"),

    OXYGEN(Element.OXYGEN, "O"),
    OXYGEN_DELTA1(Element.OXYGEN, "OD1"),
    OXYGEN_DELTA2(Element.OXYGEN, "OD2"),
    OXYGEN_EPSILON1(Element.OXYGEN, "OE1"),
    OXYGEN_EPSILON2(Element.OXYGEN, "OE2"),
    OXYGEN_GAMMA(Element.OXYGEN, "OG"),
    OXYGEN_GAMMA1(Element.OXYGEN, "OG1"),
    OXYGEN_HYDROXYL(Element.OXYGEN, "OH"),

    SULPHUR_GAMMA(Element.SULPHUR, "SG"),
    SULPHUR_DELTA(Element.SULPHUR, "SD");

    //--------------------------------------------------------------------------
    /**
     * Element to which atom belongs.
     */
    private Element elementType;
    /**
     * Abbreviation, e.g. PDB short-name of atom.
     */
    private String shortName;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param element
     *        - Element object holding the element information to which atom
     *          belongs.
     * @param abbreviation
     *        - String object holding the short name of the atom, e.g. PDB
     *          short name
     */
    AtomType(final Element element, final String abbreviation) {
        this.elementType = element;
        this.shortName = abbreviation;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the associated element type of this atom type.
     * @return Element object representing element type.
     */
    public Element getElement() {
        return this.elementType;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the short name of a PDB atom of this atom type.
     * @return String object holding the short name.
     */
    public String getAbbreviation() {
        return this.shortName;
    }
}

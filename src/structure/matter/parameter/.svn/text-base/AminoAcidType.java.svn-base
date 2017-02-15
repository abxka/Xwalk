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
 * AminoAcid Types that are supported by Xwalk.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public enum AminoAcidType {

    // hydrophobic amino acids
    GLYCINE ('G',"GLY"),
    ALANINE ('A', "ALA"),
    VALINE ('V', "VAL"),
    LEUCINE ('L', "LEU"),
    ISOLEUCINE ('I', "ILE"),
    PHENYLALANINE ('F', "PHE"),
    PROLINE ('P', "PRO"),
    TRYPTOPHANE ('W', "TRP"),
    // amino acids containing a nitrogen atom
    ASPARAGINE ('N', "ASN"),
    GLUTAMINE ('Q', "GLN"),
    // amino acids containing a sulphur atom
    METHIONINE ('M', "MET"),
    CYSTEINE ('C', "CYS"),
    // amino acids containing a sulphur atom
    THREONINE ('T', "THR"),
    TYROSINE ('Y', "TYR"),
    SERINE ('S', "SER"),
    // positively charged amino acids
    ARGININE ('R', "ARG"),
    LYSINE ('K', "LYS"),
    HISTIDINE ('H', "HIS"),
    // negatively charged amino acids
    ASPARTIC_ACID ('D', "ASP"),
    GLUTAMIC_ACID ('E', "GLU"),

    // special type amino acids
    ASPARAGINE_ACID ('B', "ASX"),
    GLUTAMINE_ACID ('Z', "GLX");

    /**
     * One letter code of the amino acid.
     */
    private final char oneLetterCode;
    /**
     * Three letter code of the amino acid.
     */
    private final String threeLetterCode;

    /**
     * Constructor.
     * @param oneLetterCode
     *        - Sting object representing the one letter code of the amino acid.
     * @param threeLetterCode
     *        - Sting object representing the three letter code of the amino
     *          acid.
     */
    AminoAcidType(final char oneLetterCode, final String threeLetterCode) {
        this.oneLetterCode = oneLetterCode;
        this.threeLetterCode = threeLetterCode;
    }

    /**
     * Returns the one letter code of an amino acid.
     * @return character representing the one letter code.
     */
    public char getOneLetterCode() {
        return this.oneLetterCode;
    }

    /**
     * Returns the three letter code of an amino acid.
     * @return Sting object representing the three letter code.
     */
    public String getThreeLetterCode() {
        return this.threeLetterCode;
    }
}

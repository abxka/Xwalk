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

package xwalk.constants;

import java.text.DecimalFormat;
import java.text.NumberFormat;

import structure.matter.parameter.AminoAcidType;

/**
 * Class holding generic constant values related to Xwalk execution.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 *
 */
public class Constants {
    /**
     * Constructor.
     */
    protected Constants() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    // CONSTANT NUMERIC VALUES
    //--------------------------------------------------------------------------
    /**
     * Maximum dimension of protein complex after which local grid
     * calculation is used.
     */
    public static final int MAX_PROTEIN_DIMENSION = 150;
    //--------------------------------------------------------------------------
    /**
     * Default cross-linker length.
     */
    public static final float DEFAULT_CROSS_LINKER_LENGTH = 34;
    //--------------------------------------------------------------------------
    /**
     * Default size of cross-linker ends.
     */
    public static final float CROSS_LINKER_END_SIZE = 10;
    //--------------------------------------------------------------------------
    /**
     * Minimum length of peptide in order to be detected by xQuest.
     */
    public static final int MIN_PEPTIDE_LENGTH = 5;
    //--------------------------------------------------------------------------
    /**
     * Maximum length of peptide in order to be detected by xQuest.
     */
    public static final int MAX_PEPTIDE_LENGTH = 40;
    //--------------------------------------------------------------------------
    /**
     * Number of miscleavages allowed in a digested peptide.
     */
    public static final int MAX_NUMBER_OF_MISCLEAVAGES = 1;
    //--------------------------------------------------------------------------
    /**
     * Maximum allowed distance for calculating SASD.
     */
    public static final float MAX_SASD_DISTANCE = 100.0f;
    //--------------------------------------------------------------------------
    /**
     * Solvent radius when SASD calculations are done on the backbone
     * coordinates only.
     */
    public static final float SOLVENT_RADIUS_BACKBONE = 2.0f;
    //--------------------------------------------------------------------------
    /**
     * Default SASD value for cross-links that exceed with their SASD the user
     * given maximum distance.
     */
    public static final float NON_CONFORMING_CROSS_LINK = -1.0f;
    //--------------------------------------------------------------------------
    /**
     * Default SASD value for cross-links with a solvent inaccessible first
     * atom.
     */
    public static final float FIRST_ATOM_IS_SOLVENT_INACCESSIBLE = -2.0f;
    //--------------------------------------------------------------------------
    /**
     * Default SASD value for cross-links with a solvent inaccessible first
     * atom.
     */
    public static final float SECOND_ATOM_IS_SOLVENT_INACCESSIBLE = -3.0f;
    //--------------------------------------------------------------------------
    /**
     * Default SASD value for cross-links with solvent inaccessible atoms.
     */
    public static final float BOTH_ATOMS_ARE_SOLVENT_INACCESSIBLE = -4.0f;
    //--------------------------------------------------------------------------
    /**
     * Default SASD value for cross-linked atom that are buried and unreachable
     * to solvent despite having a positive solvent accessibility.
     */
    public static final float FIRST_ATOM_IS_BURIED = -5.0f;

    //--------------------------------------------------------------------------
    // CONSTANT STRING VALUES
    //--------------------------------------------------------------------------
    /**
     * Regular expression for checking peptide sequences for being
     * cross-linkable. Criteria for cross-linkable are:
     * <ol>
     *     <li>Tryptic peptide, i.e. C-terminus must be either arginine or
     *         lysine.
     *     </li>
     *     <li>One central lysine residue which represents the cross-linked
     *         amino acid.
     *     </li>
     *     <li>Up to one mis-cleavage is allowed. Here the miscleavage is
     *         positioned prior to the cross-linked lysine residue.
     *     </li>
     * </ol>
     */
    public static final String CROSS_LINKABLE_PEPTIDE_SEQUENCE_EXPRESSION1 =
         "[^"
        + AminoAcidType.LYSINE.getOneLetterCode()
        + AminoAcidType.ARGININE.getOneLetterCode()
        + "]*"
        + "["
        + AminoAcidType.LYSINE.getOneLetterCode()
        + AminoAcidType.ARGININE.getOneLetterCode()
        + "]{0,1}"
        + "[^"
        + AminoAcidType.LYSINE.getOneLetterCode()
        + AminoAcidType.ARGININE.getOneLetterCode()
        + "]*"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + "{1}"
       + "[^"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]*"
       + "["
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]$";
    //--------------------------------------------------------------------------
    /**
     * Regular expression for checking peptide sequences for being
     * cross-linkable. Criteria for cross-linkable are:
     * <ol>
     *     <li>Tryptic peptide, i.e. C-terminus must be either arginine or
     *         lysine.
     *     </li>
     *     <li>One central lysine residue which represents the cross-linked
     *         amino acid.
     *     </li>
     *     <li>Up to one mis-cleavage is allowed. Here the miscleavage is
     *         positioned post to the cross-linked lysine residue.
     *     </li>
     * </ol>
     */
    public static final String CROSS_LINKABLE_PEPTIDE_SEQUENCE_EXPRESSION2 =
         "[^"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]*"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + "{1}"
       + "[^"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]*"
       + "["
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]{0,1}"
       + "[^"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]*"
       + "["
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]$";
    //--------------------------------------------------------------------------
    // CONSTANT ENUM SETS
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // OTHER CONSTANT VALUES
    //--------------------------------------------------------------------------
    /**
     * Distances are given up to a single digit after comma.
     */
    public static final NumberFormat DISTANCE_DEC_FORMAT =
                                                       new DecimalFormat("0.0");
    /**
     * Distances are given up to a single digit after comma.
     */
    public static final NumberFormat PROBABILITY_DEC_FORMAT =
                                                    new DecimalFormat("0.0000");

}

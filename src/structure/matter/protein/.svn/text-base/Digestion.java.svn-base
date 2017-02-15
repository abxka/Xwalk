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

package structure.matter.protein;

import java.util.ArrayList;

import structure.matter.MatterUtilities;
import structure.matter.parameter.AminoAcidType;
import xwalk.constants.Constants;

/**
 * Class for digesting PolyPeptide object according to a particular protease.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 *
 */
public class Digestion {

    /**
     * Constructor.
     */
    protected Digestion() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------

    /**
     * PolyPeptide digestion according to a trypsin.
     * Following info was extracted from: <br>
     * <a href="http://www.expasy.ch/tools/peptidecutter/
     *          peptidecutter_special_enzymes.html">here</a>. <br>
     * <pre>
     *                  Cleavage
     *                     ||
     *  Pn-----P4-P3-P2-P1-||-P1'-P2'-P3'-P4'------Pm
     *   |     |  |  |  |  ||   |   |   |        |
     *  Sn-----S4-S3-S2-S1-||-S1'-S2'-S3'-S4'------Sm
     *                     ||
     *  </pre>
     *  Px = Peptide position <br>
     *  Sx = Protease position <br>
     *  <br>
     *  For Trypsin:
     *  <ol>
     *     <li> Cleavage preferable at Arg, Lys at P1,
     *          but neighbouring AA have large impact on digest, in
     *          particular
     *     </li>
     *     <li> Pro at P1' has negative influence
     *     </li>
     *     <br>
     *     Additionally following exception might occur:
     *     <li> Arg, Lys at P1' induces inhibition
     *     </li>
     *     <li> Pro usually blocks the action when found in position P1',
     *          but not when Lys is in position P1 and Trp is in position P2 at
     *          the same time. This blocking of cleavage exerted by Pro in
     *          position P1' is also negligible when Arg is in position P1 and
     *          Met is in position P2 at the same time
     *     </li>
     *     <li> Lys is found in position P1 the following situation considerably
     *          block the action of trypsin:
     *          <ul>
     *              <li>Asp in position P2 and Asp in position P1'
     *              </li>
     *              <li>Cys in position P2 and Asp in position P1'
     *              </li>
     *              <li>Cys in position P2 and His in position P1'
     *              </li>
     *              <li>Cys in position P2 and Tyr in position P1'
     *              </li>
     *         </ul>
     *     <li>Arg is in P1 and the following situations are found:
     *         <ul>
     *              <li>Arg in position P2 and His in position P1'
     *              </li>
     *              <li>Cys in position P2 and Lys in position P1'
     *              </li>
     *              <li>Arg in position P2 and Arg in positionP1'
     *              </li>
     *         </ul>
     *     </li>
     *  </ol>
     *
     * @param protein
     *        - PolyPeptide object to be digested.
     * @param useException
     *        - boolean value indicating to include exceptions to digestion.
     * @return List of PolyPeptide object being the peptides remaining after
     *         digestion.
     */
    public static ArrayList < PolyPeptide > partialTrypticDigest(
                                                     final PolyPeptide protein,
                                                     final boolean useException
                                                  ) {

        ArrayList < PolyPeptide > fullDigest = Digestion.fullTrypticDigest(
                                                                    protein,
                                                                    useException
                                                                          );

        ArrayList < PolyPeptide > digestWithCrossLink =
                                                new ArrayList < PolyPeptide >();
        for (int i = 0; i < fullDigest.size(); i++) {
            PolyPeptide peptide = new PolyPeptide(new ArrayList <AminoAcid>());
            peptide.addAll(fullDigest.get(i));
            String pepSequence = peptide.toStringOneLetterCode();
            // if peptide is longer than maxLength right from the start then
            // just skip it.
            if (pepSequence.length() > Constants.MAX_PEPTIDE_LENGTH) {
                continue;
            }

            // if peptide conforms to a cross-linkable peptide sequence
            // right from the start then added to the peptide list right from
            // the beginnin.
            if (pepSequence.matches(
                    Constants.CROSS_LINKABLE_PEPTIDE_SEQUENCE_EXPRESSION1
                                   )
                ||
                pepSequence.matches(
                    Constants.CROSS_LINKABLE_PEPTIDE_SEQUENCE_EXPRESSION2
                                   )
               ) {
                if (pepSequence.length() >= Constants.MIN_PEPTIDE_LENGTH) {
                    digestWithCrossLink.add(peptide);
                    peptide = new PolyPeptide(peptide);
                }
            }

            // if peptide is not already larger then maxLength, then check
            // whether peptide sequence can be grown continuing to conform
            // a cross-linkable peptide and without being larger then maxLength
            for (int j = i + 1; j < fullDigest.size(); j++) {
                PolyPeptide peptide2 = fullDigest.get(j);

                pepSequence += peptide2.toStringOneLetterCode();
                peptide.addAll(peptide2);

                if (pepSequence.length() > Constants.MAX_PEPTIDE_LENGTH) {
                    break;
                }
                if (pepSequence.length() < Constants.MIN_PEPTIDE_LENGTH) {
                    continue;
                }

                if (pepSequence.matches(
                        Constants.CROSS_LINKABLE_PEPTIDE_SEQUENCE_EXPRESSION1
                                       )
                    ||
                    pepSequence.matches(
                        Constants.CROSS_LINKABLE_PEPTIDE_SEQUENCE_EXPRESSION2
                                       )
                   ) {
                    digestWithCrossLink.add(peptide);
                    peptide = new PolyPeptide(peptide);
                }
            }
        }
        MatterUtilities.sort(digestWithCrossLink);
    return digestWithCrossLink;
    }
    //--------------------------------------------------------------------------
    /**
     * PolyPeptide digestion according to a trypsin.
     * Following info was extracted from: <br>
     * <a href="http://www.expasy.ch/tools/peptidecutter/
     *          peptidecutter_special_enzymes.html">here</a>. <br>
     * <pre>
     *                  Cleavage
     *                     ||
     *  Pn-----P4-P3-P2-P1-||-P1'-P2'-P3'-P4'------Pm
     *   |     |  |  |  |  ||   |   |   |        |
     *  Sn-----S4-S3-S2-S1-||-S1'-S2'-S3'-S4'------Sm
     *                     ||
     *  </pre>
     *  Px = Peptide position <br>
     *  Sx = Protease position <br>
     *  <br>
     *  For Trypsin:
     *  <ol>
     *     <li> Cleavage preferable at Arg, Lys at P1,
     *          but neighbouring AA have large impact on digest, in
     *          particular
     *     </li>
     *     <li> Pro at P1' has negative influence
     *     </li>
     *     <br>
     *     Additionally following exception might occur:
     *     <li> Arg, Lys at P1' induces inhibition
     *     </li>
     *     <li> Pro usually blocks the action when found in position P1',
     *          but not when Lys is in position P1 and Trp is in position P2 at
     *          the same time. This blocking of cleavage exerted by Pro in
     *          position P1' is also negligible when Arg is in position P1 and
     *          Met is in position P2 at the same time
     *     </li>
     *     <li> Lys is found in position P1 the following situation considerably
     *          block the action of trypsin:
     *          <ul>
     *              <li>Asp in position P2 and Asp in position P1'
     *              </li>
     *              <li>Cys in position P2 and Asp in position P1'
     *              </li>
     *              <li>Cys in position P2 and His in position P1'
     *              </li>
     *              <li>Cys in position P2 and Tyr in position P1'
     *              </li>
     *         </ul>
     *     <li>Arg is in P1 and the following situations are found:
     *         <ul>
     *              <li>Arg in position P2 and His in position P1'
     *              </li>
     *              <li>Cys in position P2 and Lys in position P1'
     *              </li>
     *              <li>Arg in position P2 and Arg in positionP1'
     *              </li>
     *         </ul>
     *     </li>
     *  </ol>
     *
     * Performs a full tryptic digestion on protein.
     * @param protein
     *        PolyPeptide object representing a protein to be digested.
     * @param useException
     *        - boolean value indicating to include exceptions to digestion.
     * @return List of digestion product peptides.
     */
    public static final ArrayList < PolyPeptide > fullTrypticDigest(
                                                      final PolyPeptide protein,
                                                      final boolean useException
                                                        ) {
        ArrayList < PolyPeptide > trypticPeptides =
                                                new ArrayList < PolyPeptide >();

        ArrayList < AminoAcid > peptide = new ArrayList < AminoAcid >();

        for (int i = 0; i < protein.size(); i++) {
            peptide.add(protein.get(i));

            AminoAcid p1 = protein.get(i);
            AminoAcidType p1type = p1.getType();
            AminoAcid p1p = null;
            AminoAcidType p1pType = null;
            AminoAcid p2 = null;
            AminoAcidType p2type = null;

            if (i + 1 < protein.size()) {
                p1p = protein.get(i + 1);
                p1pType = p1p.getType();
            }
            if (i - 1 > 0) {
                p2 = protein.get(i - 1);
                p2type = p2.getType();
            }

            if (// check for rule 1.)
                    (p1type == AminoAcidType.LYSINE
                    ||
                    p1type == AminoAcidType.ARGININE)
                    &&
                    p1pType != AminoAcidType.PROLINE) {

                if (useException) {
                    if (!Digestion.trypticDigestByExPASyRule(p1type,
                                                            p1pType,
                                                            p2type)) {
                        continue;
                    }
                }
                trypticPeptides.add(new PolyPeptide(peptide));
                peptide = new ArrayList < AminoAcid >();

            } else if (i == protein.size() - 1) {
                trypticPeptides.add(new PolyPeptide(peptide));
            }
        }
        return trypticPeptides;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns back whether trypsin could digest according to the
     * <a href="http://www.expasy.ch/tools/peptidecutter/
     *          peptidecutter_special_enzymes.html">ExPASy exception rules</a>,
     * given a tripeptide sequence.
     * <pre>
     *                  Cleavage
     *                     ||
     *  Pn-----P4-P3-P2-P1-||-P1'-P2'-P3'-P4'------Pm
     *   |     |  |  |  |  ||   |   |   |        |
     *  Sn-----S4-S3-S2-S1-||-S1'-S2'-S3'-S4'------Sm
     *                     ||
     *  </pre>
     *  Px = Peptide position <br>
     *  Sx = Protease position <br>
     *  <br>
     *  The exception rules are the following:
     *  <ol>
     *     <li> Arg, Lys at P1' induces inhibition
     *     </li>
     *     <li> Pro usually blocks the action when found in position P1',
     *          but not when Lys is in position P1 and Trp is in position P2 at
     *          the same time. This blocking of cleavage exerted by Pro in
     *          position P1' is also negligible when Arg is in position P1 and
     *          Met is in position P2 at the same time
     *     </li>
     *     <li> Lys is found in position P1 the following situation considerably
     *          block the action of trypsin:
     *          <ul>
     *              <li>Asp in position P2 and Asp in position P1'
     *              </li>
     *              <li>Cys in position P2 and Asp in position P1'
     *              </li>
     *              <li>Cys in position P2 and His in position P1'
     *              </li>
     *              <li>Cys in position P2 and Tyr in position P1'
     *              </li>
     *         </ul>
     *     <li>Arg is in P1 and the following situations are found:
     *         <ul>
     *              <li>Arg in position P2 and His in position P1'
     *              </li>
     *              <li>Cys in position P2 and Lys in position P1'
     *              </li>
     *              <li>Arg in position P2 and Arg in positionP1'
     *              </li>
     *         </ul>
     *     </li>
     *  </ol>
     *
     * @param p1type
     *        AminoAcidType of the position p1.
     * @param p1pType
     *        AminoAcidType of the position p1-prime.
     * @param p2type
     *        AminoAcidType of the position p2.
     * @return {@code TRUE} if digestion is possible at position p1,
     *         {@code FALSE} otherwise.
     */
    private static boolean trypticDigestByExPASyRule(final AminoAcidType p1type,
                                              final AminoAcidType p1pType,
                                              final AminoAcidType p2type) {
        // check for rule 2.)
        if (p1pType != null) {
            if (p1type == AminoAcidType.LYSINE
                &&
                p2type == AminoAcidType.TRYPTOPHANE) {
                return false;
            }
            if (p1type == AminoAcidType.ARGININE
                &&
                p2type == AminoAcidType.METHIONINE) {
                return false;
            }
        }

        // check for rule 3.)
        if (p1pType != null) {
            if (p1pType == AminoAcidType.LYSINE
                ||
                p1pType == AminoAcidType.ARGININE) {
                return true;
            }
        }
        // check for rule 5.)
        if (p1type == AminoAcidType.LYSINE) {

            if (p2type == AminoAcidType.ASPARTIC_ACID
                &&
                p1pType == AminoAcidType.ASPARTIC_ACID
            ) {
                return true;
            }
            if (p2type == AminoAcidType.CYSTEINE
                &&
               (p1pType == AminoAcidType.ASPARTIC_ACID
                ||
                p1pType == AminoAcidType.HISTIDINE
                ||
                p1pType == AminoAcidType.TYROSINE)
            ) {
                return true;
            }
        }
        // check for rule 6.)
        if (p1type == AminoAcidType.ARGININE) {

            if (p2type == AminoAcidType.ARGININE
                &&
                p1pType == AminoAcidType.HISTIDINE
            ) {
                return true;
            }
            if (p2type == AminoAcidType.CYSTEINE
                &&
                p1pType == AminoAcidType.LYSINE
            ) {
                return true;
            }
            if (p2type == AminoAcidType.ARGININE
                &&
                p1pType == AminoAcidType.ARGININE
            ) {
                return true;
            }
        }
        return false;
    }
    //--------------------------------------------------------------------------
}

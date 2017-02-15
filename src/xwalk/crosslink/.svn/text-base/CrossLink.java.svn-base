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

package xwalk.crosslink;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.TreeSet;

import structure.constants.Constants;
import structure.constants.Constants.BondTypes;
import structure.grid.Path;
import structure.math.Mathematics;
import structure.matter.Atom;
import structure.matter.Bond;
import structure.matter.parameter.ParameterReader;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import xwalk.crosslink.CrossLinkParameter.Parameter;


/**
 * Class for representing Cross-Link objects.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class CrossLink extends Bond {

    //--------------------------------------------------------------------------
    /**
     * Distance that the cross-link spans in sequence space.
     * Default value is -1.
     */
    private int seqDist = -1;
    /**
     * Distance that the cross-link spans in Euclidean space.
     * Default value is -1.
     */
    private float eucDist = -1.0f;
    /**
     * Distance that the cross-link spans in Solvent-Path distance space.
     * Default value is -0.5.
     */
    private float solventPathDistance = -0.5f;

    /**
     * Probability of finding a cross-link with this Euclidean distance in a
     * cross-linking experiment. The probability is based on observed
     * cross-link distances in the literature and in the Aebersold lab.
     * Default value is 0.
     */
    private float eucDistProbability = -1.0f;
    /**
     * Probability of finding a cross-link with this SAS distance in a
     * cross-linking experiment. The probability is based on observed
     * cross-link distances in the literature and in the Aebersold lab.
     * Default value is 0.
     */
    private float sasdDistProbability = -1.0f;
    /**
     * Stores the information whether a probability calculation has been
     * requested. As a consequence probabilities will be printed out
     * in the toString() method.
     * Default value is false.
     */
    private boolean doProbability = false;
    /**
     * First protein atom that is connected by the cross-link.
     */
    private Atom preAtom;
    /**
     * Second protein atom that is connected by the cross-link.
     */
    private Atom postAtom;
    /**
     * Second protein atom that is connected by the cross-link.
     */
    private Path solventDistancePath;
    /**
     * Ranking index of this cross-link within a list of cross-links.
     * Default value is -1.
     */
    private int index = -1;
    /**
     * String object holding the path to the PDB file in which the cross-link
     * has been found.
     * Default value is an empty string.
     */
    private String filePath = "";
    /**
     * Peptide sequence of first protein atom that is connected by the
     * cross-link.
     */
    private PolyPeptide preAtomPeptide;
    /**
     * Peptide sequence of second protein atom that is connected by the
     * cross-link.
     */
    private PolyPeptide postAtomPeptide;
    /**
     * Hash value.
     */
    //private int hashValue;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param atom1
     *        - First protein atom to be connected by the virtual cross-linker.
     * @param atom2
     *        - Second protein atom to be connected by the virtual cross-linker.
     */
    public CrossLink(final Atom atom1, final Atom atom2) {
        this(atom1,
             atom2,
             Math.abs(atom1.getRank() - atom2.getRank()),
             Mathematics.distance(atom1.getXYZ(), atom2.getXYZ())
            );
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor. No Euclidean distance calculations are performed.
     * @param atom1
     *        - First protein atom to be connected by the virtual cross-linker.
     * @param atom2
     *        - Second protein atom to be connected by the virtual cross-linker.
     * @param euclideanDistance
     *        - float value representing the Euclidean distance between both
     *          atoms.
     * @param sequenceDistance
     *        - integer value representing the distance in sequence space of
     *          both atoms.
     */
    public CrossLink(
                     final Atom atom1,
                     final Atom atom2,
                     final int sequenceDistance,
                     final float euclideanDistance
                    ) {
        super(atom1, atom2, BondTypes.CROSS_LINK);

        ArrayList < Atom > list = new ArrayList < Atom >();
        list.add(atom1);
        list.add(atom2);
        // sorting atom pair by chain id.
        if (CrossLinkParameter.getParameter(
                                        Parameter.DISTANCE_FILE_PATH).equals("")
                                            ) {
            Collections.sort(list, new Comparator < Atom >() {
                                      public int compare(final Atom atom1,
                                                         final Atom atom2) {
                                      String chainId1 = atom1.getChainId() + "";
                                      String chainId2 = atom2.getChainId() + "";
                                      return chainId1.compareTo(chainId2);
                                      }
                                  }
                            );
        }
        this.preAtom = list.get(0);
        this.postAtom = list.get(1);

        this.seqDist = sequenceDistance;
        this.eucDist = euclideanDistance;
        /*
        this.hashValue = this.getPreAtom().getResidueName().hashCode()
                       + this.getPreAtom().getResidueNumber()
                       + this.getPreAtom().getChainId()
                       + this.getPreAtom().getName().hashCode()
                       + this.getPreAtom().getAlternativeLocation()
                       + this.getPostAtom().getResidueName().hashCode()
                       + this.getPostAtom().getResidueNumber()
                       + this.getPostAtom().getChainId()
                       + this.getPostAtom().getName().hashCode()
                       + this.getPostAtom().getAlternativeLocation();
        */
    }
    //--------------------------------------------------------------------------

    /**
     * Sets the Solvent-Path distance. In order to calculate the Solvent-Path
     * distance please see in the Class CrossLinkList.
     * @param dist
     *        - Distance in the Solvent-Path space.
     */
    public final void setSolventPathDistance(final float dist) {
        this.solventPathDistance = dist;
    }
    //--------------------------------------------------------------------------

    /**
     * Sets the grid path of the Solvent-Path distance.
     * @param path -
     *        Path object holding the list of GridCell object between source
     *        and target grid cells within a Grid object.
     */
    public final void setPath(final Path path) {
        this.solventDistancePath = path;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the distance in sequence space.
     * @return integer number representing the distance in sequence space.
     */
    public final int getSequenceDistance() {
        return this.seqDist;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the distance in Euclidean space.
     * @return float number representing the distance in Euclidean space.
     */
    public final float getEuclideanDistance() {
        return Float.parseFloat(
            xwalk.constants.Constants.DISTANCE_DEC_FORMAT.format(this.eucDist));
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the distance in Solvent-Path distance space.
     * @return float number representing the distance in Solvent-Path distance
     *         space.
     */
    public final float getSolventPathDistance() {
        return Float.parseFloat(
                xwalk.constants.Constants.DISTANCE_DEC_FORMAT.format(
                                                        this.solventPathDistance
                                                                    ));
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the probability of this cross-link to be observed with its
     * SAS distance in real experiments.
     * @return float number representing the probability.
     */
    public final float getSolventPathDistanceProbability() {
        return Float.parseFloat(
                xwalk.constants.Constants.PROBABILITY_DEC_FORMAT.format(
                                                        this.sasdDistProbability
                                                                       ));
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the probability of this cross-link to be observed with its
     * Euclidean distance in real experiments.
     * @return float number representing the probability.
     */
    public final float getEuclideanDistanceProbability() {
        return Float.parseFloat(
                xwalk.constants.Constants.PROBABILITY_DEC_FORMAT.format(
                                                         this.eucDistProbability
                                                                       ));
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the grid path of the Solvent-Path distance.
     * @return Path object holding the list of GridCell object between source
     *         and target grid cells within a Grid object.
     */
    public final Path getPath() {
        return this.solventDistancePath;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether a second cross-link has the same residue name and residue
     * number.
     * @param crossLink
     *        - CrossLink object to be compared to this CrossLink object.
     * @return {@code TRUE} if both CrossLink object are equal in homology,
     * {@code FALSE} otherwise.
     */
    public final boolean equalsInHomolog(final CrossLink crossLink) {
        Atom preAtom1 = this.getPreAtom();
        Atom postAtom1 = this.getPostAtom();
        Atom preAtom2 = crossLink.getPreAtom();
        Atom postAtom2 = crossLink.getPostAtom();

        String residueId1 = preAtom1.getResidueName() + "#"
                            + preAtom1.getResidueNumber();
        String residueId2 = postAtom1.getResidueName() + "#"
                            + postAtom1.getResidueNumber();
        String residueId3 = preAtom2.getResidueName() + "#"
                            + preAtom2.getResidueNumber();
        String residueId4 = postAtom2.getResidueName() + "#"
                            + postAtom2.getResidueNumber();
        if ((residueId1.equals(residueId3) && residueId2.equals(residueId4))
            ||
            (residueId2.equals(residueId3) && residueId1.equals(residueId4))) {
            return true;
        }
        return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether a second cross-link has the same atom identifier, i.e.
     * residue name, residue number, chain Id and atom name.
     * @param crossLink
     *        - CrossLink object to be compared to this CrossLink object.
     * @return {@code TRUE} if both CrossLink object are equal in all atom
     *         identifier, {@code FALSE} otherwise.
     */
    public final boolean equals(final CrossLink crossLink) {
        Atom preAtom1 = this.getPreAtom();
        Atom postAtom1 = this.getPostAtom();
        Atom preAtom2 = crossLink.getPreAtom();
        Atom postAtom2 = crossLink.getPostAtom();

        if (this.equalsInHomolog(crossLink)) {
            String residueId1 = preAtom1.getChainId() + "#"
                                + preAtom1.getName().trim();
            String residueId2 = postAtom1.getChainId() + "#"
                                + postAtom1.getName().trim();
            String residueId3 = preAtom2.getChainId() + "#"
                                + preAtom2.getName().trim();
            String residueId4 = postAtom2.getChainId() + "#"
                                + postAtom2.getName().trim();

            if ((residueId1.equals(residueId3) && residueId2.equals(residueId4))
               ||
               (residueId2.equals(residueId3) && residueId1.equals(residueId4)))
            {
                return true;
            }
        }
        return false;
    }
    //-------------------------------------------------------------------------
    /*
     * Returns the hash value of this cross-link, which incorporates information
     * from its residue names, residue numbers, chainID and atom names.
     * @return integer variable representing the hash value of this CrossLink
     *         object.
     */
//    public final int hashCode() {
//        return this.hashValue;
//    }
    //--------------------------------------------------------------------------
    /**
     * Sets the ranking index of this cross-link.
     * @param rank
     *        - integer value representing the ranking index.
     */
    public final void setIndex(final int rank) {
        this.index = rank;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the ranking index of this cross-link.
     * @return Integer value representing the ranking index.
     */
    public final int getIndex() {
        return index;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a String representation of this cross-link in distance file
     * format.
     * @return String object holding the representation of this cross-link in
     * distance file format.
     */
    public final String toString() {
        String atomId1 = AminoAcid.getAminoAcidId(preAtom)
                         + "-"
                         + preAtom.getName().trim();
        String atomId2 = AminoAcid.getAminoAcidId(postAtom)
                         + "-"
                         + postAtom.getName().trim();

        if (preAtom.getAlternativeLocation() != ' ') {
            atomId1 += "-" + preAtom.getAlternativeLocation();
        }

        if (postAtom.getAlternativeLocation() != ' ') {
            atomId2 += "-" + postAtom.getAlternativeLocation();
        }

        StringBuffer output = new StringBuffer();
        int maxPeptideLength = xwalk.constants.Constants.MAX_PEPTIDE_LENGTH;
        int minPeptideLength = xwalk.constants.Constants.MIN_PEPTIDE_LENGTH;
        boolean outputPeptide = false;
        if (preAtomPeptide != null && postAtomPeptide != null) {
            outputPeptide = preAtomPeptide.size() <= maxPeptideLength
                            &&
                            preAtomPeptide.size() >= minPeptideLength
                            &&
                            postAtomPeptide.size() <= maxPeptideLength
                            &&
                            postAtomPeptide.size() >= minPeptideLength;
        }

        output.append(this.filePath + "\t" + atomId1 + "\t" + atomId2 + "\t"
                    + this.seqDist + "\t"
                    + this.getEuclideanDistance() + "\t");
        if (this.getSolventPathDistance() > -0.1
                ||
            this.getSolventPathDistance() < -0.9 ) {
            output.append(this.getSolventPathDistance() + "\t");
        } else {
            output.append("-\t");
        }
        if (this.doProbability) {
            output.append(this.getEuclideanDistanceProbability() + "\t");
            if (this.getSolventPathDistance() > -0.1
                    ||
                this.getSolventPathDistance() < -0.9 ) {
                output.append(this.getSolventPathDistanceProbability() + "\t");
            } else {
                output.append("-\t");
            }
        } else {
            output.append("-\t-\t");
        }

        if (outputPeptide) {
            output.append(this.preAtomPeptide.toStringOneLetterCode() + "-"
                       +  this.postAtomPeptide.toStringOneLetterCode());
        } else {
            output.append("-");
        }
        output.append(Constants.LINE_SEPERATOR);
    return output.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the path to the file in which this cross-link has been found.
     * @param fileName
     *        - String object holding the path to the file.
     */
    public final void setFileName(final String fileName) {
        this.filePath = fileName;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the path to the file in which this cross-link has been found.
     * @return String object holding the path to the file.
     */
    public final String getFileName() {
        return filePath;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets both peptide object that are interconnected by this cross-linker.
     * @param prePeptide
     *        - PolyPeptide object of the preAtom.
     * @param postPeptide
     *        - PolyPeptide object of the postAtom.
     * @return {@code TRUE} if both peptides contain both cross-linked atoms,
     *         {@code FALSE} otherwise.
     */
    public final boolean setPeptides(final PolyPeptide prePeptide,
                                     final PolyPeptide postPeptide) {
        for (AminoAcid aa : prePeptide) {
            if (aa.getAllAtoms().contains(this.preAtom)) {
                this.preAtomPeptide = prePeptide;
            }
            if (aa.getAllAtoms().contains(this.postAtom)) {
                this.postAtomPeptide = prePeptide;
            }
        }
        // here the order must be inverse, i.e. first check for post than for
        // pre atom to allow for "self-cross-links". Self-cross-links are
        // cross-links between two identical atom, which is physically
        // incorrect but might be requested over a distance file.
        for (AminoAcid aa : postPeptide) {
            if (aa.getAllAtoms().contains(this.postAtom)) {
                this.postAtomPeptide = postPeptide;
            }
            if (aa.getAllAtoms().contains(this.preAtom)) {
                this.preAtomPeptide = postPeptide;
            }
        }
        if (this.preAtomPeptide != null && this.postAtomPeptide != null) {
            return true;
        }
        return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the peptide object to which preAtom belongs.
     * @return PolyPeptide object of preAtom.
     */
    public final PolyPeptide getPrePeptide() {
        return this.preAtomPeptide;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the peptide object to which postAtom belongs.
     * @return PolyPeptide object of postAtom.
     */
    public final PolyPeptide getPostPeptide() {
        return this.postAtomPeptide;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the probability of observing this cross-link with its Euclidean
     * distance.
     */
    public final void setEucProbability() {
        float preProb = -1;
        TreeSet<Float> sortedBins = new TreeSet<Float>(
                   ParameterReader.getEuclideanDistanceProbabilitySet().keySet()
                                                        );
        for (float bin : sortedBins) {
            float prob =
                  ParameterReader.getEuclideanDistanceProbabilitySet().get(bin);
            if (bin > this.eucDist) {
                this.eucDistProbability = preProb;
                break;
            }
            preProb = prob;
        }
        if (this.eucDistProbability == -1 && this.eucDist != -1) {
            this.eucDistProbability = 0;
        }
        this.doProbability = true;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the probability of observing this cross-link with its SAS distance.
     */
    public final void setSASDprobability() {
        float preProb = -1;
        TreeSet<Float> sortedBins = new TreeSet<Float>(
                         ParameterReader.getSASdistanceProbabilitySet().keySet()
                                                        );
        for (float bin : sortedBins) {
            float prob =
                        ParameterReader.getSASdistanceProbabilitySet().get(bin);
            if (bin > this.solventPathDistance) {
                this.sasdDistProbability = preProb;
                break;
            }
            preProb = prob;
        }
        if (this.sasdDistProbability == -1 && this.solventPathDistance != -1) {
            this.sasdDistProbability = 0;
        }
        this.doProbability = true;
    }
}

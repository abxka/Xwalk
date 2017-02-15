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
 */
package sequence.alignment;

import java.util.Hashtable;

import structure.exceptions.FileFormatException;
import structure.io.ReadFile;


/**
 * Class for handling protein sequence alignments as produced by the
 * NEEDLE or WATER algorithms in the
 * <a href="http://emboss.sourceforge.net">EMBOSS</a> package.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.5
 */
public class Alignment {
    //--------------------------------------------------------------------------
    // CLASS OBJECTS
    //--------------------------------------------------------------------------
    /**
     * Sequence with gaps of 1st protein in the alignment.
     */
    private String sequence1st = "";
    //--------------------------------------------------------------------------
    /**
     * Sequence with gaps of 2nd protein in the alignment.
     */
    private String sequence2nd = "";
    //--------------------------------------------------------------------------
    /**
     * Name of 1st protein as found in the alignment file.
     */
    private String sequence1stName = "";
    //--------------------------------------------------------------------------
    /**
     * Name of 2nd protein as found in the alignment file.
     */
    private String sequence2ndName = "";
    //--------------------------------------------------------------------------
    /**
     * Starting position of the 1st protein in the alignment.
     */
    private int startPosition1st = Integer.MIN_VALUE;
    //--------------------------------------------------------------------------
    /**
     * Starting position of the 2nd protein in the alignment.
     */
    private int startPosition2nd = Integer.MIN_VALUE;
    //--------------------------------------------------------------------------
    /**
     * Alignment characters found between both protein sequences in the
     * alignment file.
     */
    private String alignmentSequence = "";
    //--------------------------------------------------------------------------
    /**
     * Gap open penalty score used for the alignment.
     */
    private double gapOpenPenalty;
    //--------------------------------------------------------------------------
    /**
     * Gap extension penalty score used for the alignment.
     */
    private double gapExtensionPenalty;
    //--------------------------------------------------------------------------
    /**
     * Alignment length.
     */
    private int length;
    //--------------------------------------------------------------------------
    /**
     * Number of identical amino acids found in the alignment between both
     * proteins.
     */
    private int identity;
    //--------------------------------------------------------------------------
    /**
     * Number of similar amino acids found in the alignment between both
     * proteins.
     */
    private int similarity;
    /**
     * Number of gaps found in the alignment between both proteins.
     */
    private int gaps;
    //--------------------------------------------------------------------------
    /**
     * Alignment score between both proteins.
     */
    private double score;
    //--------------------------------------------------------------------------
    /**
     * Similarity matrix used for the alignment.
     */
    private String matrix;
    //--------------------------------------------------------------------------
    /**
     * Stores the residue numbers of the 2nd protein that correspond to the
     * residue number in the 1st protein.
     */
    private Hashtable<Integer, Integer> secondVsFirst;
    //--------------------------------------------------------------------------
    /**
     * Stores the residue numbers of the 1st protein that correspond to the
     * residue number in the 2nd protein.
     */
    private Hashtable<Integer, Integer> firstVsSecond;
    //--------------------------------------------------------------------------
    // CONSTRUCTORS
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param alignmentFile
     *        ReadFile object holding the content of a protein sequence
     *        alignment file.
     * @throws FileFormatException if an error occurs during parsing the
     *         alignment file.
     */
    public Alignment(final ReadFile alignmentFile) throws FileFormatException {
        this.read(alignmentFile);
    }
    //--------------------------------------------------------------------------
    // CLASS METHODS
    //--------------------------------------------------------------------------
    /**
     * Reads in an alignment file and stores all information within this class.
     * @param alignmentFile
     *        ReadFile object holding the content of a protein sequence
     *        alignment file.
     * @throws FileFormatException if an error occurs during parsing.
     */
    private void read(final ReadFile alignmentFile) throws FileFormatException {
        for (String line : alignmentFile) {
            try {
                String[] lineArray = line.split(" +");
                if (line.startsWith("# 1:")) {
                    this.sequence1stName = lineArray[2];
                } else if (line.startsWith("# 2:")) {
                    this.sequence1stName = lineArray[2];
                } else if (line.startsWith("# Matrix:")) {
                    this.matrix += lineArray[2];
                } else if (line.startsWith("# Gap_penalty:")) {
                    this.gapOpenPenalty = Double.parseDouble(lineArray[2]);
                } else if (line.startsWith("# Extend_penalty:")) {
                    this.gapExtensionPenalty = Double.parseDouble(lineArray[2]);
                } else if (line.startsWith("# Length:")) {
                    this.length = Integer.parseInt(lineArray[2]);
                } else if (line.startsWith("# Identity:")) {
                    String[] fraction = lineArray[2].split("/");
                    this.identity = Integer.parseInt(fraction[0]);
                } else if (line.startsWith("# Similarity:")) {
                    String[] fraction = lineArray[2].split("/");
                    this.similarity = Integer.parseInt(fraction[0]);
                } else if (line.startsWith("# Gaps:")) {
                    String[] fraction = lineArray[2].split("/");
                    this.gaps = Integer.parseInt(fraction[0]);
                } else if (line.startsWith("# Score:")) {
                    this.score = Double.parseDouble(lineArray[2]);
                } else if (line.startsWith(this.sequence1stName)) {
                    this.sequence1st += lineArray[2];
                    if (this.startPosition1st == Integer.MIN_VALUE) {
                        this.startPosition1st = Integer.parseInt(lineArray[1]);
                    }
                } else if (line.startsWith(this.sequence2ndName)) {
                    this.sequence2nd += lineArray[2];
                    if (this.startPosition2nd == Integer.MIN_VALUE) {
                        this.startPosition2nd = Integer.parseInt(lineArray[1]);
                    }
                } else if (line.contains("[|:.]") && !sequence1st.equals("")) {
                    this.alignmentSequence += line.substring(21, 71);
                }
            } catch (Exception e) {
                throw new FileFormatException("ERROR while reading alignment "
                        + "file parameters: " + e);
            }
        }
        if (this.sequence1st.length() != this.sequence2nd.length()) {
            throw new FileFormatException("ERROR while reading alignment file. "
                                        + "Aligned sequences have different "
                                        + "lenghts: "
                                        + this.sequence1st.length()
                                        + " vs "
                                        + this.sequence2nd.length());
        }
        this.setAminoAcidNumbers();
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the mutual corresponding amino acid numbers. {@code NULL} will be
     * assigned to amino acids numbers that fall within gaps.
     */
    private void setAminoAcidNumbers() {
        int n1 = -1;
        int n2 = -1;
        for (int i = 0; i < this.sequence1st.length(); i++) {
            char aa1 = this.sequence1st.charAt(i);
            char aa2 = this.sequence2nd.charAt(i);
            if (aa1 == '-') {
                n1++;
            }
            if (aa2 == '-') {
                n2++;
            }
            if (aa1 != '-' && aa2 != '-') {
                firstVsSecond.put(this.startPosition1st + n1,
                                  this.startPosition2nd + n2);
                secondVsFirst.put(this.startPosition2nd + n2,
                                  this.startPosition1st + n1);
            } else if (aa1 == '-') {
                secondVsFirst.put(this.startPosition2nd + n2,
                                  null);
            } else if (aa2 == '-') {
                firstVsSecond.put(this.startPosition1st + n1,
                                  null);
            }
        }
    }
    //--------------------------------------------------------------------------
    /** Returns the sequence of the 1st protein with gaps.
     * @return String object holding the sequence of the 1st protein with gaps.
     */
    public final String getSequence1st() {
        return sequence1st;
    }
    //--------------------------------------------------------------------------
    /** Returns the sequence of the 2nd protein with gaps.
     * @return String object holding the sequence of the 2nd protein with gaps.
     */
    public final String getSequence2nd() {
        return sequence2nd;
    }
    //--------------------------------------------------------------------------
    /** Returns the name of the 1st protein as it appears in the alignment file.
     * @return String object holding the name of the 1st protein.
     */
    public final String getSequence1stName() {
        return sequence1stName;
    }
    //--------------------------------------------------------------------------
    /** Returns the name of the 2nd protein as it appears in the alignment file.
     * @return String object holding the name of the 2nd protein.
     */
    public final String getSequence2ndName() {
        return sequence2ndName;
    }
    //--------------------------------------------------------------------------
    /** Returns the alignment characters found between both protein sequences in
     * the alignment file.
     * @return String object holding the alignment characters.
     */
    public final String getAlignmentSequence() {
        return alignmentSequence;
    }
    //--------------------------------------------------------------------------
    /** Returns the gap open penalty score.
     * @return double value representing the gap open penalty score.
     */
    public final double getGapOpenPenalty() {
        return gapOpenPenalty;
    }
    //--------------------------------------------------------------------------
    /** Returns the gap extension penalty score.
     * @return double value representing the gap extension penalty score.
     */
    public final double getGapExtensionPenalty() {
        return gapExtensionPenalty;
    }
    //--------------------------------------------------------------------------
    /** Returns the length of the alignment.
     * @return double value representing the length of the alignment.
     */
    public final int getLength() {
        return length;
    }
    //--------------------------------------------------------------------------
    /** Returns the number of identical amino acids within the alignment.
     * @return integer number representing the number of identical amino acids.
     */
    public final int getIdentity() {
        return identity;
    }
    //--------------------------------------------------------------------------
    /** Returns the number of similar amino acids within the alignment.
     * @return integer number representing the number of similar amino acids.
     */
    public final int getSimilarity() {
        return similarity;
    }
    //--------------------------------------------------------------------------
    /** Returns the number of gaps within the alignment.
     * @return integer number representing the number of gaps.
     */
    public final int getGaps() {
        return gaps;
    }
    //--------------------------------------------------------------------------
    /** Returns the alignment score.
     * @return double value representing the alignment score.
     */
    public final double getScore() {
        return score;
    }
    //--------------------------------------------------------------------------
    /** Returns the name of the similarity matrix used for the alignment.
     * @return String object holding the name of the similarity matrix.
     */
    public final String getMatrix() {
        return matrix;
    }
    //--------------------------------------------------------------------------
    /** Returns the amino acid number of the 2nd protein that corresponds to
     * an amino acid number from the 1st protein.
     * @param residueNumberOf1stProtein
     *        integer number of the amino acid in the 1st protein.
     * @return integer number representing the correspond amino acid from the
     *         2nd protein. Returns {@code -1}, if residueNumberOf1stProtein
     *         corresponds to a gap or is outside the alignment region.
     */
    public final int getResidueNumberOf2ndProtein(
                                          final int residueNumberOf1stProtein) {
        if (firstVsSecond.get(residueNumberOf1stProtein) == null) {
            return -1;
        }
        return firstVsSecond.get(residueNumberOf1stProtein);
    }
    //--------------------------------------------------------------------------
    /** Returns the amino acid number of the 1st protein that corresponds to
     * an amino acid number from the 2nd protein.
     * @param residueNumberOf2ndProtein
     *        integer number of the amino acid in the 2nd protein.
     * @return integer number representing the correspond amino acid from the
     *         1st protein. Returns {@code -1}, if residueNumberOf1stProtein
     *         corresponds to a gap or is outside the alignment region.
     */
    public final int getResidueNumberOf1stProtein(
                                          final int residueNumberOf2ndProtein) {
        if (secondVsFirst.get(residueNumberOf2ndProtein) == null) {
            return -1;
        }
        return secondVsFirst.get(residueNumberOf2ndProtein);
    }
}

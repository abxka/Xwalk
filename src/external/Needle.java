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

package external;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;

import sequence.alignment.Alignment;
import structure.constants.Constants;
import structure.exceptions.FileFormatException;
import structure.io.ReadFile;

/**
 * Class making use of functionality provided by Simon Hubbard's NACCESS
 * application.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Needle {

    /**
     * Path to the NEEDLE application.
     */
    private File needle;
    //--------------------------------------------------------------------------
    /**
     * Gap open penalty score. Default = 10.
     */
    private double gapOpenPenalty = 10;
    //--------------------------------------------------------------------------
    /**
     * Gap extension penalty score. Default = 0.5.
     */
    private double gapExtensionPenalty = 0.5;
    //--------------------------------------------------------------------------
    /**
     * Alignment file.
     */
    private Alignment alignment;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param needleApp
     *        String holding the path to the NEEDLE application.
     * @throws FileNotFoundException
     *         if NEEDLE application can not be found at given path or is not
     *         executable.
     */
    public Needle(final File needleApp) throws FileNotFoundException {
        if (!needleApp.exists()) {
            throw new FileNotFoundException("ERROR: NEEDLE application does "
                                          + "not reside at "
                                          + needleApp.getAbsolutePath());
        } else if (!needleApp.canExecute()) {
            throw new FileNotFoundException("ERROR: NEEDLE application is "
                                          + "not executable!");
        } else {
            this.needle = needleApp;
        }
    }
    //--------------------------------------------------------------------------
    /** Returns the gap open penalty score.
     * @return double value representing the gap open penalty score.
     * @see #setGapOpenPenalty(double)
     */
    public final double getGapOpenPenalty() {
        return gapOpenPenalty;
    }
    //--------------------------------------------------------------------------
    /** Sets the gap open penalty score. Default = 10.0.
     * @param gapOpenPenaltyScore
     *        double value representing the new gap open penalty score.
     * @see #getGapOpenPenalty()
     */
    public final void setGapOpenPenalty(final double gapOpenPenaltyScore) {
        this.gapOpenPenalty = gapOpenPenaltyScore;
    }
    //--------------------------------------------------------------------------
    /** Returns the gap extension penalty score.
     * @return double value representing the gap extension penalty score.
     * @see #setGapExtensionPenalty(double)
     */
    public final double getGapExtensionPenalty() {
        return gapExtensionPenalty;
    }
    //--------------------------------------------------------------------------
    /** Sets the gap extension penalty score. Default = 0.5.
     * @param gapExtensionPenaltyScore
     *        double value representing the new gap extension penalty score.
     * @see #getGapExtensionPenalty()
     */
    public final void setGapExtensionPenalty(
                                        final double gapExtensionPenaltyScore) {
        this.gapExtensionPenalty = gapExtensionPenaltyScore;
    }
    //--------------------------------------------------------------------------
    /**
     * Runs NACCESS and stores the residue and atom accessibility files in
     * a ReadFile object.
     * @param fasta1
     *        File object holding a text file of a protein sequence in FASTA
     *        format.
     * @param fasta2
     *        File object holding a text file of a protein sequence in FASTA
     *        format.
     * @return {@code TRUE} if execution of NACCESS was successful,
     *         {@code FALSE} otherwise.
     * @throws FileFormatException if an error occurs during the execution or
     *         parsing of the alignment file.
     */
    public final boolean run(final File fasta1, final File fasta2)
                                                    throws FileFormatException {
        String command = needle + " " + fasta1.getAbsolutePath()
                                + " " + fasta2.getAbsolutePath()
                                + " -gapopen " + this.gapOpenPenalty
                                + " -gapextend " + this.gapExtensionPenalty
                                + " stdout";
        String[] alignFileContent = ExternCommand.execute(command, false).split(
                                                      Constants.LINE_SEPERATOR);
        ArrayList<String> alignFile = new ArrayList<String>();
        for (String line : alignFileContent) {
            alignFile.add(line);
        }
        this.alignment = new Alignment(new ReadFile(alignFile));
    return true;
    }
    //--------------------------------------------------------------------------
    /** Returns the previously calculated alignment between two FASTA files.
     * @return Alignment object holding the alignment.
     * @see #run(File, File)
     */
    public final Alignment getAlignment() {
        return alignment;
    }
}

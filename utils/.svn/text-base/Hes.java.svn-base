
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

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;

import mm.hydrophobicity.Hydrophobicity;

import structure.constants.Constants;
import structure.io.Commandline;
import structure.io.pdb.PDBreader;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.hetgroups.SmallMolecule;
import structure.matter.protein.PolyPeptideList;

/**
 * Class for calculating the Hydrophobicity Environment Scores (HES) for a
 * protein structure in PDB format.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Hes {
    //--------------------------------------------------------------------------
    /**
     * Path to the PDB formatted file of the protein complex.
     */
    private File pdbFile;
    /**
     * Path to a PDB file, holding sample points coordinates.
     */
    private File sampleFile;
    //--------------------------------------------------------------------------
    /**
     * Reads all parameter from the commandline.
     * @param args
     *        String array holding all commandline arguments.
     */
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param args
     *        String array holding all commandline arguments.
     */
    public Hes(final String[] args) {
        this.readCommandline(args);
    }

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param args
     *        String array holding all commandline arguments.
     */
    private void readCommandline(final String[] args) {
        String nL = Constants.LINE_SEPERATOR;

        //-----------user information-------------------------------------------
        if (args.length == 0) {
            System.out.println(nL
                             + "java " + this.getClass().getName() + " -help"
                             + nL);
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-help", false).equals("EXISTS")) {
            System.err.println(nL
                           + "Information:"
                           + nL
                           + "\tThis program calculates Hydrophobicity "
                           + "Environment Scores (HES) for a protein "
                           + "structure, where the underlying atomic XlogP "
                           + "values are listed as occupancy values and the "
                           + "HES values are given in the temperature factor "
                           + "column. In addition basic statistics will be "
                           + "printed out as REMARK entries."
                           + nL
                           + nL
                           + "Usage: "
                           + nL
                           + "\tjava " + this.getClass().getName()
                           + " -in 1b14.pdb"
                           + nL
                           + nL
                           + "Parameters:"
                           + nL
                           + "\t-in <path>\tany structure file in PDB format. "
                           + "Only ATOM entries will be read in (required)."
                           + nL
                           + "\t-sample <path>\tstructure file in PDB format "
                           + "on which HES values will be mapped. "
                           + "If not set, HES values will be calculated on "
                           + "protein atoms. ATOM and HETATM entries will be "
                           + "read in (optional)."
                           + nL
                           + nL
                           );
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-in", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-in\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            String file = Commandline.get(args, "-in", true);
            this.pdbFile = new File(file);

            if (!this.pdbFile.exists()) {
                System.err.print(nL
                              + "Couldn't open PDB infile \""
                              + this.pdbFile.getAbsolutePath()
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (!Commandline.get(args, "-sample", true).equals("ERROR")) {
            String file = Commandline.get(args, "-sample", true);
            this.sampleFile = new File(file);

            if (!this.sampleFile.exists()) {
                System.err.print(nL
                              + "Couldn't open sample PDB file \""
                              + this.sampleFile.getAbsolutePath()
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
    }
    /**
     * Returns REMARK lines with some statistics on HES.
     * @param atoms
     *        AtomList object with HES values set.
     * @return String object holding PDB-like REMARK lines with statistics on
     *         the non-interface part of the protein surface.
     */
    private String getRemarks(final AtomList atoms) {
        String nL = Constants.LINE_SEPERATOR;
        NumberFormat dec = xwalk.constants.Constants.DISTANCE_DEC_FORMAT;
        StringBuffer output = new StringBuffer();

        float hesSum = 0;
        for (Atom atom : atoms) {
            hesSum += atom.getHes();
        }

        output.append("HEADER   " + this.pdbFile.getName()
                    + nL);
        output.append("REMARK   ATOM COUNT: "
                    + atoms.size() + nL);
        output.append("REMARK   HES SUM: "
                    + dec.format(hesSum) + nL);
        output.append("REMARK   AVERAGE: "
                    + dec.format(hesSum / atoms.size())
                    + nL);
        output.append("REMARK" + nL);
        output.append("REMARK   Occupancy Column: atomic XlogP value" + nL);
        output.append("REMARK   B-factor Column: HES value" + nL);
        output.append("REMARK" + nL);
    return output.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Main method. Reads in a PDB file, whose path is given as a first argument
     * on the commandline and calculates binding interfaces to all protein
     * chains given in the PDB file.
     * @param args
     *        Array of Strings holding all commandline arguments.
     */
    public static void main(final String[] args) {
        Hes hes = new Hes(args);

        //----------------------------------------------------------------------
        // Read in PDB files
        //----------------------------------------------------------------------
        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(hes.pdbFile.getAbsolutePath());

        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

        PolyPeptideList proteinComplex =
                                readers.get(0).getEntireProteinComplex().get(0);

        AtomList sampleCoords;
        if (hes.sampleFile == null) {
            sampleCoords = proteinComplex.getAllAtoms();
        } else {
            try {
                readers = PDBreader.createPDBreaders(
                                              hes.sampleFile.getAbsolutePath());
            } catch (Exception e) {
                System.err.println(e.getMessage());
            }
            sampleCoords = readers.get(0).getAllAtoms();
            if (readers.get(0).getAllSmallMolecules().size() > 0) {
                for (SmallMolecule molecule
                               : readers.get(0).getAllSmallMolecules().get(0)) {
                    sampleCoords.addAll(molecule.getAllAtoms());
                }
            }
        }

        Hydrophobicity hydrophobicity = new Hydrophobicity(proteinComplex);
        hydrophobicity.mapHydrophobicity(sampleCoords);
        for (int i = 0; i < sampleCoords.size(); i++) {
            sampleCoords.get(i).setTemperatureFactor(
                                                  sampleCoords.get(i).getHes());
            if (hes.sampleFile == null) {
                sampleCoords.get(i).setOccupancy(
                                proteinComplex.getAllAtoms().get(i).getXlogP());
            } else {
                sampleCoords.get(i).setOccupancy(0);
            }
        }
        //------------------------------------------------------------------
        // Output HES statistics and coordinates with HES in temperature factor
        // columns and XlogP values in the occupancy column.
        //------------------------------------------------------------------
        String output = hes.getRemarks(sampleCoords);
        System.out.print(output + sampleCoords.toString());
    }
}

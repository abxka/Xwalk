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
import java.util.ArrayList;

import external.Naccess;


import structure.io.Commandline;
import structure.io.pdb.PDBreader;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptideList;
import structure.matter.protein.Surface;

/**
 * Class holding a main method to calculate the number of atoms with a certain
 * radius around each amino acid.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Compactness {
    /**
     * Empty Constructor.
     */
    protected Compactness() {
    }
    //--------------------------------------------------------------------------
    /**
     * Path to the PDB formatted file of the protein complex.
     */
    private File pdbFile;
    /**
     * Path to the naccess program.
     */
    private File naccessPath;
    /**
     * Radius to determine environment. Default 5.0 Angstroem.
     */
    private double radius = 5;
    //--------------------------------------------------------------------------
    /**
     * Reads all parameter from the commandline.
     * @param args
     *        String array holding all commandline arguments.
     */
    private void readCommandline(final String[] args) {
        String nL = structure.constants.Constants.LINE_SEPERATOR;

        //-----------user information-------------------------------------------
        if (args.length == 0) {
            System.out.println(nL
                             + "java " + Compactness.class.getName() + " -help"
                             + nL);
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-help", false).equals("EXISTS")) {
            System.err.println(nL
                           + "USAGE:" + nL
                           + "\tjava " + Compactness.class.getName()
                           + " -in 1b14.pdb" + nL
                           + nL
                           + "INFORMATION:" + nL
                           + "\tThis program calculates for every residue "
                           + "in a protein complex the number of atoms it is "
                           + "surrounded with."
                           + nL
                           + "\tNOTE: During the calculations temporary files "
                           + "starting with proteinA, proteinB and proteinAB "
                           + "and with suffixes .asa, .rsa and .log will be "
                           + "generated, which can be deleted manually "
                           + "afterwards."
                           + nL
                           + nL
                           + "PARAMETERS:"
                           + nL
                           + "\t-in <path>\tany structure file in PDB format "
                           + "(required)."
                           + nL
                           + "\t-naccess <string>\tPath to the naccess "
                           + "executable, which will be used to focus "
                           + "calculations only on solvent accessible amino "
                           + "acids. (optional)."
                           + nL
                           + "\t-radius <double>\tRadius of environment. "
                           + "Default are 5.0 Angstroem (optional)."
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

            this.pdbFile = new File(Commandline.get(args, "-in", true));

            if (!this.pdbFile.exists()) {
                System.err.print(nL
                              + "Couldn't open infile \""
                              + this.pdbFile.getAbsolutePath()
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------
        if (!Commandline.get(args, "-naccess", true).equals("ERROR")) {
            this.naccessPath =
                              new File(Commandline.get(args, "-naccess", true));
            if (!this.naccessPath.exists()) {
                System.err.print("ERROR: " + this.naccessPath + " is not "
                               + "executable" + nL);
                System.exit(1);
            }
        }
        //----------------------------
        if (!Commandline.get(args, "-radius", true).equals("ERROR")) {
            this.radius = Double.parseDouble(
                                        Commandline.get(args, "-radius", true));
        }
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
        Compactness compactness = new Compactness();
        compactness.readCommandline(args);

        //----------------------------------------------------------------------
        // Read in PDB files
        //----------------------------------------------------------------------
        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(
                                         compactness.pdbFile.getAbsolutePath());

        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

        PolyPeptideList proteinComplex =
                                readers.get(0).getEntireProteinComplex().get(0);

        ArrayList<AminoAcid> sample = proteinComplex.getAllAminoAcids();
        if (compactness.naccessPath != null) {
            try {
                Surface surface = new Surface(proteinComplex,
                                    new Naccess(
                                    compactness.naccessPath.getAbsolutePath()));
                sample = surface.getSurface();
            } catch (Exception e) {
                System.err.println("ERROR during NACCESS execution: " + e);
            }
        }
        AtomList sampleAtoms = new AtomList();
        for (AminoAcid aa : sample) {
            for (Atom atom : aa.getAllAtoms()) {
                AtomList env = MatterUtilities.getEnvironment(
                                    atom,
                                    proteinComplex.getAllAtoms(),
                                    compactness.radius);
                atom.setTemperatureFactor(env.size());
                sampleAtoms.add(atom);
            }
        }
        System.out.print(sampleAtoms);
    }
}

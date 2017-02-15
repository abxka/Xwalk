
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

import external.APBS;
import external.PDB2PQR;
import external.SIMS;

import mm.constants.Constants.ForceField;
import mm.electrostatics.ElectroStaticsUtilities;

import structure.constants.Constants;
import structure.io.Commandline;
import structure.io.pdb.PDBreader;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.hetgroups.SmallMolecule;
import structure.matter.protein.DotSurface;
import structure.matter.protein.PolyPeptideList;

/**
 * Class for calculating the electrostatic potential on a protein structure
 * in PDB format using PDB2PQR and APBS software.
 * @author Abdullah Kahraman
 * @version 0.5
 * @since 0.5
 */
public class Electro {
    //--------------------------------------------------------------------------
    // OBJECTS
    //--------------------------------------------------------------------------
    /**
     * Path to the PDB formatted file of the protein complex.
     */
    private File inFile;
    /**
     * Path to a PDB file, holding sample points coordinates.
     */
    private File sampleFile;
    /**
     * Path to a PQR file of the protein complex.
     */
    private File pqrFile;
    /**
     * Path to the APBS executable.
     */
    private File apbs;
    /**
     * Path to the multivalue executable.
     */
    private File multivalue;
    /**
     * Path to the PDB2PQR executable.
     */
    private File pdb2pqr;
    /**
     * Path to the SIMS executable.
     */
    private File sims;
    /**
     * Solve the linear Poisson-Boltzmann equation? Default {@code TRUE}.
     */
    private boolean lpbe = true;
    //--------------------------------------------------------------------------
    // CLASS METHODS
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param args
     *        String array holding all commandline arguments.
     */
    public Electro(final String[] args) {
        this.readCommandline(args);
    }
    //--------------------------------------------------------------------------
    /**
     * Reads all parameter from the commandline.
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
                           + "Information:" + nL
                           + "This program calculates the electrostatic "
                           + "potential of the protein on a dot surface "
                           + "representation of the same protein and maps them "
                           + "back on the protein atoms. The potential values "
                           + "will be written into the temperature factor "
                           + "column, while atomic partial charges from PARSE "
                           + "will be written into the occupancy column." + nL
                           + nL
                           + "Usage:" + nL
                           + "java " + this.getClass().getName()
                           + " -in 1b14.pdb "
                           + "-apbs /Applications/apbs-1.3/bin/apbs "
                           + "-multivalue /Applications/apbs-1.3/share/tools/"
                           + "mesh/multivalue "
                           + "-pdb2pqr /Applications/pdb2pqr-1.8/pdb2pqr.py "
                           + "-sims /Applications/SIMS/Sims_ex " + nL
                           + nL
                           + "Parameters:" + nL
                           + "\t-in <path>\tAny structure file in PDB format "
                           + "(required)." + nL
                           + "\t-sample <path>\tstructure file in PDB format "
                           + "on which HES values will be mapped. "
                           + "If not set, HES values will be calculated on "
                           + "protein atoms. ATOM and HETATM entries will be "
                           + "read in (optional)." + nL
                           + "\t-pqr <path>\tA PQR file of the input file, "
                           + "which will be used for the electrostatic "
                           + "calculation. PDB2PQR calculation will be skipped "
                           + "(optional)." + nL
                           + "\t-apbs <path>\tFull path to the APBS "
                           + "application (required). You can download it from "
                           + "http://apbs.sourceforge.net." + nL
                           + "\t-multivalue <path>\tFull path to APBS's "
                           + "multivalue application (required)." + nL
                           + "\t-pdb2pqr <path>\tFull path to the PDB2PQR "
                           + "application (required). Please ensure that "
                           + "PROPKA is compiled. You can download it from "
                           + "http://www.poissonboltzmann.org/pdb2pqr." + nL
                           + "\t-sims <path>\tFull path to the SIMS "
                           + "application (required if -sample is not set). "
                           + "You can download it from "
                           + "http://hekto.med.unc.edu:8080/HERMANS/software/"
                           + "SIMS/SIMS.html." + nL
                           );
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-in", true).equals("ERROR")) {
            System.err.println(nL + "ERROR: Could not find -in parameter. "
                                  + "Please check the help page with -help "
                                  + "!!!" + nL);
            System.exit(1);
        } else {
            String file = Commandline.get(args, "-in", true);
            this.inFile = new File(file);
            if (!this.inFile.exists()) {
                System.err.print("ERROR: " + file + " does not exist!" + nL);
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
        //----------------------------------------------------------------------
        if (!Commandline.get(args, "-pqr", true).equals("ERROR")) {
            String file = Commandline.get(args, "-pqr", true);
            this.pqrFile = new File(file);

            if (!this.pqrFile.exists()) {
                System.err.print(nL
                              + "Couldn't open PQR file \""
                              + this.pqrFile.getAbsolutePath()
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (Commandline.get(args, "-apbs", true).equals("ERROR")) {
            System.err.println(nL + "ERROR: Could not find -apbs parameter. "
                                  + "Please check the help page with -help "
                                  + "!!!" + nL);
            System.exit(1);
        } else {
            String file = Commandline.get(args, "-apbs", true);
            this.apbs = new File(file);

            if (!this.apbs.exists()) {
                System.err.print("ERROR: " + file + " does not exist!" + nL);
                System.exit(1);
            }
            if (!this.apbs.canExecute()) {
                System.err.print("ERROR: Can't execute " + file + "!" + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (Commandline.get(args, "-multivalue", true).equals("ERROR")) {
            System.err.println(nL + "ERROR: Could not find -multivalue "
                                  + "parameter. Please check the help page "
                                  + "with -help !!!" + nL);
            System.exit(1);
        } else {
            String file = Commandline.get(args, "-multivalue", true);
            this.multivalue = new File(file);

            if (!this.multivalue.exists()) {
                System.err.print("ERROR: " + file + " does not exist!" + nL);
                System.exit(1);
            }
            if (!this.multivalue.canExecute()) {
                System.err.print("ERROR: Can't execute " + file + "!" + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (this.pqrFile == null) {
            if (Commandline.get(args, "-pdb2pqr", true).equals("ERROR")) {
                System.err.println(nL + "ERROR: Could not find -pdb2pqr "
                                      + "parameter. Please check the help page "
                                      + "with -help !!!" + nL);
                System.exit(1);
            } else {
                String file = Commandline.get(args, "-pdb2pqr", true);
                this.pdb2pqr = new File(file);

                if (!this.pdb2pqr.exists()) {
                    System.err.print("ERROR: " + file + " does not exist!"
                                     + nL);
                    System.exit(1);
                }
                if (!this.pdb2pqr.canExecute()) {
                    System.err.print("ERROR: Can't execute " + file + "!" + nL);
                    System.exit(1);
                }
            }
        }
        //----------------------------------------------------------------------
        // consider SIMS calculations only if no sample file has been provided.
        if (this.sampleFile == null) {
            if (Commandline.get(args, "-sims", true).equals("ERROR")) {
                System.err.println(nL
                                  + "ERROR: Could not find -sims parameter. "
                                  + "Please check the help page with -help "
                                  + "!!!" + nL);
                System.exit(1);
            } else {
                String file = Commandline.get(args, "-sims", true);
                this.sims = new File(file);

                if (!this.sims.exists()) {
                    System.err.print("ERROR: " + file + " does not exist!"
                                    + nL);
                    System.exit(1);
                }
                if (!this.sims.canExecute()) {
                    System.err.print("ERROR: Can't execute " + file + "!" + nL);
                    System.exit(1);
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns REMARK lines with statistics on the non-interface part of a
     * protein surface.
     * @param sampleCoords
     *        AtomList object holding the coordinates of sample points for
     *        which electrostatic potentials were calculated.
     * @return String object holding PDB-like REMARK lines with statistics on
     *         the non-interface part of the protein surface.
     */
    private String getRemarks(final AtomList sampleCoords) {
        String nL = Constants.LINE_SEPERATOR;
        NumberFormat dec = xwalk.constants.Constants.DISTANCE_DEC_FORMAT;
        StringBuffer output = new StringBuffer();

        float potentialSum = 0;
        for (Atom atom : sampleCoords) {
            potentialSum += atom.getPotential();
        }

        output.append("HEADER   " + this.inFile.getName()
                    + nL);
        output.append("REMARK   ATOM COUNT: "
                    + sampleCoords.size() + nL);
        output.append("REMARK   ELECTROPOT SUM: "
                    + dec.format(potentialSum) + nL);
        output.append("REMARK   AVERAGE: "
                    + dec.format(potentialSum / sampleCoords.size())
                    + nL);
        output.append("REMARK" + nL);
        output.append("REMARK   Occupancy Column: PARSE partial charge" + nL);
        output.append("REMARK   B-factor Column: ELECTROPOT value" + nL);
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
        Electro electro = new Electro(args);
        //----------------------------------------------------------------------
        // Read in PDB files
        //----------------------------------------------------------------------
        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(
                                              electro.inFile.getAbsolutePath());

        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

        PolyPeptideList proteinComplex =
                                readers.get(0).getEntireProteinComplex().get(0);

        if (readers.get(0).getAllSmallMolecules().size() > 0) {
            proteinComplex.addSmallMolecules(
                                  readers.get(0).getAllSmallMolecules().get(0));
        }

        AtomList sampleCoords;
        if (electro.sampleFile == null) {
            sampleCoords =  new DotSurface(proteinComplex,
                                           new SIMS(electro.sims));
        } else {
            try {
                readers = PDBreader.createPDBreaders(
                                          electro.sampleFile.getAbsolutePath());
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

        if (electro.pqrFile == null) {
            //------------------------------------------------------------------
            // Run PDB2PQR on protein complex
            //------------------------------------------------------------------

            String pqrParameter = "--chain --nodebump --with-ph=7";
            electro.pqrFile = PDB2PQR.run(electro.pdb2pqr,
                                       electro.inFile,
                                       ForceField.PARSE,
                                       pqrParameter,
                                       true);
        }
//        try {
//            PDB2PQR.adjust(pqrFile, true);
//        } catch (Exception e) {
//            System.err.println("ERROR while trying to adjust PQR file: " + e);
//        }
        ElectroStaticsUtilities.assignPartialCharges(proteinComplex,
                                                     electro.pqrFile);
        String pqrFileName =
           electro.pqrFile.getName().substring(0,
                                    electro.pqrFile.getName().lastIndexOf("."));
        File apbsParameterFile = APBS.createInputFile(
                        electro.pqrFile,
                        pqrFileName,
                        sampleCoords,
                        electro.lpbe,
                        mm.constants.Constants.PROTEIN_DIELECTRIC_CONSTANT,
                        mm.constants.Constants.SOLVENT_DIELECTRIC_CONSTANT,
                        mm.constants.Constants.TWENTYFIVE_CELSIUS_IN_FAHRENHEIT
                                                       );
        APBS.run(electro.apbs, apbsParameterFile);
        ElectroStaticsUtilities.assignElectrostaticPotentials(
                                                  sampleCoords,
                                                  new File(pqrFileName + ".dx"),
                                                  electro.multivalue);
        AtomList hetgroupAtoms = new AtomList();
        if (electro.sampleFile == null) {
            ElectroStaticsUtilities.mapElectrostaticPotentials(
                                                   proteinComplex.getAllAtoms(),
                                                   sampleCoords);
            sampleCoords = proteinComplex.getAllAtoms();
            for (SmallMolecule hetgroup : proteinComplex.getSmallMolecules()) {
                sampleCoords.addAll(hetgroup.getAllAtoms());
                hetgroupAtoms.addAll(hetgroup.getAllAtoms());
            }
        }

        for (int i = 0; i < sampleCoords.size(); i++) {
            sampleCoords.get(i).setTemperatureFactor(
                                            sampleCoords.get(i).getPotential());
            if (electro.sampleFile == null) {
                if (i < proteinComplex.getAllAtoms().size()) {
                    sampleCoords.get(i).setOccupancy(
                        proteinComplex.getAllAtoms().get(i).getPartialCharge());
                } else {
                       sampleCoords.get(i).setOccupancy(
                       hetgroupAtoms.get(i - proteinComplex.getAllAtoms().size()
                                                       ).getPartialCharge());
                }
            } else {
                sampleCoords.get(i).setOccupancy(0);
            }
        }
        //------------------------------------------------------------------
        // Output EP statistics and coordinates with EP in temperature factor
        // columns and partial charges in the occupancy column.
        //------------------------------------------------------------------
        String output = electro.getRemarks(sampleCoords);
        System.out.print(output + sampleCoords.toString());
    }
}

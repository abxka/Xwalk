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
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;

import external.ExternCommand;
import external.Naccess;

import mm.evolution.Consurf;
import mm.hydrophobicity.Hydrophobicity;

import structure.constants.Constants;
import structure.io.Commandline;
import structure.io.ReadFile;
import structure.io.pdb.PDBreader;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;
import structure.sas.BindingInterface;

/**
 * Class holding a main method to calculate the binding interfaces for
 * a protein complex.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Interface {
    /**
     * Empty Constructor.
     */
    protected Interface() {
    }
    //--------------------------------------------------------------------------
    /**
     * Path to the PDB formatted file of the protein complex.
     */
    private String pdbFile;
    /**
     * Do XlogP hydrophobicity calculation.
     */
    private boolean doHes = false;
    /**
     * Do evolutionary conservation estimation with ConSurf.
     */
    private boolean doConsurf = false;
    /**
     * Path to the text file with list of ConSurf files to be used in the
     * conservation calculation.
     */
    private String consurfFile;
    /**
     * PDB ID of the protein complex.
     */
    private String pdbId;
    /**
     * Path to the naccess program.
     */
    private String naccessPath;
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
                             + "java " + Interface.class.getName() + " -help"
                             + nL);
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-help", false).equals("EXISTS")) {
            System.err.println(nL
                           + "USAGE:" + nL
                           + "\tjava " + Interface.class.getName()
                           + " -in 1b14.pdb" + nL
                           + nL
                           + "INFORMATION:" + nL
                           + "\tThis program determines all binary interfaces "
                           + "in a protein complex. An interface is defined by "
                           + "those amino acids that change their absolute "
                           + "total solvent accessiblilty by at least 5% upon "
                           + "complex formation. In addition hydrophobic "
                           + "environment scores can be calculated and "
                           + "the evolutionary conservation can be determined "
                           + "and mapped on interface amino acids." + nL
                           + nL
                           + "\tNOTE: During the calculations temporary files "
                           + "starting with proteinA, proteinB and proteinAB "
                           + "and with suffixes .asa, .rsa and .log will be "
                           + "generated, which can be deleted manually "
                           + "afterwards." + nL
                           + nL
                           + "PARAMETERS:" + nL
                           + "\t-in <path>\tany structure file in PDB format "
                           + "(required)." + nL
                           + "\t-naccess <string>\tPath to the naccess "
                           + "executable, which will be used to calculate "
                           + "the buried surface area at the interfaces and "
                           + "the average HES and conservation grades for "
                           + "non-interface surface amino acids (required)."
                           + nL
                           + "\t-hes [switch]\tCalculates hydrophobic "
                           + "enviroment scores for each interface atom and "
                           + "assigns them as occupancy values (optional)." + nL
                           + "\t-consurf [switch]\tDownloads evolutionary "
                           + "conservation files from the ConSurf server and "
                           + "assings them to each interface amino acid's "
                           + "temperature factor (optional)." + nL
                           + "\t-consurf <path>\ttext file listing ConSurf "
                           + "grade files to all protein components of the "
                           + "structure file (optional)." + nL
                           + "\t-id [string]\tThe PDB Id of the structure "
                           + "file, in the case that the structure file is not "
                           + "named by its PDB Id and a ConSurf file should "
                           + "be downloaded." + nL
                           );
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-in", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-in\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            this.pdbFile = Commandline.get(args, "-in", true);

            if (!ReadFile.exists(this.pdbFile)) {
                System.err.print(nL
                              + "Couldn't open infile \"" + this.pdbFile
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------
        if (Commandline.get(args, "-naccess", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter"
                             + "\"-naccess\" !!!" + nL);
            System.exit(1);
        } else {
            this.naccessPath = Commandline.get(args, "-naccess", true);
            if (!ExternCommand.exists(this.naccessPath)) {
                System.err.print("ERROR: " + this.naccessPath + " is not "
                               + "executable" + nL);
                System.exit(1);
            }
        }
        //----------------------------
        if (Commandline.get(args, "-hes", false).equals("EXISTS")) {
            doHes = true;
        }
        //----------------------------
        if (Commandline.get(args, "-consurf", false).equals("EXISTS")) {
            doConsurf = true;

            // check whether path has been given by the user
            this.consurfFile = Commandline.get(args, "-consurf", true);
            if (this.consurfFile.startsWith("-")) {
                this.consurfFile = null;
            }
        }
        //----------------------------
        if (!Commandline.get(args, "-id", true).equals("ERROR")) {
            this.pdbId = Commandline.get(args, "-id", true);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Determines and sets the conservation grade for each protein in a protein
     * complex using either the ConSurf-DB or local copies of ConSurf
     * conservation grade files.
     * @param proteinComplex
     *        Protein complex to which conservation will be determined.
     * @throws IOException when an error occurs while reading in ConSurf files
     *         from the local drive or from the ConSurf web database.
     */
    private void setConservation(final PolyPeptideList proteinComplex)
                                                            throws IOException {

        String fileName = new File(this.pdbFile).getName().replaceAll(
                                                                     "\\..*", ""
                                                                     );
        String pdbIdent = null;
        if (fileName.matches("[1-9][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z]")
            &&
            fileName.length() == 4) {
            pdbIdent = fileName;
        } else if (this.pdbId != null) {
            pdbIdent = this.pdbId;
        }
        if (pdbIdent != null) {
            for (PolyPeptide protein : proteinComplex) {
                try {
                    Consurf.setConservation(protein, pdbIdent);
                } catch (IOException e) {
                    throw new IOException("ERROR: Could not download ConSurf "
                                        + "grade file for PDB Id " + fileName
                                        + ". " + e.getMessage());
                }
            }
        } else if (this.consurfFile != null) {
            ReadFile read = new ReadFile(this.consurfFile);
            for (PolyPeptide protein : proteinComplex) {
                char proteinChainId = protein.get(0).getAtom(0).getChainId();
                boolean found = false;

                Consurf consurf = null;
                for (String line : read) {
                     consurf = new Consurf(line);
                     if (consurf.getChainId() == proteinChainId) {
                         found = true;
                         break;
                     }
                }
                if (!found) {
                    System.err.println("ERROR: ConSurf file for chain "
                                      + proteinChainId + " is missing.");
                    System.exit(0);
                }
                consurf.assignConservation(protein);
            }
        } else {
            System.err.println("Please specify a ConSurf grade file");
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all BindingInterfaces found between two proteins in the protein
     * complex.
     * @param proteinComplex
     *        Protein complex from which binding sites will be extracted.
     * @param naccess
     *        Naccess object to extract SASA for a protein complex and its
     *        components.
     * @return List of BindingInterface object found in the protein complex.
     */
    private ArrayList<BindingInterface> getInterfaces(
                                           final PolyPeptideList proteinComplex,
                                           final Naccess naccess) {
        ArrayList<BindingInterface> complexInterfaces =
                                              new ArrayList<BindingInterface>();
        for (int j = 0; j < proteinComplex.size(); j++) {
            for (int k = j + 1; k < proteinComplex.size(); k++) {
                BindingInterface bi = new BindingInterface(
                                                          proteinComplex.get(j),
                                                          proteinComplex.get(k),
                                                          naccess);
                if (bi.getInterface().get(0).size() > 0) {
                    complexInterfaces.add(bi);
                }
            }
        }
        return complexInterfaces;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines all amino acids on the surface of a protein complex using
     * the third party external application NACCESS and a relative ASA of 5%.
     * See Janin, Miller, Chothia (1988), JMB.
     * @param proteinComplex
     *        ProteinComplex object
     * @return List of AminoAcid objects where each amino acid is found at the
     *         surface as given by its relative ASA.
     * @throws IOException if an error occurs while executing NACCESS.
     */
    private ArrayList<AminoAcid> getSurfaceAminoAcids(
                                            final PolyPeptideList proteinComplex
                                                     ) throws IOException {
        ArrayList<AminoAcid> surfaceAminoAcids = new ArrayList<AminoAcid>();

        if (this.naccessPath != null) {
            Naccess naccess = new Naccess(this.naccessPath);
            naccess.run(this.pdbFile);
            naccess.setSolventAccessibility(proteinComplex.getAllAminoAcids());

            for (AminoAcid aa : proteinComplex.getAllAminoAcids()) {
                if (aa.getRelativeSas()
                                      > Constants.MINIMUM_REL_SASA_ON_SURFACE) {
                    surfaceAminoAcids.add(aa);
                }
            }
        }
        return surfaceAminoAcids;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns REMARK lines with statistics on the non-interface part of a
     * protein surface.
     * @param nonInterfaceSurfaceAminoAcids
     *        List of AminoAcids that form the non-interface surface of a
     *        protein complex.
     * @return String object holding PDB-like REMARK lines with statistics on
     *         the non-interface part of the protein surface.
     */
    private String getRemarksOnNonInterface(
                     final ArrayList<AminoAcid> nonInterfaceSurfaceAminoAcids) {
        String nL = Constants.LINE_SEPERATOR;
        NumberFormat dec = xwalk.constants.Constants.DISTANCE_DEC_FORMAT;
        StringBuffer output = new StringBuffer();

        int nonInterfaceAaCount = 0;
        int nonInterfaceAtomCount = 0;
        float nonInterfaceSumHes = 0;
        int nonInterfaceSumConservation = 0;
        int nonInterfaceSumTypicInterfaceAA = 0;
        int nonInterfaceSumTypicNonInterfaceAA = 0;

        int i = 0;
        if (nonInterfaceSurfaceAminoAcids != null) {
            for (AminoAcid aa : nonInterfaceSurfaceAminoAcids) {
                nonInterfaceAaCount++;
                nonInterfaceSumConservation += aa.getConservationGrade();
                if (MatterUtilities.isInterfaceTypic(aa)) {
                    nonInterfaceSumTypicInterfaceAA++;
                }
                if (MatterUtilities.isNonInterfaceTypic(aa)) {
                    nonInterfaceSumTypicNonInterfaceAA++;
                }
                for (Atom atom : aa.getAllAtoms()) {
                    nonInterfaceAtomCount++;
                        nonInterfaceSumHes += atom.getHes();
                }
            }
        }

        output.append("HEADER   " + new File(this.pdbFile).getName()
                    + nL);
        output.append("REMARK   0  AMINO ACID COUNT AT NON-INTERFACE: "
                    + nonInterfaceAaCount + nL);
        output.append("REMARK   0  ATOM COUNT AT NON-INTERFACE: "
                    + nonInterfaceAtomCount + nL);
        output.append("REMARK   0  INTERFACE TYPIC RESIDUES AT "
                    + "NON-INTERFACE: " + nonInterfaceSumTypicInterfaceAA + nL);
        output.append("REMARK   0  NON-INTERFACE TYPIC RESIDUES AT "
                    + "NON-INTERFACE: "
                    + nonInterfaceSumTypicNonInterfaceAA + nL);
        output.append("REMARK   0  AVERATE INTERFACE TYPIC RESIDUES AT "
                    + "NON-INTERFACE: "
                    + dec.format((float) nonInterfaceSumTypicInterfaceAA
                                 /
                                 nonInterfaceAaCount) + nL);
        output.append("REMARK   0  AVERAGE NON-INTERFACE TYPIC RESIDUES AT "
                    + "NON-INTERFACE: "
                    + dec.format((float) nonInterfaceSumTypicNonInterfaceAA
                                 /
                                 nonInterfaceAaCount) + nL);
        output.append("REMARK   0  HES SUM AT NON-INTERFACE: "
                    + dec.format(nonInterfaceSumHes) + nL);
        output.append("REMARK   0  CONSERVATION SUM AT NON-INTERFACE: "
                    + nonInterfaceSumConservation + nL);
        output.append("REMARK   0  AVERAGE HES AT NON-INTERFACE: "
                    + dec.format(nonInterfaceSumHes / nonInterfaceAtomCount)
                    + nL);
        output.append("REMARK   0  AVERAGE CONSERVATION AT NON-INTERFACE: "
                    + dec.format((float) nonInterfaceSumConservation
                                 /
                                 nonInterfaceAaCount) + nL);
        output.append("REMARK" + nL);
        output.append("REMARK   Occupancy Column: HES value" + nL);
        output.append("REMARK   B-factor Column: Conservation grade" + nL);
        output.append("REMARK" + nL);
    return output.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns REMARK lines with statistics on the interfaces of a protein
     * complexes.
     * @param complexInterfaces
     *        List of BindingInterface objects found in a protein complex.
     * @return String object holding PDB-like REMARK lines with statistics on
     *         the protein interfaces.
     */
    private String getRemarksOnInterface(
                          final ArrayList<BindingInterface> complexInterfaces) {

        StringBuffer output = new StringBuffer();
        String nL = Constants.LINE_SEPERATOR;
        NumberFormat dec = xwalk.constants.Constants.DISTANCE_DEC_FORMAT;

        int j = 0;
        int k = 1;
        for (BindingInterface bi : complexInterfaces) {
            int interfaceAaCount = 0;
            int interfaceAtomCount = 0;
            float interfaceSumHes = 0;
            int interfaceSumConservation = 0;
            int interfaceTypicInterfaceAA = 0;
            int interfaceTypicNonInterfaceAA = 0;
            AtomList atoms = new AtomList();

            for (ArrayList<AminoAcid> halfInterface : bi.getInterface()) {
                for (AminoAcid aa : halfInterface) {
                    interfaceAaCount++;
                    interfaceSumConservation += aa.getConservationGrade();
                    if (MatterUtilities.isInterfaceTypic(aa)) {
                        interfaceTypicInterfaceAA++;
                    }
                    if (MatterUtilities.isNonInterfaceTypic(aa)) {
                        interfaceTypicNonInterfaceAA++;
                    }
                    for (Atom atom : aa.getAllAtoms()) {
                        interfaceAtomCount++;
                        if (this.doHes) {
                            atom.setOccupancy(atom.getHes());
                            interfaceSumHes += atom.getOccupancy();
                        }
                        if (this.doConsurf) {
                            atom.setTemperatureFactor(
                                                       aa.getConservationGrade()
                                                     );
                        }
                        atoms.add(atom);
                    }
                }
            }
            float bsa = -1;
            if (this.naccessPath != null) {
                bsa = bi.getBSA();
            }

            output.append("REMARK   " + k + "  TITLE " + " "
                       + bi.getInterface().get(0).get(0).getAtom(0).getChainId()
                       + "-"
                       + bi.getInterface().get(1).get(0).getAtom(0).getChainId()
                       + nL);
            output.append("REMARK   " + k + "  AMINO ACID COUNT AT INTERFACE: "
                        + interfaceAaCount + nL);
            output.append("REMARK   " + k + "  ATOM COUNT AT INTERFACE: "
                        + interfaceAtomCount + nL);
            output.append("REMARK   " + k + "  BURIED SURFACE AREA AT "
                        + "INTERFACE: "
                        + dec.format(bsa) + nL);
            output.append("REMARK   " + k + "  INTERFACE TYPIC RESIDUES AT "
                        + "INTERFACE: " + interfaceTypicInterfaceAA + nL);
            output.append("REMARK   " + k + "  NON-INTERFACE TYPIC RESIDUES AT "
                        + "INTERFACE: " + interfaceTypicNonInterfaceAA + nL);
            output.append("REMARK   " + k + "  AVERAGE INTERFACE TYPIC "
                        + "RESIDUES AT INTERFACE: "
                        + dec.format((float) interfaceTypicInterfaceAA
                                     /
                                     interfaceAaCount) + nL);
        output.append("REMARK   " + k + "  AVERAGE NON-INTERFACE TYPIC "
                        + "RESIDUES AT INTERFACE: "
                        + dec.format((float) interfaceTypicNonInterfaceAA
                                     /
                                     interfaceAaCount) + nL);
            output.append("REMARK   " + k + "  HES SUM AT INTERFACE: "
                        + dec.format(interfaceSumHes) + nL);
            output.append("REMARK   " + k + "  CONSERVATION SUM AT INTERFACE: "
                        + interfaceSumConservation + nL);
            output.append("REMARK   " + k + "  AVERAGE HES AT INTERFACE: "
                        + dec.format(interfaceSumHes / interfaceAtomCount)
                        + nL);
            output.append("REMARK   " + k + "  AVERAGE CONSERVATION AT "
                        + "INTERFACE: "
                        + dec.format((float) interfaceSumConservation
                                     /
                                     interfaceAaCount) + nL);
            output.append("REMARK" + nL);
            k++;
        }
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
        Interface interfase = new Interface();
        interfase.readCommandline(args);

        //----------------------------------------------------------------------
        // Read in PDB files
        //----------------------------------------------------------------------
        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(interfase.pdbFile);

        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

        PolyPeptideList proteinComplex =
                                readers.get(0).getEntireProteinComplex().get(0);

        if (interfase.doConsurf) {
            //------------------------------------------------------------------
            // Determine conservation grade for all amino acids in the complex.
            //------------------------------------------------------------------
            try {
                interfase.setConservation(proteinComplex);
            } catch (IOException e) {
                System.err.println("ERROR: Problems while assigning "
                                 + "conservation grades: " + e.getMessage());
            }
        }

        //----------------------------------------------------------------------
        // Calculate interface
        //----------------------------------------------------------------------
        Naccess naccess = null;
        try {
            naccess = new Naccess(interfase.naccessPath);
        } catch (IOException e) {
            System.err.println("ERROR while executing NACCESS: "
                             + e.getMessage());
            System.exit(-1);
        }
        ArrayList<BindingInterface> complexInterfaces =
                                                interfase.getInterfaces(
                                                                proteinComplex,
                                                                naccess
                                                                       );
        ArrayList<AminoAcid> interfaceAminoAcids =
                                                     new ArrayList<AminoAcid>();
        for (BindingInterface bi : complexInterfaces) {
            if (bi.getInterface().size() >= 2) {
                interfaceAminoAcids.addAll(bi.getInterface().get(0));
                interfaceAminoAcids.addAll(bi.getInterface().get(1));
            }
        }
        //----------------------------------------------------------------------
        // Calculate non-interface surface amino acids with NACCESS
        //----------------------------------------------------------------------

        ArrayList<AminoAcid> nonInterfaceSurfaceAminoAcids =
                                                     new ArrayList<AminoAcid>();

        if (interfase.naccessPath != null) {
            ArrayList<AminoAcid> surfaceAminoAcids = null;
            try {
                surfaceAminoAcids =
                                 interfase.getSurfaceAminoAcids(proteinComplex);
            } catch (IOException e) {
                System.err.println("ERROR: while executing NACCESS: "
                                 + e.getMessage());
            }

            for (AminoAcid surfaceAA : surfaceAminoAcids) {
                boolean found = false;
                for (AminoAcid interfaceAA : interfaceAminoAcids) {
                    if (interfaceAA == surfaceAA) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    nonInterfaceSurfaceAminoAcids.add(surfaceAA);
                }
            }
        }
        //----------------------------------------------------------------------
        // Extract Cartesian coordinates of interface atoms
        //----------------------------------------------------------------------
        AtomList interfaceCoords = new AtomList();
        for (BindingInterface bi : complexInterfaces) {
            for (ArrayList<AminoAcid> halfInterface : bi.getInterface()) {
                for (AminoAcid aa : halfInterface) {
                    for (Atom atom : aa.getAllAtoms()) {
                        interfaceCoords.add(atom);
                    }
                }
            }
        }
        if (interfase.doHes) {
            //------------------------------------------------------------------
            // Calculate HES for each interface atom
            //------------------------------------------------------------------
            Hydrophobicity hes = new Hydrophobicity(proteinComplex);
            hes.mapHydrophobicity(interfaceCoords);

            if (interfase.naccessPath != null) {
                AtomList nonInterfaceSurfaceCoords = new AtomList();
                for (AminoAcid nonInterfaceSurfaceAA
                                              : nonInterfaceSurfaceAminoAcids) {
                    for (Atom atom : nonInterfaceSurfaceAA.getAllAtoms()) {
                        nonInterfaceSurfaceCoords.add(atom);
                    }
                }
                hes.mapHydrophobicity(nonInterfaceSurfaceCoords);
            }
        }

        //------------------------------------------------------------------
        // Output interface with HES and conservation grades listed in the
        // occupancy and temperature factor columns.
        //------------------------------------------------------------------

        String output = interfase.getRemarksOnNonInterface(
                                                   nonInterfaceSurfaceAminoAcids
                                                          );
        output += interfase.getRemarksOnInterface(complexInterfaces);
        if (interfase.doHes || interfase.doConsurf) {
            for (BindingInterface complexInterface : complexInterfaces) {
                for (ArrayList <AminoAcid> interfaceResidues
                                            : complexInterface.getInterface()) {
                    for (AminoAcid interfaceResidue : interfaceResidues) {
                        for (Atom atom
                                             : interfaceResidue.getAllAtoms()) {
                            if (interfase.doHes) {
                                atom.setOccupancy(atom.getHes());
                            }
                            if (interfase.doConsurf) {
                                atom.setTemperatureFactor(
                                         interfaceResidue.getConservationGrade()
                                                         );
                            }
                        }
                    }
                }
            }
        }
        System.out.print(output);
        for (BindingInterface bindingSite : complexInterfaces) {
            System.out.print(bindingSite);
        }
    }
}

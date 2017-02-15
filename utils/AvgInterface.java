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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import external.Naccess;

import structure.constants.Constants;
import structure.io.pdb.PDBreader;
import structure.matter.Atom;
import structure.matter.MatterUtilities;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;
import structure.sas.BindingInterface;

/**
 * Class holding a main method to calculate the average binding interface for
 * a protein complex have different conformational topologies.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class AvgInterface {

    /**
     * List of Atom objects forming all binding site atoms.
     */
    private ArrayList<AminoAcid> allInterfacesAminoAcids =
                                                     new ArrayList<AminoAcid>();
    /**
     * Number of occurrence of an amino acid at an interface. Need to use a
     * String representation of an amino acid as otherwise identical amino
     * acids from different conformations wont be recognized.
     */
    private Hashtable<String, Integer> interfaceCount =
                                            new Hashtable<String, Integer>();

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param readers
     *        List of PDBreader objects holding each the content of a single PDB
     *        file. The PDB files might be compressed into a tar or gzip file.
     * @param naccess
     *        Naccess object to extract SASA for a protein complex and its
     *        components.
     */
    public AvgInterface(final ArrayList<PDBreader> readers,
                        final Naccess naccess) {
        this.setUnionOfBindingInterfaces(readers, naccess);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a string representation of an amino acid.
     * @param aa
     *        AminoAcid object
     * @return String object holding the string representation of the amino
     *         acid.
     */
    public static String getResidueId(final AminoAcid aa) {
        return "#"
             + aa.getAtom(0).getResidueName()
             + aa.getAtom(0).getResidueNumber()
             + aa.getAtom(0).getChainId()
             + "#";
    }
    //--------------------------------------------------------------------------
    /**
     * Determines all binding interfaces found in a list of PDB files.
     * @param readers
     *        List of PDBreader objects holding each the content of a single PDB
     *        file. The PDB files might be compressed into a tar or gzip file.
     * @param naccess
     *        Naccess object to extract SASA for a protein complex and its
     *        components.
     */
    private void setUnionOfBindingInterfaces(
                                             final ArrayList<PDBreader> readers,
                                             final Naccess naccess
                                            ) {

        String allInterfacesIds = "";


        for (PDBreader file : readers) {
            ArrayList<PolyPeptideList> fileMolecules =
                                                file.getEntireProteinComplex();
            for (int i = 0; i < fileMolecules.size(); i++) {
                PolyPeptideList proteinComplex = fileMolecules.get(i);

                ArrayList<AminoAcid> complexInterfaceAminoAcids =
                                                     new ArrayList<AminoAcid>();
                String complexInterfaceAminoAcidIds = "";

                for (int j = 0; j < proteinComplex.size(); j++) {
                    for (int k = j + 1; k < proteinComplex.size(); k++) {

                        //------------------------------------------------------
                        // get binding interface between both proteins
                        //------------------------------------------------------
                        BindingInterface bi = new BindingInterface(
                                                          proteinComplex.get(j),
                                                          proteinComplex.get(k),
                                                          naccess
                                                                  );
                        ArrayList<ArrayList<AminoAcid>> complexInterface =
                                                              bi.getInterface();

                        //------------------------------------------------------
                        // calculate average number of occurrence of an amino
                        // acid at the binding interface
                        //------------------------------------------------------
                        for (ArrayList<AminoAcid> interfaceHalf
                                                           : complexInterface) {
                            for (AminoAcid aa : interfaceHalf) {
                                String id = AvgInterface.getResidueId(aa);

                                if (allInterfacesIds.indexOf(id) == -1
                                    &&
                                    complexInterfaceAminoAcidIds.indexOf(id)
                                                                        == -1) {
                                    this.interfaceCount.put(id, 1);

                                    complexInterfaceAminoAcids.add(aa);
                                    complexInterfaceAminoAcidIds += id;
                                } else if (complexInterfaceAminoAcidIds.indexOf(
                                                                              id
                                                                               )
                                                                        == -1) {
                                    this.interfaceCount.put(
                                                 id,
                                                 this.interfaceCount.get(id) + 1
                                                               );

                                    complexInterfaceAminoAcids.add(aa);
                                    complexInterfaceAminoAcidIds += id;
                                }
                            }
                        }
                    }
                }
                allInterfacesIds += complexInterfaceAminoAcidIds;
                this.allInterfacesAminoAcids.addAll(complexInterfaceAminoAcids);
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a list of all binding site atoms found in the union of the PDB
     * files.
     * @return List of AminoAcid objects that are found at the interface of the
     *         dimer.
     */
    public final ArrayList<AminoAcid> getInterfacesAminoAcids() {
        return this.allInterfacesAminoAcids;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the number of occurrences for each interfacial amino acid.
     * @return Hashtable with keys being the amino acid identifier, where the
     *         identifier consist of chainId, residueNumber and residueName.
     */
    public final Hashtable<String, Integer> getInterfacesCounts() {
        return this.interfaceCount;
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
        String nL = Constants.LINE_SEPERATOR;
        //----------------------------------------------------------------------
        // output small USAGE information if no commandline argument is given.
        //----------------------------------------------------------------------
        if (args.length < 3) {
            System.err.print(nL
                          + "Usage: " + nL
                          + nL
                          +  "java " + AvgInterface.class.getName()
                          + " decoys.tar.gz map.pdb "
                          + "/Application/naccess/naccess" + nL
                          + nL);
            System.exit(1);
        }

        //----------------------------------------------------------------------
        // read tar.gz file with all complexes to be analysed for amino acid
        // occurrences at binding interfaces.
        //----------------------------------------------------------------------
        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(args[0]);

        } catch (Exception e) {
            System.err.print(e.getMessage() + nL);
        }

        //----------------------------------------------------------------------
        // determine number of complexes in the tar.gz file
        //----------------------------------------------------------------------
        int complexCount = 0;
        for (PDBreader files : readers) {
            for (int i = 0; i < files.getEntireProteinComplex().size(); i++) {
                complexCount++;
             }
        }

        //----------------------------------------------------------------------
        // initiate NACCESS
        //----------------------------------------------------------------------
        Naccess naccess = null;
        try {
            naccess = new Naccess(args[2]);
        } catch (IOException e) {
            System.err.println("ERROR: Problems encounted while executing "
                             + "NACCESS " + e.getMessage());
            System.exit(-1);
        }

        //----------------------------------------------------------------------
        // calculate occurrence numbers for each amino acid at an interface
        //----------------------------------------------------------------------
        AvgInterface avgInterface = new AvgInterface(readers, naccess);
        ArrayList<AminoAcid> allInterfacesAminoAcids =
                                         avgInterface.getInterfacesAminoAcids();
        Hashtable<String, Integer> interfaceCounts =
                                             avgInterface.getInterfacesCounts();

        //----------------------------------------------------------------------
        // read map protein
        //----------------------------------------------------------------------
        PDBreader map = null;
        try {
            map = new PDBreader(args[1]);
        } catch (Exception e) {
            System.err.print(e.getMessage() + nL);
        }
        //----------------------------------------------------------------------
        // initialize occupancy and temperature factor values to 0.0.
        //----------------------------------------------------------------------
        ArrayList<PolyPeptideList> complexes = map.getEntireProteinComplex();
        for (PolyPeptideList complex : complexes) {
            for (PolyPeptide protein : complex) {
                for (AminoAcid aa : protein) {
                    for (Atom atom : aa.getAllAtoms()) {
                        atom.setOccupancy(0.0f);
                        atom.setTemperatureFactor(0.0f);
                    }
                }
            }
        }

        //----------------------------------------------------------------------
        // assign residues of map protein the occurrence numbers.
        //----------------------------------------------------------------------
        for (PolyPeptideList complex : complexes) {
            for (PolyPeptide protein : complex) {
                for (AminoAcid aa1 : protein) {
                    String id = AvgInterface.getResidueId(aa1);
                    for (AminoAcid aa2 : allInterfacesAminoAcids) {

                        if (MatterUtilities.equalsResidue(aa1.getAtom(0),
                                                          aa2.getAtom(0))) {

                            int count = interfaceCounts.get(id);
                            for (Atom atom : aa1.getAllAtoms()) {
                                atom.setOccupancy(count);
                                atom.setTemperatureFactor((float) count
                                                          /
                                                          complexCount);
                            }
                        }
                    }
                }
            }
            System.out.print(complex + "END" + nL);
        }
    }
}

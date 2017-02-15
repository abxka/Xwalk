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

package structure.sas;

import java.util.ArrayList;
import java.util.Hashtable;

import external.Naccess;

import structure.constants.Constants;
import structure.io.pdb.PDBwriter;
import structure.matter.Atom;
import structure.matter.MatterUtilities;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;


/**
 * Class for handling Interfaces between protein molecules.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class BindingInterface {

    //--------------------------------------------------------------------------
    /**
     * The binding interface consists of two list of amino acids, where each
     * list comes from one protein partner.
     */
    private ArrayList<ArrayList<AminoAcid>> bindingInterface =
                                         new ArrayList<ArrayList<AminoAcid>>(2);
    //--------------------------------------------------------------------------
    /**
     * Buried surface area of this interface.
     */
    private float bsa;
    //--------------------------------------------------------------------------
    /**
     * Constructor to calculate the amino acids at the binding interface
     * between two protein structures.
     * @param proteinA
     *        First protein of two to which interface is to be calculated.
     * @param proteinB
     *        Second protein of two to which interface is to be calculated.
     * @param naccess
     *        Naccess object to extract SASA for a protein complex and its
     *        components.
     */
    public BindingInterface(final PolyPeptide proteinA,
                            final PolyPeptide proteinB,
                            final Naccess naccess) {

        //----------------------------------------------------------------------
        // initiate bindingInterface object.
        //----------------------------------------------------------------------
        ArrayList<AminoAcid> interfaceHalf1 = new ArrayList<AminoAcid>();
        ArrayList<AminoAcid> interfaceHalf2 = new ArrayList<AminoAcid>();
        this.bindingInterface.add(interfaceHalf1);
        this.bindingInterface.add(interfaceHalf2);

        this.setInterface(proteinA, proteinB, naccess);
    }

    //--------------------------------------------------------------------------
    /**
     * Determines all amino acids at the binding interface of two proteins,
     * where the interface is defined by a absolute change in the total SASA of
     * an amino acid by > 0.1 Angstroem upon protein complex formation.
     *
     * @param proteinA
     *        First protein of two to which interface is to be calculated.
     * @param proteinB
     *        Second protein of two to which interface is to be calculated.
     * @param naccess
     *        Naccess object to extract SASA for a protein complex and its
     *        components.
     */
    private void setInterface(final PolyPeptide proteinA,
                              final PolyPeptide proteinB,
                              final Naccess naccess) {

        //----------------------------------------------------------------------
        // Extract amino acids on the protein surface
        //----------------------------------------------------------------------
        String tempFileName = "proteinA.pdb";
        PDBwriter write = new PDBwriter();
        write.setFile(tempFileName);
        write.write(proteinA);
        naccess.run(tempFileName);
        naccess.setSolventAccessibility(proteinA);
        float sasaA = naccess.getTotalSolventAccessibility();

        tempFileName = "proteinB.pdb";
        write.setFile(tempFileName);
        write.write(proteinB);
        naccess.run(tempFileName);
        naccess.setSolventAccessibility(proteinB);
        float sasaB = naccess.getTotalSolventAccessibility();

        tempFileName = "proteinAB.pdb";
        write.setFile(tempFileName);
        write.write(proteinA);
        write.setFile(tempFileName, true);
        write.write(proteinB);
        naccess.run(tempFileName);
        ArrayList<AminoAcid> aminoAcidsAB = new ArrayList<AminoAcid>();
        aminoAcidsAB.addAll(proteinA.copy());
        aminoAcidsAB.addAll(proteinB.copy());
        PolyPeptide proteinAB = new PolyPeptide(aminoAcidsAB);
        naccess.setSolventAccessibility(proteinAB);
        float sasaAB = naccess.getTotalSolventAccessibility();

        //----------------------------------------------------------------------
        // Calculate buried surface area from SASA of single proteins & complex
        //----------------------------------------------------------------------
        this.bsa = sasaA + sasaB - sasaAB;
        //----------------------------------------------------------------------
        // Put surface of protein A+B into a hash
        //----------------------------------------------------------------------
        Hashtable<Character, PolyPeptide> surfaces =
                                        new Hashtable<Character, PolyPeptide>();
        char chainIdA = proteinA.get(0).getAtom(0).getChainId();
        char chainIdB = proteinB.get(0).getAtom(0).getChainId();
        surfaces.put(chainIdA, proteinA);
        surfaces.put(chainIdB, proteinB);

        char minChainId = Character.toChars(Math.min(chainIdA, chainIdB))[0];

        //----------------------------------------------------------------------
        // Determine binding interface
        //----------------------------------------------------------------------
        for (AminoAcid aa1 : proteinAB) {
            char chainId = aa1.getAtom(0).getChainId();
            for (AminoAcid aa2 : surfaces.get(chainId)) {
                if (MatterUtilities.equalsResidue(aa1.getAtom(0),
                                                  aa2.getAtom(0))) {
                    float sasDiff = aa1.getTotalSas() - aa2.getTotalSas();
                    if (Math.abs(sasDiff)
                        >
                        Constants.MINUMUM_SASA_DIFFERECE_FOR_INTERFACE) {

                        if (minChainId == chainId) {
                            this.bindingInterface.get(0).add(aa1);
                        } else {
                            this.bindingInterface.get(1).add(aa2);
                        }
                    }
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the binding interface between protein1 and protein2.
     * @return List of AminoAcid objects that are found at the interface of the
     *         dimer.
     */
    public final ArrayList<ArrayList<AminoAcid>> getInterface() {
        return this.bindingInterface;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the buried buried surface area at this interface.
     * @return float value representing the buried surface area.
     */
    public final float getBSA() {
        return bsa;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the atom coordinates of this binding interface in PDB format.
     * @return String representation of the binding interface.
     */
    public final String toString() {
        StringBuffer output = new StringBuffer();
        for (AminoAcid aa : this.bindingInterface.get(0)) {
            for (Atom atom : aa.getAllAtoms()) {
                output.append(atom.toString());
            }
        }
        for (AminoAcid aa : this.bindingInterface.get(1)) {
            for (Atom atom : aa.getAllAtoms()) {
                output.append(atom.toString());
            }
        }
        return output.toString();
    }
    //--------------------------------------------------------------------------
}

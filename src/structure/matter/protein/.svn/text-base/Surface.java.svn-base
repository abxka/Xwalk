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

import structure.io.pdb.PDBwriter;
import external.Naccess;

/**
 * Class for calculating and storing amino acids that are located on the surface
 * of protein molecules.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Surface {
    //--------------------------------------------------------------------------
    /**
     * Minimum relative SASA of an amino acid in order to be considered as
     * a surface residue. See Miller, S., Janin, J., Lesk, A. M. & Chothia, C.
     * Interior and surface of monomeric proteins. J Mol Biol 196, 641-656
     * (1987).
     */
    private static final int MIN_RELATIVE_SOLVENT_ACCESSIBILITY = 5;
    //--------------------------------------------------------------------------
    /**
     * List of amino acids that are forming the surface of a protein.
     */
    private ArrayList<AminoAcid> surfaceAminoAcids = new ArrayList<AminoAcid>();
    //--------------------------------------------------------------------------
    /**
     * Total SASA of the protein molecule.
     */
    private double sasa;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param complex
     *        List of proteins to which the surface will be calculated.
     * @param naccess
     *        The SASA will be calculated using NACCESS.
     */
    public Surface(final PolyPeptideList complex, final Naccess naccess) {
        String tempFileName = "temp.pdb";
        PDBwriter write = new PDBwriter();
        write.setFile(tempFileName);
        write.write(complex);
        naccess.run(tempFileName);

        ArrayList<AminoAcid> allAA = new ArrayList<AminoAcid>();
        for (AminoAcid aa : complex.getAllAminoAcids()) {
            allAA.add(aa.copy());
        }

        naccess.setSolventAccessibility(allAA);
        this.sasa = naccess.getTotalSolventAccessibility();

        for (AminoAcid aa : allAA) {
            if (aa.getRelativeSas()
                >
                Surface.MIN_RELATIVE_SOLVENT_ACCESSIBILITY) {
                this.surfaceAminoAcids.add(aa);
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the amino acids that are located at the surface of a protein
     * molecule.
     * @return List of surface amino acids.
     */
    public final ArrayList<AminoAcid> getSurface() {
        return this.surfaceAminoAcids;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the total SASA of the protein complex.
     * @return double value representing the total SASA.
     */
    public final double getTotalSasa() {
        return this.sasa;
    }
}

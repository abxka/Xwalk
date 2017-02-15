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
import java.io.IOException;
import java.util.ArrayList;

import structure.io.ReadFile;
import structure.matter.Atom;
import structure.matter.MatterUtilities;
import structure.matter.protein.AminoAcid;

/**
 * Class making use of functionality provided by Simon Hubbard's NACCESS
 * application.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Naccess {

    /**
     * Path to the NACCESS application.
     */
    private String naccess;
    //--------------------------------------------------------------------------
    /**
     * Residue solvent accessibility file.
     */
    private ReadFile rsaFile;
    //--------------------------------------------------------------------------
    /**
     * Atom solvent accessibility file.
     */
    private ReadFile asaFile;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param naccessPath
     *        String holding the path to the NACCESS application.
     * @throws FileNotFoundException
     *         if NACCESS application can not be found at given path.
     */
    public Naccess(final String naccessPath) throws FileNotFoundException {
        if (!ReadFile.exists(naccessPath)) {
            throw new FileNotFoundException("ERROR: NACCESS application does "
                                          + "not reside at " + naccessPath);
        } else {
            this.naccess = naccessPath;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Runs NACCESS and stores the residue and atom accessibility files in
     * a ReadFile object.
     * @param pdbFileName
     *        String holding the path to the PDB file of the protein of which
     *        the total SAS is to be calculated.
     * @return {@code TRUE} if execution of NACCESS was successful,
     *         {@code FALSE} otherwise.
     */
    public final boolean run(final String pdbFileName) {
        String fileNameWithoutSuffix =
                        new File(pdbFileName).getName().replaceAll("\\..*", "");
        ExternCommand.execute(this.naccess + " " + pdbFileName, false);
        try {
            if (ReadFile.exists(fileNameWithoutSuffix + ".rsa")) {
                this.rsaFile = new ReadFile(fileNameWithoutSuffix + ".rsa");
            }
            if (ReadFile.exists(fileNameWithoutSuffix + ".asa")) {
                this.asaFile = new ReadFile(fileNameWithoutSuffix + ".asa");
            }
        } catch (IOException e) {
            return false;
        }
    return true;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the total solvent accessible surface area of the protein, which
     * is found in the TOTAL line of NACCESS's *.rsa files.
     * @return float value, representing the total SASA. If no information on
     *         the total SASA can be found in the .rsa file, then -1 will be
     *         returned.
     */
    public final float getTotalSolventAccessibility() {
        for (String line : this.rsaFile) {
            if (line.startsWith("TOTAL")) {
                String[] array = line.split("\\s+");
                return Float.parseFloat(array[1]);
            }
        }
        return -1;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the solvent accessibility of each amino acid.
     * @param aminoAcids
     *        List of AminoAcid objects for which the SASA is to be extracted.
     */
    @SuppressWarnings("unchecked")
    public final void setSolventAccessibility(
                                           final ArrayList<AminoAcid> aminoAcids
                                             ) {
        ArrayList<AminoAcid> copy = (ArrayList<AminoAcid>) aminoAcids.clone();

        for (String line : this.rsaFile) {
            if (line.startsWith("RES")) {
                Atom dummy = new Atom();
                dummy.setResidueName(line.substring(4, 7).trim());
                dummy.setResidueNumber(Integer.parseInt(
                                            line.substring(9, 13).trim()
                                                       ));
                dummy.setChainId(line.substring(8, 9).charAt(0));
                float totalSasa = Float.parseFloat(
                                                   line.substring(13, 22).trim()
                                                  );
                float relSasa = Float.parseFloat(
                                                   line.substring(22, 28).trim()
                                                  );
                for (AminoAcid aa : copy) {
                    if (MatterUtilities.equalsResidue(dummy, aa.getAtom(0))) {
                        aa.setRelativeSas(relSasa);
                        aa.setTotalSas(totalSasa);
                        copy.remove(aa);
                        break;
                    }
                }
            }
        }
    }
}

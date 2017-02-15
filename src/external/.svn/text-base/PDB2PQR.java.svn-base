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

import mm.constants.Constants.ForceField;

/**
 * Class provids functionalities to generate PDB2PQR files.
 * @author Abdullah Kahraman
 * @version 0.5
 * @since 0.5
 */
public class PDB2PQR {

    /**
     * Constructor with prevention against calls from subclass.
     */
    protected PDB2PQR() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Executes APBS on the command line.
     * @param pdb2pqr
     *        File object of the PDB2PQR executable.
     * @param pdbFile
     *        File object of the protein structure in PDB format.
     * @param forceField
     *        {@link ForceField} enum type of the force field from which partial
     *        charges and atom radii will be extracted.
     *        format.
     * @param addParam
     *        String with additional parameters that can be supplied before
     *        executing PDB2PQR.
     * @param verbose
     *        Switch for outputting additional information while PDB2PQR
     *        is progressing.
     * @return File object of the PQR file.
     */
    public static final File run(final File pdb2pqr,
                                 final File pdbFile,
                                 final ForceField forceField,
                                 final String addParam,
                                 final boolean verbose) {
        File pqrFile =
                new File(pdbFile.getName().replaceAll(".pdb$", ".pqr"));
        String version = ExternCommand.execute(pdb2pqr.getAbsolutePath()
                                             + " --version", false);
        if (!forceField.toString().equals("ROSETTA")
            ||
            (forceField.toString().equals("ROSETTA")
             &&
             (version.indexOf("1.5") != -1
              ||
              version.indexOf("1.6") != -1))) {
            ExternCommand.execute(pdb2pqr.getAbsolutePath()
                                + " --ff=" + forceField.toString()
                                + " " + addParam + " " + pdbFile
                                + " " + pqrFile, verbose);
        } else {
            String pdb2pqrDir = pdb2pqr.getParent();
            ExternCommand.execute(pdb2pqr.getAbsolutePath()
                                + " --userff=" + pdb2pqrDir
                                + "/dat/ROSETTA.DAT --usernames=" + pdb2pqrDir
                                + "/dat/ROSETTA.names " + addParam + " -v "
                                + pdbFile + " " + pqrFile,
                                verbose);
        }
        return pqrFile;
    }
}

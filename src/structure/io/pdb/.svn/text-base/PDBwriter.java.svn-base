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

package structure.io.pdb;

import java.util.Locale;

import structure.constants.Constants;
import structure.io.WriteFile;
import structure.matter.Atom;
import structure.matter.Molecule;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;


/**
 * Class for writing PDB files on the hard drive.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class PDBwriter extends WriteFile {

    /**
     * Operating system dependent newLine character.
     */
    private final String nL = Constants.LINE_SEPERATOR;

    //--------------------------------------------------------------------------
    /**
     * Writes a single ATOM or HETATM line into a PDB file.
     * @param atom
     *        - Atom object holding the information of the atom to be
     *          written out.
     * @return {@code TRUE} if writing process was successful; {@code FALSE}
     *          otherwise.
     * @see #write(Molecule)
     * @see #write(PolyPeptide)
     * @see #write(PolyPeptideList)
     */
    public final boolean write(final Atom atom) {
        StringBuffer output = new StringBuffer();
        output.append(PDBwriter.pdbFormat(atom) + nL);
        output.append("END" + nL);
        return super.write(output.toString());
    }
    //--------------------------------------------------------------------------

    /**
     * Writes all ATOM or HETATM lines of an entire molecule into a PDB file.
     * @param molecule
     *        - Molecule object.
     * @return {@code TRUE} if writing process was successful; {@code FALSE}
     *          otherwise.
     * @see #write(Atom)
     * @see #write(PolyPeptide)
     * @see #write(PolyPeptideList)
     */
    public final boolean write(final Molecule molecule) {
        StringBuffer output = new StringBuffer();
        for (Atom atom : molecule.getAllAtoms()) {
            output.append(PDBwriter.pdbFormat(atom) + nL);
        }
        output.append("END" + nL);
        return super.write(output.toString());
    }
    //--------------------------------------------------------------------------

    /**
     * Writes all ATOM or HETATM lines of a poly-peptide into a PDB file.
     * @param polyPeptide
     *        - PolyPeptide object having a list of {@code AminoAcid} objects,
     *          which in turns has a list of {@code Atom} objects
     * @return {@code TRUE} if writing process was successful; {@code FALSE}
     *         otherwise.
     * @see #write(Atom)
     * @see #write(Molecule)
     * @see #write(PolyPeptideList)
     */
    public final boolean write(final PolyPeptide polyPeptide) {
        StringBuffer output = new StringBuffer();
        for (AminoAcid aminoAcid : polyPeptide) {
            for (Atom atom : aminoAcid.getAllAtoms()) {
                output.append(PDBwriter.pdbFormat(atom) + nL);
            }
        }
        output.append("END" + nL);
        return super.write(output.toString());
    }
    //--------------------------------------------------------------------------

    /**
     * Writes all ATOM or HETATM lines of an entire molecular complex into a PDB
     * file.
     * @param complex
     *        - MolecularComplex object
     * @return {@code TRUE} if writing process was successful; {@code FALSE}
     *         otherwise.
     * @see #write(Atom)
     * @see #write(Molecule)
     * @see #write(PolyPeptide)
     */
    public final boolean write(final PolyPeptideList complex) {
        StringBuffer output = new StringBuffer();
        for (PolyPeptide polyPeptide : complex) {
             for (AminoAcid aminoAcid : polyPeptide) {
                  for (Atom atom : aminoAcid.getAllAtoms()) {
                       output.append(PDBwriter.pdbFormat(atom) + nL);
                }
            }
            output.append("TER" + nL);
        }
        output.append("END" + nL);
        return super.write(output.toString());
    }
    //--------------------------------------------------------------------------

    /**
     * Returns a {@code String} object holding the atom information in PDB
     * format.
     * @param atom
     *        Atom object holding the information of the atom to be written out.
     * @return String object holding the atom information in PDB format.
     * @see #pdbFormat(Atom, double, double)
     */
    public static String pdbFormat(final Atom atom) {
        return PDBwriter.pdbFormat(atom,
                                   atom.getOccupancy(),
                                   atom.getTemperatureFactor()
                                  );
    }
    //--------------------------------------------------------------------------
    /**
     * A more ugly way of writing out the information of an {@code Atom} object
     * in PDB format as described in
     * <a href = "http://www.wwpdb.org/documentation/format32/sect9.html">
     * here </a>.
     * Might need a more intelligent rewrite.
     * @param atom
     *        Atom object holding the information of the atom to be written out.
     * @param occupancy
     *        - Double value representing the occupancy value.
     * @param temperatureFactor
     *        - Double value representing the temperature factor.
     * @return String object holding the atom information in PDB format.
     * @see    #pdbFormat(Atom)
     */
    public static String pdbFormat(final Atom atom,
                                   final double occupancy,
                                   final double temperatureFactor) {
        Locale.setDefault(Locale.US);
        StringBuffer output = new StringBuffer();
        for (int j = atom.getFlag().length(); j < 6; j++) {
            output.append(" ");
        }
        output.append(atom.getFlag());
        for (int j = Integer.toString(atom.getSerialNumber()).length();
             j < 5;
             j++
            ) {
            output.append(" ");
        }
        output.append(atom.getSerialNumber());
        output.append(" ");
        for (int j = atom.getName().length(); j < 4; j++) {
            output.append(" ");
        }
        output.append(atom.getName());
        output.append(atom.getAlternativeLocation());
        for (int j = atom.getResidueName().length(); j < 3; j++) {
            output.append(" ");
        }
        output.append(atom.getResidueName());
        output.append(" ");
        output.append(atom.getChainId());
        for (int j = Integer.toString(atom.getResidueNumber()).length();
             j < 4;
             j++
            ) {
            output.append(" ");
        }
        output.append(atom.getResidueNumber());
        output.append(atom.getICode());
        output.append("   ");
        for (int j = Constants.CARTESIAN_DEC_FORMAT.format(
                                                        atom.getXYZ().getX()
                                                          ).length();
             j < 8;
             j++
            ) {
            output.append(" ");
        }
        output.append(Constants.CARTESIAN_DEC_FORMAT.format(
                                                        atom.getXYZ().getX()
                                                           ));
        for (int j = Constants.CARTESIAN_DEC_FORMAT.format(
                                                        atom.getXYZ().getY()
                                                          ).length();
             j < 8;
             j++
            ) {
            output.append(" ");
        }
        output.append(Constants.CARTESIAN_DEC_FORMAT.format(
                                                        atom.getXYZ().getY()
                                                           ));
        for (int j = Constants.CARTESIAN_DEC_FORMAT.format(
                                                        atom.getXYZ().getZ()
                                                          ).length();
             j < 8;
             j++
            ) {
            output.append(" ");
        }
        output.append(Constants.CARTESIAN_DEC_FORMAT.format(
                                                        atom.getXYZ().getZ()
                                                           ));
        for (int j = Constants.OCCUPANCY_TEMP_DEC_FORMAT.format(
                                                                occupancy
                                                               ).length();
             j < 6;
             j++) {
            output.append(" ");
        }
        output.append(Constants.OCCUPANCY_TEMP_DEC_FORMAT.format(occupancy));

        for (int j = Constants.OCCUPANCY_TEMP_DEC_FORMAT.format(
                                                               temperatureFactor
                                                               ).length();
             j < 6;
             j++) {
            output.append(" ");
        }
        output.append(Constants.OCCUPANCY_TEMP_DEC_FORMAT.format(
                                                               temperatureFactor
                                                                )
                     );
        output.append("         ");
        for (int j = atom.getElement().getSymbol().length(); j < 2; j++) {
            output.append(" ");
        }
        output.append(atom.getElement().getSymbol());
        output.append(" ");
    return output.toString();
    }
}

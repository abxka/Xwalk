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

package mm.electrostatics;

import java.io.File;
import java.io.IOException;

import external.ExternCommand;
import external.PDB2PQR;

import mm.constants.Constants.ForceField;
import structure.constants.Constants;
import structure.io.ReadFile;
import structure.io.WriteFile;
import structure.math.Mathematics;
import structure.math.Point3d;
import structure.math.Point3i;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.hetgroups.SmallMolecule;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;

/**
 * A generic class that holds various methods to operate on classes defined in
 * the matter package.
 * @author Abdullah Kahraman
 * @version 0.5
 * @since 0.5
 */
public class ElectroStaticsUtilities {

    /**
     * Constructor with prevention against calls from subclass.
     */
    protected ElectroStaticsUtilities() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the grid size for the finite difference method in electrostatics
     * potential calculations. Dimension calculations were extracted from
     * the psize.py script within the PDB2PQR software package.
     * @param grid
     *        Point3d object holding in the X,Y,Z coordinates the length of the
     *        grid in X,Y,Z axis direction.
     * @param fineGridSpacing
     *        double value by which the grid sizes will be divided.
     * @return A Point3i object that holds the new grid dimensions
     */
    public static Point3i getGridDim(final Point3d grid,
                                     final double fineGridSpacing) {
        double[] grid3d = grid.convert();
        int[] gridDim = new int[3];
        for (int i = 0; i < 3; i++){
            int tmp = (int) (grid3d[i] / fineGridSpacing + 0.5);
            gridDim[i] = 32 * ((int)((tmp-1) / 32 + 0.5)) + 1;
            if (gridDim[i] < 33){
                gridDim[i] = 33;
            }
        }
        return new Point3i(gridDim[0], gridDim[1], gridDim[2]);
    }
    //--------------------------------------------------------------------------
    /**
     * Protonates a protein complex by running PDB2PQR.
     * @param pdb2pqr
     *        String object holding the path to the pdb2pqr executable.
     * @param proteinComplex
     *        PolyPeptideList object holding the coordinates for the proteins to
     *        be protonated.
     * @param forceField
     *        {@link ForceField} enum type.
     * @param pH
     *        float variable holding the pH value at which the protonation
     *        shall be calculated.
     */
    public final void protonate(final File pdb2pqr,
                                PolyPeptideList proteinComplex,
                                final ForceField forceField,
                                final float pH) {

        // first write out protein coordinates on which PDB2PQR will run
        File pdbFile = new File("proteinBeforePH123.pdb");
        WriteFile write = new WriteFile();
        write.setFile(pdbFile.getAbsolutePath());
        write.write(proteinComplex.toString());

        // set pqr file path
        String pqrFilePath = "proteinPH123.pqr";

        // run PDB2PQR and write out output into file in Trash directory
        String param = "--chain --nodebump --noopt --with-ph=" + pH + "";
        PDB2PQR.run(pdb2pqr, pdbFile, forceField,
                    param, true);

        // read in PQR file
        AtomList pqrAtoms = new AtomList();
        try {
            PQRreader pqrFile = new PQRreader(pqrFilePath);
            pqrAtoms = pqrFile.getAllAtoms();
        } catch (Exception e) {
            System.err.println("ERROR while reading in PQR file "
                             + pqrFilePath + "." + e);
        }

        // assign to protein atoms radii and charges as calculated by PDB2PQR
        AtomList addedAtoms = new AtomList();
        AtomList removedAtoms = new AtomList();
        int i = -1;
        for (PolyPeptide protein : proteinComplex) {
            for (AminoAcid aa : protein) {
                for (Atom pdbAtom : aa.getAllAtoms()) {
                    i++;
                    Atom pqrAtom = pqrAtoms.get(i);
                    // if coordinate is not found, most likely reason is that
                    // new hydrogen was added or old hydrogen was removed.
                    if (!pdbAtom.getXYZ().equals(pqrAtom.getXYZ())) {
                        // If new hydrogen was added, then this coordinate
                        // should not be found in the PDB file.
                        boolean foundPDB = false;
                        for (Atom pdbAtom2 : proteinComplex.getAllAtoms()) {
                            if (pqrAtom.getXYZ().equals(pdbAtom2)) {
                                // however if it is found, this simply means
                                // that PDB2PQR has placed coordinate at a
                                // different location in the file.
                                foundPDB = true;
                                break;
                            }
                        }
                        boolean foundPQR = false;
                        for (Atom pqrAtom2 : pqrAtoms) {
                            if (pdbAtom.getXYZ().equals(pqrAtom2.getXYZ())) {
                                foundPQR = true;
                                break;
                            }
                        }

                        // pqrAtom is new to the file, so it must be a new
                        // hydrogen
                        if (!foundPDB) {
                            aa.add(pqrAtom);
                            addedAtoms.add(pqrAtom);
                        }

                        if (!foundPQR){
                            aa.remove(pdbAtom);
                            removedAtoms.add(pdbAtom);
                        }
                    }
                }
            }
        }

        // Redo serial numbering for hydrogens
        int maxSerial = 0;
        for (Atom atom : proteinComplex.getAllAtoms()) {
            maxSerial = Math.max(atom.getSerialNumber(), maxSerial);
        }
        for (Atom atom : proteinComplex.getAllAtoms().getHydrogenAtoms()) {
            atom.setSerialNumber(maxSerial++);
        }

        if (addedAtoms.size() != 0) {
            System.err.println("Following hydrogen were added to protein "
                             + "structure while pH dependent protontation");
            for (Atom atom : addedAtoms) {
                System.err.print("\t" + atom);
            }
            if (removedAtoms.size() != 0) {
                System.err.println("Following hydrogen were removed from the "
                                 + "protein structure while pH dependent"
                                 + "protontation");
                for (Atom atom : removedAtoms) {
                    System.err.print("\t" + atom);
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Assigns the electrostatic potentials that were calculated using the
     * {@link #run(String, String)} method to an AtomList object. The AtomList
     * object should be ideally those atoms that were also used to focus
     * the electrostatic potential calculations, see {@link #createInputFile(
     * String, String, AtomList, boolean, double, double, double)}
     * @param samplePoints
     *        AtomList object on which the electrostatic potentials will be
     *        mapped.
     * @param dxFile
     *        File object of the APBS .dx output file.
     * @param multivalue
     *        File object of the multivalue application, which is part of the
     *        APBS software package.
     */
    public static void assignElectrostaticPotentials(
                                                  final AtomList samplePoints,
                                                  final File dxFile,
                                                  final File multivalue) {

        File sampleFile = new File("temp2Bdeleted.xyz");
        File outfile = new File("temp2Bdeleted.txt");

        StringBuffer coordiantes = new StringBuffer();
        for (Atom atom : samplePoints) {
            coordiantes.append(atom.getXYZ().getX() + ", "
                             + atom.getXYZ().getY() + ", "
                             + atom.getXYZ().getZ() + Constants.LINE_SEPERATOR);
        }

        WriteFile write = new WriteFile();
        write.setFile(sampleFile.getName());
        write.write(coordiantes.toString());

        // execute the multivalue application
        ExternCommand.execute(multivalue.getAbsolutePath() + " "
                            + sampleFile.getName() + " "
                            + dxFile.getAbsolutePath() + " "
                            + outfile.getAbsolutePath(), true);
        // read in the multivalue outfile and assign the potential values to
        // the sample points
        try {
            ReadFile read = new ReadFile(outfile.getAbsolutePath());
            for (int i = 0; i < samplePoints.size(); i++) {
                float pot = Float.parseFloat(read.get(i).split(",")[3]);

                // potential are given in kT/e. To convert to kcal/mol multiply
                // (kT/e)*0.592 see http://dx.doi.org/10.1006/jmbi.1999.3157
                // caption of Table 2
                pot *= 0.592;

                samplePoints.get(i).setPotential(pot);
            }
        } catch (Exception e) {
            System.err.println("ERROR while reading the multivalue output file "
                             + outfile + " and trying to assign the "
                             + "electrostatic potential values to the sample "
                             + "point coordinates: " + e);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Map the electrostatic potentials from a DotSurface object to a
     * AtomList object. All dot points within 5 Angstroem to an atom are
     * determined and their average potential assigned to a protein atom.
     * @param atomList
     *        AtomList object on which the electrostatic potentials will be
     *        mapped.
     * @param sampleCoords
     *        AtomList object from which the electrostatic potentials will be
     *        extracted.
     */
    public static void mapElectrostaticPotentials(AtomList atomList,
                                                  final AtomList sampleCoords) {
        for (Atom atom : atomList) {
            int n = 0;
            atom.setPotential(0.0f);
            for (Atom dot : sampleCoords) {
                double dist = Mathematics.distance(atom.getXYZ(), dot.getXYZ());
                if (dist < 5.0) {
                    atom.setPotential(atom.getPotential() + dot.getPotential());
                    n++;
                }
            }
            if (n > 0) {
                atom.setPotential(atom.getPotential() / n);
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Assigns partial charges to atoms based on a previous PDB2PQR run, see
     * {@link #run(String, String, ForceField, String, boolean)}.
     * @param complex
     *        Protein complex on which PDB2PQR was executed.
     * @param pqr
     *        File object of the PQR file.
     */
    public static void assignPartialCharges(final PolyPeptideList complex,
                                            final File pqr) {
        try {
            PQRreader pqrReader = new PQRreader(pqr.getAbsolutePath());
            AtomList allAtoms = complex.getAllAtoms();
            for (SmallMolecule hetgroup : complex.getSmallMolecules()) {
                allAtoms.addAll(hetgroup.getAllAtoms());
            }
            for (Atom pdbAtom : allAtoms) {
                for (Atom pqrAtom : pqrReader.getAllAtoms()) {
                    if (MatterUtilities.equalsType(pdbAtom, pqrAtom)
                        &&
                        MatterUtilities.equalsResidue(pdbAtom, pqrAtom)) {
                        pdbAtom.setPartialCharge(pqrAtom.getPartialCharge());
                    }
                }
            }
        } catch (Exception e) {
            System.err.println("ERROR while reading the PQR file "
                             + pqr.getAbsolutePath()
                             + " and trying to extract partial "
                             + "charges: " + e);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Adjusts a PDB2PQR file such that APBS can read it without error.
     * @param pqrFilePath
     *        String holding the absolute path to the protein structure in PQR
     *        format.
     * @param removeHETATM
     *        Switch for removing HETATM groups.
     * @throws IOException if error occurs while reading the PQR file.
     * @return
     *        String object holding the path to the adjusted PQR file.
     */
    public static String adjustPQRfile(final String pqrFilePath,
                                       final boolean removeHETATM)
                                                            throws IOException {
        ReadFile read = new ReadFile(pqrFilePath);
        StringBuffer newFileContent = new StringBuffer();
        for (String line : read) {
            StringBuffer newLine = new StringBuffer();

            if (line.startsWith("HETATM") && removeHETATM) {
                continue;
            }
            if (line.startsWith("ATOM") || line.startsWith("HETATM")) {
                String flag = line.substring(0,6);
                newLine.append(flag);
                String serial = line.substring(6,11);
                newLine.append(serial);
                newLine.append(" ");
                String name = line.substring(12,16);
                name = name.replaceAll(" ", "");
                newLine.append(name);

                char altLoc = line.charAt(16);
                // if(name.charAt(3) != ' ' && altLoc != ' ') newLine+=" ";
                newLine.append(" ");
                //newLine += altLoc;

                String resName = line.substring(17,20);
                if (altLoc != ' ' && resName.charAt(0) != ' ') {
                    newLine.append(" ");
                }
                newLine.append(resName);
                newLine.append(" ");

                char chainId = line.charAt(21);
                newLine.append(chainId);

                String resNo = line.substring(22,26);
                if (chainId != ' ' && resNo.charAt(0) != ' ') {
                    newLine.append(" ");
                }
                newLine.append(resNo);

                char iCode = line.charAt(26);
                if (resNo.charAt(3) != ' ' && iCode != ' ') {
                    newLine.append(" ");
                }
                newLine.append(iCode);
                newLine.append(" ");
                newLine.append(" ");
                newLine.append(" ");
                String x = line.substring(30,38);
                newLine.append(x);

                String y = line.substring(38,46);
                if (x.charAt(7) != ' ' && y.charAt(0) != ' ') {
                    newLine.append(" ");
                }
                newLine.append(y);

                String z = line.substring(46,54);
                if (y.charAt(7) != ' ' && z.charAt(0) != ' ') {
                    newLine.append(" ");
                }
                newLine.append(z);

                String occ = line.substring(54,62);
                if (z.charAt(7) != ' ' && occ.charAt(0) != ' ') {
                    newLine.append(" ");
                }
                newLine.append(occ);

                String temp = line.substring(62,69);
                if (occ.charAt(7) != ' ' && temp.charAt(0) != ' ') {
                    newLine.append(" ");
                }
                newLine.append(temp);

                newLine.append("\n");
                newFileContent.append(newLine);
            } else {
                newFileContent.append(line);
            }
        }
        String newFile = pqrFilePath.substring(0, pqrFilePath.lastIndexOf("."))
                       + "_adjusted.pqr";
        WriteFile write = new WriteFile();
        write.setFile(newFile);
        write.write(newFileContent.toString());
    return newFile;
    }
}

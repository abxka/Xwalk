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

import structure.constants.Constants;
import structure.io.ReadFile;
import structure.io.WriteFile;
import structure.math.Point3f;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.hetgroups.SmallMolecule;
import structure.matter.protein.PolyPeptideList;

/**
 * Class providing functionalities to generate dot surfaces with
 * <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1180969">SIMS</a>.
 * @author Abdullah Kahraman
 * @version 0.5
 * @since 0.5
 */
public class SIMS {

    /**
     * Path to the SIMS executable.
     */
    private static File sims;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param sims
     *        File object of the SIMS executable.
     */
    public SIMS(final File sims) {
        this.sims = sims;
    }
    //--------------------------------------------------------------------------
    /**
     * Creates a parameter file for a SIMS computation. Returns the name of the
     * input file, which is by default input.sims.
     * @return String object holding the name of the SIMS input file.
     */
    public final String createParameterFile() {
        String parameterFile = "input.sims";

        StringBuffer param = new StringBuffer();
        String nL = Constants.LINE_SEPERATOR;
        param.append("!commad file for  SIMS" + nL);
        param.append("$dot_density=" + Constants.DOT_DENSITY + nL);
        param.append("$probe_rad=" + Constants.SOLVENT_RADIUS + nL);
        param.append("$smooth_rad=" + Constants.SMOOTHING_RADIUS + nL);
        param.append("$dot_file\n");
        param.append("$pdb_rext\n");
        param.append("!end\n");

        WriteFile write = new WriteFile();
        write.setFile(parameterFile);
        write.write(param.toString());

        return parameterFile;
    }
    //--------------------------------------------------------------------------
    /**
     * Executes SIMS on the command line.
     * @param complex
     *        PolyPeptideList object holding the coordinates of the protein
     *        complex.
     */
    public final void run(final PolyPeptideList complex) {

        WriteFile write = new WriteFile();
        write.setFile("molec.cor");

        AtomList atoms = complex.getAllAtoms();
        for (SmallMolecule hetgroup : complex.getSmallMolecules()) {
            atoms.addAll(hetgroup.getAllAtoms());
        }
        for (Atom atom : atoms) {
            atom.setOccupancy(atom.getVanDerWaalsRadius());
        }
        write.write(atoms.toString());

        ExternCommand.execute(SIMS.sims.getAbsolutePath(), true);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the dot surface of a previously executed SIMS run. Will read
     * all coordinates from the SIMS_dot.xyz file and convert them to a PDB
     * format.
     * @return AtomList object holding all molecular surface dots.
     */
    public final AtomList getDotSurface() {
        AtomList dotSurface = new AtomList();

        ReadFile simsOutput = null;
        try {
            simsOutput = new ReadFile("SIMS_dot.xyz");
        } catch (Exception e) {
            System.out.println("ERROR while reading in SIMS output file"
                             + " SIMS_dot.txt. " + e);
        }
        for (String line : simsOutput) {
            if (line.startsWith("DOT")) {
                //String[] lineArray = line.split("\\s+");
                //String ename = lineArray[0];
                //int id = Integer.parseInt(lineArray[1]);
                //int iat = Integer.parseInt(lineArray[2]);
                // 1,2,3 - contact, saddle, reentrant concave
                int type = Integer.parseInt(line.substring(15, 19).trim());
                float area = Float.parseFloat(line.substring(19,26).trim());
                float x = Float.parseFloat(line.substring(27,35).trim());
                float y = Float.parseFloat(line.substring(35,43).trim());
                float z = Float.parseFloat(line.substring(43,51).trim());
                //double nvx = Double.parseDouble(lineArray[8]);
                //double nvy = Double.parseDouble(lineArray[9]);
                //double nvz = Double.parseDouble(lineArray[10]);

                Atom dot = new Atom();
                dot.setName("C");
                dot.setResidueName("SUR");
                dot.setResidueNumber(type);
                dot.setChainId('Z');
                dot.setXYZ(new Point3f(x, y, z));
                dot.setOccupancy(area);
                dot.setVanDerWaalsRadius(0.1f);
                dot.setSerialNumber(dotSurface.size() + 1);
                // add only contact dots, so in order to decrease complexity of
                // surface
                if (type == 1) {
                    dotSurface.add(dot);
                }
            }
        }
        return dotSurface;
    }
}

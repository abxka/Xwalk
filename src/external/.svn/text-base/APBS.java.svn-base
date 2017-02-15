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
import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.util.Locale;

import mm.electrostatics.ElectroStaticsUtilities;
import mm.electrostatics.PQRreader;

import structure.constants.Constants;
import structure.io.WriteFile;
import structure.math.Point3d;
import structure.math.Point3f;
import structure.math.Point3i;
import structure.math.pdb.Transformation;
import structure.matter.AtomList;


/**
 * Class providing functionalities to generate APBS parameter files,
 * executing APBS calculations and mapping electrostatic potentials on a list
 * of coordinates.
 * @author Abdullah Kahraman
 * @version 0.5
 * @since 0.5
 */
public class APBS {

    /**
     * Constructor with prevention against calls from subclass.
     */
    protected APBS() {
        throw new UnsupportedOperationException();
    }
   //--------------------------------------------------------------------------
    /**
     * Executes APBS on the command line.
     * @param apbsPath
     *        File object of the APBS application.
     * @param apbsParameterFile
     *        File object of the APBS parameter input file,
     *        see {@link #createInputFile(String, String, AtomList, boolean,
     *        double, double, double)}
     */
    public static final void run(final File apbsPath,
                                 final File apbsParameterFile) {
        ExternCommand.execute(apbsPath + " " + apbsParameterFile, true);
    }
    //--------------------------------------------------------------------------
    /**
     * Creates an APBS parameter file.
     * @param proteinPqr
     *        File object of the PQR file.
     * @param apbsOutputFile
     *        String object holding the name of the APBS file, without the
     *        .dx suffix.
     * @param samplePoints
     *        AtomList object holding coordinates on which the electrostatic
     *        potential calculations will be focused.
     * @param lpbe
     *        Boolean variable switch for solving the linear Poisson-Boltzmann
     *        equation, instead of the nonlinear PB equation.
     * @param pdie
     *        Double variable holding the protein dielectric constant.
     * @param sdie
     *        Double variable holding the solvent dielectric constant.
     * @param temperature
     *        Double variable holding the temperature.
     * @return File object of the parameter file, which is
     *         inferred from the APBS .dx file name, but with the suffix .in.
     */
    public static final File createInputFile(
                                    final File proteinPqr,
                                    final String apbsOutputFile,
                                    final AtomList samplePoints,
                                    final boolean lpbe,
                                    final double pdie,
                                    final double sdie,
                                    final double temperature) {

        double coarseGridFactor = 1.7;
        int fineGridAddition = 30;
        double fineGridSpacing = 1;

        String file = apbsOutputFile.replaceAll(".dx$", "");

        //----------------------------------------------------------------------
        // Read in PQR files
        //----------------------------------------------------------------------
        AtomList protein = new AtomList();
        try {
           protein =
                  new PQRreader(proteinPqr.getAbsolutePath()).getAllAtoms();

        } catch (Exception e) {
            System.err.println("ERROR while reading input file " + e);
        }

        Point3f proteinCentre = Transformation.centerOfGeometry(protein);
        Point3f sampleCentre = Transformation.centerOfGeometry(samplePoints);

        Point3f protDim = Transformation.dimension(protein);
        Point3f sampleDim = Transformation.dimension(samplePoints);


        Point3d protCoarseGrid = new Point3d(
                           protDim.getX() * coarseGridFactor + fineGridAddition,
                           protDim.getY() * coarseGridFactor + fineGridAddition,
                           protDim.getZ() * coarseGridFactor + fineGridAddition
                                  );

        Point3d protFineGrid = new Point3d(protDim.getX() + fineGridAddition,
                                           protDim.getY() + fineGridAddition,
                                           protDim.getZ() + fineGridAddition);

        Point3d sampleFineGrid = new Point3d(
                                           sampleDim.getX() + fineGridAddition,
                                           sampleDim.getY() + fineGridAddition,
                                           sampleDim.getZ() + fineGridAddition);

        Point3i protGridDim = ElectroStaticsUtilities.getGridDim(
                                                               protFineGrid,
                                                               fineGridSpacing);

        Locale.setDefault(Locale.US);
        NumberFormat decFormat = new DecimalFormat("0.000");

        String nL = Constants.LINE_SEPERATOR;
        String fileContent =
             "#" + nL
           + "#" + nL
           + "read" + nL
           + "    mol pqr " + proteinPqr.getAbsolutePath() + "  # read molecule 1" + nL
           + "end" + nL
           + "elec" + nL
           + "    mg-auto" + nL
           + "    dime   " + protGridDim.getI() + " " + protGridDim.getJ() + " " + protGridDim.getK() + "  # number of grid points" + nL
           + "    cglen  " + decFormat.format(protCoarseGrid.getX()) + " " + decFormat.format(protCoarseGrid.getY()) + " " + decFormat.format(protCoarseGrid.getZ()) + " # coarse mesh lengths (A)" + nL
           + "    cgcent " + decFormat.format(proteinCentre.getX()) + " " + decFormat.format(proteinCentre.getY()) + " " + decFormat.format(proteinCentre.getZ()) + " # coarse mesh centre of geometry" + nL
           + "    fglen   " + decFormat.format(sampleFineGrid.getX()) + " " + decFormat.format(sampleFineGrid.getY()) + " " + decFormat.format(sampleFineGrid.getZ()) + " # fine mesh lengths (A)" + nL
           + "    fgcent " + decFormat.format(sampleCentre.getX()) + " " + decFormat.format(sampleCentre.getY()) + " " + decFormat.format(sampleCentre.getZ()) + " # fine grid/ligand centre of geometry" + nL;
           if (lpbe){
               fileContent += "    lpbe                       # solve the linear PBE with lpbe" + nL;
           }
           else{
               fileContent += "    npbe                       # solve the full nonlinear PBE with npbe" + nL;
           }
           fileContent +=
             "    bcfl sdh                    # Boundary condition flag\n"
           + "                                # zero => Zero\n"
           + "                                # sdh => Single DH sphere\n"
           + "                                # mdh => Multiple DH spheres\n"
           + "                                # focus => Focusing\n"
           + "                                #\n"
           + "    ion  1 0.100000 2.000000    # Counterion declaration:\n"
           + "    ion -1 0.100000 2.000000    # ion <charge> <conc (M)> <radius>" + nL
           + "    ion  2 0.000000 2.000000    # ion <charge> <conc (M)> <radius>" + nL
           + "    ion -2 0.000000 2.000000    # ion <charge> <conc (M)> <radius>" + nL
           + "    pdie " + pdie + "               # Solute dielectric: default=4" + nL
           + "    sdie " + sdie + "              # Solvent dielectric: default=78 (vacuum is 1)" + nL
           + "    chgm spl2                   # Charges mapping on grid method" + nL
           + "                                # spl0 is: the charge is mapped onto the nearest-neighbor grid points." + nL
           + "                                # spl2 is: the charge is mapped onto the nearest- and next-nearest-neighbor grid points." + nL
           + "                                # spl4 is: similar to spl2, except the charge/multipole is additionally mapped to include next-next-nearest neighbors." + nL
           + "    mol 1                       # which molecule to use" + nL
           + "    srfm smol                   # Surface calculation method" + nL
           + "                                # mol => Mol surface for epsilon; inflated VdW for kappa; no smoothing" + nL
           + "                                # smol => As 0 with harmoic average smoothing" + nL
           + "                                # spl2 => Cubic spline " + nL
           + "    srad 1.400000               # Solvent radius (1.4 for water)" + nL
           + "    swin 0.3                    # Surface cubic spline window .. default 0.3" + nL
           + "    temp " + temperature +  "             # System temperature (298.15 default)" + nL
           + "    sdens 10.000000             # Specify the number of grid points per square-angstrom to use in Vacc object. Ignored when srad is 0.0 (see srad) or srfm is spl2 (see srfm). There is a direct correlation between the value used for the Vacc sphere density, the accuracy of the Vacc object, and the APBS calculation time. APBS default value is 10.0." + nL
           + "    gamma 0.105                 # Surface tension parameter for apolar forces (in kJ/mol/A^2)" + nL
           + "                                # only used for force calculations, so we don't care, but" + nL
           + "                                # it's always required, and 0.105 is the default" + nL
           + "#    calcenergy no               # Energy I/O to stdout" + nL
           + "#                                #  no        => don't write out energy" + nL
           + "#                                #  total     => write out total energy" + nL
           + "#                                #  comps     => write out total energy and all components" + nL
           + "#    calcforce no                # Atomic forces I/O (to stdout)" + nL
           + "#                                #  no        => don't write out forces" + nL
           + "#                                #  total     => write out net forces on molecule" + nL
           + "#                                #  comps     => write out atom-level forces" + nL
           + "    write pot dx " + file + "     # What to write .. this says write the potential in dx" + nL
           + "                                 # format to a file." + nL
           + "end" + nL
           + "quit" + nL;

           WriteFile write = new WriteFile();
           write.setFile(file + ".in");
           write.write(fileContent);

       return new File(file + ".in");
    }

    //-------------------------------------------------------------------------
    /**
     * Main method. Writes an APBS parameter file given a PQR file of a protein
     * and a PDB file with a list of atom coordinates on which the electrostatic
     * potentials will be calculated.
     * @param args
     *        Array of Strings holding all commandline arguments.
     */
    public static void main(final String[] args) {
        if (args.length == 0
            ||
            args[0].equals("-help")
            ||
            args[0].equals("-h")) {
            System.err.println("\nUsage: java " + APBS.class.getName()
                             + " protein.pqr sample.pdb\n");
            System.exit(0);
        }

//        APBS.createInputFile(args[0], new DotSurface(), true, 4.0, 78.0, 298.15);
    }
}

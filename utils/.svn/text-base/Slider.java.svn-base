/*
 * (C) 2011 Abdullah Kahraman
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
import java.util.Hashtable;

import structure.constants.Constants;
import structure.io.Commandline;
import structure.io.ReadFile;
import structure.io.WriteFile;
import structure.io.pdb.PDBreader;
import structure.math.Mathematics;
import structure.math.Point3f;
import structure.math.pdb.Transformation;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;
import xwalk.crosslink.CrossLinkList;
import xwalk.crosslink.CrossLinkUtilities;
import xwalk.io.DistanceReader;

/**
 * Class holding a main method to slide two protein towards each other given
 * a list of distance constraints between both proteins.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Slider {
    /**
     * Empty Constructor.
     */
    protected Slider() {
    }

    //--------------------------------------------------------------------------

    /**
     * Path to the first PDB file.
     */
    private String refFilePath;
    /**
     * Path to the second PDB file.
     */
    private String mobFilePath;
    /**
     * Path to the distance file that holds the distance constraints.
     */
    private String distFilePath;
    /**
     * verbose information to print out.
     */
    private boolean verbose = false;
    /**
     * maximum temperature in simulated annealing.
     */
    private double maxTemperature = 10;
    /**
     * maximum MC iteration cycles.
     */
    private int maxIterationCycles = 1000;
    /**
     * last accepted distance sum.
     */
    private double lastAcceptedDistanceSum = Double.MAX_VALUE;
    /**
     * lowest accepted distance sum.
     */
    private double lowestDistanceSum = Double.MAX_VALUE;
    /**
     * conformation with lowest distance sum.
     */
    private PolyPeptideList proteinMobLowest = null;

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
            System.out.println(nL + nL
                            + "java " + Slider.class.getName() + " -help"
                            + nL);
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-help", false).equals("EXISTS")) {
            System.err.println(nL
                           + "\"Slider\" slides two proteins towards "
                           + "each other while minimizing a list of "
                           + "distance constraints."
                           + nL
                           + nL
                           + "Usage:"
                           + nL
                           + "java " + Slider.class.getName()
                           + " -ref 1brsA.pdb -mob 1brsD.pdb -dst 1brs.dist"
                           + nL
                           + nL
                           + "Parameters:"
                           + nL
                           + "\t-ref\t<path>\tFirst protein, which coordinates "
                           + "will be kept fiexed (required)."
                           + nL
                           + "\t-mob\t<path>\tSecond protein, which "
                           + "coordinates will be moved to fullfill the "
                           + " distance constraints (required)."
                           + nL
                           + "\t-dist\t<path>\tDistance file holding at least "
                           + "the first 4 columns of the Xwalk output format. "
                           + "The file will be used to extract the indices and "
                           + "the residue pairs for the distance calculation "
                           + " (required)"
                           + nL
                           + "\t-temp\t[double]\tMaximum temperature used in "
                           + "simulated annealing (default: " + maxTemperature
                           + ", optional)"
                           + nL
                           + "\t-cycles\t[int]\tNumber of Monte Carlo cycles "
                           + "(default: " + maxIterationCycles + ", optional)"
                           + nL
                           + "\t-v\t[switch]\tPrint out verbose information on "
                           + "the STDERR channel (optional)"
                           + nL
                           + nL
                    );
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-ref", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-ref\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            this.refFilePath = Commandline.get(args, "-ref", true);

            if (!ReadFile.exists(this.refFilePath)) {
                System.err.print(nL
                              + "Couldn't open file \"" + this.refFilePath
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (Commandline.get(args, "-mob", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-mob\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            this.mobFilePath = Commandline.get(args, "-mob", true);

            if (!ReadFile.exists(this.mobFilePath)) {
                System.err.print(nL
                              + "Couldn't open file \"" + this.mobFilePath
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (Commandline.get(args, "-dist", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-dist\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            this.distFilePath = Commandline.get(args, "-dist", true);

            if (!ReadFile.exists(this.distFilePath)) {
                System.err.print(nL
                              + "Couldn't open file \"" + this.distFilePath
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (Commandline.get(args, "-v", false).equals("EXISTS")) {
            this.verbose = true;
        }
        //----------------------------------------------------------------------
        if (!Commandline.get(args, "-temp", true).equals("ERROR")) {
            this.maxTemperature = Double.parseDouble(
                                     Commandline.get(args, "-temp", true)
                                                 );
        }
        //----------------------------------------------------------------------
        if (!Commandline.get(args, "-cycles", true).equals("ERROR")) {
            this.maxIterationCycles = Integer.parseInt(
                                     Commandline.get(args, "-cycles", true)
                                                      );
        }
    }
    //--------------------------------------------------------------------------
    private static Point3f getRandomTranslationVector() {
        Point3f translationVector = new Point3f(
                                             (float) (1 - (Math.random() * 1)),
                                             (float) (1 - (Math.random() * 2)),
                                             (float) (1 - (Math.random() * 2)));
        return translationVector;
    }
    //--------------------------------------------------------------------------
    private static double[][] getRandomRatationMatrix(final boolean verbose){
        // rotate around each axis by +/- 5 degree.
        double phi = (Math.PI / 18) - ((Math.PI / 9) * Math.random());
        double theta = (Math.PI / 18) - ((Math.PI / 9) * Math.random());
        double teta = (Math.PI / 18) - ((Math.PI / 9) * Math.random());
        if (verbose) {
            System.err.println("ROTATION: "
            + "PHI="
            + Constants.CARTESIAN_DEC_FORMAT.format(phi)
            + ", "
            + "THETA="
            + Constants.CARTESIAN_DEC_FORMAT.format(theta)
            + ", "
            + "TETA="
            + Constants.CARTESIAN_DEC_FORMAT.format(teta)
                              );
         }

        double[][] rotationMatrix = Mathematics.getEulerRotationMatrix(phi,
                                                                       theta,
                                                                       teta);
        return rotationMatrix;
    }

    //--------------------------------------------------------------------------
    private static double getAvgerageDistance(
                                          final PolyPeptideList proteinRef,
                                          final PolyPeptideList proteinMobCopy,
                                          final CrossLinkList constraintsList) {
        // calculate distance sum
        PolyPeptideList complexCopy = new PolyPeptideList();
        complexCopy.addAll(proteinRef);
        complexCopy.addAll(proteinMobCopy);
        Hashtable < Atom, AtomList > relevantAtomPairs = null;
        try {
            relevantAtomPairs = CrossLinkUtilities.extractRelevantPairs(
                                                             complexCopy,
                                                             constraintsList
                                                                       );
        } catch (IOException e) {
            System.err.println("ERROR while attempting to extract "
                             + "cross-linked residues from the input PDB "
                             + "files: " + e);
            System.exit(1);
        }

        // redundant cross-links are those that are found between
        // alternative locations, in which case the shortest cross-link
        // should only be considered.
        Hashtable<String, Double> redundant = new Hashtable<String, Double>();
        for (Atom atom1 : relevantAtomPairs.keySet()) {
            String atom1Id = atom1.getResidueName()
                           + atom1.getResidueNumber()
                           + atom1.getChainId();
            for (Atom atom2 : relevantAtomPairs.get(atom1)) {
                String atom2Id = atom2.getResidueName()
                               + atom2.getResidueNumber()
                               + atom2.getChainId();
                double currentDist = Mathematics.distance(atom1.getXYZ(),
                                                          atom2.getXYZ());
                if (redundant.get(atom1Id + atom2Id) != null) {
                    double prevDist = redundant.get(atom1Id + atom2Id);
                    if (currentDist < prevDist) {
                        redundant.put(atom1Id + atom2Id, currentDist);
                    }
                } else if (redundant.get(atom2Id + atom1Id) != null) {
                    double prevDist = redundant.get(atom2Id + atom1Id);
                    if (currentDist < prevDist) {
                        redundant.put(atom2Id + atom1Id, currentDist);
                    }
                } else {
                    redundant.put(atom1Id + atom2Id, currentDist);
                }
            }
        }
        double avgDist = 0;
        int n = 0;
//        for (String id : redundant.keySet()) {
//            double dist = redundant.get(id);
//            if (dist > 30) {
//                avgDist += dist;
//                n++;
//            }
//        }
//        if (n==0) {
            for (String id : redundant.keySet()) {
                double dist = redundant.get(id);
                avgDist += dist;
                n++;
            }
            avgDist /= n;
  //      }

    return avgDist;
    }

    //--------------------------------------------------------------------------

    private static double getBoltzmannProbability(
                                           final double lastAcceptedDistanceSum,
                                           final double distSum,
                                           final double temperature) {
        // check whether distance is smaller, in which case do the move
        // on the original mobile protein
        double boltzmannFactor = (lastAcceptedDistanceSum - distSum)
                                 /
                                 temperature;
        double probability = Math.exp(boltzmannFactor);
    return probability;
    }

    //--------------------------------------------------------------------------
    private static boolean doTransformation(
                                           final double lastAcceptedDistanceSum,
                                           final double distSum,
                                           final double temperature,
                                           final boolean verbose) {

        double probability = Slider.getBoltzmannProbability(
                                                        lastAcceptedDistanceSum,
                                                        distSum,
                                                        temperature);

        boolean doMove = false;
        if (probability >= 1) {
            doMove = true;
            if (verbose) {
                System.err.println("ACCEPT: Better score: "
                          + Constants.CARTESIAN_DEC_FORMAT.format(distSum));
            }
        } else if (probability >= Math.random()) {
            doMove = true;
            if (verbose) {
                System.err.println("ACCEPT: By thermal probability: "
                        + Constants.CARTESIAN_DEC_FORMAT.format(probability)
                        + "\t"
                        + Constants.CARTESIAN_DEC_FORMAT.format(distSum));
            }
        } else {
            doMove = false;
            if (verbose) {
                System.err.println("FAILED: "
                        + Constants.CARTESIAN_DEC_FORMAT.format(probability)
                        + "\t"
                        + Constants.CARTESIAN_DEC_FORMAT.format(distSum));
            }
        }
        return doMove;
    }

    //--------------------------------------------------------------------------
    private PolyPeptideList doTranslation(final PolyPeptideList proteinRef,
                                          final PolyPeptideList proteinMob,
                                          final CrossLinkList constraintsList) {
    
        Point3f translationVector = new Point3f(-0.1f, 0.0f, 0.0f);
        if (this.verbose) {
            System.err.println("TRANSLATION: "
            + "X="
            + Constants.CARTESIAN_DEC_FORMAT.format(translationVector.getX())
            + ", "
            + "Y="
            + Constants.CARTESIAN_DEC_FORMAT.format(translationVector.getY())
            + ", "
            + "Z="
            + Constants.CARTESIAN_DEC_FORMAT.format(translationVector.getZ())
                              );
         }
        // create copy of mobile protein to test move
        PolyPeptideList proteinMobCopy = new PolyPeptideList();
        for (PolyPeptide protein : proteinMob) {
            proteinMobCopy.add(protein.copy());
        }

        PolyPeptideList proteinMobMin = new PolyPeptideList();
        double min = Double.MAX_VALUE;
        int cycles = -1;
        while (cycles++ < this.maxIterationCycles) {
            Transformation.move(proteinMobCopy.getAllAtoms(),
                                translationVector);

            double distSum = Slider.getAvgerageDistance(proteinRef,
                                                        proteinMobCopy,
                                                        constraintsList);

            if (distSum < min) {
                if (verbose) {
                    System.err.println("ACCEPT: Shorter distance at "
                           + Constants.CARTESIAN_DEC_FORMAT.format(cycles * 0.1)
                           + " A x-axis translation: "
                           + Constants.CARTESIAN_DEC_FORMAT.format(distSum));
                }

                min = distSum;
                // create copy of mobile protein to test move
                proteinMobMin.clear();
                for (PolyPeptide protein : proteinMobCopy) {
                    proteinMobMin.add(protein.copy());
                }
            }
            else {
                if (verbose) {
                    System.err.println("FAILED: Longer distance at "
                           + Constants.CARTESIAN_DEC_FORMAT.format(cycles * 0.1)
                           + " A x-axis translation: "
                           + Constants.CARTESIAN_DEC_FORMAT.format(distSum)
                           + " vs "
                           + Constants.CARTESIAN_DEC_FORMAT.format(min));
                }
            }
        }
        this.lowestDistanceSum = min;
        return proteinMobMin;
    }
    //--------------------------------------------------------------------------
    private void doMC(final PolyPeptideList proteinRef,
                      final PolyPeptideList proteinMob,
                      final CrossLinkList constraintsList,
                      int maxIterationCycles,
                      final double temperature) {
        // continue move attempts for a maximum number of move attempts or
        // until move resulted in a sufficient conformation
        Point3f translationVector = new Point3f(0.1f, 0.0f, 0.0f);
        if (this.verbose) {
            System.err.println("TRANSLATION: "
            + "X="
            + Constants.CARTESIAN_DEC_FORMAT.format(translationVector.getX())
            + ", "
            + "Y="
            + Constants.CARTESIAN_DEC_FORMAT.format(translationVector.getY())
            + ", "
            + "Z="
            + Constants.CARTESIAN_DEC_FORMAT.format(translationVector.getZ())
                              );
         }
        // create copy of mobile protein to test move
        PolyPeptideList proteinMobCopy = new PolyPeptideList();
        for (PolyPeptide protein : proteinMob) {
            proteinMobCopy.add(protein.copy());
        }
        while (maxIterationCycles-- > 0) {
            //----------------------RANDOM TRANSLATION--------------------------
            // generate a random translation vector with -1 to 1 coordinate
            // values.
//            Point3f translationVector = Slider.getRandomTranslationVector();
            // set y and z to 0.0 to allow translation only in x direction
//            translationVector = translationVector.add(0.0f,
//                                                     -translationVector.getY(),
//                                                    -translationVector.getZ());


            // do the random move
            Transformation.move(proteinMobCopy.getAllAtoms(),
                                translationVector);

            double distSum = Slider.getAvgerageDistance(proteinRef,
                                                        proteinMobCopy,
                                                        constraintsList);


            boolean doMove = Slider.doTransformation(
                                                   this.lastAcceptedDistanceSum,
                                                   distSum,
                                                   temperature / 100,
                                                   this.verbose);

            if (doMove) {
                Transformation.move(proteinMob.getAllAtoms(),
                                    translationVector);
                this.lastAcceptedDistanceSum = distSum;
                this.proteinMobLowest = proteinMobCopy;

                if (distSum < this.lowestDistanceSum) {
                    this.proteinMobLowest = proteinMobCopy;
                    this.lowestDistanceSum = distSum;
                }
            }

            //----------------------RANDOM ROTATION-----------------------------
            // do the random move
            // do rotation with Euler angles for which we first need to
            // translate the protein to the coordinate center.
/*            double[][] rotationMatrix = Mathematics.getEulerRotationMatrix(0, Math.PI/180, 0);
                    //Slider.getRandomRatationMatrix(
                    //                                                this.verbose
                    //                                                  );

            // create copy of mobile protein to test move
            proteinMobCopy = null;
            proteinMobCopy = new PolyPeptideList();
            for (PolyPeptide protein : proteinMob) {
                proteinMobCopy.add(protein.copy());
            }

            Transformation.rotateAtOrigin(proteinMobCopy.getAllAtoms(),
                                          rotationMatrix);

            distSum = Slider.getAvgerageDistance(proteinRef,
                                        proteinMobCopy,
                                        constraintsList);

            boolean doRotation = Slider.doTransformation(
                                                   this.lastAcceptedDistanceSum,
                                                   distSum,
                                                   temperature,
                                                   this.verbose);
            if (doRotation) {
                Transformation.rotateAtOrigin(proteinMob.getAllAtoms(),
                                              rotationMatrix);
                this.lastAcceptedDistanceSum = distSum;
            }

            if (distSum < this.lowestDistanceSum) {
                this.proteinMobLowest = proteinMobCopy;
                this.lowestDistanceSum = distSum;
            }
*/        }
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
        Slider slider = new Slider();
        slider.readCommandline(args);

        //----------------------------------------------------------------------
        // Read in PDB files
        //----------------------------------------------------------------------
        PDBreader readersRef = null;
        PDBreader readersMob = null;
        CrossLinkList constraintsList = null;
        try {
            readersRef = PDBreader.createPDBreaders(slider.refFilePath).get(0);
            readersMob = PDBreader.createPDBreaders(slider.mobFilePath).get(0);
            constraintsList = DistanceReader.getCrossLinks(slider.distFilePath,
                                                           false,
                                                           true,
                                                           true);
        } catch (Exception e) {
            System.err.println("ERROR while reading in input files: " + e);
            System.exit(1);
        }

        PolyPeptideList proteinRef =
                                    readersRef.getEntireProteinComplex().get(0);
        PolyPeptideList proteinMob =
                                    readersMob.getEntireProteinComplex().get(0);

        Point3f translationVector = new Point3f(-1.0f, 0.0f, 0.0f);
        if (slider.verbose) {
            System.err.println("TRANSLATION: "
            + "X="
            + Constants.CARTESIAN_DEC_FORMAT.format(translationVector.getX())
            + ", "
            + "Y="
            + Constants.CARTESIAN_DEC_FORMAT.format(translationVector.getY())
            + ", "
            + "Z="
            + Constants.CARTESIAN_DEC_FORMAT.format(translationVector.getZ())
                              );
         }
        double phi = 0;
        double theta = 1;
        double teta = 0;
        double[][] rotationMatrix = Mathematics.getEulerRotationMatrix(
                                                                  phi,
                                                                  theta,
                                                                  teta);

        for (int i = 0; i<=50 ; i++) {

            // create copy of mobile protein to test move
            PolyPeptideList proteinMobCopy = new PolyPeptideList();
            for (PolyPeptide protein : proteinMob) {
                proteinMobCopy.add(protein.copy());
            }
            Transformation.move(proteinMobCopy.getAllAtoms(),
                                translationVector.add((float)i, 0.0f, 0.0f));

            double min = Double.MAX_VALUE;
            PolyPeptideList proteinMobMin = new PolyPeptideList();
            for (double j = 0; j <= 2 * Math.PI; j += Math.PI / 180) {
            
                Transformation.rotateCenterOfMassAtOrigin(
                                                   proteinMobCopy.getAllAtoms(),
                                                   rotationMatrix
                                                     );
                double distSum = Slider.getAvgerageDistance(proteinRef,
                                                            proteinMobCopy,
                                                            constraintsList);

                if (distSum < min) {
                    if (slider.verbose) {
                        System.err.println("ACCEPT: Shorter distance at " + i
                                         + " A x-axis translation; "
                                         + "ROTATION: "
                                         + "PHI="
                                         + Constants.CARTESIAN_DEC_FORMAT.format(phi)
                                         + ", "
                                         + "THETA="
                                         + Constants.CARTESIAN_DEC_FORMAT.format(theta)
                                         + ", "
                                         + "TETA="
                                         + Constants.CARTESIAN_DEC_FORMAT.format(teta)
                                         + "\t"
                                         + Constants.CARTESIAN_DEC_FORMAT.format(distSum));
                    }
                    min = distSum;
                    // create copy of mobile protein to test move
                    proteinMobMin.clear();
                    for (PolyPeptide protein : proteinMobCopy) {
                        proteinMobMin.add(protein.copy());
                    }
                }
                else {
                    if (slider.verbose) {
                        System.err.println("FAILED: Longer distance at " + i
                                + " A x-axis translation; "
                                + "ROTATION: "
                                + "PHI="
                                + Constants.CARTESIAN_DEC_FORMAT.format(phi)
                                + ", "
                                + "THETA="
                                + Constants.CARTESIAN_DEC_FORMAT.format(theta)
                                + ", "
                                + "TETA="
                                + Constants.CARTESIAN_DEC_FORMAT.format(teta)
                                + "\t"
                                + Constants.CARTESIAN_DEC_FORMAT.format(distSum));
                    }
                }
            }
            WriteFile write = new WriteFile();
            write.setFile(i + "th.pdb");
            write.write("REMARK Distance sum is "
                    + Constants.CARTESIAN_DEC_FORMAT.format(min)
                    + Constants.LINE_SEPERATOR
                    + proteinMobMin + "\n");
        }
    }
}

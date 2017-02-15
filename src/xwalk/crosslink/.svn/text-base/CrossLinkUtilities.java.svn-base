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

package xwalk.crosslink;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.TreeSet;
import java.util.zip.DataFormatException;

import structure.constants.Constants;
import structure.exceptions.FileFormatException;
import structure.grid.AtomGrid;
import structure.grid.GridCell;
import structure.grid.GridUtilities;
import structure.grid.Path;
import structure.io.pdb.PDBreader;
import structure.math.Mathematics;
import structure.math.algorithms.BreadthFirstSearch;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.parameter.AtomType;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.Digestion;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;

import xwalk.io.DistanceReader;
import xwalk.math.SolventPathDistance;
import xwalk.crosslink.CrossLinkParameter.Parameter;

/**
 * Class that holds all relevant methods to calculate virtual cross-links on PDB
 * protein structures.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public final class CrossLinkUtilities {

    /**
     * Constructor.
     */
    protected CrossLinkUtilities() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a list of potential virtual cross-links.
     * @param complexes
     *        - List of protein complex objects holding all protein complex
     *          molecules for which virtual cross-links should be calculated.
     *        necessary for the virtual cross-link calculation.
     * @throws IOException if an error occurred while reading the infile.
     * @return CrossLinkList object that holds all virtual cross-links on a
     *         protein complex.
     */
    public static CrossLinkList getVirtualCrossLinks(
                                  final ArrayList < PolyPeptideList > complexes
                                                    )
                                                  throws IOException {
        //--------------------------------
        // read in cross-links from distance file
        CrossLinkList distXlList = null;
        if (!CrossLinkParameter.getParameter(
                                     Parameter.DISTANCE_FILE_PATH).equals("")) {
            String fileName = CrossLinkParameter.getParameter(
                                                    Parameter.DISTANCE_FILE_PATH
                                                    );
            boolean onlyIntra = Boolean.parseBoolean(
                                       CrossLinkParameter.getParameter(
                                            Parameter.DO_INTRAMOLECULAR_DISTANCE
                                                                      )
                                                    );
            boolean onlyInter = Boolean.parseBoolean(
                                CrossLinkParameter.getParameter(
                                            Parameter.DO_INTERMOLECULAR_DISTANCE
                                                               )
                                                    );
            distXlList = DistanceReader.getCrossLinks(fileName,
                                                      onlyIntra,
                                                      onlyInter,
                                                      true);
            if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                      Parameter.DO_PROBABILITY
                                                                    ))) {
                for (CrossLink xl : distXlList) {
                    xl.setEucProbability();
                    xl.setSASDprobability();
                }
            }
            if (distXlList.size() == 0) {
                System.err.println("WARNING: No suitable cross-links found in"
                                 + " the distance file \"" + fileName + "\"");
                return new CrossLinkList();
            }
        }
        //--------------------------------
        // digest protein
        ArrayList < PolyPeptideList > allDigest = null;
        if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                 Parameter.DO_TRYPSIN_DIGEST
                                                       )
                                )
        ) {
            allDigest = CrossLinkUtilities.digestProteinComplex(complexes);
        }

        //--------------------------------
        // calculate distances

        CrossLinkList allCrossLinkList = new CrossLinkList();
        // find and create virtual cross-links on the protein complexes.
        for (int i = 0; i < complexes.size(); i++) {
            PolyPeptideList complex = complexes.get(i);
            PolyPeptideList digest = null;

            if (allDigest != null) {
                digest = allDigest.get(i);
            }

            //---------------------------------
            // First find cross-links based on Euclidean distance.
            CrossLinkList crossLinkList =
                            CrossLinkUtilities.crossLinkByEuclideanDistance(
                                                                     complex,
                                                                     distXlList,
                                                                     digest
                                                                           );
            //---------------------------------
            // remove redundant cross-links if the complex should be labeled by
            // the user as homomeric.
            if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                          Parameter.IS_HOMOMERIC
                                                           ))) {
                CrossLinkUtilities.removeRedundanciesInHomomers(crossLinkList);
            }

            //---------------------------------
            // If requested by the user check further for Solvent Path distance.
            if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                           ))) {
                if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                               ))) {
                    System.err.print("Checking \"" + crossLinkList.size() + "\""
                                   + " Euclidean distances for becoming solvent"
                                   + " accessible surface distances.\n");
                }
                CrossLinkUtilities.calculatesSolventPathDistance(complex,
                                                                 crossLinkList
                                                                );
            }

            // if no distance file has been given, than output only
            // cross-links that can be conform, distance-wise and SAS-wise.
            if (distXlList == null) {
                CrossLinkUtilities.cleanCrossLinkList(crossLinkList);
            }

            //---------------------------------
            // set indices of cross-links
            CrossLinkUtilities.setCrossLinkIndicesAndFileName(crossLinkList);
            allCrossLinkList.addAll(crossLinkList);
        }

        return allCrossLinkList;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a list of potential virtual mono-links.
     * @param complexes
     *        - List of protein complex objects holding all protein complex
     *          molecules for which virtual cross-links should be calculated.
     * @throws IOException if an error occurred while reading the infile.
     * @return MonoLinkList object that holds all virtual mono-links on a
     *         protein complex.
     */
    public static MonoLinkList getVirtualMonoLinks(
                                  final ArrayList < PolyPeptideList > complexes
                                                  )
                                                  throws IOException {

        // read in mono-links from distance file
        MonoLinkList distMlList = null;
        if (!CrossLinkParameter.getParameter(
                                    Parameter.DISTANCE_FILE_PATH).equals("")) {
            String fileName = CrossLinkParameter.getParameter(
                                                    Parameter.DISTANCE_FILE_PATH
                                                    );

            distMlList = DistanceReader.getMonoLinks(fileName);
        }

        MonoLinkList allMonoLinkList = new MonoLinkList();
        // find and create virtual mono-links on the protein complexes.
        for (int i = 0; i < complexes.size(); i++) {
            PolyPeptideList complex = complexes.get(i);

            float gridCellSize = Float.parseFloat(
                                            CrossLinkParameter.getParameter(
                                                        Parameter.GRID_CELL_SIZE
                                                                           ));
            // assume that cross-linker requires around 10 Angstroem of space.
            // This number is purely random.
            float maxDist = xwalk.constants.Constants.CROSS_LINKER_END_SIZE;

            MonoLinkList monoLinkList = new MonoLinkList();

            if (distMlList != null) {
                // clone distMlList
                MonoLinkList distMlListClone = new MonoLinkList();
                for (MonoLink monoLink : distMlList) {
                    distMlListClone.add(monoLink.copy());
                }

                for (MonoLink monoLink : distMlListClone) {
                    ArrayList<AtomList> candidates =
                        CrossLinkUtilities.proteinAtomsMatchingXLatom(
                                                                    complex,
                                                                    monoLink
                                                                     );
                    for (ArrayList<Atom> monoLinkAtoms : candidates) {
                        for (Atom monoLinkAtom : monoLinkAtoms) {
                            // if a candidate has already been found to be
                            // solvent accessible than stop here and continue
                            // with next mono link.
                            if (!monoLink.isSolventAccessible()) {
                                monoLink.set(monoLinkAtom);
                                AtomGrid grid = new AtomGrid(
                                                          complex.getAllAtoms(),
                                                          monoLink,
                                                          maxDist,
                                                          gridCellSize);

                                monoLink.setSolventAccessibility(
                                      GridUtilities.isAccessible(monoLink,
                                                                 grid));
                            }
                            if (monoLinkList.get(monoLink) == null) {
                                if (!Boolean.parseBoolean(
                                        CrossLinkParameter.getParameter(
                                                      Parameter.DO_KEEP_FILENAME
                                                                       ))) {
                                    monoLink.setFileName(complex.getName());
                                }
                                monoLinkList.add(monoLink);
                            }
                        }
                    }
                }
            } else {
                TreeSet < AtomList > relevantAtoms1 =
                                       CrossLinkUtilities.findAllRelevantAtoms1(
                                                                       complex);
                TreeSet < AtomList > relevantAtoms2 =
                                       CrossLinkUtilities.findAllRelevantAtoms2(
                                                                       complex);

                // Create TreeSet for unique list of relevant atoms.
                HashSet < Atom > relevantAtoms = new HashSet < Atom >();

                for (AtomList list : relevantAtoms1) {
                    relevantAtoms.addAll(list);
                }
                for (AtomList list : relevantAtoms2) {
                    relevantAtoms.addAll(list);
                }

                for (Atom atom : relevantAtoms) {
                    MonoLink monoLink = new MonoLink();
                    monoLink.set(atom);

                    AtomGrid grid = new AtomGrid(complex.getAllAtoms(),
                                                 monoLink,
                                                 maxDist,
                                                 gridCellSize);
                    monoLink.setSolventAccessibility(
                                      GridUtilities.isAccessible(monoLink,
                                                                 grid));
                    if (monoLink.isSolventAccessible()) {
                        monoLink.setFileName(complex.getName());
                        monoLinkList.add(monoLink);
                    }
                }
                monoLinkList.setIndices();
            }
            allMonoLinkList.addAll(monoLinkList);
        }
    return allMonoLinkList;
    }
    //--------------------------------------------------------------------------
    /**
     * Digest a list of protein complex objects.
     * @param complexes
     *        List of Protein complex object to be all digested.
     * @return List of PolyPeptides that are formed by digestion.
     */
    public static ArrayList <PolyPeptideList> digestProteinComplex(
            final ArrayList < PolyPeptideList > complexes
        ) {
        ArrayList < PolyPeptideList > allDigest =
                                            new ArrayList < PolyPeptideList >();
        for (PolyPeptideList complex : complexes) {

            PolyPeptideList digest = CrossLinkUtilities.trypsinate(complex,
                                          Boolean.parseBoolean(
                                             CrossLinkParameter.getParameter(
                                                        Parameter.DO_EXPASY_RULE
                                                                            )
                                                              ));
            allDigest.add(digest);
            if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                           ))) {
                // output digested peptides
                for (int i = 0; i < digest.size(); i++) {
                    System.err.println(i + 1 + ". "
                                     + digest.get(i).toStringOneLetterCode());
                }
            }
        }
        return allDigest;
    }

    //--------------------------------------------------------------------------
    /**
     * CrossLinks all atoms in a protein complex that either have an Euclidean
     * distance smaller than a user set -max value or if a distance file is
     * given calculates the Euclidean distance between all amino acid pairs and
     * sets SASD = -2 of those that have a Euclidean distance smaller then
     * maxDist, otherwise leaves SASD = -1.
     * @param complex
     *        Protein complex object.
     * @param distanceFileCrossLinks
     *        List of CrossLink objects extracted from a distance file. If no
     *        distance file has been read, just submit a {@code NULL}.
     * @param digest
     *        - List of PolyPeptides that are formed by digestion. If no
     *          digested peptides exist, than submit {@code NULL}.
     * @return List of CrossLink object that all have a Euclidean distance
     *         smaller then the user set maxDist. If a distance file is given
     *         than also CrossLink objects are created for amino acid pairs that
     *         exceed the maxDist, which however have an SASD == -2.0.
     * @throws IOException if an error occurred while reading the distance file.
     */
    public static CrossLinkList crossLinkByEuclideanDistance(
                                     final PolyPeptideList complex,
                                     final CrossLinkList distanceFileCrossLinks,
                                     final PolyPeptideList digest)
                                            throws IOException {
        Hashtable < Atom, AtomList > relevantAtomPairs =
            new Hashtable < Atom, AtomList >();

        if (distanceFileCrossLinks == null) {
            // find all atom pairs in the complex that fulfill the user set
            // criteria and have a Euclidean distance < max.
            relevantAtomPairs = CrossLinkUtilities.findRelevantPairs(complex);
        } else {
            relevantAtomPairs = CrossLinkUtilities.extractRelevantPairs(
                                                        complex,
                                                        distanceFileCrossLinks);
        }


        // create CrossLinks object from all relevant atom pairs.
        CrossLinkList crossLinks = new CrossLinkList();
        for (Atom atom1 : relevantAtomPairs.keySet()) {
            PolyPeptide atom1TrypticPeptide = null;
            for (Atom atom2 : relevantAtomPairs.get(atom1)) {

                boolean onlyIntra = Boolean.parseBoolean(
                                      CrossLinkParameter.getParameter(
                                         Parameter.DO_INTRAMOLECULAR_DISTANCE));
                boolean onlyInter = Boolean.parseBoolean(
                                      CrossLinkParameter.getParameter(
                                         Parameter.DO_INTERMOLECULAR_DISTANCE));
                if (onlyIntra) {
                    if (atom1.getChainId() != atom2.getChainId()) {
                        continue;
                    }
                }
                if (onlyInter) {
                    if (atom1.getChainId() == atom2.getChainId()) {
                        continue;
                    }
                }
                PolyPeptide atom2TrypticPeptide = null;
                boolean conforming = false;

                float errorRange = 0;
                if (Boolean.parseBoolean(
                       CrossLinkParameter.getParameter(Parameter.DO_BFACTOR))) {
                   errorRange += Constants.getCoordinateUncertainty(atom1)
                                 +
                                 Constants.getCoordinateUncertainty(atom2);
                }

                float dist = Mathematics.distance(atom1.getXYZ(),
                                                  atom2.getXYZ());
                if (dist <= Float.parseFloat(CrossLinkParameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                     )
                                              ) + errorRange) {
                    conforming = true;
                }

                CrossLink xl = new CrossLink(atom1, atom2);

                // set file name of each cross-link to the file name of
                // its protein complex
                if (distanceFileCrossLinks == null
                    ||
                    !Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                     Parameter.DO_KEEP_FILENAME
                                                               ))) {
                    xl.setFileName(complex.getName());
                }

                if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                        Parameter.DO_PROBABILITY
                                                               ))) {
                    xl.setEucProbability();
                }
                if (!conforming) {
                    xl.setSolventPathDistance(
                            xwalk.constants.Constants.NON_CONFORMING_CROSS_LINK
                                             );
                }

                if (digest == null) {
                    crossLinks.add(xl);
                } else {
                    for (PolyPeptide peptide : digest) {
                        for (AminoAcid aa : peptide) {
                            // cross-link can only be at the central mis-cleaved
                            // site not at the C-term
                            if (aa.getAllAtoms().contains(atom1)
                                &&
                                !peptide.get(
                                        peptide.size() - 1
                                           ).getAllAtoms().contains(atom1)) {
                                atom1TrypticPeptide = peptide;
                                // Added these lines to allow for
                                // self cross-links within a peptide.
                                if (aa.getAllAtoms().contains(atom2)
                                    &&
                                    !peptide.get(
                                                peptide.size() - 1
                                                ).getAllAtoms().contains(
                                                                         atom2
                                                                        )) {
                                    atom2TrypticPeptide = peptide;
                                }
                                break;
                            }
                            if (aa.getAllAtoms().contains(atom2)
                                &&
                                !peptide.get(
                                             peptide.size() - 1
                                            ).getAllAtoms().contains(atom2)) {
                                atom2TrypticPeptide = peptide;
                                break;
                            }
                        }
                    }
                    if (atom1TrypticPeptide != null
                        &&
                        atom2TrypticPeptide != null) {

                        xl.setPeptides(atom1TrypticPeptide,
                                       atom2TrypticPeptide);
                        crossLinks.add(xl);
                    }
                }
            }
        }

        //---------------------------------
        // sort list of cross-links by Euclidean distance to ensure
        // consistency in distance list. The distance list can otherwise deviate
        // slightly from run to run depending on the order of the crossLinks
        // in the crossLink list.
        crossLinks.sort();
        return crossLinks;
    }

    //--------------------------------------------------------------------------
    /**
     * Extracts all necessary atom coordinates from the PDB file as defined in
     * the CrossLinkParameters.
     * @return List of PolyPeptideList objects, each holding all protein
     *         coordinates labeled as ATOM up to an END flag or end of file.
     * @throws IOException if input file could not be read.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         <a href="http://www.wwpdb.org/documentation/format32/sect9.html">
     *         PDB standards</a>.
     * @throws DataFormatException if .gz format is unknown.
     */
    public static ArrayList < PolyPeptideList > getComplexesCoordinates()
                                                   throws IOException,
                                                          FileFormatException,
                                                          DataFormatException {

        ArrayList < PDBreader > pdbReaders =
            PDBreader.createPDBreaders(CrossLinkParameter.getParameter(
                                                           Parameter.INFILE_PATH
                                                             ));

        ArrayList < PolyPeptideList > proteinComplexes =
              CrossLinkUtilities.extractProteinComplexes(pdbReaders);

        // assign vdW radius to protein atoms
        for (PolyPeptideList polyPeptideList : proteinComplexes) {
             CrossLinkUtilities.setRadius(polyPeptideList);
        }

    return proteinComplexes;
    }
    //--------------------------------------------------------------------------
    /**
     * Extracts user set chain and alternative location based PDBcomplex objects
     * from PDBreader objects.
     * @param pdbReaders
     *        - List of PDBreader objects each holding the content of a single
     *          PDB file.
     * @return List of PolyPeptideList objects holding only coordinates
     *         of atoms with user defined chain and alternative location
     *         information.
     */
    public static ArrayList < PolyPeptideList > extractProteinComplexes(
                                       final ArrayList < PDBreader > pdbReaders
                                                                       ) {
        ArrayList < PolyPeptideList > proteinComplexes =
            new ArrayList < PolyPeptideList >();

        /* All coordinates in a PDB file should be read in. The -c1 and -c2
         * flags should not limit the coordinate retrieval but only
         * the amino acid selection for calculating the SASD.
        if (parameter.getParameter(
                Parameter.CHAIN_ID1).equals(
                        parameter.getParameter(Parameter.CHAIN_ID2)
                                           )
                                   &&
                 parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1).equals(
                 parameter.getParameter(Parameter.ALTERNATIVE_LOCATION2)
                                   )
            ) {
            // as the restrictions are equal for both ends of the virtual
            // cross-links, it is unimportant which Id informations are taken.
            for (PDBreader reader : pdbReaders) {
                proteinComplexes.addAll(reader.getProteinComplex(
                         parameter.getParameter(Parameter.CHAIN_ID1),
                         parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1)
                                                         )
                                );
            }
        } else {

            for (PDBreader reader : pdbReaders) {
                proteinComplexes.addAll(reader.getProteinComplex(
                         parameter.getParameter(Parameter.CHAIN_ID1)
                       + parameter.getParameter(Parameter.CHAIN_ID2),
                         parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1)
                       + parameter.getParameter(Parameter.ALTERNATIVE_LOCATION2)
                                                         )
                                );
            }
        }
*/
        for (PDBreader reader : pdbReaders) {
            proteinComplexes.addAll(reader.getProteinComplex(
                                                         Constants.ALPHANUMERIC,
                                                         Constants.ALPHANUMERIC
                                                            ));
        }

        if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                 Parameter.DO_BACKBONE_READ))) {
            // replace full atom coordinates of protein with backbone and
            // beta-carbon only coordinates
            for (PolyPeptideList proteinComplex : proteinComplexes) {
                for (PolyPeptide protein : proteinComplex) {
                    for (int i = 0; i < protein.size(); i++) {
                        AtomList backbone = new AtomList();
                        for (Atom atom : protein.get(i).getAllAtoms()) {
                            if (atom.getType() == AtomType.CARBON_ALPHA
                               ||
                               atom.getType() == AtomType.CARBON_BETA
                               ||
                               atom.getType() == AtomType.NITROGEN
                               ||
                               atom.getType() == AtomType.OXYGEN) {
                                backbone.add(atom);
                            }
                        }
                        AminoAcid bb = new AminoAcid(backbone);
                        protein.set(i, bb);
                    }
                }
            }
        }
        return proteinComplexes;
    }

    //--------------------------------------------------------------------------
    /**
     * Digest all protein components of a protein complex.
     * @param proteinComplex
     *        - PolyPeptideList objects to be digested.
     * @param useExpasyRules
     *        - boolean value indicating to use
     *          <a href="http://www.expasy.ch/tools/peptidecutter/
     *          peptidecutter_special_enzymes.html">ExPASy Exception rules</a>
                for digestion.
     * @return new PolyPeptideList object with PolyPeptide objects holding the
     *         digested peptide segments.
     */
    public static PolyPeptideList trypsinate(
                           final PolyPeptideList proteinComplex,
                           final boolean useExpasyRules
                                                     ) {
        PolyPeptideList digestedComplex = new PolyPeptideList();
        for (PolyPeptide protein : proteinComplex) {
            ArrayList < PolyPeptide > digest =
                                     Digestion.partialTrypticDigest(
                                                                  protein,
                                                                  useExpasyRules
                                                                   );
            digestedComplex.addAll(digest);
            digestedComplex.setName(proteinComplex.getName());
        }
        return digestedComplex;
    }

    //--------------------------------------------------------------------------
    /**
     * Sets the van der Waals radius of the protein complex atoms appropriately,
     * either only to SURFNET radii or if SAS calculation should be carried out
     * to SURFNET radius + Solvent molecule radius.
     * @param complex -
     *        Protein complex object.
     * @throws IOException if an error occurs while reading the parameter file.
     */
    private static void setRadius(final PolyPeptideList complex)
                                                            throws IOException {

        // Finally set atom radii to SURFNET ones.
        complex.setAtomRadii();

        // If cross-links should be excluded by SAS, then van der Waals radii
        // must be increased by the radius of the solvent molecule.
        /*        if (
                Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                       )
                                )
           &&
           Boolean.parseBoolean(parameter.getParameter(Parameter.DO_SAS))) {
*/

            float solventRadius = Float.parseFloat(
                                             CrossLinkParameter.getParameter(
                                                        Parameter.SOLVENT_RADIUS
                                                                            )
                                                     );
            for (PolyPeptide protein : complex) {
                for (AminoAcid aa : protein) {
                    for (Atom atom : aa.getAllAtoms()) {
                        atom.setVanDerWaalsRadius(atom.getVanDerWaalsRadius()
                                                + solventRadius);
                    }
                }
            }
//        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns pairs of atoms that conform to the atom and amino acid
     * identifiers set by the user and which have a Euclidean distance smaller
     * then the user set maxDist.
     * @param complex -
     *        Protein complex object.
     * @return Hashtable of atom pairs that conform to the user set identifiers.
     */
    private static Hashtable < Atom, AtomList > findRelevantPairs(
                                             final PolyPeptideList complex
                                                                 ) {
        ArrayList <TreeSet < AtomList >> relevantAtoms =
                                       new ArrayList < TreeSet <AtomList >>();
        relevantAtoms.add(CrossLinkUtilities.findAllRelevantAtoms1(complex));
        relevantAtoms.add(CrossLinkUtilities.findAllRelevantAtoms2(complex));

        Hashtable < Atom, AtomList > pairs =
                    CrossLinkUtilities.createPairsBetweenRelevantAtoms(
                                                           relevantAtoms.get(0),
                                                           relevantAtoms.get(1)
                                                                      );

        pairs = CrossLinkUtilities.fixIntraInterSelection(pairs);
        return pairs;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a non-redundant/unique list of all atoms found in the
     * crossLinks object.
     * @param crossLinks
     *        - List of CrossLink objects extracted from the distance file.
     * @return List of non-redundant/unique cross-linked atoms.
     */
    private static AtomList uniqueXlAtomList(final CrossLinkList crossLinks) {
        // get all cross-linked atoms
        AtomList uniqueCrossLinkAtoms = new AtomList();
        for (CrossLink xl : crossLinks) {
             Atom preAtom = xl.getPreAtom();
             Atom postAtom = xl.getPostAtom();
             if (!uniqueCrossLinkAtoms.contains(preAtom)) {
                 uniqueCrossLinkAtoms.add(preAtom);
             }
             if (!uniqueCrossLinkAtoms.contains(postAtom)) {
                 uniqueCrossLinkAtoms.add(postAtom);
             }
        }
        return uniqueCrossLinkAtoms;
    }

    //--------------------------------------------------------------------------
    /**
     * Returns all atoms from the complex that match a cross-linked atom.
     * @param complex -
     *        Protein complex object.
     * @param crossLinkAtom -
     *        Atom object found to be cross-linked in the complex.
     * @return List of atoms in the complex matching the cross-linked atom.
     */
    private static ArrayList<AtomList> proteinAtomsMatchingXLatom(
                                             final PolyPeptideList complex,
                                             final Atom crossLinkAtom
                                                                   ) {
        ArrayList<AtomList> candidate = new ArrayList<AtomList>();
        for (AminoAcid aa1 : complex.getAllAminoAcids()) {
            // if atom2 could originate from residue aa1 than continue
            if (aa1.getType().getThreeLetterCode().equals(
                                                  crossLinkAtom.getResidueName()
                                                         )
                &&
                aa1.getNumber() == crossLinkAtom.getResidueNumber()) {
                // if user hasn't specified any atom name in the distance
                // file than simply consider all atoms from the residue
                if (crossLinkAtom.getName().equals("")) {
                    if (crossLinkAtom.getChainId() == ' ') {
                       candidate.add(aa1.getAllAtoms());
                    } else if (
                       crossLinkAtom.getChainId() == aa1.getAtom(0).getChainId()
                            ) {
                         candidate.add(aa1.getAllAtoms());
                    }
                } else {
                    for (Atom atom : aa1.getAllAtoms()) {

                        if (atom.getName().trim().equals(
                                                  crossLinkAtom.getName().trim()
                                                         )) {
                            AtomList list = new AtomList();
                            if (crossLinkAtom.getChainId() == ' ') {
                                list.add(atom);
                                candidate.add(list);
                            } else if (
                                 crossLinkAtom.getChainId() == atom.getChainId()
                                     ) {
                                if (crossLinkAtom.getAlternativeLocation()
                                                                       == ' ') {
                                    list.add(atom);
                                    candidate.add(list);
                                } else {
                                    if (crossLinkAtom.getAlternativeLocation()
                                        ==
                                        atom.getAlternativeLocation()) {
                                        list.add(atom);
                                        candidate.add(list);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return candidate;
    }
    //--------------------------------------------------------------------------
    /**
     * Extracts all relevant pairs of atoms from a user given distance file.
     * @param complex -
     *        Protein complex object.
     * @param crossLinks
     *        - List of CrossLink objects extracted from the distance file.
     * @return Hashtable of atom pairs that conform to the user set identifiers.
     * @throws IOException if input file could not be read.
     */
    public static Hashtable < Atom, AtomList > extractRelevantPairs(
                                             final PolyPeptideList complex,
                                             final CrossLinkList crossLinks)
                                                            throws IOException {

        // The idea is: 1.) get a non-redundant list of XL atoms form the
        // distance file. 2.) find all possible atoms in the complex that could
        // match to the list of non-redundant XL atoms in (1.) and place them
        // in a hash. 3.) Go through the list of cross-links in the distance
        // file, get the possible atoms from the hash and determine the pair of
        // atoms that are closest according to their Euclidean distance. 4.)
        // Return back the list of the closest atom pairs.

        Hashtable < Atom, AtomList > pairs = new Hashtable < Atom, AtomList >();

        AtomList uniqueCrossLinkAtoms =
                                CrossLinkUtilities.uniqueXlAtomList(crossLinks);

        Hashtable<Atom, ArrayList<AtomList>> xlAtomMatchingComplexAtoms =
                                     new Hashtable<Atom, ArrayList<AtomList>>();
        for (Atom xlAtom : uniqueCrossLinkAtoms) {
            ArrayList<AtomList> candidates =
                            CrossLinkUtilities.proteinAtomsMatchingXLatom(
                                                                        complex,
                                                                        xlAtom
                                                                         );
            if (candidates.size() != 0) {
                xlAtomMatchingComplexAtoms.put(xlAtom, candidates);
            }
        }
        // Assign atoms from complex to the cross-links found in the distance
        // file
        for (CrossLink xl : crossLinks) {
             Atom preAtom = xl.getPreAtom();
             Atom postAtom = xl.getPostAtom();

             // find both cross-link atoms in the xlAtomMatchingComplexAtoms
             // object
             boolean foundPreAtom = false;
             boolean foundPostAtom = false;
             for (Atom xlAtom : xlAtomMatchingComplexAtoms.keySet()) {
                 if (xlAtom.equals(preAtom)) {
                     preAtom = xlAtom;
                     foundPreAtom = true;
                 }
                 if (xlAtom.equals(postAtom)) {
                     postAtom = xlAtom;
                     foundPostAtom = true;
                 }
                 if (foundPreAtom && foundPostAtom) {
                     break;
                 }
             }
             // if either of both could not be found, then either of both
             // has no coordinates in the protein complex.
             if (!foundPreAtom || !foundPostAtom) {
                 // Output a WARNING message that atoms in distance files could
                 // not be found in the input file
                 System.err.print("WARNING: ");
                 if (!foundPreAtom && foundPostAtom) {
                     System.err.print("1st atom ");
                 } else if (foundPreAtom && !foundPostAtom) {
                     System.err.print("2nd atom ");
                 } else {
                     System.err.print("Both atoms ");
                 }
                 System.err.print("could not be found : " + xl.getIndex() + "\t"
                                + xl.toString());
                 continue;
             }

             ArrayList<AtomList> preAtomCandidates =
                                        xlAtomMatchingComplexAtoms.get(preAtom);
             ArrayList<AtomList> postAtomCandidates =
                                       xlAtomMatchingComplexAtoms.get(postAtom);

             for (AtomList preAtoms : preAtomCandidates) {
                 for (AtomList postAtoms : postAtomCandidates) {
                     float minDist = Integer.MAX_VALUE;
                     Atom minPreAtom = null;
                     Atom minPostAtom = null;
                     for (Atom preAtomCandidate : preAtoms) {
                         for (Atom postAtomCandidate : postAtoms) {
                             float dist = Mathematics.distance(
                                                  preAtomCandidate.getXYZ(),
                                                  postAtomCandidate.getXYZ()
                                                               );
                             if (dist < minDist) {
                                 minDist = dist;
                                 minPreAtom = preAtomCandidate;
                                 minPostAtom = postAtomCandidate;
                             }
                         }
                     }
                     AtomList pair = pairs.get(minPreAtom);
                     if (pair == null) {
                         pair = new AtomList();
                         pair.add(minPostAtom);
                         pairs.put(minPreAtom, pair);
                     } else {
                         pair.add(minPostAtom);
                         pairs.put(minPreAtom, pair);
                     }
                 }
             }
        }
        return pairs;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the indices of the cross-linked objects either iteratively or if
     * a distance file is set by the user according to the indices in the
     * distance file.
     * @param crossLinkList
     *        - List of cross-link object.
     * @throws IOException if an error occurs while reading a distance file.
     */
    private static void setCrossLinkIndicesAndFileName(
                                               final CrossLinkList crossLinkList
                                                      ) throws IOException {
        if (CrossLinkParameter.getParameter(
                                     Parameter.DISTANCE_FILE_PATH).equals("")) {
            for (int i = 0; i < crossLinkList.size(); i++) {
                crossLinkList.get(i).setIndex(i + 1);
            }
        } else {
            String fileName = CrossLinkParameter.getParameter(
                                                    Parameter.DISTANCE_FILE_PATH
                                                    );
                boolean onlyIntra = Boolean.parseBoolean(
                                      CrossLinkParameter.getParameter(
                                            Parameter.DO_INTRAMOLECULAR_DISTANCE
                                                                     ));
                boolean onlyInter = Boolean.parseBoolean(
                                      CrossLinkParameter.getParameter(
                                            Parameter.DO_INTERMOLECULAR_DISTANCE
                                                                     ));

            CrossLinkList distanceFileCrossLinks = DistanceReader.getCrossLinks(
                                                                      fileName,
                                                                      onlyIntra,
                                                                      onlyInter,
                                                                      false);

            // assign indices of distance file to newly found cross-links.
            if (distanceFileCrossLinks.size() != 0) {
                for (CrossLink dxl : distanceFileCrossLinks) {
                    for (CrossLink xl : crossLinkList) {
                        int reverse = -1;
                        // first adapt the order of the cross-link atoms in
                        // the vXL and the XL in the distance file.
                        if (dxl.getPreAtom().getResidueName().equals(
                                                xl.getPreAtom().getResidueName()
                                                                    )
                           &&
                           dxl.getPostAtom().getResidueName().equals(
                                               xl.getPostAtom().getResidueName()
                                                                    )
                           &&
                           dxl.getPreAtom().getResidueNumber()
                               ==
                           xl.getPreAtom().getResidueNumber()
                           &&
                           dxl.getPostAtom().getResidueNumber()
                               ==
                           xl.getPostAtom().getResidueNumber()
                        ) {
                            reverse = 0;
                        } else if (dxl.getPreAtom().getResidueName().equals(
                                               xl.getPostAtom().getResidueName()
                                                                           )
                                &&
                                dxl.getPostAtom().getResidueName().equals(
                                                xl.getPreAtom().getResidueName()
                                                                         )
                                &&
                                dxl.getPreAtom().getResidueNumber()
                                    ==
                                xl.getPostAtom().getResidueNumber()
                                &&
                                dxl.getPostAtom().getResidueNumber()
                                    ==
                                xl.getPreAtom().getResidueNumber()
                        ) {
                            reverse = 1;
                        } else {
                            reverse = -1;
                        }
                        Atom preAtom = null;
                        Atom postAtom = null;
                        switch (reverse) {
                            case  0 :
                                preAtom = dxl.getPreAtom();
                                postAtom = dxl.getPostAtom();
                                break;
                            case  1 :
                                preAtom = dxl.getPostAtom();
                                postAtom = dxl.getPreAtom();
                                break;
                            default :
                                continue;
                        }
                        // Now find vXL in the distance file and assign index
                        // in the distance file to vXL.
                        if (preAtom.getChainId() != ' ') {
                            if (preAtom.getChainId()
                                !=
                                xl.getPreAtom().getChainId()) {
                                continue;
                            }
                        }
                        if (postAtom.getChainId() != ' ') {
                            if (postAtom.getChainId()
                                !=
                                xl.getPostAtom().getChainId()) {
                                continue;
                            }
                        }
                        if (!preAtom.getName().equals("")) {
                            if (!preAtom.getName().trim().equals(
                                               xl.getPreAtom().getName().trim())
                                                         ) {
                                continue;
                            }
                        }
                        if (!postAtom.getName().trim().equals("")) {
                            if (!postAtom.getName().equals(
                                              xl.getPostAtom().getName().trim())
                                                          ) {
                                continue;
                            }
                        }
                        if (preAtom.getResidueNumber()
                            !=
                            xl.getPreAtom().getResidueNumber()) {
                            continue;
                        }
                        if (postAtom.getResidueNumber()
                            !=
                            xl.getPostAtom().getResidueNumber()) {
                            continue;
                        }
                        if (!preAtom.getResidueName().equalsIgnoreCase(
                                                xl.getPreAtom().getResidueName()
                                                                      )) {
                            continue;
                        }
                        if (!postAtom.getResidueName().equalsIgnoreCase(
                                xl.getPostAtom().getResidueName()
                                                      )) {
                            continue;
                        }

                        xl.setIndex(dxl.getIndex());
                        if (Boolean.parseBoolean(
                                                CrossLinkParameter.getParameter(
                                                      Parameter.DO_KEEP_FILENAME
                                                                       ))) {
                            xl.setFileName(dxl.getFileName());
                        }
                    }
                }
            }
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Returns all atoms that conform to the identifier of the first atom as set
     * by the user.
     * @param complex -
     *        Protein complex object.
     * @return TreeSet of AtomList objects that hold all atoms of amino acids
     *         that conform to the user set identifiers.
     */
    private static TreeSet < AtomList > findAllRelevantAtoms1(
                                             final PolyPeptideList complex) {

        TreeSet < AtomList > candidates1 = new TreeSet < AtomList >(
                new Comparator<AtomList>() {
                    public int compare(final AtomList a1, final AtomList a2) {
                        for (Atom atom1 : a1) {
                            for (Atom atom2 : a2) {
                                if (atom1.equals(atom2)) {
                                    return 0;
                                }
                            }
                        }
                        return 1;
                    }
                });

        for (PolyPeptide protein : complex) {
            for (AminoAcid residue : protein) {
                if (CrossLinkUtilities.isAminoAcid1Relevant(residue)) {
                    if (CrossLinkParameter.getParameter(
                                                Parameter.ATOM_TYPE1).equals("")
                                              ) {
                        candidates1.add(residue.getAllAtoms());
                    } else {
                        AtomList list = new AtomList();
                        for (Atom atom : residue.getAllAtoms()) {
                            if (CrossLinkUtilities.isAtomRelevant1(atom)) {
                                list.add(atom);
                            }
                        }
                        if (list.size() != 0) {
                            candidates1.add(list);
                        } else {
                            System.err.println("WARNING: "
                                 + CrossLinkParameter.getParameter(
                                                           Parameter.INFILE_PATH
                                                                  )
                                 + "\tAtom \"-a1 "
                                 + CrossLinkParameter.getParameter(
                                                            Parameter.ATOM_TYPE1
                                                   ).replaceAll("#", " ").trim()
                                 + "\" not found in residue "
                                 + residue.getAtom(0).getResidueName()
                                 + residue.getAtom(0).getResidueNumber()
                                 + residue.getAtom(0).getChainId());
                        }
                    }
                }
            }
        }
        return candidates1;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns all atoms that conform to the identifier of the second atom as
     * set by the user.
     * @param complex
     *        - Protein complex object
     * @return TreeSet of AtomList objects that hold all atoms of amino acids
     *         that conform to the user set identifiers.
     */
    private static TreeSet < AtomList > findAllRelevantAtoms2(
                                             final PolyPeptideList complex
                                                               ) {

        // TreeSet to sort atoms in a cross-link.
        TreeSet < AtomList > candidates2 = new TreeSet < AtomList >(
                new Comparator<AtomList>() {
                    public int compare(final AtomList a1, final AtomList a2) {
                        for (Atom atom1 : a1) {
                            for (Atom atom2 : a2) {
                                if (atom1.equals(atom2)) {
                                    return 0;
                                }
                            }
                        }
                        return 1;
                    }
                });

        StringBuffer dataNotFoundMessage = new StringBuffer();

        for (PolyPeptide protein : complex) {
            for (AminoAcid residue : protein) {
                if (CrossLinkUtilities.isAminoAcid2Relevant(residue)) {
                    if (CrossLinkParameter.getParameter(
                                               Parameter.ATOM_TYPE2).equals(""))
                    {
                        candidates2.add(residue.getAllAtoms());
                    } else {
                        AtomList list = new AtomList();
                        for (Atom atom : residue.getAllAtoms()) {
                            if (CrossLinkUtilities.isAtomRelevant2(atom)) {
                                list.add(atom);
                            }
                        }
                        if (list.size() != 0) {
                            candidates2.add(list);
                        } else {
                            dataNotFoundMessage.append("WARNING: "
                                 + CrossLinkParameter.getParameter(
                                                           Parameter.INFILE_PATH
                                                                  )
                                 + "\tAtom \"-a2 "
                                 + CrossLinkParameter.getParameter(
                                                     Parameter.ATOM_TYPE2
                                                   ).replaceAll("#", " ").trim()
                                 + "\" not found in residue "
                                 + residue.getAtom(0).getResidueName()
                                 + residue.getAtom(0).getResidueNumber()
                                 + residue.getAtom(0).getChainId()
                                 + Constants.LINE_SEPERATOR);
                        }
                    }
                }
            }
        }

    return candidates2;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether an AminoAcid object conforms to the name, number and chain
     * Id of the first amino acid as set by the user.
     * @param acid
     *        - Amino acid object
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links.
     * @return {@code TRUE} if atom corresponds to the first amino acid
     *         identifier as set by the user, {@code FALSE} otherwise.
     */
    private static boolean isAminoAcid1Relevant(final AminoAcid acid) {
        boolean residueNumberFound = CrossLinkParameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NAME1
                                                           ).equals("")
                                     &&
                                     CrossLinkParameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER1
                                                           ).indexOf(
                                              "#"
                                            + acid.getAtom(0).getResidueNumber()
                                            + "#"
                                                                    ) != -1;
        boolean residueNameAndNumberFound = !CrossLinkParameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NAME1
                                                                   ).equals("")
                                            &&
                                         !CrossLinkParameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NUMBER1
                                                                ).equals("-999")
                                            &&
                                            CrossLinkParameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NAME1
                                                                ).indexOf(
                                                "#"
                                              + acid.getAtom(0).getResidueName()
                                              + "#"
                                                                         ) != -1
                                            &&
                                            CrossLinkParameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER1
                                                                 ).indexOf(
                                              "#"
                                            + acid.getAtom(0).getResidueNumber()
                                            + "#"
                                                                        ) != -1;

        boolean residueNameFound = !CrossLinkParameter.getParameter(
                                           Parameter.AMINO_ACID_RESIDUE_NAME1
                                                          ).equals("")
                                    &&
                                    CrossLinkParameter.getParameter(
                                           Parameter.AMINO_ACID_RESIDUE_NUMBER1
                                                          ).equals("-999")
                                    &&
                                    CrossLinkParameter.getParameter(
                                           Parameter.AMINO_ACID_RESIDUE_NAME1
                                                          ).indexOf(
                                                "#"
                                              + acid.getAtom(0).getResidueName()
                                              + "#"
                                                                   ) != -1;

        boolean chainFound = CrossLinkParameter.getParameter(
                                                    Parameter.CHAIN_ID1
                                                    ).indexOf(
                           Character.toString(acid.getAtom(0).getChainId())
                                                             ) != -1;
        if (residueNumberFound || residueNameAndNumberFound
                || residueNameFound) {
            if (chainFound) {
                return true;
            }
        }
        return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether an AminoAcid object conforms to the name, number and chain
     * Id of the second amino acid as set by the user.
     * @param acid
     *        - Amino acid object
     * @return {@code TRUE} if atom corresponds to the second amino acid
     *         identifier as set by the user, {@code FALSE} otherwise.
     */
    private static boolean isAminoAcid2Relevant(
                                              final AminoAcid acid
                                               ) {
        boolean residueNumberFound = CrossLinkParameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NAME2
                                                           ).equals("")
                                     &&
                                     CrossLinkParameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER2
                                                           ).indexOf(
                                              "#"
                                            + acid.getAtom(0).getResidueNumber()
                                            + "#"
                                                                     ) != -1;

        boolean residueNameAndNumberFound = !CrossLinkParameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NAME2
                                                                   ).equals("")
                                            &&
                                         !CrossLinkParameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NUMBER2
                                                                ).equals("-999")
                                            &&
                                            CrossLinkParameter.getParameter(
                                              Parameter.AMINO_ACID_RESIDUE_NAME2
                                                                ).indexOf(
                                                "#"
                                              + acid.getAtom(0).getResidueName()
                                              + "#"
                                                                         ) != -1
                                            &&
                                            CrossLinkParameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER2
                                                                ).indexOf(
                                              "#"
                                            + acid.getAtom(0).getResidueNumber()
                                            + "#"
                                                                        ) != -1;

        boolean residueNameFound = !CrossLinkParameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NAME2
                                                          ).equals("")
                                    &&
                                    CrossLinkParameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER2
                                                          ).equals("-999")
                                    &&
                                    CrossLinkParameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NAME2
                                                          ).indexOf(
                                                "#"
                                              + acid.getAtom(0).getResidueName()
                                              + "#"
                                                                   ) != -1;
        boolean chainFound = CrossLinkParameter.getParameter(
                                                    Parameter.CHAIN_ID2
                                                   ).indexOf(
                                Character.toString(acid.getAtom(0).getChainId())
                                                            ) != -1;

        if (residueNumberFound || residueNameAndNumberFound
            || residueNameFound) {
            if (chainFound) {
                return true;
            }
        }
        return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether an Atom object conforms to the name and alternative
     * location of the first atom as set by the user.
     * @param atom
     *        - Atom object
     * @return {@code TRUE} if atom corresponds to the first atom identifier as
     *         set by the user, {@code FALSE} otherwise.
     */
    private static boolean isAtomRelevant1(final Atom atom) {
        if (!CrossLinkParameter.getParameter(Parameter.ATOM_TYPE1).equals("")) {
            if (CrossLinkParameter.getParameter(Parameter.ATOM_TYPE1).indexOf(
                                                           "#"
                                                         + atom.getName().trim()
                                                         + "#"
                                                                    ) != -1
                &&
                CrossLinkParameter.getParameter(
                        Parameter.ALTERNATIVE_LOCATION1).indexOf(
                               Character.toString(atom.getAlternativeLocation())
                                                                ) != -1) {
                return true;
            }
            return false;
        } else {
            return true;
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether an Atom object conforms to the name and alternative
     * location of the second atom as set by the user.
     * @param atom
     *        - Atom object
     * @return {@code TRUE} if atom corresponds to the second atom identifier as
     *         set by the user, {@code FALSE} otherwise.
     */
    private static boolean isAtomRelevant2(final Atom atom) {
        if (!CrossLinkParameter.getParameter(Parameter.ATOM_TYPE2).equals("")) {
            if (CrossLinkParameter.getParameter(Parameter.ATOM_TYPE2).indexOf(
                                                           "#"
                                                         + atom.getName().trim()
                                                         + "#") != -1
                &&
                CrossLinkParameter.getParameter(
                        Parameter.ALTERNATIVE_LOCATION2).indexOf(
                               Character.toString(atom.getAlternativeLocation())
                                                                ) != -1) {
                return true;
            }
            return false;
        } else {
            return true;
        }
    }

    //--------------------------------------------------------------------------

    /**
     * Returns those pairs of potential cross-linkable atoms that have an
     * Euclidean distance smaller than the user set maximum distance.
     * @param candidates1
     *        - Set of AtomList objects that fulfill all criteria set by the
     *          user on the commandline for the first cross-linked atoms.
     * @param candidates2
     *        - Set of AtomList objects that fulfill all criteria set by the
     *          user on the commandline for the second cross-linked atoms.
     * @return Hashtable of potential cross-linkable atoms.
     */
    private static Hashtable < Atom, AtomList > createPairsBetweenRelevantAtoms(
                                       final TreeSet < AtomList > candidates1,
                                       final TreeSet < AtomList > candidates2) {
        Hashtable < Atom, AtomList > pairs = new Hashtable < Atom, AtomList >();
        for (AtomList list1 : candidates1) {
            for (AtomList list2 : candidates2) {
                // An amino acid can not be self-cross-linked.
                if (MatterUtilities.equalsResidue(list1.get(0), list2.get(0))) {
                    continue;
                } else {
                    AtomList minimumDistanceAtomPair =
                            MatterUtilities.getClosestAtomPair(list1, list2);

                    float dist = Mathematics.distance(
                                    minimumDistanceAtomPair.get(0).getXYZ(),
                                    minimumDistanceAtomPair.get(1).getXYZ()
                                                      );

                    float errorRange = 0;
                    if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                       Parameter.DO_BFACTOR))) {
                       errorRange += Constants.getCoordinateUncertainty(
                                                 minimumDistanceAtomPair.get(0))
                                     +
                                     Constants.getCoordinateUncertainty(
                                                minimumDistanceAtomPair.get(1));
                   }

                    if (dist > Float.parseFloat(CrossLinkParameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                           )
                                                 ) + errorRange
                       ) {
                        continue;
                    } else {
                        // The distance between these two amino acids would
                        // indicate that both could be cross-linked, at least in
                        // terms of their Euclidean distance.
                        // At a later stage, these potential candidates should
                        // be checked further to have a SASD that
                        // conforms to the length of the cross-linker.
                        Atom atom0 = minimumDistanceAtomPair.get(0);
                        Atom atom1 = minimumDistanceAtomPair.get(1);
                        AtomList associate0 = pairs.get(atom0);
                        AtomList associate1 = pairs.get(atom1);
                        if (associate0 == null && associate1 == null) {
                            associate0 = new AtomList();
                            associate0.add(atom1);
                            pairs.put(atom0, associate0);
                        } else {
                            if (associate0 == null) {
                                if (!associate1.contains(atom0)) {
                                    associate1.add(atom0);
                                    pairs.put(atom1, associate1);
                                }
                            } else if (associate1 == null) {
                                if (!associate0.contains(atom1)) {
                                    associate0.add(atom1);
                                    pairs.put(atom0, associate0);
                                }
                            }
                        }
                    }
                }
            }
        }
        return pairs;
    }

    //--------------------------------------------------------------------------

    /**
     * Returns pairs of atoms that depending on the user set parameter, contain
     * only intra, inter or all potential cross-links.
     * @param pairs
     *        - Hashtable of all pairs of atom that conform to the atom and
     *          amino acid identifiers as set by the user.
     * @return Hashtable of intra-, inter or intra/inter atom pairs.
     */
    private static Hashtable < Atom, AtomList > fixIntraInterSelection(
                                        final Hashtable < Atom, AtomList > pairs
                                                                      ) {
        Hashtable < Atom, AtomList > newPairs =
                                            new Hashtable < Atom, AtomList >();

        for (Atom atom1 : pairs.keySet()) {
            for (Atom atom2 : pairs.get(atom1)) {
                if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                           Parameter.DO_INTRAMOLECULAR_DISTANCE)
                                                               )
                                         &&
                   !Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                           Parameter.DO_INTERMOLECULAR_DISTANCE)
                                                               )
                                        ) {
                    if (atom1.getChainId() != atom2.getChainId()) {
                        continue;
                    }
                } else if (!Boolean.parseBoolean(
                                        CrossLinkParameter.getParameter(
                                            Parameter.DO_INTRAMOLECULAR_DISTANCE
                                                                       )
                                                )
                           &&
                           Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                            Parameter.DO_INTERMOLECULAR_DISTANCE
                                                                      )
                                               )
                          ) {
                    if (atom1.getChainId() == atom2.getChainId()) {
                        continue;
                    }
                }
                AtomList atomList = newPairs.get(atom1);
                if (atomList == null) {
                    atomList = new AtomList();
                    atomList.add(atom2);
                } else {
                    atomList.add(atom2);
                }
                newPairs.put(atom1, atomList);
            }
        }
    return newPairs;
    }

    //--------------------------------------------------------------------------

    /**
     * Removes redundant atom pairs within homologous structures.
     * @param crossLinkList
     *        - CrossLinkList object holding all user set parameters for
     *          calculating cross-links
     */
    private static void removeRedundanciesInHomomers(
                                               final CrossLinkList crossLinkList
                                                    ) {

        // get all redundant cross links
        Hashtable < CrossLink, ArrayList < CrossLink > >
                                                   redundantCrossLinksCandidates
                      = new Hashtable < CrossLink, ArrayList < CrossLink > >();
        ArrayList <CrossLink> homologList = new ArrayList <CrossLink> ();
        for (int i = 0; i < crossLinkList.size(); i++) {
            CrossLink crossLink1 = crossLinkList.get(i);
            if (!homologList.contains(crossLink1)) {
                for (int j = i + 1; j < crossLinkList.size(); j++) {
                    CrossLink crossLink2 = crossLinkList.get(j);
                    if (crossLink1.equalsInHomolog(crossLink2)) {
                        ArrayList < CrossLink > list =
                            redundantCrossLinksCandidates.get(crossLink1);
                        if (list == null) {
                          list = new ArrayList < CrossLink >();
                        }
                        if (!list.contains(crossLink1)) {
                            homologList.add(crossLink2);
                            list.add(crossLink2);
                            redundantCrossLinksCandidates.put(crossLink1, list);
                        }
                    }
                }
            }
        }

        // remove all redundant cross links except of the one with the lowest
        // Euclidean/SolventPath distance.
        for (CrossLink crossLink1 : redundantCrossLinksCandidates.keySet()) {
            ArrayList <CrossLink> toBremoved = new ArrayList<CrossLink>();
            CrossLink minXL = crossLink1;
            float minDist =
                (minXL.getSolventPathDistance() < 0.0 ?
                 minXL.getEuclideanDistance() : minXL.getSolventPathDistance());
            for (CrossLink crossLink2 : redundantCrossLinksCandidates.get(
                                                                      crossLink1
                                                                          )) {
                float dist2 =
                    (crossLink2.getSolventPathDistance() < 0.0 ?
                     crossLink2.getEuclideanDistance()
                                         : crossLink2.getSolventPathDistance());

                if (minDist <= dist2) {
                    toBremoved.add(crossLink2);
                } else {
                    toBremoved.add(minXL);
                    minDist = dist2;
                    minXL = crossLink2;
                }
            }
            crossLinkList.removeAll(toBremoved);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Removes side chains from cross-linked amino acids.
     * @param complex
     *      - PolyPeptideList object holding all atoms of the protein.
     * @param crossLinksByEuclideanDistance
     *      - List of CrossLinks object that have all a Euclidean distance
     *        smaller then a user set maxDist value between cross-linked atoms.
     */
    private static void removeSideChainsFromCrossLinkedResidues(
                            final PolyPeptideList complex,
                            final CrossLinkList crossLinksByEuclideanDistance) {
        StringBuffer crossLinkedAtoms = new StringBuffer();
        for (CrossLink xl : crossLinksByEuclideanDistance) {
            crossLinkedAtoms.append("#"
                                  + AminoAcid.getAminoAcidId(xl.getPreAtom())
                                  + "#"
                                  + AminoAcid.getAminoAcidId(xl.getPostAtom())
                                  + "#");
        }
        for (PolyPeptide protein : complex) {
            for (AminoAcid aa : protein) {
                if (crossLinkedAtoms.indexOf(
                                         AminoAcid.getAminoAcidId(aa.getAtom(0))
                                            ) != -1) {
                        aa.removeSideChain();
                }
            }
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Creates CrossLink objects between potential cross-linkable atom pairs
     * that are closer than the user set maximum distance. The cross-links
     * objects will be generated on a local grid.
     * @param complex
     *      - PolyPeptideList object holding all atoms of the protein.
     * @param crossLinksByEuclideanDistance
     *      - List of CrossLinks object that have all a Euclidean distance
     *        smaller then a user set maxDist value between cross-linked atoms.
     */
    private static void calculatesSolventPathDistance(
                               final PolyPeptideList complex,
                               final CrossLinkList crossLinksByEuclideanDistance
                                                          ) {

        float gridCellSize = Float.parseFloat(CrossLinkParameter.getParameter(
                                                        Parameter.GRID_CELL_SIZE
                                                                       ));

        float maxDist = Float.parseFloat(CrossLinkParameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                        ));
        //---------------------------------
        // If requested remove side chain of cross-linked amino acids.
        // Not necessary however, if backbone read was set too.
        if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                Parameter.DO_REMOVE_SIDECHAINS))
                                                                ) {
                CrossLinkUtilities.removeSideChainsFromCrossLinkedResidues(
                                                   complex,
                                                   crossLinksByEuclideanDistance
                                                                          );
        }


        Hashtable <Atom, AtomList> pairs =
                                         crossLinksByEuclideanDistance.toHash();

        for (Atom atom : pairs.keySet()) {
            AtomList pairedAtoms = pairs.get(atom);

            String pairedAminoAcidId = "#" + AminoAcid.getAminoAcidId(atom)
                                           + atom.getName() + "#";
            for (Atom pairedAtom : pairedAtoms) {
                pairedAminoAcidId += "#" + AminoAcid.getAminoAcidId(pairedAtom)
                                         + pairedAtom.getName() + "#";
            }
            AtomList nonXLedAtoms = new AtomList();
            for (Atom complexAtom : complex.getAllAtoms()) {
                String complexAminoAcidId = "#" + AminoAcid.getAminoAcidId(
                                                                     complexAtom
                                                                          )
                                                + complexAtom.getName() + "#";
                if (pairedAminoAcidId.indexOf(complexAminoAcidId) == -1) {
                    nonXLedAtoms.add(complexAtom);
                }
            }

            AtomGrid grid = new AtomGrid(nonXLedAtoms,
                                         atom,
                                         maxDist,
                                         gridCellSize);

            ArrayList <Path> paths =
                CrossLinkUtilities.calculateShortestPathThroughSolvent(
                                                                     grid,
                                                                     atom,
                                                                     pairedAtoms
                                                                      );
            if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                      Parameter.DO_GRID_OUTPUT))
                                                           ) {
                System.out.println("HEADER " + atom.getResidueName().trim()
                                       + "-" + atom.getResidueNumber()
                                       + "-" + atom.getChainId()
                                       + "-" + atom.getName().trim()
                                 + Constants.LINE_SEPERATOR
                                 + grid.toString()
                                 + "TER");
            }
            for (int i = 0; i < pairedAtoms.size(); i++) {
                CrossLink crossLink = crossLinksByEuclideanDistance.get(
                                                              atom,
                                                              pairedAtoms.get(i)
                                                                       );
                // if the paths contain any grid cell, it means that a SASD path
                // could be calculation. If the paths are empty, it means
                // that the path calculation had to be stopped prematurely due
                // to solvent inaccessibility.
                float dist = 0;
                if (paths.size() > 0) {
                    dist = SolventPathDistance.extractTargetDistances(
                                                                    paths.get(i)
                                                                     );
                } else {
                    dist = Float.MAX_VALUE;
                }

                maxDist = Float.parseFloat(CrossLinkParameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                 ));
                float errorRange = 0;
                if (Boolean.parseBoolean(
                       CrossLinkParameter.getParameter(Parameter.DO_BFACTOR))) {
                    errorRange += Constants.getCoordinateUncertainty(atom)
                                  +
                                  Constants.getCoordinateUncertainty(
                                                            pairedAtoms.get(i));
                }

                // Conforming distance found!
                if (dist <= maxDist + errorRange) {
                    // for very short distances, the SASD can be shorter than
                    // the Euclidean distance due to the grid-ification of
                    // distance space. In such cases, set the SASD to the
                    // Euclidean distance.
                    if (dist < crossLink.getEuclideanDistance() && dist > 0) {
                        crossLink.setSolventPathDistance(
                                                crossLink.getEuclideanDistance()
                                                        );
                    } else {
                        crossLink.setSolventPathDistance(dist);
                    }
                    crossLink.setPath(paths.get(i));

                    if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                        Parameter.DO_PROBABILITY
                                                                   ))) {
                        crossLink.setSASDprobability();
                    }
                } else {
                    // if its not due to solvent inaccessibility, than the
                    // distance is simply to large
                    if (Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                    Parameter.DO_VERBOSE_OUTPUT)
                                                             )) {
                        System.err.println("Following XL exceed the maximum "
                                         + "distance of "
                                         + (maxDist + errorRange)
                                         + " with a distance of " + dist
                                         + Constants.LINE_SEPERATOR
                                         + crossLink);
                    }
                    crossLink.setSolventPathDistance(
                             xwalk.constants.Constants.NON_CONFORMING_CROSS_LINK
                                                    );
                }
            }
            // run garbage collector to free up space in particular from
            // AtomGrid object from the previous run.
            System.gc();
        }
        //---------------------------------
        // sort list of cross-links by SASD.
        crossLinksByEuclideanDistance.sort();

    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether both amino acid atoms in the CrossLink object are solvent
     * accessible. If not, it assigns the crossLink object one of the following
     * three SASD values:
     * {@link xwalk.constants.Constants.BOTH_ATOMS_ARE_SOLVENT_INACCESSIBLE}
     * {@link xwalk.constants.Constants.FIRST_ATOM_IS_SOLVENT_INACCESSIBLE}
     * {@link xwalk.constants.Constants.SECOND_ATOM_IS_SOLVENT_INACCESSIBLE}
     * @param crossLink
     *        - CrossLink object to checked for solvent accessibility.
     * @param grid
     *      - AtomGrid object holding the grid, which was used to calculate
     *        the shortest path for the crossLink object.
     * @return {@code TRUE} if both atoms of the CrossLink object are solvent
     *         accessible, {@code FALSE} otherwise
     */
    /*
    private static boolean isSolventAccessible(final CrossLink crossLink,
                                               final AtomGrid grid) {

        Atom atom1 = crossLink.getPreAtom();
        Atom atom2 = crossLink.getPostAtom();

        if (!GridUtilities.isAccessible(atom1, grid)
            &&
            !GridUtilities.isAccessible(atom2, grid)) {
            crossLink.setSolventPathDistance(
                   xwalk.constants.Constants.BOTH_ATOMS_ARE_SOLVENT_INACCESSIBLE
                                            );
            return false;
        } else if (!GridUtilities.isAccessible(atom1, grid)) {
                crossLink.setSolventPathDistance(
                    xwalk.constants.Constants.FIRST_ATOM_IS_SOLVENT_INACCESSIBLE
                                                );
                return false;
        } else if (!GridUtilities.isAccessible(atom2, grid)) {
                crossLink.setSolventPathDistance(
           xwalk.constants.Constants.SECOND_ATOM_IS_SOLVENT_INACCESSIBLE
                                                );
            return false;
        }
        return true;
    }
    */
    //--------------------------------------------------------------------------
    /**
     * Calculates solvent path distances using a local grid.
     * @param grid
     *      - Local AtomGrid object build around one potential cross-linkable
     *        atom.
     * @param atom1
     *        - First protein atom to be connected by the virtual cross-linker.
     * @param atoms2
     *        - List of atoms to be cross-linked to atom1, if distance is
     *          shorter than maxDist.
     * @return List of Path objects that form the shortest path between atom1
     *         and atom2 objects. Those paths that exceed -max with their length
     *         have only a single grid cell stored, namely the grid cell of the
     *         target cell.
     */
    private static ArrayList < Path > calculateShortestPathThroughSolvent(
                                              final AtomGrid grid,
                                              final Atom atom1,
                                              final AtomList atoms2
                                                            ) {
        // as soon as one of atom2 is solvent accessible calculate
        // shortest path.
        SolventPathDistance solvDist  = new SolventPathDistance(
                                                               atom1,
                                                               atoms2,
                                                               grid);
        // check which atoms are accessible
        ArrayList < Path > paths  = new ArrayList < Path >();
        boolean atom1isAccessible = false;
        if (GridUtilities.isAccessible(atom1, grid)) {
            atom1isAccessible = true;
        }
        boolean[] atoms2areAccessible = new boolean[atoms2.size()];
        boolean atom2isAccessible = false;
        for (int i = 0; i < atoms2.size(); i++) {
            if (GridUtilities.isAccessible(atoms2.get(i), grid)) {
                atoms2areAccessible[i] = true;
                atom2isAccessible = true;
            }
        }
        // only continue if atom1 is accessible
        if (atom1isAccessible && atom2isAccessible) {
            paths = solvDist.getShortestPath(
                                       Float.parseFloat(
                                               CrossLinkParameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                       )));
        }
        boolean successful = paths.size() == 0 ? false : true;
        for (int i = 0; i < atoms2.size(); i++) {
            // if distance calculations were performed but an empty path array
            // was returned, than the first atom is buried.
            if (atom1isAccessible && atom2isAccessible && !successful) {
                Path path = new Path();
                GridCell dummy = new GridCell(atom1.getXYZ(),
                                              GridCell.getSize());
                dummy.setDistance(
                              xwalk.constants.Constants.FIRST_ATOM_IS_BURIED
                                 );
                path.add(dummy);
                paths.add(i, path);
            } else if (!atom1isAccessible && atoms2areAccessible[i]) {
                Path path = new Path();
                GridCell dummy = new GridCell(atom1.getXYZ(),
                                              GridCell.getSize());
                dummy.setDistance(
                    xwalk.constants.Constants.FIRST_ATOM_IS_SOLVENT_INACCESSIBLE
                                 );
                path.add(dummy);
                paths.add(i, path);
            } else if (atom1isAccessible
                       &&
                       !atoms2areAccessible[i]
                       &&
                       atom2isAccessible) {
                paths.get(i).get(
                        BreadthFirstSearch.CELL_NO_OF_TARGET_CELL_IN_PATH
                                 ).setDistance(
                   xwalk.constants.Constants.SECOND_ATOM_IS_SOLVENT_INACCESSIBLE
                                 );
            } else if (atom1isAccessible
                    &&
                    !atoms2areAccessible[i]
                    &&
                    !successful) {
                Path path = new Path();
                GridCell dummy = new GridCell(atom1.getXYZ(),
                                              GridCell.getSize());
                dummy.setDistance(
                   xwalk.constants.Constants.SECOND_ATOM_IS_SOLVENT_INACCESSIBLE
                              );
                path.add(dummy);
                paths.add(i, path);
            } else if (!atom1isAccessible && !atoms2areAccessible[i]) {
                // if both are inaccessible
                Path path = new Path();
                GridCell dummy = new GridCell(atom1.getXYZ(),
                                              GridCell.getSize());
                dummy.setDistance(
                   xwalk.constants.Constants.BOTH_ATOMS_ARE_SOLVENT_INACCESSIBLE
                                 );
                path.add(dummy);
                paths.add(i, path);
            } else if (atom1isAccessible
                       &&
                       atoms2areAccessible[i]
                       &&
                       successful
                       &&
                       paths.get(i).get(
                               BreadthFirstSearch.CELL_NO_OF_TARGET_CELL_IN_PATH
                                        ).getDistance()
                                           == Constants.DEFAULT_GRID_DISTANCE) {
                // if non of the above apply, than the distance exceeds max
                paths.get(i).get(
                        BreadthFirstSearch.CELL_NO_OF_TARGET_CELL_IN_PATH
                                 ).setDistance(
                             xwalk.constants.Constants.NON_CONFORMING_CROSS_LINK
                                 );
            }
        }
    return paths;
    }
    //--------------------------------------------------------------------------
    /**
     * Removes all crossLink object in a CrossLinkList that have a SASD < -0.6,
     * which causes the removal of any non-conforming crossLink object from the
     * CrossLinkList.
     * @param list
     *        - List of crossLink object, which will be scanned for SASD < -0.6.
     */
    private static void cleanCrossLinkList (CrossLinkList list) {
        CrossLinkList toBremoved = new CrossLinkList();
        for (CrossLink xl : list) {
            if (xl.getSolventPathDistance() < -0.6) {
                toBremoved.add(xl);
            }
        }
        list.removeAll(toBremoved);
    }
    //--------------------------------------------------------------------------


}

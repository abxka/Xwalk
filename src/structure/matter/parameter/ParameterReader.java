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

package structure.matter.parameter;

import java.io.File;
import java.io.IOException;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.io.ReadFile;


/**
 * Class for reading in and handling parameter files, e.g. atom radius parameter
 * files.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class ParameterReader {

    /**
     * Holds van der Waals radii from the various supported force fields.
     */
    private static Hashtable < Element, Float > vdWradii =
                                           new Hashtable < Element, Float >();
    /**
     * Holds atomic XlogP values for all standard PDB amino acid atom types.
     */
    private static Hashtable <AminoAcidType, Hashtable <AtomType, Float >> xlogP
               = new Hashtable <AminoAcidType, Hashtable <AtomType, Float >>();

    /**
     * Holds the probabilities for Euclidean cross-link distances.
     */
    private static Hashtable < Float, Float > euc_probablities =
                                           new Hashtable < Float, Float >();

    /**
     * Holds the probabilities for SAS cross-link distances.
     */
    private static Hashtable < Float, Float > sasd_probablities =
                                           new Hashtable < Float, Float >();

    /**
     * Location of parameter molecular mechanics parameter files.
     */
    private static final String PARAMETER_MM_DIR = "mm"
                                                 + File.separatorChar
                                                 + "parameter"
                                                 + File.separatorChar;

    /**
     * Location of parameter distance probability files.
     */
    private static final String PARAMETER_DP_DIR = "xwalk"
                                                 + File.separatorChar
                                                 + "parameter"
                                                 + File.separatorChar;
    /**
     * MMRR94 radius identifier.
     */
    private static final String MMFF94_FILENAME = "MMFF94_radii.txt";
    /**
     * PARSE radius identifier.
     */
    private static final String PARSE_FILENAME = "PARSE_radii.txt";
    /**
     * SURFNET radius identifier.
     */
    private static final String SURFNET_FILENAME = "SURFNET_radii.txt";
    /**
     * RASMOL radius identifier.
     */
    private static final String RASMOL_FILENAME = "RASMOL_radii.txt";
    /**
     * CHARMM radius identifier.
     */
    private static final String CHARMM_FILENAME = "CHARMM_radii.txt";
    /**
     * XLOGP radius identifier.
     */
    private static final String XLOGP_FILENAME = "atomic_XlogP.txt";
    /**
     * Solvent Accessible Surface Distance identifier.
     */
    private static final String SASD_PROB_FILENAME = "sasd_prob.txt";
    /**
     * Euclidean Distance identifier.
     */
    private static final String EUC_PROB_FILENAME = "euc_prob.txt";

    /**
     * Column number of element name in two column parameter files. Should be 0.
     */
    private static int elementNameColumn = 0;
    /**
     * Column number of van der Waals radius in two column parameter files.
     * Should be 1.
     */
    private static int vdWradiusColumn = 1;
    /**
     * Column number of amino acid name in three column parameter files.
     * Should be 0.
     */
    private static int aminoAcidNameColumn = 0;
    /**
     * Column number of atom name in three column parameter files. Should be 1.
     */
    private static int atomNameColumn = 1;
    /**
     * Column number of xlogP value in a three column parameter files.
     * Should be 2.
     */
    private static int xlogPcolumn = 2;
    /**
     * Column number of distance bin in two column parameter files. Should be 0.
     */
    private static int distanceBinColumn = 0;
    /**
     * Column number of probability in two column parameter files.
     * Should be 1.
     */
    private static int probabilityColumn = 1;

    //--------------------------------------------------------------------------
    /**
     * Sets the parameter.
     * @param parameter
     *        - One of the supported ParameterSet object in Xwalk.
     * @throws IOException if an error occurs while reading the parameter file.
     */
    public static void setParameterReader(
                                         final Constants.ParameterSets parameter
                                         )
                                         throws IOException {
        if (parameter == Constants.ParameterSets.RASMOL) {
            ParameterReader.readParameterSet(ParameterReader.PARAMETER_MM_DIR
                                + ParameterReader.RASMOL_FILENAME);
        }
        if (parameter == Constants.ParameterSets.SURFNET) {
            ParameterReader.readParameterSet(ParameterReader.PARAMETER_MM_DIR
                                + ParameterReader.SURFNET_FILENAME);
        }
        if (parameter == Constants.ParameterSets.MMFF94) {
            ParameterReader.readParameterSet(ParameterReader.PARAMETER_MM_DIR
                                + ParameterReader.MMFF94_FILENAME);
        }
        if (parameter == Constants.ParameterSets.PARSE) {
            ParameterReader.readParameterSet(ParameterReader.PARAMETER_MM_DIR
                                + ParameterReader.PARSE_FILENAME);
        }
        if (parameter == Constants.ParameterSets.CHARMM) {
            ParameterReader.readParameterSet(ParameterReader.PARAMETER_MM_DIR
                                + ParameterReader.CHARMM_FILENAME);
        }
        if (parameter == Constants.ParameterSets.XLOGP) {
            ParameterReader.readParameterSet(ParameterReader.PARAMETER_MM_DIR
                                + ParameterReader.XLOGP_FILENAME);
        }
        if (parameter == Constants.ParameterSets.SASD_PROB) {
            ParameterReader.readParameterSet(ParameterReader.PARAMETER_DP_DIR
                                + ParameterReader.SASD_PROB_FILENAME);
        }
        if (parameter == Constants.ParameterSets.EUC_PROB) {
            ParameterReader.readParameterSet(ParameterReader.PARAMETER_DP_DIR
                                + ParameterReader.EUC_PROB_FILENAME);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Reads in the information from a parameter file into a Hashtable.
     * @param parameterFileName
     *        - String object holding the path to the desired parameter file.
     * @throws IOException if an error occurs while reading the parameter file.
     */
    private static void readParameterSet(final String parameterFileName)
                                                            throws IOException {
        ReadFile read = new ReadFile(parameterFileName);
        for (String line : read) {
            if (!line.startsWith("#")) {
                String[] column = line.split("\t");
                if (column.length == 2) {
                    if (parameterFileName.indexOf(
                                               ParameterReader.EUC_PROB_FILENAME
                                                 ) != -1) {
                        float distBin = Float.parseFloat(
                                       column[ParameterReader.distanceBinColumn]
                                                           );
                        float prob = Float.parseFloat(
                                       column[ParameterReader.probabilityColumn]
                                                        );
                        ParameterReader.euc_probablities.put(distBin, prob);
                    } else if (parameterFileName.indexOf(
                                              ParameterReader.SASD_PROB_FILENAME
                                                        ) != -1) {
                        float distBin = Float.parseFloat(
                                       column[ParameterReader.distanceBinColumn]
                                                        );
                        float prob = Float.parseFloat(
                                       column[ParameterReader.probabilityColumn]
                                                     );
                        ParameterReader.sasd_probablities.put(distBin, prob);
                    } else {
                        String elementName =
                                      column[ParameterReader.elementNameColumn];
                        float radius = Float.parseFloat(
                                         column[ParameterReader.vdWradiusColumn]
                                                      );
                        for (Element e : Element.values()) {
                            if (e.toString().equals(elementName)) {
                                ParameterReader.vdWradii.put(e, radius);
                            }
                        }
                    }
                } else if (column.length == 3) {
                    String aminoAcidName =
                                    column[ParameterReader.aminoAcidNameColumn];
                    String atomName = column[ParameterReader.atomNameColumn];
                    float xlogPvalue = Float.parseFloat(
                                             column[ParameterReader.xlogPcolumn]
                                                          );

                    for (AminoAcidType aat : AminoAcidType.values()) {
                        if (aat.getThreeLetterCode().equals(aminoAcidName)) {
                            for (AtomType at : AtomType.values()) {
                                if (at.getAbbreviation().equals(atomName)) {
                                    if (ParameterReader.xlogP.get(aat)
                                        ==
                                        null) {
                                        Hashtable<AtomType, Float> atomXlogP =
                                              new Hashtable<AtomType, Float>();
                                        atomXlogP.put(at, xlogPvalue);
                                        ParameterReader.xlogP.put(aat,
                                                                  atomXlogP);
                                    } else {
                                        ParameterReader.xlogP.get(aat).put(
                                                                      at,
                                                                      xlogPvalue
                                                                          );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Returns the set of atom van der Waals radius parameter.
     * @return Hashtable with Element keys and float elements as radii.
     */
    public static final Hashtable < Element, Float >
                                                    getVdwRadiusParameterSet() {
        return ParameterReader.vdWradii;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the set of atomic XlogP values for all protein amino acid atoms.
     * @return Hashtable with XlogP values for all protein amino acid atoms.
     */
    public static final Hashtable <AminoAcidType, Hashtable <AtomType, Float >>
                                                        getXlogPparameterSet() {
        return ParameterReader.xlogP;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the set of probabilities for a certrain Euclidean distance bin.
     * @return Hashtable with distance keys and float elements as
     * probabilities.
     */
    public static final Hashtable < Float, Float >
                                          getEuclideanDistanceProbabilitySet() {
        return ParameterReader.euc_probablities;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the set of probabilities for a certrain SAS distance bin.
     * @return Hashtable with distance keys and float elements as
     * probabilities.
     */
    public static final Hashtable < Float, Float >
                                          getSASdistanceProbabilitySet() {
        return ParameterReader.sasd_probablities;
    }
}

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
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.constants.Constants.ParameterSets;
import structure.matter.parameter.ParameterReader;

import xwalk.io.CommandlineArguments;

/**
 * Class that stores and handles parameters for calculating cross-links with
 * Xwalk.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class CrossLinkParameter {

    /**
     * Hashtable that stores all parameters necessary to calculate cross-links.
     */
    private static Hashtable < Parameter, String > param =
                                          new Hashtable < Parameter, String >();
    //--------------------------------------------------------------------------
    /**
     * Supported cross-link parameters, that effect the behavior of the Xwalk
     * application.
     */
    public enum Parameter { ALTERNATIVE_LOCATION1, ALTERNATIVE_LOCATION2,
                            AMINO_ACID_RESIDUE_NUMBER1,
                            AMINO_ACID_RESIDUE_NUMBER2,
                            AMINO_ACID_RESIDUE_NAME1,
                            AMINO_ACID_RESIDUE_NAME2,
                            ATOM_TYPE1,
                            ATOM_TYPE2,
                            CHAIN_ID1,
                            CHAIN_ID2,
                            DISTANCE_FILE_PATH,
                            DO_KEEP_FILENAME,
                            DO_FORCE_OUTPUT,
                            DO_VERBOSE_OUTPUT,
                            DO_GRID_OUTPUT,
                            DO_TRYPSIN_DIGEST,
                            DO_EXPASY_RULE,
                            DO_BACKBONE_READ,
                            DO_REMOVE_SIDECHAINS,
                            DO_INTRAMOLECULAR_DISTANCE,
                            DO_INTERMOLECULAR_DISTANCE,
                            DO_MONO_CROSSLINK,
                            DO_PROBABILITY,
                            DO_SOLVENT_PATH_DISTANCE,
                            DO_BFACTOR,
                            DO_PYMOL_OUTPUT,
                            INFILE_PATH,
                            IS_HOMOMERIC,
                            GRID_CELL_SIZE,
                            MAXIMUM_DISTANCE,
                            MINIMUM_SOLVENT_ACCESSIBILITY_RATIO,
                            OUTFILE_PATH,
                            SOLVENT_RADIUS,
    };
    //--------------------------------------------------------------------------
    /**
     * Constructor that sets all relevant Xwalk parameters from the user given
     * commandline arguments.
     * @param arg
     *        - CommandlineArguments object holding all arguments that were
     *          un/set by the user.
     */
    public CrossLinkParameter(final CommandlineArguments arg) {
        this.setParameter(Parameter.ALTERNATIVE_LOCATION1,
                          arg.getAlternativeLocation1Argument());
        this.setParameter(Parameter.ALTERNATIVE_LOCATION2,
                          arg.getAlternativeLocation2Argument());
        this.setParameter(Parameter.AMINO_ACID_RESIDUE_NAME1,
                          arg.getAminoAcidName1Argument());
        this.setParameter(Parameter.AMINO_ACID_RESIDUE_NAME2,
                          arg.getAminoAcidName2Argument());
        this.setParameter(Parameter.AMINO_ACID_RESIDUE_NUMBER1,
                          arg.getAminoAcidNumber1Argument());
        this.setParameter(Parameter.AMINO_ACID_RESIDUE_NUMBER2,
                          arg.getAminoAcidNumber2Argument());
        this.setParameter(Parameter.ATOM_TYPE1, arg.getAtomType1Argument());
        this.setParameter(Parameter.ATOM_TYPE2, arg.getAtomType2Argument());
        this.setParameter(Parameter.CHAIN_ID1, arg.getChainIds1Argument());
        this.setParameter(Parameter.CHAIN_ID2, arg.getChainIds2Argument());
        this.setParameter(Parameter.DO_INTERMOLECULAR_DISTANCE,
                          Boolean.toString(arg.isInterMolecularDistanceSet()));
        this.setParameter(Parameter.DO_INTRAMOLECULAR_DISTANCE,
                           Boolean.toString(arg.isIntraMolecularDistanceSet()));
        this.setParameter(Parameter.DO_MONO_CROSSLINK,
                                    Boolean.toString(arg.isMonoCrossLinkSet()));
        this.setParameter(Parameter.DO_PYMOL_OUTPUT, Boolean.toString(
                                                       arg.isPymolOutputSet()));
        this.setParameter(Parameter.DO_SOLVENT_PATH_DISTANCE, Boolean.toString(
                                    arg.isSolventPathDistanceCalculationSet()));
        this.setParameter(Parameter.DO_PROBABILITY, Boolean.toString(
                                               arg.isProbabilityArgumentSet()));
        this.setParameter(Parameter.GRID_CELL_SIZE, Double.toString(
                                                arg.getGridCellSizeArgument()));
        this.setParameter(Parameter.INFILE_PATH, arg.getInfileArgument());
        this.setParameter(Parameter.DO_BACKBONE_READ, Boolean.toString(
                                              arg.isBackboneOnlyArgumentSet()));
        this.setParameter(Parameter.DO_REMOVE_SIDECHAINS, Boolean.toString(
                                           arg.isRemoveSideChainArgumentSet()));
        this.setParameter(Parameter.DISTANCE_FILE_PATH,
                                               arg.getDistanceInfileArgument());
        this.setParameter(Parameter.DO_KEEP_FILENAME, Boolean.toString(
                                                  arg.isKeepNameArgumentSet()));
        this.setParameter(Parameter.OUTFILE_PATH, arg.getOutfileArgument());
        this.setParameter(Parameter.IS_HOMOMERIC, Boolean.toString(
                                                         arg.isHomomericSet()));
        this.setParameter(Parameter.MAXIMUM_DISTANCE, Double.toString(
                                             arg.getMaximumDistanceArgument()));
        this.setParameter(Parameter.SOLVENT_RADIUS, Double.toString(
                                               arg.getSolventRadiusArgument()));
        this.setParameter(Parameter.DO_FORCE_OUTPUT, Boolean.toString(
                                                       arg.isForceOutputSet()));
        this.setParameter(Parameter.DO_VERBOSE_OUTPUT, Boolean.toString(
                                                     arg.isVerboseOutputSet()));
        this.setParameter(Parameter.DO_GRID_OUTPUT, Boolean.toString(
                                                        arg.isGridOutputSet()));
        this.setParameter(Parameter.DO_TRYPSIN_DIGEST, Boolean.toString(
                                                arg.isTrypsinateArgumentSet()));
        this.setParameter(Parameter.DO_EXPASY_RULE, Boolean.toString(
                                                    arg.isExpasyArgumentSet()));
        this.setParameter(Parameter.DO_BFACTOR, Boolean.toString(
                                                   arg.isBfactorArgumentSet()));
        try {
            ParameterReader.setParameterReader(ParameterSets.SURFNET);
            ParameterReader.setParameterReader(ParameterSets.XLOGP);
            ParameterReader.setParameterReader(ParameterSets.SASD_PROB);
            ParameterReader.setParameterReader(ParameterSets.EUC_PROB);
        } catch (IOException e) {
            System.err.println("ERROR while reading parameter files.");
        }

        if (arg.isVerboseOutputSet()) {
            this.output();
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Sets a cross-link param.
     * @param xlParameter
     *        Parameter object that holds a list of all relevant cross-link
     *        parameters.
     * @param xlValue
     *           String object that holds the value of the cross-link param
     */
    public final void setParameter(final Parameter xlParameter,
                                   final String xlValue) {
        param.put(xlParameter, xlValue);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the value of a cross-link param.
     * @param parameter
     *        Parameter object that holds a list of all relevant cross-link
     *        parameters.
     * @return String object that holds the value of the cross-link param
     */
    public static final String getParameter(final Parameter parameter) {
        return param.get(parameter);
    }
    //--------------------------------------------------------------------------
    /**
     * Outputs all param via the STDERR stream to the terminal.
     */
    public final void output() {
        ArrayList < Parameter > pars = new ArrayList < Parameter >();
        for (Parameter par : param.keySet()) {
            pars.add(par);
        }
        Collections.sort(pars, new Comparator<Parameter>() {
            public int compare(final Parameter p1, final Parameter p2) {
                return p1.toString().compareTo(p2.toString());
            }
        });

        System.err.println("List of all argument values:");
        for (Parameter par : pars) {
            System.err.println(par.toString() + ": "
                             + param.get(par));
        }
        System.err.print(Constants.LINE_SEPERATOR);
    }
}

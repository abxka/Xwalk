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

package mm.constants;

import structure.matter.parameter.AminoAcidType;
/**
 * Various constants that are used by various classes.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Constants {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected Constants() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    // CONSTANT NUMERIC VALUES
    //--------------------------------------------------------------------------
    /**
     * In Coloumb.
     */
    public static final double ELECTRON_CHARGE = 1.6021765314E-19;
    //--------------------------------------------------------------------------
    /**
     * In metric gram.
     */
    public static final double ELECTRON_MASS = 9.109382616E-28;
    //--------------------------------------------------------------------------
    /**
     * In Joules times s.
     */
    public static final double DIRACS_CONSTANT = 1.05457159682E-34;
    //--------------------------------------------------------------------------
    /**
     * XlogP value of glycine residue, which serves as a reference for
     * calculating the XlogP value of all amino acids. Default -1.128564.
     */
    public static final double GLYCINE_LOGP = -1.128564;
    //--------------------------------------------------------------------------
    /**
     * XlogP value added to organic ions like COO- in order to correct for their
     * charge and thus higher solvation tendency. Default -1.083.
     */
    public static final float ORGANIC_ION_CORRECTION_FACTOR = -1.083f;
    //--------------------------------------------------------------------------
    /**
     * XlogP value given to all metal ions. Value is more or less arbitrary, but
     * should be lower than any charged amino acid. Default -3.
     */
    public static final float METAL_XLOGP = -3;
    //--------------------------------------------------------------------------
    /**
     * Radius above which any physicochemical property calculations will be
     * omitted. Default 9.0 Angstroem.
     */
    public static final float PHYSICOCHEMICAL_INFLUENCE_RADIUS = 9.0f;
    //--------------------------------------------------------------------------
    /**
     * Distance for non-bonded contact between protein and ligand atoms
     * as defined by HBPLUS/LIGPLOT (Wallace, A., Laskowski, R. & Thornton, J.
     * LIGPLOT: a program to generate schematic diagrams of protein-ligand
     * interactions. Protein Eng 8, 127-134 (1995).). Default is 3.9.
     */
    public static final float NON_BONDED_CONTACT_DISTANCE = 3.9f;
    //--------------------------------------------------------------------------
    /**
     * Dielectric constant of a solvent. Default 78.0.
     */
    public static final double SOLVENT_DIELECTRIC_CONSTANT = 78.0;
    //--------------------------------------------------------------------------
    /**
     * Dielectric constant of a protein. Default 4.0.
     */
    public static final double PROTEIN_DIELECTRIC_CONSTANT = 4.0;
    //--------------------------------------------------------------------------
    /**
     * Conversion of 25 degrees Celsius to Fahrenheit measure. Default 298.15 K.
     */
    public static final double TWENTYFIVE_CELSIUS_IN_FAHRENHEIT = 298.15;

    //--------------------------------------------------------------------------
    // CONSTANT STRING VALUES
    //--------------------------------------------------------------------------
    /**
     * Text string identifying XlogP energy values in a hash.
     */
    public static final String XLOGP_ENERGY_PROPERTY_TAG = "XlogPenergy";
    /**
     * Text string identifying XlogP potential energy values in a hash.
     */
    public static final String XLOGP_POTENTIAL_ENERGY_PROPERTY_TAG =
                                                               "XlogPpotEnergy";

    //--------------------------------------------------------------------------
    // ENUM TYPES
    //--------------------------------------------------------------------------

    /**
     * Values taken from <a href="http://www.bioinformatics.cm-uj.krakow.pl/
     * reveal/">model</a>.
     */
    public enum PhysicoChemicalProperty {

        // hydrophobic amino acids
        GLYCINE(AminoAcidType.GLYCINE, 0.550, 0.00,  0.000),
        ALANINE(AminoAcidType.ALANINE, 0.572, 0.31, 0.570),
        VALINE(AminoAcidType.VALINE, 0.811, 1.22, 1.226),
        LEUCINE(AminoAcidType.LEUCINE ,0.783, 1.70, 2.191),
        ISOLEUCINE(AminoAcidType.ISOLEUCINE, 0.883, 1.80, 1.903),
        PHENYLALANINE(AminoAcidType.PHENYLALANINE, 0.906, 1.79, 1.970),
        PROLINE(AminoAcidType.PROLINE, 0.300, 0.72, 0.639),
        TRYPTOPHANE(AminoAcidType.TRYPTOPHANE, 0.856, 2.25, 2.892),
        // amino acids containing a nitrogen atom
        ASPARAGINE(AminoAcidType.ASPARAGINE, 0.278, -0.6, -1.961),
        GLUTAMINE(AminoAcidType.GLUTAMINE, 0.250, -0.22, -1.560),
        // amino acids containing a sulphur atom
        METHIONINE(AminoAcidType.METHIONINE, 0.828, 1.23, 1.558),
        CYSTEINE(AminoAcidType.CYSTEINE, 1.000, 1.54, 1.783),
        // amino acids containing a oxygen atom
        THREONINE(AminoAcidType.THREONINE, 0.478, 0.26, -0.583),
        TYROSINE(AminoAcidType.TYROSINE, 0.700, 0.96, 1.353),
        SERINE(AminoAcidType.SERINE, 0.422, -0.04, -1.099),
        // positively charged amino acids
        ARGININE(AminoAcidType.ARGININE, 0.272, -1.01, -0.873),
        LYSINE (AminoAcidType.LYSINE, 0.000, -0.99, -0.022),
        HISTIDINE(AminoAcidType.HISTIDINE, 0.628, 0.13, -0.009),
        // negatively charged amino acids
        GLUTAMIC_ACID(AminoAcidType.GLUTAMIC_ACID, 0.083, -0.64, -2.701),
        ASPARTIC_ACID(AminoAcidType.ASPARTIC_ACID, 0.167, -0.77, -3.102);

        /**
         * Fuzzy Oil Drop hydrophobicity value.
         */
        private double fuzzyOilDropHydrophobicity;
        /**
         * Fauchere and Pliska hydrophobicity value.
         */
        private double faucherePliskaHydrophobicity;
        /**
         * XlogP hydrophobicity value.
         */
        private double xlogPhydrophobicity;
        /**
         * Type of amino acid.
         */
        private AminoAcidType aminoAcidType;

        /**
         * Constructor.
         * @param type
         *        AminoAcidType object indicating which type of amino acid this
         *        aminoacid is.
         * @param fuzzyOilDrop
         *        double value representing the amino acids hydrophobicity
         *        value as defined by the
         *        <a href="http://www.bioinformatics.cm-uj.krakow.pl/reveal/">
         *        Fuzzy Oil Drop method</a>.
         * @param faucherePliska
         *        double value representing the amino acids hydrophobicity
         *        value as defined by Fauchere and Pliska, European Journal
         *        of Medicinal Chemistry (1983) vol. 18 (4) pp. 369-375.
         * @param xlogP
         *        double value representing the amino acids hydrophobicity
         *        value as defined by the XlogP method give in Kahraman et al.,
         *        Proteins, (2010).
         */
        PhysicoChemicalProperty(final AminoAcidType type,
                                final double fuzzyOilDrop,
                                final double faucherePliska,
                                final double xlogP) {
            this.fuzzyOilDropHydrophobicity = fuzzyOilDrop;
            this.faucherePliskaHydrophobicity = faucherePliska;
            this.xlogPhydrophobicity = xlogP;
            this.aminoAcidType = type;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * <a href="http://www.poissonboltzmann.org/pdb2pqr/
     * frequently-asked-questions/what-force-fields-or-parameter-sets-are-
     * available">Force fields</a> supported by PDB2PQR.
     */
    public enum ForceField {

        /**
         *  PARSE.
         */
        PARSE,
        /**
         *  AMBER 94.
         */
        AMBER,
        /**
         *  CHARMM 27.
         */
        CHARMM,
        /**
         *  USER-DEFINED.
         */
        ROSETTA
    }
}

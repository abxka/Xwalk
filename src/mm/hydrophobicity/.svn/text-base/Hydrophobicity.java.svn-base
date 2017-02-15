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

package mm.hydrophobicity;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.math.Mathematics;
import structure.math.Point3f;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.parameter.AminoAcidType;
import structure.matter.parameter.AtomType;
import structure.matter.parameter.ParameterReader;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;

/**
 * Class providing functionalities to calculate hydrophobic properties on
 * molecules and molecular assemblies.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Hydrophobicity {

    //--------------------------------------------------------------------------
    /**
     * PolyPeptide object holding all amino acids to which atomic XlogP
     * values will be calculated.
     */
    private PolyPeptideList polyPeptideComplex;

    //--------------------------------------------------------------------------
    /**
     * Constructor. Reads in the XlogP parameter file and assigns the values
     * to the protein atoms.
     * @param complex
     *        Protein complex object holding all amino acids to which atomic
     *        XlogP values will be calculated.
     */
    public Hydrophobicity(final PolyPeptideList complex) {
//        if (ParameterReader.getXlogPparameterSet() == null) {
            try {
                ParameterReader.setParameterReader(
                                                   Constants.ParameterSets.XLOGP
                                                  );
            } catch (IOException e) {
                System.err.print(e.getMessage() + Constants.LINE_SEPERATOR
                              + "ERROR: While reading the XlogP parameter file"
                              + Constants.LINE_SEPERATOR);
            }
//        }
        this.polyPeptideComplex = complex;
        this.setAtomicXlogP();
    }
    //--------------------------------------------------------------------------
    /**
     * Sets to all amino acids atoms in a protein their associated XlogP values.
     * @return float value representing the sum of XlogP values for the
     *         protein complex.
     */
    private float setAtomicXlogP() {

        Hashtable <AminoAcidType, Hashtable <AtomType, Float>> xlogPs =
                                         ParameterReader.getXlogPparameterSet();
        float sum = 0;
        for (PolyPeptide polyPeptide : this.polyPeptideComplex) {
            for (AminoAcid aa : polyPeptide) {
                for (Atom atom : aa.getAllAtoms()) {
                    if (!atom.getElement().getSymbol().equals("H")) {
                        Hashtable<AtomType, Float> atomicXlogPs =
                                                       xlogPs.get(aa.getType());
                        float xlogP = atomicXlogPs.get(atom.getType());
                        atom.setXlogP(xlogP);
                        sum += xlogP;
                    } else {
                        // hydrogen atoms have always a XlogP value of 0.
                        atom.setXlogP(0);
                    }
                }
            }
        }
        return sum;
    }

    //--------------------------------------------------------------------------
    /**
     * Maps the hydrophobicity property of the protein on a list of atom. Note,
     * that only atoms within
     * mm.constants.Constants.PHYSICOCHEMICAL_INFLUENCE_RADIUS will be included
     * in the potential calculation.
     * @param sampleAtoms
     *        List of atoms on which potential will be calculated.
     */
    public final void mapHydrophobicity(final AtomList sampleAtoms) {

        float max = mm.constants.Constants.PHYSICOCHEMICAL_INFLUENCE_RADIUS;

        for (Atom sampleAtom : sampleAtoms) {
            for (Atom atom : this.polyPeptideComplex.getAllAtoms()) {
                float dist = Mathematics.distance(atom.getXYZ(),
                                                  sampleAtom.getXYZ());
                if (dist <= max) {
                    float h = atom.getXlogP();
                    double s = Mathematics.sigmoidFunction(dist, max);
                    sampleAtom.setHes((float) (sampleAtom.getHes() + (h * s)));
                }
            }
        }
    }
}

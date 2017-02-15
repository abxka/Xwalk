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

package structure.matter.hetgroups;

import structure.constants.Constants;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.Molecule;

/**
 * Class representing any hetgroups/small molecules.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 *
 */
public class SmallMolecule extends Molecule {
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;

    //--------------------------------------------------------------------------
    /**
     * Constructor. Checks whether all atoms really belong to a single small
     * molecule. Furthermore assigns the amino acid to one of various amino acid
     * types.
     * @param atoms
     *        - AtomList object holding the atomic coordinates of this amino
     *        acid.
     */
    public SmallMolecule(final AtomList atoms) {
        super(atoms);
        // first check whether molecule really consists of a single amino acid.
        try {
            Atom preAtom = this.getAtom(0);
            for (Atom atom : this.getAllAtoms()) {
                if (!atom.getResidueName().equals(preAtom.getResidueName())
                    &&
                    atom.getResidueNumber() != preAtom.getResidueNumber()) {
                    throw new Exception("Data contains more than one residue "
                                      + "information.");
                }
                if (!atom.getFlag().equals("HETATM")) {
                    throw new Exception("Data has non HETATM entries.");
                }
            }
        } catch (Exception e) {
            System.err.print("ERROR: in reading in SmallMolecule information. "
                           + e.getMessage() + Constants.LINE_SEPERATOR);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Creates a copy of this SmallMolecule object.
     * @return Copy of this SmallMolecule object.
     */
    public final SmallMolecule copy() {
        AtomList atomsCopy = new AtomList();
        for (Atom atom : this.getAllAtoms()) {
            atomsCopy.add(atom.copy());
        }
        SmallMolecule copy = new SmallMolecule(atomsCopy);
    return copy;
    }
}

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

package structure.matter;

import java.util.ArrayList;

import structure.constants.Constants;
import structure.matter.parameter.Element;

/**
 * Abstract class representing chemical molecules.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public abstract class Molecule {
    /**
     * List of all bonds within the molecule.
     */
    private ArrayList < Bond > bonds = new ArrayList < Bond >();
    //--------------------------------------------------------------------------
    /**
     * List of all atoms in the molecule.
     */
    private AtomList atoms = new AtomList();
    //--------------------------------------------------------------------------
    /**
     * Total solvent accessibility of this residue.
     */
    private float totalSolventAccessibility;
    //--------------------------------------------------------------------------
    /**
     * Relative solvent accessibility of this residue.
     */
    private float relativeSolventAccessibility;
    //--------------------------------------------------------------------------
    /**
     * Connects all atoms in this molecule with Bond objects that are non-metal
     * ions and have a distance smaller than the van der Waals radii of both
     * atoms + PDB uncertainty factor.
     * @param atomList
     *        - List of Atoms that are building up the molecule.
     */
    public Molecule(final AtomList atomList) {
        this.atoms = atomList;
        for (Bond bond : MatterUtilities.calculateBonds(atomList)) {
            this.bonds.add(bond);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all bonds found within this molecule.
     * @return ArrayList object with elements of type Bond.
     */
    public final ArrayList < Bond > getBonds() {
        return this.bonds;
    }
    //--------------------------------------------------------------------------

    /**
     * Determines the element type of the atom by comparing the atom name in a
     * PDB file with the periodic table element symbol of the same atom. If the
     * atom is of type calcium Ca, than the molecule argument must either be
     * {@code null}, or has 0 or 1 element. Furthermore, if special atom types
     * such as hydrogens bound to oxygens shall be recognized, the hydrogens and
     * their covalent bound atoms must have instances within the bond parameter
     * of this molecule object. Element types that are dependent on aromaticity
     * need first an aromaticity assignment.
     */
    public final void setElements() {
        for (Atom atom : this.atoms) {
            String atomName = atom.getName().trim();
            if (atom.getFlag().equals("HETATM")) {
                // METAL ELEMENTS
                // Calcium is the only element that could truly be mixed up with
                // an other element, namely carbon. Therefore, the
                if (this != null) {
                    if (this.atoms.size() == 0 || this.atoms.size() == 1) {
                        if (atomName.equals(Element.CALCIUM.getSymbol())) {
                            atom.setElement(Element.CALCIUM);
                        }
                    }
                } else if (atomName.equals(Element.LITHIUM.getSymbol())) {
                    atom.setElement(Element.LITHIUM);
                } else if (atomName.equals(Element.BERYLLIUM.getSymbol())) {
                    atom.setElement(Element.BERYLLIUM);
                } else if (atomName.equals(Element.SODIUM.getSymbol())) {
                    atom.setElement(Element.SODIUM);
                } else if (atomName.equals(Element.MAGNESIUM.getSymbol())) {
                    atom.setElement(Element.MAGNESIUM);
                } else if (atomName.equals(Element.POTASSIUM.getSymbol())) {
                    atom.setElement(Element.POTASSIUM);
                } else if (atomName.equals(Element.ALUMINUM.getSymbol())) {
                    atom.setElement(Element.ALUMINUM);
                } else if (atomName.equals(Element.IRON.getSymbol())) {
                    atom.setElement(Element.IRON);
                } else if (atomName.equals(Element.ZINC.getSymbol())) {
                    atom.setElement(Element.ZINC);
                } else if (atomName.equals(Element.MANGANESE.getSymbol())) {
                    atom.setElement(Element.MANGANESE);
                } else if (atomName.equals(Element.LUTETIUM.getSymbol())) {
                    atom.setElement(Element.LUTETIUM);
                } else if (atomName.equals(Element.URANIUM.getSymbol())) {
                    atom.setElement(Element.URANIUM);
                } else if (atomName.equals(Element.MERCURY.getSymbol())) {
                    atom.setElement(Element.MERCURY);
                } else if (atomName.equals(Element.CADMIUM.getSymbol())) {
                    atom.setElement(Element.CADMIUM);

                // NON-METAL ELEMENTS
                } else if (atomName.equals(Element.IODINE.getSymbol())) {
                    atom.setElement(Element.IODINE);
                } else if (atomName.equals(Element.SELENIUM.getSymbol())) {
                    atom.setElement(Element.SELENIUM);
                } else if (atomName.equals(Element.FLUORINE.getSymbol())) {
                    atom.setElement(Element.FLUORINE);
                } else if (atomName.equals(Element.CHLORINE.getSymbol())) {
                    atom.setElement(Element.CHLORINE);
                } else if (atomName.equals(Element.BROMINE.getSymbol())) {
                    atom.setElement(Element.BROMINE);

                // METALLOIDS
                } else if (atomName.equals(Element.SILICON.getSymbol())) {
                    atom.setElement(Element.SILICON);
                }
            }
            // if non of the above names have been found than
            // it is likely that the atom is of organic type.
            if (atomName.startsWith(Element.CARBON.getSymbol())) {
                atom.setElement(Element.CARBON);
            } else if (atomName.startsWith(Element.NITROGEN.getSymbol())) {
                atom.setElement(Element.NITROGEN);
            } else if (atomName.startsWith(Element.OXYGEN.getSymbol())) {
                atom.setElement(Element.OXYGEN);
            } else if (atomName.startsWith(Element.SULPHUR.getSymbol())) {
                atom.setElement(Element.SULPHUR);
            } else if (atomName.startsWith(Element.HYDROGEN.getSymbol())) {
                atom.setElement(Element.HYDROGEN);
            } else if (atomName.matches("\\d+"
                                      + Element.HYDROGEN.getSymbol()
                                      + ".*")) {
                atom.setElement(Element.HYDROGEN);
            } else if (atomName.startsWith(Element.PHOSPHOR.getSymbol())) {
                atom.setElement(Element.PHOSPHOR);
            }
/*            // treat special cases of atom types
            for (Bond bond : this.getBonds()) {
                if (atom.getElement() == Element.HYDROGEN) {
                    if (bond.getPreAtom().equals(atom)
                        ||
                        bond.getPostAtom().equals(atom)) {
                        // determine covalent bound atom
                        Atom covalent;
                        if (bond.getPreAtom().equals(atom)) {
                            covalent = bond.getPostAtom();
                        } else {
                            covalent = bond.getPreAtom();
                        }
                        // assign hydrogen atomtype depending on covalent bound
                        // atom.
                        if (covalent.getName().trim().startsWith(
                                                  Element.CARBON.getSymbol()
                                                                )
                           ) {
                            atom.setElement(Element.HYDROGEN_BOUND_TO_CARBON);
                        } else if (covalent.getName().trim().startsWith(
                                                      Element.OXYGEN.getSymbol()
                                                                       )
                                   ||
                                   covalent.getName().trim().startsWith(
                                                    Element.NITROGEN.getSymbol()
                                                                       )
                                  ) {
                            atom.setElement(
                                    Element.HYDROGEN_BOUND_TO_OXYGEN_OR_NITROGEN
                                           );
                        } else if (covalent.getName().trim().startsWith(
                                                     Element.SULPHUR.getSymbol()
                                                                       )
                                  ) {
                            atom.setElement(Element.HYDROGEN_BOUND_TO_SULPHUR);
                        } else if (covalent.getElement() == Element.CARBON
                                   &&
                                   !covalent.isAromatic()) {
                            atom.setElement(Element.HYDROGEN_ALIPHATIC);
                        }
                    }
                }
            }
            if (atom.getElement() == Element.CARBON && atom.isAromatic()) {
                atom.setElement(Element.CARBON_AROMATIC);
            }
*/
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the total solvent accessibility of this amino acid.
     * @return float value representing this amino acid's total SAS or
     * {@code NULL} if not set previously via {@link #setTotalSas(float)}.
     * @see #setTotalSas(float)
     */
    public final float getTotalSas() {
        return totalSolventAccessibility;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the absolute solvent accessibility of this amino acid.
     * @param sas
     *        float value value representing this amino acid's total SAS.
     */
    public final void setTotalSas(final float sas) {
        this.totalSolventAccessibility = sas;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the relative solvent accessibility of this amino acid.
     * @return float value representing this amino acid's relative SAS or
     * {@code NULL} if not set previously via {@link #setRelativeSas(float)}.
     * @see #setRelativeSas(float)
     */
    public final float getRelativeSas() {
        return relativeSolventAccessibility;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the relative solvent accessibility of this amino acid.
     * @param sas
     *        float value value representing this amino acid's relative SAS.
     */
    public final void setRelativeSas(final float sas) {
        this.relativeSolventAccessibility = sas;
    }
    //--------------------------------------------------------------------------

    /**
     * Gets a specific atom from this molecule.
     * @param i
     *        - Index of atom object to be retrieved.
     * @return Atom object with index i
     */
    public final Atom getAtom(final int i) {
        return this.atoms.get(i);
    }
    //--------------------------------------------------------------------------
    /**
     * Gets all atoms in this molecule.
     * @return AtomList object holding all atom objects.
     */
    public final AtomList getAllAtoms() {
        return this.atoms;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all atom and bond related information of the molecule in PDB
     * format.
     * @return String object holding the text information of this molecule in
     *         PDB format.
     */
    public final String toString() {
        StringBuffer buffer = new StringBuffer();
        buffer.append(this.atoms.toString());
//        buffer.append(this.bonds.toString());
        buffer.append("END" + Constants.LINE_SEPERATOR);
    return buffer.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Removes atom from this molecule object. The atom will be searched by
     * object equality and removed from the list of atoms. At the same time
     * bonds which where formed between atom and other atoms of this molecule
     * object will be removed.
     * @param atom
     *        - Atom to be removed.
     */
    public final void remove(final Atom atom) {
        AtomList toBremoved1 = new AtomList();
        for (Atom a : this.atoms) {
            if (a.equals(atom)) {
                toBremoved1.add(a);
                break;
            }
        }
        this.atoms.removeAll(toBremoved1);

        ArrayList<Bond> toBremoved2 = new ArrayList<Bond>();
        for (Bond bond : this.bonds) {
            if (bond.isInBond(atom)) {
                toBremoved2.add(bond);
            }
        }
        this.bonds.removeAll(toBremoved2);
    }
    //--------------------------------------------------------------------------

}

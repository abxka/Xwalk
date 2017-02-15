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
import java.util.Collections;
import java.util.Comparator;

import structure.constants.Constants;
import structure.constants.Constants.BondTypes;
import structure.constants.Constants.ElementTypes;
import structure.math.Mathematics;
import structure.math.Point3f;
import structure.matter.parameter.AminoAcidType;
import structure.matter.parameter.Element;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;


/**
 * A generic class that holds various methods to operate on classes defined in
 * the matter package.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public abstract class MatterUtilities {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected MatterUtilities() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------

    /**
     * Assess whether two atoms are of the same type.
     * @param atom1
     *        - First Atom object.
     * @param atom2
     *        - Second Atom object to assess the equivalence with.
     * @return {@code TRUE} only if both atoms have the same atom name and
     *         residue name, {@code FALSE} otherwise.
     */
    public static boolean equalsType(final Atom atom1, final Atom atom2) {
        if (atom1.getName().trim().equals(atom2.getName().trim())
            &&
            atom1.getResidueName().trim().equals(atom2.getResidueName().trim())
           ) {
               return true;
        }
    return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Assess whether two atoms are in the same residue.
     * @param atom1
     *        - First Atom object.
     * @param atom2
     *        - Second Atom object to assess the equivalence with.
     * @return {@code TRUE} only if both atoms have the same residue name,
     *         residue number and chain Id, {@code FALSE} otherwise.
     */
    public static boolean equalsResidue(final Atom atom1, final Atom atom2) {
        if (atom1.getResidueName().trim().equals(atom2.getResidueName().trim())
            &&
           (atom1.getResidueNumber() == atom2.getResidueNumber())
            &&
           (atom1.getChainId() == atom2.getChainId())) {
            return true;
        }
    return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Assess whether two atoms occupy the same point in Cartesian space.
     * @param atom1
     *        - First Atom object.
     * @param atom2
     *        - Second Atom object to assess the equivalence with.
     * @return {@code TRUE} only if both atoms have the same Cartesian
     *         coordinates, {@code FALSE} otherwise.
     */
    public static boolean equalsPosition(final Atom atom1, final Atom atom2) {
        return atom1.getXYZ().equals(atom2.getXYZ());
    }
    //--------------------------------------------------------------------------

    /**
     * Returns those two atoms that are closest in two AtomList objects.
     * @param list1
     *        - AtomList object holding a first list of atom coordinates.
     * @param list2
     *        - AtomList object holding a second list of atom coordinates.
     * @return AtomList object holding the coordinates of the two closest atoms
     *         in both atom lists.
     */
     public static AtomList getClosestAtomPair(final AtomList list1,
                                               final AtomList list2) {
        float minDist = Integer.MAX_VALUE;
        AtomList minList = new AtomList();
        minList.add(list1.get(0));
        minList.add(list2.get(0));

        for (Atom atom1 : list1) {
            for (Atom atom2 : list2) {
                float dist = Mathematics.distance(atom1.getXYZ(),
                                                  atom2.getXYZ()
                                                 );
                if (dist < minDist) {
                    minDist = dist;
                    minList.set(0, atom1);
                    minList.set(1, atom2);
                }
            }
        }
    return minList;
    }

    //--------------------------------------------------------------------------

    /**
     * Calculates the maximum Cartesian coordinates of an atom list.
     * @param coords
     *        - AtomList object holding the Cartesian coordinates of a list of
     *          atoms.
     * @return Point3d object holding the maximum coordinates in each Cartesian
     *         dimension in its XYZ fields.
     */
    public static Point3f getMaximumCooridnate(final AtomList coords) {
        float maxX = Integer.MIN_VALUE;
        float maxY = Integer.MIN_VALUE;
        float maxZ = Integer.MIN_VALUE;
        for (Atom atom : coords) {
             Point3f xyz = atom.getXYZ();
             float r = atom.getVanDerWaalsRadius();
             maxX = Math.max(maxX, xyz.getX() + r);
             maxY = Math.max(maxY, xyz.getY() + r);
             maxZ = Math.max(maxZ, xyz.getZ() + r);
        }
        return new Point3f(maxX, maxY, maxZ);
    }
    //--------------------------------------------------------------------------

    /**
     * Calculates the minimum Cartesian coordinates of an atom list.
     * @param coords
     *        AtomList object holding the Cartesian coordinates of a list of
     *        atoms.
     * @return Point3f object holding the minimum coordinates in each Cartesian
     *         dimension in its XYZ fields.
     */
    public static Point3f getMinimumCooridnate(final AtomList coords) {
        float minX = Integer.MAX_VALUE;
        float minY = Integer.MAX_VALUE;
        float minZ = Integer.MAX_VALUE;
        for (Atom atom : coords) {
             Point3f xyz = atom.getXYZ();
             float r = atom.getVanDerWaalsRadius();
             minX = Math.min(minX, xyz.getX() - r);
             minY = Math.min(minY, xyz.getY() - r);
             minZ = Math.min(minZ, xyz.getZ() - r);
        }
        return new Point3f(minX, minY, minZ);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the dimension of a protein complex.
     * @param complex
     *        - PolyPeptideList object to which dimension should be calculated.
     * @return float value representing the dimension of the protein complex.
     */
    public static float getDimension(final PolyPeptideList complex) {
        AtomList allAtoms = complex.getAllAtoms();
        Point3f max = MatterUtilities.getMaximumCooridnate(allAtoms);
        Point3f min = MatterUtilities.getMinimumCooridnate(allAtoms);
        Point3f diff = new Point3f(max.getX() - min.getX(),
                                   max.getY() - min.getY(),
                                   max.getZ() - min.getZ());
        float sum = (float) (Math.pow(diff.getX(), 2)
                           + Math.pow(diff.getY(), 2)
                           + Math.pow(diff.getZ(), 2));

        float dim = (float) Math.sqrt(sum);
    return dim;
    }
    //--------------------------------------------------------------------------

    /**
     * Connects a list of atoms depending on the distance and their
     * van der Waals radii + PDB uncertainty factor.
     * @param atoms
     *        - AtomList object holding the coordinates of a PDB molecule.
     * @return Bond array holding all potential covalently bound atoms.
     */
    public static ArrayList < Bond > calculateBonds(final AtomList atoms) {
        // Covalent bonds in organic molecules should have a distance less than
        // 1.54 Angstroem between both atom centers.
        // (see http://en.wikipedia.org/wiki/Bond_length).
        // To account for resolution error distance is taken to be less than
        // 1.8 (0.28 Angstroem estimated standard error for X-ray structures)
        ArrayList < Bond > bonds = new ArrayList < Bond >();
        for (Atom atom1 : atoms) {
            for (Atom atom2 : atoms) {
                if (atom1 != atom2) {
                    if (atom1.getElement().getType() == ElementTypes.METAL
                        ||
                        atom2.getElement().getType() == ElementTypes.METAL) {
                        continue;
                    }
                    float dist = Mathematics.distance(atom1.getXYZ(),
                                                       atom2.getXYZ());
                    float maxDist;
                    if (atom1.getElement() == Element.HYDROGEN
                        ||
                        atom2.getElement() == Element.HYDROGEN) {
                        maxDist = Constants.BOND_LENGTH_TO_HYDROGEN;
                    } else {
                        maxDist = (
                                   atom1.getVanDerWaalsRadius()
                                   +
                                   atom2.getVanDerWaalsRadius()
                                  )
                                + (
                                   Constants.COORDINATE_UNCERTAINTY * 2
                                  );
                    }
                    if (dist <= maxDist) {
                       bonds.add(new Bond(atom1, atom2, BondTypes.SINGLE_BOND));
                    }
                }
            }
        }
        return bonds;
    }
    //--------------------------------------------------------------------------
    /**
     * Sort a list of PolyPeptide objects according to their sequence length.
     * @param peptides
     *        - List of PolyPeptide object to be sorted.
     */
    public static void sort(final ArrayList < PolyPeptide > peptides) {
        Collections.sort(peptides, new Comparator<PolyPeptide>() {
            public int compare(final PolyPeptide p1, final PolyPeptide p2) {
                if (p1.size() < p2.size()) {
                    return 1;
                }
                if (p1.size() > p2.size()) {
                    return -1;
                }
                return 0;
            } });
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    /**
     * Checks whether an amino acid is among those that are were significantly
     * often observed to be at protein-protein interfaces.
     * @param aa
     *        {@link AminoAcid} object.
     * @return {@code TRUE} if amino acid is ILE, HIS, CYS, PHE, MET, TYR or
     *         TRP
     */
    public static boolean isInterfaceTypic(final AminoAcid aa) {
        if (aa.getType().equals(AminoAcidType.ISOLEUCINE)
         || aa.getType().equals(AminoAcidType.HISTIDINE)
         || aa.getType().equals(AminoAcidType.CYSTEINE)
         || aa.getType().equals(AminoAcidType.PHENYLALANINE)
         || aa.getType().equals(AminoAcidType.METHIONINE)
         || aa.getType().equals(AminoAcidType.TYROSINE)
         || aa.getType().equals(AminoAcidType.TRYPTOPHANE)) {
            return true;
        }
    return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether an amino acid is among those that are were not
     * significantly often observed to be at protein-protein interfaces.
     * @param aa
     *        {@link AminoAcid} object.
     * @return {@code TRUE} if amino acid is LYS, GLU, ASP, SER, THR, ALA, PRO
     *         or GLN
     */
    public static boolean isNonInterfaceTypic(final AminoAcid aa) {
        if (aa.getType().equals(AminoAcidType.LYSINE)
         || aa.getType().equals(AminoAcidType.GLUTAMIC_ACID)
         || aa.getType().equals(AminoAcidType.ASPARTIC_ACID)
         || aa.getType().equals(AminoAcidType.SERINE)
         || aa.getType().equals(AminoAcidType.THREONINE)
         || aa.getType().equals(AminoAcidType.ALANINE)
         || aa.getType().equals(AminoAcidType.PROLINE)
         || aa.getType().equals(AminoAcidType.GLUTAMINE)
        ) {
        return true;
        }
    return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns for an atom all atoms from a atomList that are within a certain
     * radius threshold.
     * @param atom
     *        Atom object around which the environment should be determined.
     * @param atomList
     *        AtomList object constituting the potential environment.
     * @param radius
     *        Radius of the environment.
     * @return AtomList object representing the environment.
     */
    public static AtomList getEnvironment(final Atom atom,
                                          final AtomList atomList,
                                          final double radius) {
        AtomList proximity = new AtomList();
        for (Atom atom2 : atomList) {
            if (Mathematics.distance(atom.getXYZ(), atom2.getXYZ()) < radius) {
                proximity.add(atom2);
            }
        }
        return proximity;
    }
}

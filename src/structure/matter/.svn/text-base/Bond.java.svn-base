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

import structure.constants.Constants.BondTypes;

/**
 * Class handles chemical bonds between Atom objects.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Bond {
    /**
     * Atom at pre position within the bond.
     */
    private Atom preAtom;
    /**
     * Atom at post position within the bond.
     */
    private Atom postAtom;
    /**
     * Type of bond.
     */
    private BondTypes bondType;
    /**
     * Constructor.
     * @param pre
     *        - Atom object of first atom.
     * @param post
     *        - Atom object of second atom.
     * @param type
     *        - BondType of bond.
     */
    public Bond(final Atom pre, final Atom post, final BondTypes type) {
        this.preAtom = pre;
        this.postAtom = post;
        this.bondType = type;
    }
    /**
     * Returns the type of this bond.
     * @return BondType object
     */
    public final BondTypes getBondType() {
        return this.bondType;
    }

    /**
     * Returns the first atom in the bond.
     * @return Atom object.
     */
    public final Atom getPreAtom() {
        return this.preAtom;
    }
    /**
     * Returns the second atom in the bond.
     * @return Atom object.
     */
    public final Atom getPostAtom() {
        return this.postAtom;
    }
    /**
     * Checks whether {@code atom} is part of the bond by assessing the
     * identity of the atom via its properties and not by reference.
     * @param atom
     *        - Atom object to be checked for existent in this bond.
     * @return {@code TRUE} if atom is part of this bond, {@code FALSE}
     *         otherwise.
     */
    public final boolean isInBond(final Atom atom) {
        if (preAtom.equals(atom) || postAtom.equals(atom)) {
            return true;
        }
        return false;
    }
}

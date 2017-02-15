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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import structure.matter.Atom;
import structure.matter.AtomList;

import xwalk.math.DistanceComparator;

/**
 * Container for a set of cross-link objects, where the order is determined by
 * the Solvent-Path distance or if not existent by the Euclidean distance.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class CrossLinkList extends ArrayList < CrossLink > {

    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;

    //--------------------------------------------------------------------------
    /**
     * Constructor, which will use the DistanceComparator class for initiating
     * this CrossLinkList.
     */
    public CrossLinkList() {
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a String representation of all cross-links in this container.
     * @return String object holding the representation of all cross-links in
     * this container.
     */
    public final String toString() {
        StringBuffer output = new StringBuffer();
        for (CrossLink crossLink : this) {
            output.append(crossLink.getIndex() + "\t" + crossLink.toString());
        }
        return output.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Sorts the cross-links by Solvent-Path distance, or if not existent
     * by Euclidean distance. Sets in addition the index value for each
     * cross-link value.
     */
    public final void sort() {
        Collections.sort(this, new DistanceComparator());
        for (int i = 0; i < this.size(); i++) {
             this.get(i).setIndex(i + 1);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Puts all cross-links in this list to a Hashtable.
     * @return Hashtable object.
     */
    public final Hashtable < Atom, AtomList > toHash() {
        Hashtable < Atom, AtomList > hash = new Hashtable < Atom, AtomList >();
        for (CrossLink crossLink : this) {
            Atom atomPre  = crossLink.getPreAtom();
            Atom atomPost = crossLink.getPostAtom();

            if (hash.get(atomPre) != null) {
                hash.get(atomPre).add(atomPost);
            } else if (hash.get(atomPost) != null) {
                hash.get(atomPost).add(atomPre);
            } else {
                AtomList paired = new AtomList();
                paired.add(atomPost);
                hash.put(atomPre, paired);
            }
        }
        return hash;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the cross-link which contains both atoms.
     * @param atom1
     *        - First Atom object to be checked in this list of cross-links.
     * @param atom2
     *        - Second Atom object to be checked in this list of cross-links.
     * @return CrossLink object that holds both atoms, NULL otherwise.
     */
    public final CrossLink get(final Atom atom1, final Atom atom2) {
        for (CrossLink crossLink : this) {
            if (crossLink.getPreAtom().equals(atom1)) {
                if (crossLink.getPostAtom().equals(atom2)) {
                    return crossLink;
                }
            }
            if (crossLink.getPreAtom().equals(atom2)) {
                if (crossLink.getPostAtom().equals(atom1)) {
                    return crossLink;
                }
            }
        }
        return null;
    }
}


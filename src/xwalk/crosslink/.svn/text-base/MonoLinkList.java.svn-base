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

import java.util.Comparator;
import java.util.TreeSet;

import structure.matter.Atom;

/**
 * Container for a set of mono-link objects, where the order is determined by
 * the size of the solvent accessible surface area.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class MonoLinkList extends TreeSet < MonoLink > {
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    //--------------------------------------------------------------------------
    /**
     * Default constructor, with Comparator initiation that sorts mono-links
     * first according to the file name that they have been found, then
     * according to their chain ID, then to their residue number.
     */
    public MonoLinkList() {
        super(new Comparator<MonoLink>() {
            /**
             * Comparator methods, which sorts mono-links in this container
             * first according to their chain ID, then to their residue number.
             * @param monoLink1 - 1st MonoLink object.
             * @param monoLink2 - 2nd MonoLink object.
             * @return Returns 1 if sasa1 is smaller than sasa2, -1 if sasa1 is
             *         larger than sasa2 and 0 if sasa1=sasa2.
             */
            public final int compare(final MonoLink monoLink1,
                                     final MonoLink monoLink2) {

                String fileName1 = monoLink1.getFileName();
                String fileName2 = monoLink2.getFileName();
                int compare = fileName1.compareTo(fileName2);
                if (compare == 0) {
                    String chainId1 = monoLink1.getChainId() + "";
                    String chainId2 = monoLink2.getChainId() + "";
                    compare = chainId1.compareTo(chainId2);
                    if (compare == 0) {
                        if (monoLink1.getResidueNumber()
                            >
                            monoLink2.getResidueNumber()) {
                            return 1;
                        } else if (monoLink1.getResidueNumber()
                                   <
                                   monoLink2.getResidueNumber()) {
                            return -1;
                        } else {
                            return 0;
                        }
                    } else {
                        return compare;
                    }
                } else {
                    return compare;
                }
            }
        });
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a String representation of all cross-links in this container.
     * @return String object holding the representation of all cross-links in
     * this container.
     */
    public final String toString() {
        StringBuffer output = new StringBuffer();
        for (MonoLink monoLink : this) {
            output.append(monoLink.getIndex() + "\t" + monoLink.toString());
        }
        return output.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the index value for each mono-link, where the index is determined
     * by the mono-links' decreasing solvent accessible surface area.
     */
    public final void setIndices() {
        int i = 1;
        for (MonoLink monoLink : this) {
             monoLink.setIndex(i++);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the mono-link of an atom.
     * @param atom
     *        - Atom object to be checked in this list of mono-links.
     * @return MonoLink object of the atom, NULL otherwise.
     */
    public final MonoLink get(final Atom atom) {
        for (MonoLink monoLink : this) {
            if (monoLink.equals(atom)) {
                return monoLink;
            }
        }
        return null;
    }
    //--------------------------------------------------------------------------
}


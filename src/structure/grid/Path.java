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

package structure.grid;

import java.util.ArrayList;

import structure.matter.Atom;


/**
 * Class that stores a list of GridCells objects as a Path object.
 * @author abdullah
 * @version 0.1
 * @since 0.1
 */
public class Path extends ArrayList < GridCell > {

    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    //--------------------------------------------------------------------------
    /**
     * Returns the GridCell objects of the Path in PDB format.
     * @param residueNumber
     *        - integer value representing the residue number which is to
     *          be assigned to the atoms forming the path through the grid.
     * @return String object holding the path's GridCell objects in PDB format.
     */
    public final String toString(final int residueNumber) {
        StringBuffer buffer = new StringBuffer();
        int i = 1;
        for (GridCell cell : this) {
             Atom dummy = cell.toAtom();
             dummy.setSerialNumber(i++);
             dummy.setResidueName("PTH");
             dummy.setResidueNumber(residueNumber);
             dummy.setName("C");
             buffer.append(dummy.toString());
        }
        return buffer.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the GridCell objects of the Path in PDB format with residue
     * number = 1.
     * @return String object holding the path's GridCell objects in PDB format.
     */
    public final String toString() {
        return this.toString(1);
    }
}

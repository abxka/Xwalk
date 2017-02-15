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

import structure.matter.AtomList;
import structure.matter.Molecule;
import xwalk.crosslink.CrossLinkerType;

/**
 * Class representing Cross-linker molecules.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class CrossLinker extends Molecule {
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    /**
     * Object to store the type of this CrossLinker object.
     */
    private CrossLinkerType type;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param xlType
     *        - Type of cross-linker as defined by the CrossLinkerType enum.
     */
    public CrossLinker(final CrossLinkerType xlType) {
        super(new AtomList());
        this.type = xlType;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the type of this cross-linker.
     * @return CrossLinkerType object.
     */
    public final CrossLinkerType getType() {
        return type;
    }
}

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

package structure.matter.protein;

import structure.matter.AtomList;
import external.SIMS;

/**
 * Class for calculating the dot surface of a protein using SIMS.
 * @author Abdullah Kahraman
 * @version 0.5
 * @since 0.5
 */
public class DotSurface extends AtomList {
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param complex
     *        List of proteins to which the surface will be calculated.
     * @param sims
     *        The SASA will be calculated using NACCESS.
     */
    public DotSurface(final PolyPeptideList complex, final SIMS sims) {
        sims.createParameterFile();
        sims.run(complex);
        this.addAll(sims.getDotSurface());
    }
}

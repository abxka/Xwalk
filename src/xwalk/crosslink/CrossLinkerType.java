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

/**
 * List of CrossLinker types available to Xwalk.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 *
 */
public enum CrossLinkerType {

    /**
     * Available cross linker types.
     */
    DSS(25, 6);

    //--------------------------------------------------------------------------
    /**
     * Length of the cross-linker.
     */
    private double length;

    //--------------------------------------------------------------------------

    /**
     * Diameter of the cross-linker.
     */
    private double diameter;

    //--------------------------------------------------------------------------

    /**
     * Constructor.
     * @param xlLength
     *        - Length of the cross-linker
     * @param xlDiameter
     *        - Diameter of the cross-linker
     */
    CrossLinkerType(final double xlLength, final double xlDiameter) {
        this.length = xlLength;
        this.diameter = xlDiameter;
    }

    //--------------------------------------------------------------------------

    /**
     * Returns the length of the cross-linker.
     * @return double number representing the length of the cross-linker
     */
    public double getLength() {
        return this.length;
    }

    //--------------------------------------------------------------------------

    /**
     * Returns the diameter of the cross-linker.
     * @return double number representing the diameter of the cross-linker
     */
    public double getdiameter() {
        return this.diameter;
    }
}

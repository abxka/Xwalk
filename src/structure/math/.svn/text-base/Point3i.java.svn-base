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

package structure.math;

import structure.constants.Constants;

/**
 * Simple class for handling Cartesian coordinates in 3 dimensional space with
 * integer numbers.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 * @see Point3d
 */
public class Point3i {
    /**
     * Index I.
     */
    private int iIndex;
    //--------------------------------------------------------------------------
    /**
     * Index J.
     */
    private int jIndex;
    //--------------------------------------------------------------------------
    /**
     * Index K.
     */
    private int kIndex;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param i
     *        - Index i.
     * @param j
     *        - Index j.
     * @param k
     *        - Index k.
     */
    public Point3i(final int i, final int j, final int k) {
        this(i, j, k, 0);
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param i
     *        - Index i.
     * @param j
     *        - Index j.
     * @param k
     *        - Index k.
     * @param radius
     *        - Radius of this point.
     */
    public Point3i(final int i, final int j, final int k, final double radius) {
        this.iIndex = i;
        this.jIndex = j;
        this.kIndex = k;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns index i of this point.
     * @return integer number representing the index i of this
     *         point.
     */
    public final int getI() {
        return iIndex;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns index j of this point.
     * @return integer number representing the index j of this
     *         point.
     */
    public final int getJ() {
        return jIndex;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns index k of this point.
     * @return integer number representing the index k of this
     *         point.
     */
    public final int getK() {
        return kIndex;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether two indices have the same index values.
     * @param  point3iObject
     *              - Point3i object
     * @return {@code TRUE} if both points have the same index values,
     *         {@code FALSE} otherwise.
     */
    public final boolean equals(final Point3i point3iObject) {
        Point3i point3i = (Point3i) point3iObject;
        if (this.iIndex == point3i.getI()
            &&
            this.jIndex == point3i.getJ()
            &&
            this.kIndex == point3i.getK()) {
            return true;
        }
    return false;
    }
    //-------------------------------------------------------------------------
    /*
     * Returns the hash value of this point, which corresponds to the
     * hashCode of the toString() method.
     * @return integer variable representing the hash value of this Point3i
     *         object.
     */
//    public final int hashCode() {
//        return this.toString().hashCode();
//    }
    //--------------------------------------------------------------------------
    /**
     * Creates a copy of this Point3i object.
     * @return Copy of this Point3i object
     */
    public final Point3i copy() {
        return new Point3i(this.getI(),
                           this.getJ(),
                           this.getK());
    }
    //--------------------------------------------------------------------------
    /**
     * Returns an integer array of three elements holding all three elements
     * of this object.
     * @return integer array of three elements.
     */
    public final int[] convert() {
        int[] point = {this.getI(),
                       this.getJ(),
                       this.getK()};
        return point;
    }

    //--------------------------------------------------------------------------
    /**
     * Creates a copy of this Point3i object.
     * @return Copy of this Point3i object
     */
    public final String toString() {
        return this.getI() + ", " + this.getJ() + ", " + this.getK() + ", "
             + Constants.LINE_SEPERATOR;
    }

}

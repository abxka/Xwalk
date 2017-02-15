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

import java.text.DecimalFormat;
import java.text.NumberFormat;

import structure.constants.Constants;


/**
 * Simple class for handling Cartesian coordinates in 3 dimensional space with
 * double numbers.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 1.0
 * @see Point3i
 */
public class Point3d {
    /**
     * Cartesian X-coordinates.
     */
    private double xCoordinate;
    //--------------------------------------------------------------------------
    /**
     * Cartesian Y-coordinates.
     */
    private double yCoordinate;
    //--------------------------------------------------------------------------
    /**
     * Cartesian Z-coordinates.
     */
    private double zCoordinate;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param x
     *        - Cartesian X-coordinates.
     * @param y
     *        - Cartesian Y-coordinates.
     * @param z
     *        - Cartesian Z-coordinates.
     */
    public Point3d(final double x, final double y, final double z) {
        this(x, y, z, 0);
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param x
     *        - Cartesian X-coordinates.
     * @param y
     *        - Cartesian Y-coordinates.
     * @param z
     *        - Cartesian Z-coordinates.
     * @param radius
     *        - Radius of this point.
     */
    public Point3d(final double x, final double y, final double z,
                   final double radius) {
        this.xCoordinate = x;
        this.yCoordinate = y;
        this.zCoordinate = z;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the Cartesian X-coordinate of this point.
     * @return double number representing the Cartesian X-coordinate of this
     *         point.
     */
    public final double getX() {
        return xCoordinate;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the Cartesian Y-coordinate of this point.
     * @return double number representing the Cartesian Y-coordinate of this
     *         point.
     */
    public final double getY() {
        return yCoordinate;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the Cartesian Z-coordinate of this point.
     * @return double number representing the Cartesian Z-coordinate of this
     *         point.
     */
    public final double getZ() {
        return zCoordinate;
    }
    //--------------------------------------------------------------------------
    /**
     * Performs addition to the Cartesian coordinates of this point.
     * @param x
     *        - double to add to X-coordinate
     * @param y
     *        - double to add to Y-coordinate
     * @param z
     *        - double to add to Z-coordinate
     * @return Point3d object with new Cartesian coordinates
     */
    public final Point3d add(final double x, final double y, final double z) {
        return new Point3d(this.getX() + x,
                           this.getY() + y,
                           this.getZ() + z);
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether two points have the same Cartesian coordinates, within
     * an error of 0.001 Angstroem.
     * @param  point3dObject
     *              - Point3d object
     * @return {@code TRUE} if both points have the same coordinates,
     *         {@code FALSE} otherwise.
     */
    public final boolean equals(final Point3d point3dObject) {
        final double errorMargin = 0.001;
        Point3d point3d = (Point3d) point3dObject;

        if (Math.abs(this.xCoordinate - point3d.getX()) <= errorMargin
            &&
            Math.abs(this.yCoordinate - point3d.getY()) <= errorMargin
            &&
            Math.abs(this.zCoordinate - point3d.getZ()) <= errorMargin) {
              return true;
        }
    return false;
    }
    //-------------------------------------------------------------------------
    /*
     * Returns the hash value of this point, which corresponds to the
     * hashCode of the toString() method.
     * @return integer variable representing the hash value of this Point3d
     *         object.
     */
//    public final int hashCode() {
//        return this.toString().hashCode();
//    }
    //--------------------------------------------------------------------------
    /**
     * Creates a copy of this Point3d object.
     * @return Copy of this Point3d object
     */
    public final Point3d copy() {
        return new Point3d(this.getX(),
                           this.getY(),
                           this.getZ());
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a double array of three elements holding all three elements
     * of this object.
     * @return double array of three elements.
     */
    public final double[] convert() {
        double[] point = {this.getX(),
                          this.getY(),
                          this.getZ()};
        return point;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a String representation of this point.
     * @return String object holding the comma separated X,Y,Z coordinates of
     * this point.
     */
    public final String toString() {
        NumberFormat decFormat = new DecimalFormat("0.00");
        return "X:" + decFormat.format(this.getX()) + ", "
             + "Y:" + decFormat.format(this.getY()) + ", "
             + "Z:" + decFormat.format(this.getZ())
             + Constants.LINE_SEPERATOR;
    }
}

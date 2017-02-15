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

package structure.math.pdb;

import structure.math.Point3f;
import structure.matter.Atom;
import structure.matter.AtomList;

/**
 * Class holding various functions for transforming, i.e. translating and
 * rotating protein atoms.
 * @author Abdullah Kahraman
 * @version 0.1
 * @version 0.1
 */
public class Transformation {
    //--------------------------------------------------------------------------
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected Transformation() {
        throw new UnsupportedOperationException();
    }
    /**
     * Returns the dimension of the AtomList object on the X,Y,Z axis.
     * @param atomList
     *           - AtomList object
     * @return A Point3f object that holds the dimensions in X,Y,Z axis
     *         in the X,Y,Z fields of the Point3f object.
     * @see #min(AtomList)
     * @see #max(AtomList)
     */
    public static Point3f dimension(final AtomList atomList) {
        Point3f min = Transformation.min(atomList);
        Point3f max = Transformation.max(atomList);
        Point3f dim = new Point3f(max.getX() - min.getX(),
                                  max.getY() - min.getY(),
                                  max.getZ() - min.getZ());
    return dim;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the minimum X,Y,Z coordinates within a list of atom coordinates.
     * @param atomList
     *           - AtomList object
     * @return A Point3f object that holds the Cartesian coordinates of the
     *         minimum
     *         X,Y,Z coordinates.
     * @see #max(AtomList)
     * @see #dimension(AtomList)
     */
    public static Point3f min(final AtomList atomList) {
        float xMin = Integer.MAX_VALUE;
        float yMin = Integer.MAX_VALUE;
        float zMin = Integer.MAX_VALUE;
        for (Atom atom : atomList) {
             float x = atom.getXYZ().getX();
             float y = atom.getXYZ().getY();
             float z = atom.getXYZ().getZ();
             float r = atom.getVanDerWaalsRadius();
             xMin = Math.min(xMin, x - r);
             yMin = Math.min(yMin, y - r);
             zMin = Math.min(zMin, z - r);
        }
        return new Point3f(xMin, yMin, zMin);
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the maximum X,Y,Z coordinates within a list of atom coordinates.
     * @param atomList
     *           - AtomList object
     * @return A Point3f object that holds the Cartesian coordinates of the
     *         maximum
     *         X,Y,Z coordinates.
     * @see #min(AtomList)
     * @see #dimension(AtomList)
     */
    public static Point3f max(final AtomList atomList) {
        float xMax = Integer.MIN_VALUE;
        float yMax = Integer.MIN_VALUE;
        float zMax = Integer.MIN_VALUE;
        for (Atom atom : atomList) {
             float x = atom.getXYZ().getX();
             float y = atom.getXYZ().getY();
             float z = atom.getXYZ().getZ();
            float r = atom.getVanDerWaalsRadius();
            xMax = Math.max(xMax, x + r);
            yMax = Math.max(yMax, y + r);
            zMax = Math.max(zMax, z + r);
        }
        return new Point3f(xMax, yMax, zMax);
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the center of geometry of any AtomList object.
     * @param atomList
     *        - AtomList object
     * @return Point3f object that holds the Cartesian coordinates of the
     *         center of geometry of the AtomList object.
     * @see #centerOfMass(AtomList)
     */
    public static Point3f centerOfGeometry(final AtomList atomList) {
        Point3f min = Transformation.min(atomList);
        Point3f max = Transformation.max(atomList);
        return new Point3f(min.getX() + (
                                         (max.getX() - min.getX())
                                                     /
                                                     2),
                           min.getY() + (
                                         (max.getY() - min.getY())
                                                     /
                                                     2),
                           min.getZ() + (
                                         (max.getZ() - min.getZ())
                                                     /
                                                     2)
                                        );
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the center of mass of any AtomList object.
     * @param atomList
     *           - AtomList object
     * @return A Point3f object that holds the Cartesian coordinates of the
     *         center of mass of the AtomList object.
     * @see #centerOfGeometry(AtomList)
     */
    public static Point3f centerOfMass(final AtomList atomList) {
        float[] center = new float[3];
        for (Atom atom : atomList) {
             center[0] += atom.getWeight() * atom.getXYZ().getX();
             center[1] += atom.getWeight() * atom.getXYZ().getY();
             center[2] += atom.getWeight() * atom.getXYZ().getZ();
        }
        center[0] /= atomList.size();
        center[1] /= atomList.size();
        center[2] /= atomList.size();
    return new Point3f(center[0], center[1], center[2]);
    }
    //--------------------------------------------------------------------------

    /**
     * Moves all atoms in an AtomList object by X,Y,Z given by the newPosition
     * Point3f object.
     * @param atomList
     *           - AtomList object
     * @param extent
     *           - Point3f object
     */
    public static void move(final AtomList atomList,
                            final Point3f extent) {
        for (Atom atom : atomList) {
             atom.setXYZ(new Point3f(
                                  atom.getXYZ().getX() + extent.getX(),
                                  atom.getXYZ().getY() + extent.getY(),
                                  atom.getXYZ().getZ() + extent.getZ()
                                        )
                            );
        }
    }

    //--------------------------------------------------------------------------

    /**
     * Rotates all atoms in an AtomList object by a rotation matrix.
     * @param atomList
     *           - AtomList object
     * @param rotationMatrix
     *           - double[][] rotation matrix
     * @see Mathematics.getEulerRotationMatrix(double, double, double).
     */
    public static void rotate(final AtomList atomList,
                       final double[][] rotationMatrix) {

        for (Atom atom : atomList) {
            double[] dxyz = {atom.getXYZ().getX(),
                             atom.getXYZ().getY(),
                             atom.getXYZ().getZ()};

            double[] rotatedXYZ =
                    structure.math.Mathematics.multiplyMatrixWithVector(
                                                       rotationMatrix, dxyz
                                                       );

            atom.setXYZ(new Point3f(Float.parseFloat(
                    structure.constants.Constants.CARTESIAN_DEC_FORMAT.format(
                                                                rotatedXYZ[0])),
                                    Float.parseFloat(
                    structure.constants.Constants.CARTESIAN_DEC_FORMAT.format(
                                                                rotatedXYZ[1])),
                                    Float.parseFloat(
                    structure.constants.Constants.CARTESIAN_DEC_FORMAT.format(
                                                                rotatedXYZ[2])))
                                                                              );
        }
    } // end of method turnDependOnInertia()

    //--------------------------------------------------------------------------
    /**
     * Rotates all atoms in an AtomList object by a rotation matrix only
     * after translating the AtomList to the coordinate origin such that
     * the coordinate origin aligns with the AtomList's center of mass.
     * @param atomList
     *           - AtomList object
     * @param rotationMatrix
     *           - double[][] rotation matrix
     * @see Mathematics.getEulerRotationMatrix(double, double, double).
     * @see Transformation.rotate(AtomList, double[][])
     */
    public static void rotateCenterOfMassAtOrigin(
                                              final AtomList atomList,
                                              final double[][] rotationMatrix) {

        Point3f com = Transformation.centerOfMass(atomList);
        Point3f negCom = new Point3f(-com.getX(), -com.getY(), -com.getZ());
        Transformation.move(atomList, negCom);
        Transformation.rotate(atomList, rotationMatrix);
        Point3f comDiff = Transformation.centerOfMass(atomList);
        com = com.add(-comDiff.getX(), -comDiff.getY(), -comDiff.getZ());
        Transformation.move(atomList, com);
    }
    //--------------------------------------------------------------------------
    /**
     * Rotates all atoms in an AtomList object by a rotation matrix only
     * after translating the AtomList to the coordinate origin such that
     * the coordinate origin aligns with the AtomList's center of geometry.
     * @param atomList
     *           - AtomList object
     * @param rotationMatrix
     *           - double[][] rotation matrix
     * @see Mathematics.getEulerRotationMatrix(double, double, double).
     * @see Transformation.rotate(AtomList, double[][])
     */
    public static void rotateCenterOfGeometryAtOrigin(
                                              final AtomList atomList,
                                              final double[][] rotationMatrix) {

        Point3f cog = Transformation.centerOfGeometry(atomList);
        Point3f negCog = new Point3f(-cog.getX(), -cog.getY(), -cog.getZ());
        Transformation.move(atomList, negCog);
        Transformation.rotate(atomList, rotationMatrix);
        Point3f cogDiff = Transformation.centerOfGeometry(atomList);
        cog = cog.add(-cogDiff.getX(), -cogDiff.getY(), -cogDiff.getZ());
        Transformation.move(atomList, cog);
    }
}

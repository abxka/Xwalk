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

import java.util.Comparator;
import java.util.Vector;
import java.util.Collections;
import java.text.DecimalFormat;
import java.util.Locale;

import structure.constants.Constants;
import structure.exceptions.ConversionException;


/**
 * Class holding various useful mathematical functions.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Mathematics {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected Mathematics() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a value between 0 and 1 for {@code x} according to a bell-curved
     * sigmoid function with {@code x=0} having value 1 and {@code x=maxX}
     * having value 0.
     * @param x
     *        - X-value to which Y-value according to sigmoid function will be
     *          calculated.
     * @param maxX
     *        - Normalization factor. Should be the maximum X-value for which
     *          meaningful Y-values should be calculated.
     * @return For x <= maxX the sigmoid function value will be returned,
     *         otherwise 0.
     */
    public static double sigmoidFunction(final double x, final double maxX) {
        double relDist = x
                         /
                         maxX;
        if (relDist <= 1) {
            return 1
                     - 0.5 * (7 * Math.pow(relDist, 2)
                              - 9 * Math.pow(relDist, 4)
                              + 5 * Math.pow(relDist, 6)
                              - 1 * Math.pow(relDist, 8)
                             );
        } else {
            return 0;
        }
    }

    //--------------------------------------------------------------------------

    /**
     * Returns the minimum value from a list of numbers.
     * @param values
     *        - Vector with a list of Number values
     * @return double with minimum value of the elements in values.
     * @see #max(Vector)
     */
    public static double min(final Vector < Number > values) {
        double min = values.get(0).doubleValue();
        for (int i = 0; i < values.size() - 1; i++) {
             min = Math.min(min, values.get(i).doubleValue());
        }
    return min;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the maximum value from a list of numbers.
     * @param values
     *        - Vector with a list of Number values
     * @return double with maximum value of the elements in values.
     * @see #min(Vector)
     */
    public static double max(final Vector < Number > values) {
        double max = values.get(0).doubleValue();
        for (int i = 0; i < values.size() - 1; i++) {
             max = Math.max(max, values.get(i).doubleValue());
        }
    return max;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the average value from a list of numbers.
     * @param values
     *        - Vector with a list of Number values
     * @return double with average value of the elements in values.
     * @see #median(Vector)
     */
    public static double average(final Vector < Number > values) {
        double avg = 0.0;
        for (int i = 0; i < values.size(); i++) {
             avg += values.get(i).doubleValue();
        }
        avg = avg
              /
              values.size();
    return avg;
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates the median from a Vector of Number objects.
     * @param values
     *        - Vector of Number objects.
     * @return median as a double value.
     * @see #average(Vector)
     */
    public static double median(final Vector < Number > values) {
        double median = 0.0;
        Collections.sort(values, new Comparator < Number >() {
            public int compare(final Number t1, final Number t2) {
                Number p1 = ((Number) t1);
                Number p2 = ((Number) t2);
                // sort such that largest value comes first
                if (p1.doubleValue() < p2.doubleValue()) {
                    return 1;
                } else if (p1.doubleValue() > p2.doubleValue()) {
                    return -1;
                } else {
                    return 0;
                }
            }
        });
        // the median is the value at the middle of a sorted list
        int arraySize = values.size();
        if (arraySize % 2 == 0) {
           median = (values.get(arraySize
                                /
                                2 - 1).doubleValue()
                   + values.get(arraySize
                                /
                                2).doubleValue())
                     / 2;
        } else {
           median = (values.get(arraySize
                                /
                                2
                               )
                    ).doubleValue();
        }
    return median;
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates the variance of a Vector of Number objects.
     * @param values
     *        - Vector of Number objects.
     * @return variance as a double value.
     * @see #standardDeviation(Vector)
     */
    public static double variance(final Vector < Number > values) {
        double avg = 0.0;
        double div = 0.0;
        avg = average(values);
        for (int i = 0; i < values.size(); i++) {
            div += Math.pow(values.get(i).doubleValue() - avg, 2);
        }
        if (values.size() == 1) {
            div = 0;
        } else {
            div = div
                  /
                  (values.size() - 1);
        }
    return div;
    }
    //--------------------------------------------------------------------------
    /** Calculates the standard deviation of a Vector of Number objects.
     * @param values
     *         - Vector of Number objects.
     * @return standard deviation as a double value.
     * @see #variance(Vector)
     */
    public static double standardDeviation(final Vector < Number > values) {
        return Math.sqrt(Mathematics.variance(values));
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates the minimum Cartesian coordinates of a list of Point3d[]
     * objects.
     * @param points
     *        - List of Point3d objects
     * @return Point3d object holding the mimumum Cartesian coordinate
     */
    public static Point3d min(final Point3d[] points) {
        double minX = Integer.MAX_VALUE;
        double minY = Integer.MAX_VALUE;
        double minZ = Integer.MAX_VALUE;
        for (Point3d point : points) {
            minX = Math.min(minX, point.getX());
            minY = Math.min(minY, point.getY());
            minZ = Math.min(minZ, point.getZ());
        }
        return new Point3d(minX, minY, minZ);
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates the maximum Cartesian coordinates of a list of Point3d[]
     * objects.
     * @param points
     *        - List of Point3d objects
     * @return Point3d object holding the maximum Cartesian coordinate
     */
    public static Point3d max(final Point3d[] points) {
        double maxX = -Integer.MAX_VALUE;
        double maxY = -Integer.MAX_VALUE;
        double maxZ = -Integer.MAX_VALUE;
        for (Point3d point : points) {
            maxX = Math.min(maxX, point.getX());
            maxY = Math.min(maxY, point.getY());
            maxZ = Math.min(maxZ, point.getZ());
        }
        return new Point3d(maxX, maxY, maxZ);
    }
    //--------------------------------------------------------------------------
    /**
     * Converts Cartesian coordinates into spherical coordinates.
     * @param x
     *        - X-coordinate of a point
     * @param y
     *        - Y-coordinate of a point
     * @param z
     *        - Z-coordinate of a point
     * @return double array with 1st element being the radius, 2nd being the
     *         angle phi and 3rd being the angle theta.
     * @throws ConversionException if an error occurs while converting
     *         coordinates.
     * @see    #rtp2xyz(double, double, double)
     */
    public static double[] xyz2rtp(final double x,
                                   final double y,
                                   final double z) throws ConversionException {

        double r = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2));
        double theta = 0;
        double phi = 0;
        if (Math.abs(r) != 0.0) {
            theta = Math.acos(z
                              /
                              r);
            if (Math.abs(x) > 0.0 && Math.abs(y) > 0.0) {
                phi = Math.atan(y
                                /
                                x);
            } else if (Math.abs(x) > 0.0 && Math.abs(y) < 0.0) {
                phi = 0;
            } else if (Math.abs(x) < 0.0 && y < 0.0) {
                phi = -Math.PI
                       /
                       2;
            } else if (Math.abs(x) < 0.0 && y > 0.0) {
                phi = Math.PI
                      /
                      2;
            }
        } else {
            // doesn't matter which value theta has, because there is no angle
            // necessary when there is no vector (caused of radius = 0)
            theta = 0;
            phi = 0;
        }
        /* phi describes the angle in the xy-plane from the x-axis, whereas
         * theta describes the angle from the z-axis to the radius vector. The
         * equation above for phi works, but has the problem that arctan can not
         * distinguish between y/x and -y/-x or -y/x and y/-x. In both cases the
         * function will return the value for y/x and -y/x which conform to
         * values between -PI/2 and PI/2. With other words, we get only values
         * in the first and fourth coordinate area. To get values in the second
         * and third coordinate, which would conform to a negative x value
         * (x < 0) we had to sum a PI to the phi value. Furthermore the phi
         * value should be used for the recalculation of the x and y, which
         * caused that we had to sum 2*PI to every negative phi value, to get a
         * identical positiv phi value.
         */
        if (x < 0.0 && Math.abs(x) > 0.0) {
            phi = phi + Math.PI;
        }
        if (phi < 0.0 && Math.abs(phi) > 0.0) {
            phi = phi + 2 * Math.PI;
        }
        // testing with recalculating x, y, z
        double[] xyz = new double[3];
        xyz[0] = r * Math.cos(phi) * Math.sin(theta);
        xyz[1] = r * Math.sin(phi) * Math.sin(theta);
        xyz[2] = r * Math.cos(theta);
        double x1 = xyz[0];
        double y1 = xyz[1];
        double z1 = xyz[2];

        final String nL = Constants.LINE_SEPERATOR;
        final double errorMargin = 0.01;

        if (Math.abs(x1 - x) > errorMargin) {
            throw new ConversionException("x:" + x + " " + x1
                                        + "\tinconsistence of value x" + nL
                                        + "y:" + y + " " + y1 + nL
                                        + "z:" + z + " " + z1 + nL
                                        + "radius:" + r + "\t"
                                        + "theta:" + theta + "\t"
                                        + "phi:" + phi + nL);
        }
        if (Math.abs(y1 - y) > errorMargin) {
            throw new ConversionException("y:" + y + " " + y1
                                        + "\tinconsistence of value y" + nL
                                        + "x:" + x + " " + x1 + nL
                                        + "z:" + z + " " + z1 + nL
                                        + "radius:" + r + "\t"
                                        + "theta:" + theta + "\t"
                                        + "phi:" + phi + nL);
        }
        if (Math.abs(z1 - z) > errorMargin) {
            throw new ConversionException("z:" + z + " " + z1
                                        + "\tinconsistence of value z" + nL
                                        + "x:" + x + " " + x1 + nL
                                        + "y:" + y + " " + y1 + nL
                                        + "radius:" + r + "\t"
                                        + "theta:" + theta + "\t"
                                        + "phi:" + phi + nL);
        }
        double[] sphericalCoord = new double[3];
        sphericalCoord[0] = r;
        sphericalCoord[1] = theta;
        sphericalCoord[2] = phi;
    return sphericalCoord;
    } // End of method xyz2rpt()
    //--------------------------------------------------------------------------
    /**
     * Calculates the Cartesian coordinates out of a spherical coordinates
     * radius, phi, theta.
     * @param  radius
     *         Radius
     * @param  phi
     *            Angle in the XY-plane
     * @param  theta
     *         Angle at the Z-axis
     * @return double array with {x,y,z} and "0.000" DecimalFormat
     * @throws ConversionException if an error occurs while converting
     *         coordinates.
     * @see    #xyz2rtp(double, double, double)
     */
    public static double[] rtp2xyz(final double radius,
                                   final double theta,
                                   final double phi)
                                                    throws ConversionException {

        double x = radius * Math.cos(phi) * Math.sin(theta);
        double y = radius * Math.sin(phi) * Math.sin(theta);
        double z = radius *                 Math.cos(theta);
        double[] rtp = xyz2rtp(x, y, z);
        double r1 = rtp[0];
        double t1 = rtp[1];
        double p1 = rtp[2];

        final String nL = Constants.LINE_SEPERATOR;
        final double errorMargin = 0.01;

        if (Math.abs(r1 - radius) > errorMargin) {
            throw new ConversionException("r:" + radius + " " + r1
                                        + "\tinconsistence of value r"
                                        + nL
                                        + "theta:" + theta + " " + t1
                                        + nL
                                        + "phi:" + phi + " " + p1 + nL
                                        + "x:" + x
                                        + nL
                                        + "y:" + y + nL + "z:" + z
                                        + nL);
//            System.exit(1);
        }
        if (Math.abs(t1 - theta) > errorMargin && Math.abs(radius) > 0.0) {
            throw new ConversionException("theta:" + theta + " " + t1
                                        + "\tinconsistence of value theta"
                                        + nL
                                        + "phi:" + phi + " " + p1
                                        + nL
                                        + "r:" + radius + " " + r1 + nL
                                        + "x:" + x
                                        + nL
                                        + "y:" + y
                                        + nL
                                        + "z:" + z
                                        + nL);
//            System.exit(1);
        }
        if (Math.abs(p1 - phi) > errorMargin
            &&
            Math.abs(theta) > 0.0
            &&
            Math.abs(radius) > 0.0) {
            throw new ConversionException("phi:" + phi + " " + p1
                                        + "\tinconsistence of value phi"
                                        + nL
                                        + "theta:" + theta + " " + t1
                                        + nL
                                        + "r:" + radius + " " + r1 + nL
                                        + "x:" + x
                                        + nL
                                        + "y:" + y
                                        + nL
                                        + "z:" + z
                                        + nL);
        }
        double[] cartCoord = new double[3];
        // when calculating XYZ values should be three digits after komma. If no
        // rounded up, calculations can be unprecise.
        Locale.setDefault(Locale.US);
        DecimalFormat decFormat = new DecimalFormat("0.000");
        cartCoord[0] = Double.parseDouble(decFormat.format(x));
        cartCoord[1] = Double.parseDouble(decFormat.format(y));
        cartCoord[2] = Double.parseDouble(decFormat.format(z));
    return cartCoord;
    } // End of method xyz2rpt()
    //--------------------------------------------------------------------------
    /**
     * Calculates the Euclidean distance between two Cartesian coordinates.
     * @param  xyz1
     *         Point3f object holding the first set of Cartesian coordinates.
     * @param  xyz2
     *         Point3f object holding the second set of Cartesian coordinates.
     * @return Euclidean distance as a double value.
     */
    public static float distance(final Point3f xyz1, final Point3f xyz2) {
        return (float) Math.sqrt(Math.pow(xyz2.getX() - xyz1.getX(), 2)
                       + Math.pow(xyz2.getY() - xyz1.getY(), 2)
                       + Math.pow(xyz2.getZ() - xyz1.getZ(), 2)
                        );
    } // End of method distance()
    //--------------------------------------------------------------------------
     //   z
     //   |
     //   |
     //   /------y
     //  /
     // x
    /**
     * Returns an Euler rotation matrix to rotate an object in 3D-space by the
     * angles phi, theta, teta.
     * Z-axis points upwards.
     * Y-axis points rightwards.
     * X-axis points from the screen towards you.
     * <a href=http://mathworld.wolfram.com/EulerAngles.html>Equations 6-14</a>.
     * @param  phi
     *         angle of rotation around z axis
     * @param  theta
     *         angle of rotation around x axis
     * @param  teta
     *         angle of rotation around z axis (again)
     * @return Two dimensional double array [3][3] holding the rotation matrix.
     */
    public static double[][] getEulerRotationMatrix(final double phi,
                                                    final double theta,
                                                    final double teta) {
        double[][] a = new double[3][3];
        a[0][0] = Math.cos(teta) * Math.cos(phi)
                - Math.cos(theta) * Math.sin(phi) * Math.sin(teta);
        a[0][1] = Math.cos(teta) * Math.sin(phi)
                + Math.cos(theta) * Math.cos(phi) * Math.sin(teta);
        a[0][2] = Math.sin(teta) * Math.sin(theta);
        a[1][0] = -Math.sin(teta) * Math.cos(phi)
                 - Math.cos(theta) * Math.sin(phi) * Math.cos(teta);
        a[1][1] = -Math.sin(teta) * Math.sin(phi)
                 + Math.cos(theta) * Math.cos(phi) * Math.cos(teta);
        a[1][2] =  Math.cos(teta) * Math.sin(theta);
        a[2][0] = Math.sin(theta) * Math.sin(phi);
        a[2][1] = -Math.sin(theta) * Math.cos(phi);
        a[2][2] = Math.cos(theta);
        return a;
    }
    //--------------------------------------------------------------------------
    /**
     * Length of an vector in a Cartesian coordinate system.
     * @param v
     *        - Array of Number objects representing the coordinates in a
     *          Cartesian coordinate system, e.g. {x,y,z}.
     * @return length of vector as a double value
     */
    public static double length(final Number[] v) {
        double sum = 0;
        for (Number n : v) {
             sum += Math.pow(n.doubleValue(), 2);
        }
    return Math.sqrt(sum);
    }
    //--------------------------------------------------------------------------
    /**
     * Sum of two vectors in a Cartesian coordinate system.
     * @param vector1
     *        - Array of Number objects representing the coordinates of the 1st
     *          vector in a Cartesian coordinate system, e.g. {x,y,z}.
     * @param vector2
     *        - Array of Number objects representing the coordinates of the 2nd
     *          vector in a Cartesian coordinate system, e.g. {x,y,z}.
     * @return Sum vector with the same number of elements as vector1/2, with
     *         elements being double values.
     * @throws ArrayIndexOutOfBoundsException if vector1 and vector2 arrays have
     *         different number of elements.
     * @see    #subtraction(Number[], Number[])
     */
    public static Number[] addition(final Number[] vector1,
                                    final Number[] vector2)
                                         throws ArrayIndexOutOfBoundsException {
        if (vector1.length != vector2.length) {
            throw new ArrayIndexOutOfBoundsException("ERROR: Can not add "
                                                   + "vectors in different "
                                                   + "dimensions."
                                                   + Constants.LINE_SEPERATOR
                                                   + "Vector1 has \""
                                                   + vector1.length + "\" "
                                                   + "elements, Vector2 has \""
                                                   + vector2.length
                                                   + "\" elements.");
        }
        Number[] sum = new Number[vector1.length];
        for (int i = 0; i < vector1.length; i++) {
            sum[i] = vector1[i].doubleValue() + vector2[i].doubleValue();
        }
    return sum;
    }
    //--------------------------------------------------------------------------
    /**
     * Difference vector between two vectors in a Cartesian coordinate system.
     * @param vector1
     *        - Array of Number objects representing the coordinates of the 1st
     *          vector in a Cartesian coordinate system, e.g. {x, y, z}.
     * @param vector2
     *        - Array of Number objects representing the coordinates of the 2nd
     *          vector in a Cartesian coordinate system, e.g. {x,y,z}.
     * @return Difference vector with the same number of elements as vector1/2,
     *         with elements being double values.
     * @throws ArrayIndexOutOfBoundsException if vector1 and vector2 arrays
     *         have different number of elements.
     * @see    #addition(Number[], Number[])
     */
    public static Number[] subtraction(final Number[] vector1,
                                       final Number[] vector2
                                      ) throws ArrayIndexOutOfBoundsException {
        if (vector1.length != vector2.length) {
            throw new ArrayIndexOutOfBoundsException("ERROR: Can not substract "
                                                   + "vectors in different "
                                                   + "dimensions."
                                                   + Constants.LINE_SEPERATOR
                                                   + "Vector1 has \""
                                                   + vector1.length + "\" "
                                                   + "elements, Vector2 has \""
                                                   + vector2.length
                                                   + "\" elements.");
        }
        Number[] diff = new Number[vector1.length];
        for (int i = 0; i < vector1.length; i++) {
             diff[i] = vector2[i].doubleValue() - vector1[i].doubleValue();
        }
    return diff;
    }
    //--------------------------------------------------------------------------
    /**
     * In/Decreases the length of a vector in a Cartesian coordinate system, by
     * multiplying {@code vector} elements with {@code factor}.
     * @param vector
     *        - Array of Number objects representing the coordinates of an
     *          vector in a Cartesian coordinate system, e.g. {x,y,z}.
     * @param factor
     *        - In/Decrease vector by this number.
     * @return Vector with in/decreased length and with elements being double
     *         values.
     * @see    #invertVector(Number[])
     * @see    #multiplyMatrixWithVector(Number[][], Number[])
     */
    public static Number[] multVectorByFactor(final Number[] vector,
                                               final double factor) {
        Number[] mult = new Number[vector.length];
        for (int i = 0; i < vector.length; i++) {
            mult[i] = factor * vector[i].doubleValue();
        }
    return mult;
    }
    //--------------------------------------------------------------------------
    /**
     * Inverting the direction of a vector in a Cartesian coordinate system, by
     * multiplying each {@code vector} element with -1.
     * @param vector
     *        - Array of Number objects representing the coordinates of an
     *          vector in an Cartesian coordinate system, e.g. {x,y,z}.
     * @return Inverted vector with elements being double values.
     * @see    #multVectorByFactor(Number[], double)
     * @see    #multiplyMatrixWithVector(Number[][], Number[])
     */
    public static Number[] invertVector(final Number[] vector) {
        return multVectorByFactor(vector, -1);
    }
    //--------------------------------------------------------------------------
    /**
     * Multiplying a vector with a matrix.
     * @param matrix
     *        - Two dimensional array of Number objects representing a matrix in
     *          a Cartesian coordinate system.
     * @param vector
     *        - Array of Number objects representing the coordinates of an
     *          vector in a Cartesian coordinate system, e.g. {x,y,z}.
     * @return Inverted vector with elements being double values.
     * @throws ArrayIndexOutOfBoundsException if the number of {@code matrix}
     *         columns is not the number of elements in {@code vector}.
     * @see    #multVectorByFactor(Number[], double)
     */
    public static double[] multiplyMatrixWithVector(final double[][] matrix,
                                                    final double[] vector)
                                         throws ArrayIndexOutOfBoundsException {
        double[] result = new double[vector.length];
        for (int row = 0; row < matrix.length; row++) {
            if (matrix[row].length != vector.length) {
                throw new ArrayIndexOutOfBoundsException("Matrix column "
                                                      + "numbers are not "
                                                      + "equal to vector row "
                                                      + "numbers"
                                                      + Constants.LINE_SEPERATOR
                                                        );
            }
            double s = 0;
            for (int column = 0; column < matrix[row].length; column++) {
                s += matrix[row][column]
                     *
                     vector[column];
            }
            result[row] = s;
        }
        return result;
    }
}

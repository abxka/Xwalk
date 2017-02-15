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

import structure.math.Point3f;
import structure.math.Point3i;
import structure.matter.Atom;

/**
 * Class for handling general grid objects.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Grid {

    /**
     * Grid cells that constitute the grid.
     */
    private GridCell[][][] gridCells;

    /**
     * Stores the maximum Cartesian coordinates of the grid.
     */
    private Point3f max;

    /**
     * Stores the minimum Cartesian coordinates of the grid.
     */
    private Point3f min;

    /**
     * Stores the number of cells in the X dimension.
     */
    private int noOfxCells = -1;
    /**
     * Stores the number of cells in the Y dimension.
     */
    private int noOfyCells = -1;
    /**
     * Stores the number of cells in the ~ dimension.
     */
    private int noOfzCells = -1;

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param minimum
     *        - Point3f object holding the minimum X,Y,Z Cartesian coordinates
     *          of this grid.
     * @param maximum
     *        - Point3f object holding the maximum X,Y,Z Cartesian coordinates
     *          of this grid.
     * @param gridCellSize
     *        - float value representing the size of all grid cells i.e. their
     *          cell edge length.
     */
    public Grid(final Point3f minimum,
                final Point3f maximum,
                final float gridCellSize) {
        this.min = minimum;
        this.max = maximum;
        this.setNumberOfCells(gridCellSize);
        this.gridCells = new GridCell[this.noOfxCells]
                                     [this.noOfyCells]
                                     [this.noOfzCells];

        for (int i = 0; i < this.noOfxCells; i++) {
            for (int j = 0; j < this.noOfyCells; j++) {
                for (int k = 0; k < this.noOfzCells; k++) {
                     float x = this.min.getX() + (i * gridCellSize)
                                                           + (gridCellSize / 2);
                     float y = this.min.getY() + (j * gridCellSize)
                                                           + (gridCellSize / 2);
                     float z = this.min.getZ() + (k * gridCellSize)
                                                           + (gridCellSize / 2);

                    this.gridCells[i][j][k] = new GridCell(new Point3f(x, y, z),
                                                           gridCellSize);
                    this.gridCells[i][j][k].setIndices(new Point3i(i, j, k));
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Returns the minimum Cartesian coordinate of this grid.
     * @return Point3f
     *        - Point3f object, which stores the minimum Cartesian coordinate.
     */
    protected final Point3f getMin() {
        return this.min;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the maximum Cartesian coordinate of this grid.
     * @return Point3f
     *        - Point3f object, which stores the maximum Cartesian coordinate.
     */
    protected final Point3f getMax() {
        return this.max;
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates the number of cells in all three XYZ Cartesian dimensions.
     * @param gridCellSize
     *        - float value representing the size of all grid cells i.e. their
     *          cell edge length.
     */
    private void setNumberOfCells(final float gridCellSize) {
       this.noOfxCells = Math.round((float) ((this.max.getX() - this.min.getX())
                                             /
                                             gridCellSize
                                            )
                                   );
       this.noOfyCells = Math.round((float) ((this.max.getY() - this.min.getY())
                                             /
                                             gridCellSize
                                            )
                                   );
       this.noOfzCells = Math.round((float) ((this.max.getZ() - this.min.getZ())
                                             /
                                             gridCellSize
                                            )
                                   );
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the number of cells for each Cartesian dimension.
     * @return Point3i object holding the number of cells for each Cartesian
     *         dimension.
     */
    public final Point3i getNumberOfCells() {
        return new Point3i(this.noOfxCells, this.noOfyCells, this.noOfzCells);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a grid cell with the grid.
     * @param i
     *        - integer value representing the index on the X dimension.
     * @param j
     *        - integer value representing the index on the Y dimension.
     * @param k
     *        - integer value representing the index on the Z dimension.
     * @return GridCell object with indices i,j,k. If indices extend over grid
     *         borders, then {@code NULL} is returned.
     */
    public final GridCell get(final int i, final int j, final int k) {
        boolean indicesLarger0 = i >= 0 && j >= 0 && k >= 0;
        boolean indicesSmallerMax = i < this.noOfxCells
                                    &&
                                    j < this.noOfyCells
                                    &&
                                    k < this.noOfzCells;
        if (indicesLarger0 && indicesSmallerMax) {
            return this.gridCells[i][j][k];
        } else {
            return null;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Resets the value and occupied status of ALL grid cells in the
     * grid.
     */
    public final void reset() {

        for (int i = 0; i < noOfxCells; i++) {
            for (int j = 0; j < noOfyCells; j++) {
                for (int k = 0; k < noOfzCells; k++) {
                    this.gridCells[i][j][k].reset();
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Resets only the value leaving occupied status as it is.
     */
    public final void resetSoft() {
        for (int i = 0; i < noOfxCells; i++) {
            for (int j = 0; j < noOfyCells; j++) {
                for (int k = 0; k < noOfzCells; k++) {
                    this.gridCells[i][j][k].resetSoft();
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the grid cell points in PDB format.
     * @return String object holding the grid in PDB format.
     */
    public final String toString() {
        StringBuffer output = new StringBuffer();
        for (int i = 0; i < noOfxCells; i++) {
             for (int j = 0; j < noOfyCells; j++) {
                  for (int k = 0; k < noOfzCells; k++) {
                       GridCell cell = this.get(i, j, k);
                       Atom dummy = cell.toAtom();
                       dummy.setName("C");
                       dummy.setResidueName("GRD");
                       dummy.setResidueNumber(1);
                       output.append(dummy.toString());
                }
            }
        }
        return output.toString();
    }
}

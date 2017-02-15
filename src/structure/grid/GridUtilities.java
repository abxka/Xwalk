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

import structure.constants.Constants;
import structure.math.Point3d;
import structure.math.Point3i;
import structure.matter.Atom;
import xwalk.crosslink.CrossLinkParameter;
import xwalk.crosslink.CrossLinkParameter.Parameter;


/**
 * Class to create and handle grid objects for AtomList objects.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 *
 */
public class GridUtilities {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected GridUtilities() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the number of grid cells that fit within the hemisphere.
     * @param radius
     *        - double value representing the radius of the hemisphere.
     * @param gridCellSize
     *        - double value representing the size of the grid cell.
     * @return integer value representing the number of grid cells.
     */
    public static int getNumberOfGridCellsFittingIntoHemisphere(
                                                       final double radius,
                                                       final double gridCellSize
                                                               ) {
        final double minRadius = 0.01;
        double numberOfCellsInAtom;
        if (radius == 0.0) {
            numberOfCellsInAtom = minRadius;
        } else {
            numberOfCellsInAtom = 0.5 + radius
                                        /
                                        gridCellSize;
        }
        return Math.round((float) numberOfCellsInAtom);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the volume information of a grid.
     * @param grid
     *        - Grid object from which the volume is calculated from.
     * @return Point3d object with volume X-coordinate = occupied volume,
     *         Y-coordinate = unoccupied volume, Z-coordinate = total volume.
     */
    public static Point3d getVolume(final Grid grid) {
        int hit = 0;
        int nonHit = 0;
        int count = 0;
        for (int i = 0; i < grid.getNumberOfCells().getI(); i++) {
            for (int j = 0; j < grid.getNumberOfCells().getJ(); j++) {
                for (int k = 0; k < grid.getNumberOfCells().getK(); k++) {
                    count++;
                    if (grid.get(i, j, k).isOccupied()) {
                        hit++;
                    } else {
                        nonHit++;
                    }
                }
            }
        }
        // double check that everything works fine.
        if (hit + nonHit != count) {
            System.err.println("WARNING: hit (" + hit + ") nonHit (" + nonHit
                               + ") != count(" + count + ")");
        }
        double gridCellSize = grid.get(0, 0, 0).getSize();
        double gridCellVolume = Math.pow(gridCellSize, 3);
    return new Point3d(gridCellVolume * hit,
                       gridCellVolume * nonHit,
                       gridCellVolume * count);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all cells that are unoccupied.
     * @param grid
     *        - GridObject from which all unoccupied cells will be extracted.
     * @return ArrayList of GridCells being unoccupied.
     */
    public static ArrayList < GridCell > getUnoccupiedCells(final Grid grid) {
        ArrayList < GridCell > candidates = new ArrayList < GridCell >();
        for (int i = 0; i < grid.getNumberOfCells().getI(); i++) {
            for (int j = 0; j < grid.getNumberOfCells().getJ(); j++) {
                for (int k = 0; k < grid.getNumberOfCells().getK(); k++) {
                    GridCell cell = grid.get(i, j, k);
                    if (!cell.isOccupied()) {
                        candidates.add(cell);
                    }
                }
            }
        }
        return candidates;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all direct neighboring grid cells.
     * @param cell
     *        - GridCell object whose neighborhood is to be determined.
     * @param grid
     *        - Grid object that holds all grid cells including the {@code cell}
     *          and all neighboring grid cells.
     * @param cellSize
     *        - integer value representing the number of shells from which
     *          neighborhood is to be extracted.
     * @return ArrayList of neighboring grid cells.
     */
    public static ArrayList < GridCell > getNeighbouringCells(
                                                            final GridCell cell,
                                                            final Grid grid,
                                                            final int cellSize
                                                             ) {
        ArrayList < GridCell > neighbours = new ArrayList < GridCell >();
        Point3i ijk = cell.getIndices();

        for (int m = -cellSize;
             m <= cellSize && ijk.getI() + m >= 0
             &&
             ijk.getI() + m < grid.getNumberOfCells().getI();
             m++
             ) {
             for (int n = -cellSize;
                  n <= cellSize && ijk.getJ() + n >= 0
                  &&
                  ijk.getJ() + n < grid.getNumberOfCells().getJ();
                  n++
                  ) {
                  for (int o = -cellSize;
                       o <= cellSize && ijk.getK() + o >= 0
                       &&
                       ijk.getK() + o < grid.getNumberOfCells().getK();
                       o++
                       ) {
                       // do not add cell itself into the list of neighbours.
                       if (m == 0 && n == 0 && o == 0) {
                           continue;
                       }
                       neighbours.add(grid.get(ijk.getI() + m,
                                               ijk.getJ() + n,
                                               ijk.getK() + o));
                }
            }
        }
        return neighbours;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether any neighbouring cells are labeled as boundary.
     * @param atom
     *        - Atom object to be checked whether it is solvent accessible.
     * @param grid
     *        - Grid object that holds all grid cells including the {@code cell}
     *          and all neighboring grid cells. Prior to the execution of this
     *          method, the Grid object must have been searched for the boundary
     *          cells with the BoundarySearch class.
     * @return {@code TRUE} if cell is accessible, {@code FALSE} otherwise.
     * @see structure.math.algorithms.BoundarySearch
     */
    public static boolean isAccessible(final Atom atom,
                                       final AtomGrid grid) {

        boolean verbose = Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                                             ));

        ArrayList<GridCell> neighboursBorder = grid.getAllGridCells(atom, 1);
        ArrayList<GridCell> neighbours = grid.getAllGridCells(atom, 0);
        ArrayList<GridCell> border = new ArrayList<GridCell>();
        for (GridCell cell1 : neighboursBorder) {
            boolean found = false;
            for (GridCell cell2 : neighbours) {
                if (cell1 == cell2) {
                    found = true;
                }
            }
            if (!found) {
                border.add(cell1);
            }
        }

        if (!verbose) {
            for (GridCell borderCell : border) {
                if (!borderCell.isOccupied()) {
                    return true;
                }
            }
            return false;
        } else {
            int n = 0;
            for (GridCell borderCell : border) {
                if (!borderCell.isOccupied()) {
                    n++;
                }
            }
            System.err.print("Following atom has " + n + " grid cells at the "
                           + "protein/solvent boundary:"
                           + Constants.LINE_SEPERATOR
                           + atom);
            if (n > 0) {
                return true;
            }
            return false;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the GridCell element that is equal to the target GridCell object.
     * @param target
     *        - GridCell object against which {@code cells} will be matched.
     * @param cells
     *        - List of GridCell objects to be matched against {@code target}.
     * @return GridCell object if found in the list of GridCell object,
     *         otherwise {@code NULL}.
     */
    public static GridCell equals(final GridCell target,
                                  final ArrayList < GridCell > cells) {
        for (GridCell cell : cells) {
            if (cell.equals(target)) {
                return cell;
            }
        }
        return null;
    }
    //--------------------------------------------------------------------------
}

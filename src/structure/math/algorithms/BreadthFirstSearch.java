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

package structure.math.algorithms;

import java.util.ArrayList;

import structure.constants.Constants;
import structure.grid.Grid;
import structure.grid.GridCell;
import structure.grid.GridUtilities;
import structure.grid.Path;
import structure.math.Mathematics;
import xwalk.crosslink.CrossLinkParameter;
import xwalk.crosslink.CrossLinkParameter.Parameter;

/**
 * Set the distances in the surroundings of a source cell with a grid.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class BreadthFirstSearch {
    /**
     * Shortest path as calculated between a source and a target cell.
     */
    private Path path = new Path();

    /**
     * Constant, indicating that the target cell is the first in the path.
     */
    public static final int CELL_NO_OF_TARGET_CELL_IN_PATH = 0;

    /**
     * Boolean indicating whether distance assignment had to end prematurely,
     * due to solvent inaccessibility of cross-link.
     */
    private boolean hasSucceeded = false;

    /**
     * List of target grid cells, which represent the end point in the
     * distance calculation.
     */
    private ArrayList < GridCell > targets;

    /**
     * Grid object in which the entire search is done.
     */
    private Grid grid;

    /**
     * float value representing the maximum distance to search for in the grid.
     */
    private float maxDist;

    /**
     * Source grid cell, which represents the starting point for the
     * distance calculation.
     */
    private GridCell source;

    /**
     * Trash for removing grid cells from the search for the shortest path.
     */
    private static ArrayList <GridCell> trash = new ArrayList<GridCell>();

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param grid
     *        - Grid object in which the entire search is done.
     * @param source
     *        - Source grid cell, which represents the starting point for the
     *          distance calculation.
     * @param targets
     *        - List of target grid cells, which represent the end point in the
     *          distance calculation
     * @param maxDist
     *        - float value representing the maximum distance to search for in
     *          the grid
     */
    public BreadthFirstSearch(
                              final Grid grid,
                              final GridCell source,
                              final ArrayList < GridCell > targets,
                              final float maxDist) {
        this.grid = grid;
        this.source = source;
        this.targets = targets;
        this.maxDist = maxDist;

        // set value of source cell to 0.0
        this.source.setDistance(0.0f);

    }
    //--------------------------------------------------------------------------
    /**
     * Perform breadth-first search on the grid to find the shortest path
     * between a single grid cell and a list of other grid cells.
     * @return List of path objects, holding each the path between the source
     *         and one target cell. If no path could be found, than each
     *         path object holds only the target cell.
     */
    public final ArrayList < Path > findShortestPath() {

        ArrayList < GridCell > actives = new ArrayList < GridCell >();

        actives.add(this.source);

        // start breadth-first search from grid cell.
        this.setDistanceRecursively(actives);

        // trace back the path
        ArrayList < Path > paths = new ArrayList < Path >();
        for (GridCell target : this.targets) {
            // first element in the path is the target cell itself.
            // Last element will be the source cell.
            this.path.add(BreadthFirstSearch.CELL_NO_OF_TARGET_CELL_IN_PATH,
                          target.copy());
            // if distance calculation succeeded and user wants a PyMOL output
            // backtrack shortest path.
            if (target.getDistance() != Constants.DEFAULT_GRID_DISTANCE
               &&
               Boolean.parseBoolean(CrossLinkParameter.getParameter(
                                                      Parameter.DO_PYMOL_OUTPUT)
                                                                    )) {
                this.backtrackPath(target);
                paths.add(this.path);
                this.path = new Path();
            } else {
                // otherwise leave path only with targetCell.
                paths.add(this.path);
                this.path = new Path();
            }
        }
        return paths;
    }
    //--------------------------------------------------------------------------

    /**
     * Sets the distances for all grid cells that lie in-between the source cell
     * and a list of target cells.
     * @param actives
     *        - List of GridCells for which neighbouring cells have to be
     *          determined and distances calculated.
     */
    private void setDistanceRecursively(final ArrayList < GridCell > actives) {
        // neighboring grid cells become new actives for the next round of
        // breadth-first-search.
        ArrayList < GridCell > newActives = new ArrayList < GridCell >();
        for (GridCell active : actives) {

             ArrayList < GridCell > neighbours =
                            GridUtilities.getNeighbouringCells(active,
                                                               this.grid,
                                                               1);

             for (GridCell neighbour : neighbours) {
                  float currentDist = Integer.MIN_VALUE;
                  float newDist = Integer.MIN_VALUE;

                  if (!neighbour.isOccupied()) {
                      currentDist = neighbour.getDistance();
                      // The distance of the neighboring grid cell is the
                      // distance of the current active cell + the distance
                      // between the active and the neighboring cell.
                      newDist = active.getDistance()
                              + (float) Mathematics.distance(active.getXYZ(),
                                                             neighbour.getXYZ()
                                                            );
                      // distance from other active cell might be shorter.
                      if (newDist < currentDist) {
                          neighbour.setDistance(newDist);
                          if (!newActives.contains(neighbour)) {
                              newActives.add(neighbour);
                          }
                      }
                  }
                  /*
                   * -------------------------------------------------------
                   * Following section has to be removed as it produces
                   * inconsistent results when running distance calculations
                   * for multiple cross-links.
                   * -------------------------------------------------------
                  // check whether we have reached the target cell already.
                  // If so, remove target from list of targets.
                  GridCell equal = GridUtilities.equals(neighbour,
                                                        this.targets);
                  if (equal != null) {
                      this.targetsFoundInSearch.add(neighbour);
                      if (this.targets.size()
                          ==
                          this.targetsFoundInSearch.size()) {
                          this.hasFinished = true;
                          // The next return command could be removed
                          // which would allow a consistent SASD calculation
                          // between Xwalk runs with explicit AA information
                          // via command line and AA provided by a distance
                          // file. However at the same time the calculation
                          // time increases by 100%. For the moment the
                          // preference is put on speed.
                          return;
                      }
                  }
                  */
             }
        }
        // Check whether it is necessary to continue distance calculation,
        // as all newActives might have already distances larger than maxDist.
        int maxDistCount = 0;
        for (GridCell neighbour : newActives) {
             double neighbourDistance = neighbour.getDistance();
             if (neighbourDistance > this.maxDist) {
                 maxDistCount++;
                 BreadthFirstSearch.trash.add(neighbour);
                 this.hasSucceeded = true;
             }
        }
        // break up recursive loop.
        if (maxDistCount == newActives.size()) {
            return;
        }

        newActives.removeAll(BreadthFirstSearch.trash);
        BreadthFirstSearch.trash.clear();

        this.setDistanceRecursively(newActives);
    }
    //--------------------------------------------------------------------------
    /**
     * Backtraces the path starting between a target cell and source cell.
     * @param target
     *        - Target grid cell, which represent the end point in the
     *        distance calculation.
     */
    private void backtrackPath(final GridCell target) {
        ArrayList < GridCell > neighbours =
                            GridUtilities.getNeighbouringCells(target,
                                                               this.grid,
                                                               1);
        float minDist = target.getDistance();
        GridCell minGridCell = target;
        for (GridCell neighbour : neighbours) {
            float dist = neighbour.getDistance();
            if (dist < minDist) {
                minDist = dist;
                minGridCell = neighbour;
            }
        }
        if (minGridCell == this.source) {
            this.path.add(minGridCell.copy());
            return;
        }
        this.path.add(minGridCell.copy());
        this.backtrackPath(minGridCell);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether the search for a target was successful.
     * @return {@code TRUE} if search found target cell, {@code FALSE}
     * otherwise.
     */
    public final boolean hasSucceeded() {
        return this.hasSucceeded;
    }
}

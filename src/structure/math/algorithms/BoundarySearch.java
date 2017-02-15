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
/*
import java.util.ArrayList;

import structure.grid.GridCell;
import structure.grid.GridUtilities;
*/
import structure.grid.Grid;

/**
 * Class which searches for GridCell objects within a Grid object that are
 * labeled as occupied and are next to GridCells that are labeled as unoccupied.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 *
 */
public class BoundarySearch {

    /**
     * Constructor with prevention against calls from subclass.
     */
    protected BoundarySearch() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Labels those cells in the grid that form the boundary between occupied
     * and unoccupied cells.
     * @param grid
     *        - Grid object holding GridCell objects that are either occupied
     *          or unoccupied.
     */
    public static final void findBoundary(final Grid grid) {

        /*
        // initialize grid by setting all cells to visited = FALSE.
        for (int i = 0; i < grid.getNumberOfCells().getI(); i++) {
            for (int j = 0; j < grid.getNumberOfCells().getJ(); j++) {
                for (int k = 0; k < grid.getNumberOfCells().getK(); k++) {

                    GridCell cell = grid.get(i, j, k);
                    ArrayList<GridCell> neighbours =
                                        GridUtilities.getNeighbouringCells(
                                                                           cell,
                                                                           grid,
                                                                           1
                                                                          );

                    for (GridCell neighbour : neighbours) {
                        if (cell.isOccupied() != neighbour.isOccupied()) {
                            GridCell occupiedCell = cell;
                            if (cell.isOccupied()) {
                                occupiedCell = neighbour;
                            }
                            occupiedCell.setBoundaryStatus();
                            break;
                        }
                    }
                }
            }
        }
        */
    }
}

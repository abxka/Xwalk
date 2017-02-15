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

package xwalk.math;

import java.util.ArrayList;

import structure.constants.Constants;
import structure.grid.AtomGrid;
import structure.grid.Grid;
import structure.grid.GridCell;
import structure.grid.Path;
import structure.math.Point3i;
import structure.math.algorithms.BreadthFirstSearch;
import structure.matter.Atom;
import structure.matter.AtomList;


/**
 * Master class for calculating and storing Solvent-Path distances.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class SolventPathDistance {
    /**
     * Source grid cell.
     */
    private GridCell sourceCell;
    /**
     * All candidate target grid cells.
     */
    private ArrayList < GridCell > targetCells;
    /**
     * Grid object which will be used to calculate the Solvent-Path distance.
     */
    private Grid grid;

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param atom1cell
     *        - GridCell object of first protein atom to be connected by the
     *          virtual cross-linker.
     * @param atom2cells
     *        - List of GridCell objects of second protein atoms to be connected
     *          by the virtual cross-linker.
     * @param grid
     *        - AtomGrid object build on protein structure.
     */
    public SolventPathDistance(final GridCell atom1cell,
                               final ArrayList < GridCell > atom2cells,
                               final Grid grid
                              ) {
        this.sourceCell = atom1cell;
        this.targetCells = atom2cells;
        this.grid = grid;
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param atom1
     *        - First protein atom to be connected by the virtual cross-linker.
     * @param atoms2
     *        - List of second protein atoms to be connected by the virtual
     *          cross-linker.
     * @param atomGrid
     *        - Grid build on protein structure.
     */
    public SolventPathDistance(final Atom atom1,
                               final AtomList atoms2,
                               final AtomGrid atomGrid) {

        Atom atom1small = atom1.copy();
        atom1small.setVanDerWaalsRadius();
        atom1small.setVanDerWaalsRadius(atom1small.getVanDerWaalsRadius()
                                        +
                                        Constants.SOLVENT_RADIUS);
        ArrayList < GridCell > atom1Cells = atomGrid.getAllGridCells(
                                                                    atom1small
                                                                    );
        for (GridCell cell1 : atom1Cells) {
            // set all cells that are occupied by this atom to unoccupied
             cell1.reset();
        }
        ArrayList < GridCell > atom2cells = new ArrayList < GridCell >();
        for (Atom atom2 : atoms2) {
            Atom atom2small = atom2.copy();
            atom2small.setVanDerWaalsRadius();
            atom2small.setVanDerWaalsRadius(atom2small.getVanDerWaalsRadius()
                                            +
                                            Constants.SOLVENT_RADIUS);
            // if atom2 is not solvent accessible than leave atom2cells list
            // empty.
            GridCell atom2cell = atomGrid.get(atom2small);
            // given a distance file, a residue pair might have a larger
            // distance then the size of the grid, in which case a NULLPOINTER
            // error would occur. Create a dummy grid cell in those cases.
            if (atom2cell == null) {
                atom2cell = new GridCell(atom2.getXYZ(),
                                         GridCell.getSize());
                atom2cell.setIndices(new Point3i(Integer.MAX_VALUE,
                                                 Integer.MAX_VALUE,
                                                 Integer.MAX_VALUE));
            }
            atom2cell.reset();
            atom2cells.add(atom2cell);
            ArrayList <GridCell> atom2Cells =
                                           atomGrid.getAllGridCells(atom2small);
            for (GridCell cell2 : atom2Cells) {
                // set all cells that are occupied by this atom to unoccupied.
                cell2.reset();
            }
        }
        this.grid = atomGrid;
        this.sourceCell = atomGrid.get(atom1);
        this.targetCells = atom2cells;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a list of Path objects, where each path corresponds to a single
     * source-target distance measure.
     * @param maxDist
     *        - float value representing the maximum allowed distance between
     *          source and target.
     * @return List of Path objects, one for each solvent path distance
     *         calculation. Paths with single cells are once with failed
     *         distance calculations. An empty path array is returned, if the
     *         sourceCell is located within a closed cavity.
     */
    public final ArrayList < Path > getShortestPath(final float maxDist) {
        // initialize distance calculation
        BreadthFirstSearch shortestPathAlgo = new BreadthFirstSearch(
                                                               this.grid,
                                                               this.sourceCell,
                                                               this.targetCells,
                                                               maxDist
                                                                    );
        ArrayList < Path > paths = shortestPathAlgo.findShortestPath();
        if (!shortestPathAlgo.hasSucceeded()) {
            return new ArrayList < Path >();
        }
        return paths;
    }
    //--------------------------------------------------------------------------
    /**
     * Extract the distance of a target cell to its source cell from a Path
     * object.
     * @param path
     *        - Path object holding the source and target grid cells and all
     *          grid cells inbetween.
     * @return float value representing the distance extracted from the target
     *         grid cell.
     */
    public static float extractTargetDistances(final Path path) {
        return path.get(
                        BreadthFirstSearch.CELL_NO_OF_TARGET_CELL_IN_PATH
                       ).getDistance();
    }
    //--------------------------------------------------------------------------
}

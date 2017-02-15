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
import structure.math.Mathematics;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;


/**
 * Class for creating and handling grid objects for atom coordinates.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class AtomGrid extends Grid {

    //-------------------------------------------------------------------------
    /**
     * AtomList object holding all atoms for which the grid should be created.
     */
    private AtomList atoms;
    //-------------------------------------------------------------------------
    /**
     * Constructor.
     * @param atomList
     *        - AtomList object on which Grid should be build upon.
     * @param gridCellSize
     *        - float value representing the cell edge length of each grid cell
     * @param offSet
     *        - float value by which grid should be in addition increased in
     *          size
     */
    public AtomGrid(final AtomList atomList,
                    final float gridCellSize,
                    final float offSet) {
        super(MatterUtilities.getMinimumCooridnate(atomList).add(-offSet,
                                                                 -offSet,
                                                                 -offSet),
              MatterUtilities.getMaximumCooridnate(atomList).add(offSet,
                                                                 offSet,
                                                                 offSet),
              gridCellSize);
        this.atoms = atomList;
        this.setOccupancy();
    }
    //-------------------------------------------------------------------------
    /**
     * Constructor.
     * @param atomList
     *        - AtomList object on which Grid should be build upon.
     * @param atom
     *        - Atom object around which the grid is to be generated.
     * @param size
     *        - maximum size of local grid
     * @param gridCellSize
     *        - float value representing the cell edge length of each grid cell
     */
    public AtomGrid(final AtomList atomList,
                    final Atom atom,
                    final float size,
                    final float gridCellSize) {
        super(atom.getXYZ().add(-size - 1, -size - 1, -size - 1),
              atom.getXYZ().add(size + 1, size + 1, size + 1),
              gridCellSize);

        this.atoms = atomList;
        this.setOccupancy();
    }
    //-------------------------------------------------------------------------
    /**
     * Sets the occupied flag for all grid cells in a grid for a list of atom
     * coordinates.
     */
    private void setOccupancy() {

        for (Atom atom : this.atoms) {
            ArrayList < GridCell > cells = this.getAllGridCells(atom);
            for (GridCell cell : cells) {
                cell.setOccupation();
            }
        }
    }

    //-------------------------------------------------------------------------
    /**
     * Returns the grid cell that is closest to an atom.
     * @param atom
     *        - Atom object to which the closest grid cell should be returned.
     * @return GridCell object that is closest to the atoms centre. If atom lies
     *         outside the grid, then {@code NULL} is returned.
     */
    public final GridCell get(final Atom atom) {

        int i = (int) ((atom.getXYZ().getX() - this.getMin().getX())
                       /
                       GridCell.getSize());

        int j = (int) ((atom.getXYZ().getY() - this.getMin().getY())
                       /
                       GridCell.getSize());

        int k = (int) ((atom.getXYZ().getZ() - this.getMin().getZ())
                       /
                       GridCell.getSize());

        if (i > 0 && i < this.getNumberOfCells().getI()
            && j > 0 && j < this.getNumberOfCells().getJ()
            && k > 0 && k < this.getNumberOfCells().getK()) {

            GridCell cell = this.get(i, j, k);

            return cell;
        } else {
            return null;
        }
    }

    //-------------------------------------------------------------------------
    /**
     * Returns all grid cell that are occupied by the entire van der Waals atom
     * shell of an atom.
     * @param atom
     *        - Atom object to which all occupied grid cells should be returned.
     * @return List of GridCell objects that are occupying the atom.
     */
    public final ArrayList < GridCell > getAllGridCells(final Atom atom) {
        return this.getAllGridCells(atom, 0);
    }
    //-------------------------------------------------------------------------
    /**
     * Returns all grid cell that are occupied by the entire van der Waals atom
     * shell of an atom.
     * @param atom
     *        - Atom object to which all occupied grid cells should be returned.
     * @param expand
     *        - integer value representing the number of shells from which
     *          neighborhood is to be extracted.
     * @return List of GridCell objects that are occupying the atom.
     */
    public final ArrayList < GridCell > getAllGridCells(final Atom atom,
                                                        final int expand) {
        ArrayList < GridCell > allCells = new ArrayList < GridCell >();

        float radius = atom.getVanDerWaalsRadius();

        GridCell centre = this.get(atom);

        // if atom should lie outside the grid, then just return an empty list
        // of GridCell objects.
        if (centre != null) {
          int expand2 = GridUtilities.getNumberOfGridCellsFittingIntoHemisphere(
                                                              radius,
                                                              GridCell.getSize()
                                                                              )
                       + expand;

           ArrayList < GridCell > neighboursCube =
                                             GridUtilities.getNeighbouringCells(
                                                                         centre,
                                                                         this,
                                                                         expand2
                                                                              );

           // check that grid cell really is located within atom vdW sphere.
           ArrayList < GridCell > neighbours = new ArrayList < GridCell >();
           for (GridCell neighbour : neighboursCube) {
               float dist = Mathematics.distance(
                                                  atom.getXYZ(),
                                                  neighbour.getXYZ()
                                                  );
               if (dist - radius - (GridCell.getSize() * expand) < 0) {
                   neighbours.add(neighbour);
               }
           }

           allCells.add(centre);
           allCells.addAll(neighbours);
        }
    return allCells;
    }
    //-------------------------------------------------------------------------
    /**
     * Checks whether the Atom object occupies any grid cell.
     * @param atom
     *        - Atom object to be checks for existence in the grid.
     * @return {@code TRUE} if atom has been found in the grid,
     *         {@code FALSE} otherwise.
     */
    public final boolean contains(final Atom atom) {
        ArrayList < GridCell > cells = this.getAllGridCells(atom);
        return cells.size() > 0 ? true : false;
    }

    //-------------------------------------------------------------------------
    /**
     * Resets the value, visited and occupied status of the grid cells that are
     * occupied by a list of atoms.
     * @param atomList
     *        - AtomList object holding the atoms, which associated grid cells
     *          should be reset.
     */
    public final void reset(final AtomList atomList) {
        for (int i = 0; i < atomList.size(); i++) {
             Atom atom = atomList.get(i);
             ArrayList < GridCell > cells = this.getAllGridCells(atom);
             for (GridCell cell : cells) {
                  cell.reset();
             }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns this grid in PDB String format with HEADER information extracted
     * from atom object.
     * @param atom
     *        - Atom object
     * @return String object holding this grid in PDB format.
     */
    public final String toString(final Atom atom) {
        return  "HEADER "
        + atom.getSerialNumber()
        + atom.getResidueName()
        + atom.getResidueNumber()
        + atom.getChainId()
        + Constants.LINE_SEPERATOR
        + this.toString()
        + "END"
        + Constants.LINE_SEPERATOR;

    }
}

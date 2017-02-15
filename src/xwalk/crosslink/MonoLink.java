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

package xwalk.crosslink;

import structure.constants.Constants;
import structure.matter.Atom;
import structure.matter.protein.AminoAcid;


/**
 * Class for representing Cross-Link objects.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class MonoLink extends Atom {

    /**
     * Ranking index of this cross-link within a list of cross-links.
     * Default value is -1.
     */
    private int index = -1;
    /**
     * String object holding the path to the PDB file in which the cross-link
     * has been found.
     * Default value is an empty string.
     */
    private String fileName = "";
    /**
     * boolean object holding information about the solvent accessibility of
     * this mono-link atom.
     * Default value is false.
     */
    private boolean isSolventAccessible = false;
    //--------------------------------------------------------------------------
    /**
     * Sets all atom related information of this mono-link to the ones of the
     * atom object.
     * @param  atom
     *        - Atom object from which all atom related information will be
     *        extracted.
     */
    public final void set(final Atom atom) {
        this.setChainId(atom.getChainId());
        this.setElement(atom.getElement());
        this.setFlag(atom.getFlag());
        this.setName(atom.getName());
        this.setXYZ(atom.getXYZ());
        this.setOccupancy(atom.getOccupancy());
        this.setRank(atom.getRank());
        this.setResidueName(atom.getResidueName());
        this.setResidueNumber(atom.getResidueNumber());
        this.setSerialNumber(atom.getSerialNumber());
        this.setTemperatureFactor(atom.getTemperatureFactor());
        this.setType(atom.getType());
        this.setVanDerWaalsRadius(atom.getVanDerWaalsRadius());
        this.setWeight(atom.getWeight());
        this.setXlogP(atom.getXlogP());
        this.setChargeState(atom.getChargeState());
        this.setICode(atom.getICode());
    }
    //--------------------------------------------------------------------------
    /**
     * Creates a copy of this MonoLink object.
     * @return MonoLink object being a copy of this monoLink object.
     */
    public final MonoLink copy() {
        Atom atom = super.copy();
        MonoLink copy = new MonoLink();
        copy.set(atom);
        copy.setFileName(this.fileName + "");
        copy.setSolventAccessibility(this.isSolventAccessible);
        copy.setIndex(this.index);
    return copy;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether a second mono-link has the same residue name and residue
     * number.
     * @param  monoLink
     *        - MonoLink object to be compared to this MonoLink object.
     * @return {@code TRUE} if both MonoLink object are equal in homology,
     * {@code FALSE} otherwise.
     */
    public final boolean equalsInHomolog(final MonoLink monoLink) {
        String residueId1 = this.getResidueName() + "#"
                            + this.getResidueNumber();
        String residueId2 = monoLink.getResidueName() + "#"
                            + monoLink.getResidueNumber();
        if (residueId1.equals(residueId2)) {
            return true;
        }
        return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether a second mono-link has the same atom identifier, i.e.
     * residue name, residue number, chain Id and atom name.
     * @param monoLink
     *        - MonoLink object to be compared to this MonoLink object.
     * @return {@code TRUE} if both MonoLink object are equal in all atom
     *         identifier, {@code FALSE} otherwise.
     */
    public final boolean equals(final MonoLink monoLink) {
        if (this.equalsInHomolog(monoLink)) {
            String residueId1 = this.getChainId() + "#"
                                + this.getName().trim();
            String residueId2 = monoLink.getChainId() + "#"
                                + monoLink.getName().trim();

            if (residueId1.equals(residueId2)) {
                return true;
            }
        }
        return false;
    }

    //--------------------------------------------------------------------------
    /*
     * Returns the hash value of this mono-link, which incorporates information
     * of its residue name, residue number, chainID and atom name.
     * @return integer variable representing the hash value of this MonoLink
     *         object.
     */
//    public int hashCode() {
//        return this.getResidueName().hashCode()
//             + this.getResidueNumber()
//            + this.getChainId()
//             + this.getName().hashCode();
//    }
    //--------------------------------------------------------------------------
    /**
     * Sets the ranking index of this mono-link.
     * @param rank
     *        - integer value representing the ranking index.
     */
    public final void setIndex(final int rank) {
        this.index = rank;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the ranking index of this mono-link.
     * @return Integer value representing the ranking index.
     */
    public final int getIndex() {
        return index;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the path to the file in which this cross-link has been found.
     * @param fileName
     *        - String object holding the path to the file.
     */
    public final void setFileName(final String fileName) {
        this.fileName = fileName;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the path to the file in which this cross-link has been found.
     * @return String object holding the path to the file.
     */
    public final String getFileName() {
        return fileName;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the solvent accessibility property of this mono-link atom.
     * @param isSolventAccessible
     *        - boolean object holding information about whether this mono-link
     *        atom is solvent accessible
     */
    public final void setSolventAccessibility(
                                           final boolean isSolventAccessible) {
        this.isSolventAccessible = isSolventAccessible;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether this mono-link is solvent accessible or not.
     * @return boolean object holding information about the solvent
     * accessibility of this mono-link atom.
     */
    public final boolean isSolventAccessible() {
        return isSolventAccessible;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a String representation of this mono-link in distance file
     * format.
     * @return String object holding the representation of this mono-link in
     * distance file format.
     */
    public String toString() {
        String atomId = AminoAcid.getAminoAcidId(this)
                        + "-"
                        + this.getName().trim();

        StringBuffer output = new StringBuffer();

        output.append(this.fileName + "\t" + atomId + "\t");
        if (this.isSolventAccessible) {
            output.append("1");
        } else {
            output.append("0");
        }
        output.append(Constants.LINE_SEPERATOR);
    return output.toString();
    }
}

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

import java.util.Comparator;

import xwalk.crosslink.CrossLink;

/**
 * Distance comparator class for CrossLink object, which uses the Solvent-Path
 * distance to compare two CrossLink objects or if not existent uses the
 * Euclidean distance.
 * @author Abdullah Kahraman
 * @version 0.1
 * @version 0.1
 */
public class DistanceComparator implements Comparator < CrossLink > {
    /**
     * Compares two cross-link objects by their Probabilities to occur in
     * experiments or if not assessed by their Solvent-Path distance or if not
     * existent Euclidean distance.
     * @param xl1
     *        - First CrossLink Object.
     * @param xl2
     *        - Second CrossLink Object.
     * @return -1 if the first has a shorter distance/higher probability/lower
     *            chainID/lower residue number.<br>
     *          0 if both objects span equal distances/probabilities/equal
     *            chainID/equal residue number.<br>
     *          1 if the second has a larger distance/lower probability/higher
     *            chainID/higher residue number
     */
    public final int compare(final CrossLink xl1, final CrossLink xl2) {
        if (xl1.getSolventPathDistanceProbability()
                <
            xl2.getSolventPathDistanceProbability()) {
            return 1;
        }
        if (xl1.getSolventPathDistanceProbability()
                >
            xl2.getSolventPathDistanceProbability()) {
            return -1;
        }
        if (xl1.getEuclideanDistanceProbability()
                <
            xl2.getEuclideanDistanceProbability()) {
            return 1;
        }
        if (xl1.getEuclideanDistanceProbability()
                >
            xl2.getEuclideanDistanceProbability()) {
            return -1;
        }

        // cross-links that do not conform to maxsub
        // will have a SASD of -1. The default
        // SASD is -0.5.
        if (xl1.getSolventPathDistance() < -0.9
                &&
           xl2.getSolventPathDistance() > -0.6) {
            return 1;
        }
        if (xl1.getSolventPathDistance() > -0.6
                &&
           xl2.getSolventPathDistance() < -0.9) {
            return -1;
        }

        if (xl1.getEuclideanDistance() < 0
                &&
           xl2.getEuclideanDistance() > 0) {
            return 1;
        }
        if (xl1.getEuclideanDistance() > 0
                &&
           xl2.getEuclideanDistance() < 0) {
            return -1;
        }

        if (xl1.getSolventPathDistance() < xl2.getSolventPathDistance()) {
            return -1;
        }
        if (xl1.getSolventPathDistance() > xl2.getSolventPathDistance()) {
            return 1;
        }
        if (xl1.getEuclideanDistance() < xl2.getEuclideanDistance()) {
           return -1;
        }
        if (xl1.getEuclideanDistance() > xl2.getEuclideanDistance()) {
            return 1;
        }

        // if distances are equal, then sort further by chain ID.
        if (xl1.getPreAtom().getChainId() < xl2.getPreAtom().getChainId()) {
            return -1;
        }
        if (xl1.getPreAtom().getChainId() > xl2.getPreAtom().getChainId()) {
            return 1;
        }
        if (xl1.getPostAtom().getChainId() < xl2.getPostAtom().getChainId()) {
            return -1;
        }
        if (xl1.getPostAtom().getChainId() > xl2.getPostAtom().getChainId()) {
            return 1;
        }

        // if still no sorting could be achieved then sort by residue number
        if (xl1.getPreAtom().getResidueNumber()
            <
            xl2.getPreAtom().getResidueNumber()) {
            return -1;
        }
        if (xl1.getPreAtom().getResidueNumber()
            >
            xl2.getPreAtom().getResidueNumber()) {
            return 1;
        }
        if (xl1.getPostAtom().getResidueNumber()
            <
            xl2.getPostAtom().getResidueNumber()) {
            return -1;
        }
        if (xl1.getPostAtom().getResidueNumber()
            >
            xl2.getPostAtom().getResidueNumber()) {
            return 1;
        }
        return 0;
    }
}

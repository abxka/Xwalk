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

package xwalk.io;

import java.io.IOException;

import structure.constants.Constants;
import structure.io.ReadFile;
import structure.matter.Atom;
import xwalk.crosslink.CrossLink;
import xwalk.crosslink.CrossLinkList;
import xwalk.crosslink.MonoLink;
import xwalk.crosslink.MonoLinkList;

/**
 * This class converts distance files into CrossLink objects.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class DistanceReader {
    /**
     * Constructor.
     */
    protected DistanceReader() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Reads all cross-links from fileName and converts these into CrossLink
     * objects. Mono-links are ignored.
     * @param fileName
     *        - String object holding the path to a distance file.
     * @param onlyIntraLinks
     *        - boolean object determining whether only intra links should be
     *          read out.
     * @param onlyInterLinks
     *        - boolean object determining whether only inter links should be
     *          read out.
     * @param verbose
     *        - boolean object determining whether information should be output
     *        while distance file is read in.
     * @return List of CrossLink objects extracted from the distance file.
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     */
    public static CrossLinkList getCrossLinks(final String fileName,
                                              final boolean onlyIntraLinks,
                                              final boolean onlyInterLinks,
                                              final boolean verbose)
                                                            throws IOException {
        CrossLinkList set = new CrossLinkList();

        ReadFile read = new ReadFile(fileName);
        for (String line : read) {
            if (!line.startsWith("#") && line.trim().length() >= 1) {
                Atom atom1 = new Atom();
                Atom atom2 = new Atom();
                String file = "";
                int index = 0;
                int seqDist = -1;
                float eucDist = -1;
                float solvDist = -1;
                try {
                    String[] array = line.trim().split("\t");
                    // if this line holds a mono-link than skip it.
                    if (array.length < 4) {
                        continue;
                    }
                    if (array.length == 4
                        &&
                        (array[3].equals("1") || array[3].equals("0"))){
                        continue;
                    }

                    index = Integer.parseInt(array[0]);
                    file = array[1];
                    String atom1info = array[2];
                    String atom2info = array[3];
                    // the support of mono links requires the distance file
                    // to be checked for mono and cross-link information

// no distance information should be read in from a distance file
// which will leave the distance information at -1, which can be used
// to recognize whether the cross-link has been found in the current structure.
/*                    if (array.length > 4) {
                        seqDist = Integer.parseInt(array[4]);
                    }
                    if (array.length > 5) {
                        eucDist = Float.parseFloat(array[5]);
                    }
                    if (array.length > 6) {
                        solvDist = Float.parseFloat(
                                                     Value.DISTANCE.getDefault()
                                                     );
                        if (!array[6].equals("-")) {
                            solvDist = Float.parseFloat(array[6]);
                        }
                    }
                    if (array.length > 7) {
                        String peptideSequence = array[7];
                    }
*/
                    array = atom1info.split("-");
                    if (array.length < 2) {
                        if (verbose) {
                            System.err.println("WARNING: First atom of "
                                             + "cross-link number " + index
                                             + " must list a "
                                             + "residue name and residue "
                                             + "number.");
                        }
                    } else {
                        atom1.setResidueName(array[0].trim());
                        atom1.setResidueNumber(
                                               Integer.parseInt(array[1].trim())
                                              );
                    }
                    if (array.length >= 3) {
                        if (array[2].trim().length() != 0) {
                            atom1.setChainId(array[2].trim().charAt(0) == '_'
                                                    ? ' ' : array[2].charAt(0));
                        }
                    }
                    if (array.length >= 4) {
                            atom1.setName(array[3].trim());
                    }
                    if (array.length >= 5) {
                        atom1.setAlternativeLocation(array[4].trim().charAt(0));
                    }


                    array = atom2info.split("-");
                    if (array.length < 2) {
                        if (verbose) {
                            System.err.println("WARNING: Second atom of"
                                             + "cross-link "
                                             + "number " + index + " must list "
                                             + "a residue name and residue "
                                             + "number.");
                        }
                    } else {
                        atom2.setResidueName(array[0].trim().toUpperCase());
                        atom2.setResidueNumber(
                                               Integer.parseInt(array[1].trim())
                                              );
                    }
                    if (array.length >= 3) {
                        if (array[2].trim().length() != 0) {
                            atom2.setChainId(array[2].trim().charAt(0) == '_'
                                                    ? ' ' : array[2].charAt(0));
                            if (onlyInterLinks) {
                               if (atom1.getChainId() == atom2.getChainId()) {
                                   continue;
                               }
                            }
                            if (onlyIntraLinks) {
                                if (atom1.getChainId() != atom2.getChainId()) {
                                    continue;
                                }
                             }
                        }
                    }
                    if (array.length >= 4) {
                        atom2.setName(array[3].trim());
                    }

                    if (array.length >= 5) {
                        atom2.setAlternativeLocation(array[4].trim().charAt(0));
                    }
                } catch (Exception e) {
                    System.err.println("WARNING: Distance file \"" + fileName
                                     + "\" does not conform to distance file "
                                     + "format" + Constants.LINE_SEPERATOR
                                     + e.getMessage());
                }

                CrossLink crossLink = new CrossLink(atom1, atom2, seqDist,
                                                    eucDist);
                crossLink.setFileName(file);
                crossLink.setSolventPathDistance(solvDist);
                crossLink.setIndex(index);

                if (set.get(atom1, atom2) == null
                    &&
                    set.get(atom2, atom1) == null) {
                        set.add(crossLink);
                } else {
                    if (verbose) {
                        System.err.print("WARNING: Following cross-link is "
                                       + "redundant and will be ignored. "
                                       + "Consider to be more specific with "
                                       + "the atom information: " + crossLink);
                    }
                }
            }
        }
    return set;
    }
    //--------------------------------------------------------------------------
    /**
     * Reads in all mono-links from fileName and converts these into a list of
     * Atom objects.
     * @param fileName
     *        - String object holding the path to a distance file.
     * @return AtomList objects with atoms extracted from the distance file.
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     */
    public static MonoLinkList getMonoLinks(final String fileName)
                                                            throws IOException {
        MonoLinkList set = new MonoLinkList();

        ReadFile read = new ReadFile(fileName);
        for (String line : read) {
            if (!line.startsWith("#") && line.trim().length() >= 1) {
                MonoLink monoLink = new MonoLink();
                String file = "";
                int index = 0;
                try {
                    String[] array = line.trim().split("\t");
                    // a mono-link can have only three columns in a distance
                    // file.
                    if (array.length > 4 || array.length < 3) {
                        continue;
                    }
                    if (array.length == 4
                        &&
                        !array[3].equals("1") && !array[3].equals("0")){
                        continue;
                    }
                    index = Integer.parseInt(array[0]);
                    monoLink.setIndex(index);
                    file = array[1];

                    monoLink.setFileName(file);

                    String atominfo = array[2];
                    array = atominfo.split("-");
                    if (array.length < 2) {
                        System.err.println("WARNING: Atom of mono-link "
                                         + "number " + index + " must list a "
                                         + "residue name and residue "
                                         + "number.");
                    } else {
                        monoLink.setResidueName(array[0].trim());
                        monoLink.setResidueNumber(
                                               Integer.parseInt(array[1].trim())
                                              );
                    }
                    if (array.length >= 3) {
                        if (array[2].trim().length() != 0) {
                            monoLink.setChainId(array[2].trim().charAt(0) == '_'
                                                    ? ' ' : array[2].charAt(0));
                        }
                    }
                    if (array.length >= 4) {
                        monoLink.setName(array[3].trim());
                    }
                } catch (Exception e) {
                    System.err.println("WARNING: Distance file \"" + fileName
                                     + "\" does not conform to distance file "
                                     + "format" + Constants.LINE_SEPERATOR
                                     + e.getMessage());
                }

                if (set.get(monoLink) == null) {
                    set.add(monoLink);
                } else {
                    System.err.print("WARNING: Following monolink is "
                                   + "redundant and will be ignored. "
                                   + "Consider to be more specific with "
                                   + "the atom information: " + monoLink);
                }
            }
        }
    return set;
    }


}

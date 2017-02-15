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

package mm.evolution;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Iterator;
import java.util.regex.Pattern;

import structure.constants.Constants;
import structure.io.ReadFile;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;

/**
 * Class for downloading and assigning evolutionary conservation grades for a
 * PDB structures.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Consurf {
    //--------------------------------------------------------------------------
    /**
     * ReadFile object of a ConSurf file.
     */
    private ReadFile consurfFile;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param consurfFilePath
     *        String object holding the path to the consurf file
     */
    public Consurf(final String consurfFilePath) {
        try {
            consurfFile = new ReadFile(consurfFilePath);
        } catch (IOException e) {
            System.err.print(e + Constants.LINE_SEPERATOR);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param br
     *        BufferedReader object of ConSurf file.
     */
    public Consurf(final BufferedReader br) {
        try {
            consurfFile = new ReadFile(br);
        } catch (IOException e) {
            System.err.print(e + Constants.LINE_SEPERATOR);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Assigns the conservation grade for all amino acids in a protein
     * molecule as given by the consurfFile.
     * @param protein
     *        Protein molecule to which ConSurf conservation grades will be
     *        mapped at.
     */
    public final void assignConservation(final PolyPeptide protein) {

        boolean dataFollows = false;
        Pattern pattern = Pattern.compile("[\t]+");
        for (Iterator<String> i = this.consurfFile.iterator(); i.hasNext();) {

            String lineBuffer = i.next().trim();

            if (lineBuffer.indexOf("POS") != -1
                &&
                lineBuffer.indexOf("RESIDUE VARIETY") != -1) {

                lineBuffer = i.next().trim();
                if (lineBuffer.indexOf("(normalized)") != -1) {
                    dataFollows = true;
                    continue;
                } else {
                    dataFollows = true;
                }
            }

            if (dataFollows && lineBuffer.length() > 0) {
                if (!lineBuffer.substring(0, 1).matches("\\d")) {
                    continue;
                }

                String[] lineArray = pattern.split(lineBuffer);

                String tmp = lineArray[2].replaceAll(" ", "");
                String[] seqPosChainId = tmp.split(":");

                char chainId = ' ';
                int seqPos = -1;

                // in some cases there is no ":" in the ATOM column, but it has
                // to be still controlled if there is a information about the AA
                // position of the conservation grade. Therefore this if clause
                // is separated from the if clause below.
                if (seqPosChainId[0].length() > 3) {
                    seqPos = Integer.parseInt(seqPosChainId[0].substring(
                                                      3,
                                                      seqPosChainId[0].length())
                                                                        );
                } else {
                    continue;
                }

                // if there is a colon in the residue column and after the colon
                // is a letter take the letter as the chain ID
                if (seqPosChainId.length > 1 && tmp.indexOf(":") != -1) {
                    chainId = seqPosChainId[1].charAt(0);
                }

                int consGrade = Integer.parseInt(lineArray[4].replaceAll("[ *]",
                                                                         ""
                                                                        ));

                for (AminoAcid aa : protein) {
                    if (aa.getAtom(0).getChainId() == chainId
                        &&
                        aa.getAtom(0).getResidueNumber() == seqPos) {
                            aa.setConservationGrade(consGrade);
                    }
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the chain ID listed in the ConSurf file.
     * @return character representing the chain ID listed in the ConsSurf file.
     */
    public final char getChainId() {

        char chainId = ' ';

        boolean dataFollows = false;
        Pattern pattern = Pattern.compile("[\t]+");
        for (Iterator<String> i = this.consurfFile.iterator(); i.hasNext();) {

            String lineBuffer = i.next().trim();

            if (lineBuffer.indexOf("POS") != -1
                &&
                lineBuffer.indexOf("RESIDUE VARIETY") != -1) {

                lineBuffer = i.next().trim();
                if (lineBuffer.indexOf("(normalized)") != -1) {
                    dataFollows = true;
                    continue;
                } else {
                    dataFollows = true;
                }
            }

            if (dataFollows && lineBuffer.length() > 0) {
                if (!lineBuffer.substring(0, 1).matches("\\d")) {
                    continue;
                }

                String[] lineArray = pattern.split(lineBuffer);

                String tmp = lineArray[2].replaceAll(" ", "");
                String[] seqPosChainId = tmp.split(":");

                // if there is a colon in the residue column and after the colon
                // is a letter take the letter as the chain ID
                if (seqPosChainId.length > 1 && tmp.indexOf(":") != -1) {
                    chainId = seqPosChainId[1].charAt(0);
                }
            }
        }
        return chainId;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines the conservation grade for all amino acids in a protein
     * molecule by downloading necessary conservation information from the
     * ConSurf database.
     * @param protein
     *        Protein molecule to which ConSurf conservation grades will be
     *        mapped at.
     * @param pdbId
     *        String object holding the PDB id to the protein object, which is
     *        necessary for downloading the protein objects associated ConSurf
     *        file to the protein.
     * @throws IOException if no conservation grade file can be found for the
     *         protein on <a href="http://consurfdb.tau.ac.il/consurf_db/">
     *         ConSurfDB</a>.
     */
    public static void setConservation(final PolyPeptide protein,
                                       final String pdbId) throws IOException {

        char chainId = protein.get(0).getAtom(0).getChainId();
        if (chainId == ' ') {
                chainId = '_';
        }
        String urlPath = "http://consurfdb.tau.ac.il/consurf_db/"
                        + pdbId.toUpperCase() + "/" + chainId
                        + "/consurf.grades";
        Consurf consurf = null;
        URL url = new URL(urlPath);
        BufferedReader br = new BufferedReader(new InputStreamReader(
                                                                url.openStream()
                                                                    ));
        consurf = new Consurf(br);
        consurf.assignConservation(protein);
    }
}

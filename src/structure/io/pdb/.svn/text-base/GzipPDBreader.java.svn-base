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

package structure.io.pdb;

import java.io.IOException;
import java.util.zip.DataFormatException;

import structure.exceptions.FileFormatException;
import structure.io.GzipFileReader;
/**
 * Class handles GNU zipped PDB files.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class GzipPDBreader {
    //--------------------------------------------------------------------------
    /**
     * Path to the gzipped PDB file.
     */
    private String fileName;
    /**
     * Constructor.
     * @param path
     *        - String object holding the path to the Gzipped PDB file.
     */
    //--------------------------------------------------------------------------
    public GzipPDBreader(final String path) {
        this.fileName = path;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a PDBreader object holding all atom coordinates of the PDB file
     * that was gzipped.
     * @return PDBreader object.
     * @throws IOException if an error occurs while reading in the gzipped file.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         <a href="http://www.wwpdb.org/documentation/format32/sect9.html">
     *         PDB standards</a>.
     * @throws DataFormatException if .gz format is unknown.
     */
    public final PDBreader getPDBreader() throws IOException,
                                                 FileFormatException,
                                                 DataFormatException {
        PDBreader reader = new PDBreader(new GzipFileReader(
                                                            this.fileName
                                                           ).getBufferedReader()
                                        );
        reader.setFileName(this.fileName);
    return reader;
    }
    //--------------------------------------------------------------------------
}

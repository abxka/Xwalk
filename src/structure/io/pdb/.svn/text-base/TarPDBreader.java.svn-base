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

package structure.io.pdb;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Hashtable;

import com.ice.tar.TarInputStream;

import structure.exceptions.FileFormatException;
import structure.io.ReadFile;
import structure.io.TarFileReader;

/**
 * Class handles Tar compressed PDB files.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class TarPDBreader {
    //--------------------------------------------------------------------------
    /**
     * Path to the tar-ed PDB file.
     */
    private String fileName = "";
    //--------------------------------------------------------------------------
    /**
     * InputStream object holding a stream to a tar-ed file.
     */
    private InputStream inputStream;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param path
     *        - String object holding the path to the Gzipped PDB file.
     */
    public TarPDBreader(final String path) {
        this.fileName = path;
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param tarStream
     *        - An InputStream
     */
    public TarPDBreader(final InputStream tarStream) {
        this.inputStream = tarStream;
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
     */
    public final ArrayList < PDBreader > getPDBreaders() throws IOException,
                                                           FileFormatException {

        ArrayList < PDBreader > pdbReaders = new ArrayList < PDBreader >();

        TarFileReader tarReader = null;
        if (this.inputStream == null) {
            tarReader = new TarFileReader(fileName);
        } else {
            tarReader = new TarFileReader(new TarInputStream(this.inputStream));
        }
        Hashtable < String, ReadFile>  tarContent = tarReader.unTar();
        for (String filePath : tarContent.keySet()) {
            PDBreader pdbReader = new PDBreader(tarContent.get(filePath));
            pdbReader.setFileName(filePath);
            pdbReaders.add(pdbReader);
        }
        return pdbReaders;
    }
    //--------------------------------------------------------------------------
}

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

package structure.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

/**
 * Class to read gzipped files.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class GzipFileReader {

    /**
     * Input stream holding the stream to the gzipped file.
     */
    private GZIPInputStream gz;

    /**
     * Constructor.
     * @param fileName
     *        - String object holding the path to the gzipped file.
     * @throws IOException if reading file into GZIPInputStream fails.
     */
    public GzipFileReader(final String fileName) throws IOException {
        this.gz = new GZIPInputStream(new FileInputStream(fileName));
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param gzInputStream
     *        - GZIPInputStream object holding the stream to the gzipped file.
     */
    public GzipFileReader(final GZIPInputStream gzInputStream) {
        this.gz = gzInputStream;
    }
    //--------------------------------------------------------------------------
    /**
     * Return the GZIPInputStream object of this reader.
     * @return GZIPInputStream object holding the stream to the gzipped file.
     */
    public final GZIPInputStream getGZIPInputStream() {
        return this.gz;
    }
    //--------------------------------------------------------------------------
    /**
     * Converts the GzipInputStream object to a BufferedReader object.
     * @return BufferedReader object.
     */
    public final BufferedReader getBufferedReader() {
        InputStreamReader in = new InputStreamReader(this.gz);
        return new BufferedReader(in);

    }
}


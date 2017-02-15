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

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import structure.constants.Constants;

import com.ice.tar.TarEntry;
import com.ice.tar.TarInputStream;

/**
 * Class for reading in tar-ed files.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class TarFileReader {

    /**
     * TarInputStream object from which the content of the tar-ed file
     * will be read from.
     */
    private TarInputStream tin;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param tarInputStream
     *        - TarInputStream object holding the content of the tar-ed file.
     */
    public TarFileReader(final TarInputStream tarInputStream) {
        tin = tarInputStream;
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param fileName
     *        - Path to the tarred file.
     * @throws FileNotFoundException if file could not be found.
     */
    public TarFileReader(final String fileName) throws FileNotFoundException {
        tin = new TarInputStream(new FileInputStream(fileName));
    }
    //--------------------------------------------------------------------------
    // The TarInputStream object itself is of no use to the current Xwalk
    // implementation. Therefore I have commented it out for the moment.
    ///**
    // * Returns the InputStream of the tar-ed file.
    // * @return TarInputStream object holding the content of the tar-ed file.
    // */
    //public final TarInputStream getTarInputStream() {
    //    return this.tin;
    //}
    //--------------------------------------------------------------------------
    /**
     * Un-tar operation.
     * @return Hashtable with content of the Tar file, where the keys are the
     *         file names and the elements are the file contents.
     * @throws IOException if an Error occurs while reading the TarInputStream
     */
    public final Hashtable < String, ReadFile> unTar() throws IOException {
        Hashtable < String, ReadFile> entries =
                          new Hashtable < String, ReadFile>();
        TarEntry te = this.tin.getNextEntry();
          while (te != null) {
              int buf;
              StringBuffer buffer = new StringBuffer();
              ArrayList < String > lines = new ArrayList < String >();
              while ((buf = this.tin.read()) != -1) {
                  buffer.append((char) buf);
                  if (buf == Constants.LINE_SEPERATOR.charAt(0)) {
                      lines.add(buffer.toString());
                      buffer = new StringBuffer();
                  }
              }
              entries.put(te.getName(), new ReadFile(lines));
              te = this.tin.getNextEntry();
          }
    return entries;
    }
    //--------------------------------------------------------------------------
}


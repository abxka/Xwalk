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

import java.util.Vector;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.File;

import structure.constants.Constants;


/**
 * Class for reading text files into a list of String objects.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class ReadFile extends ArrayList < String > {
    //--------------------------------------------------------------------------
    /**
     * default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    //--------------------------------------------------------------------------
    /**
     * Path of file to be read in.
     */
    private File file;
    //--------------------------------------------------------------------------
    /**
     * System-dependent new line carrier.
     */
    private final String nL = Constants.LINE_SEPERATOR;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param fileName
     *         - String object holding the path of the file to be read in.
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     */
    public ReadFile(final String fileName) throws IOException {
        this.read(fileName);
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param br
     *         - BufferedReader object from which to read in the file content
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     */
    public ReadFile(final BufferedReader br) throws IOException {
        this.read(br);
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param fileContent
     *         - List of String objects representing each a line in a file.
     */
    public ReadFile(final ArrayList < String > fileContent) {
        this.addAll(fileContent);
    }
    //--------------------------------------------------------------------------
    /**
     * Reads files from BufferedReader object.
     * @param br
     *        - BufferedReader object from which to read in the file content
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     * @see #read(String)
     */
    private void read(final BufferedReader br) throws IOException {
        String lineBuffer;
        while ((lineBuffer = br.readLine()) != null) {
            this.add(lineBuffer);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Reads files either from the local hard disk or from within the class
     * path.
     * @param fileName
     *        - String object holding the path of the file to be read in.
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     * @see #read(BufferedReader)
     */
    private void read(final String fileName) throws IOException {
        String lineBuffer;
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(
                                                                        fileName
                                                                              );
        BufferedReader br = null;
        if (ins == null) {
            br = new BufferedReader(new FileReader(fileName));
        } else {
            br = new BufferedReader(new InputStreamReader(ins));
        }
        while ((lineBuffer = br.readLine()) != null) {
            this.add(lineBuffer + this.nL);
        }
        br.close();
        this.file = new File(fileName);
    }
    //--------------------------------------------------------------------------
    /**
     * @param fileName
     *         - String object holding the path of the file to be read in.
     * @return {@code TRUE} if file exists, otherwise {@code FALSE}
     * @see #dirExists(String)
     */
    public static boolean exists(final String fileName) {
        try {
            BufferedReader br = new BufferedReader(
                    new FileReader(fileName));
            while ((br.readLine()) != null) { };
            br.close();
            return true;
        } catch (Exception e) {
            return false;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Checks for the existence of a directory.
     * @param dirPath
     *        - String object holding the path of a directory to check for
     *          existence.
     * @return {@code TRUE} if directory exists, otherwise {@code FALSE}
     * @see    #exists(String)
     */
    public static boolean dirExists(final String dirPath) {
        boolean exists = (new File(dirPath)).exists();
        if (exists) {
            return true;
        }
        return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the file content in a single String.
     * @return String object holding the entire content of a file
     */
    public final String toString() {
        StringBuffer fileContent = new StringBuffer();
        for (String line : this) {
            fileContent.append(line);
        }
    return fileContent.toString();
    } // End of method read()
    //--------------------------------------------------------------------------
    /**
     * Extract a specific column from the file.
     * @param columnNumber
     *        - integer number of column to be returned. 1st column is 0.
     * @param delim
     *        - String representing the delimiter with which the columns will be
     *          distinguished
     * @return Vector of strings with each element being the text in the column
     *         per line.
     * @throws ArrayIndexOutOfBoundsException if columnNumber is larger than the
     *         number of columns in the file minus 1.
     */
    public final Vector < String > getColumn(final int columnNumber,
                                             final String delim) throws
                                                ArrayIndexOutOfBoundsException {
        Vector < String > column = new Vector < String >();
        for (String line : this) {
             String[] lineArray = line.split(delim);
             if (lineArray.length <= columnNumber) {
                throw new ArrayIndexOutOfBoundsException("WARNING: Can't find "
                                                   + "column no.\""
                                                   + columnNumber + "\". "
                                                   + "Only \""
                                                   + (lineArray.length - 1)
                                                   + "\" columns exists "
                                                   + "in file \""
                                                   + this.file.getAbsolutePath()
                                                   + "\"" + this.nL);
             } else {
                column.add(lineArray[columnNumber]);
            }
        }
    return column;
    }
    //--------------------------------------------------------------------------
    /**
     * Rename a file.
     * @param newName
     *        - String object holding the new path and filename of the file.
     * @return {@code TRUE} if file renaming is successful, {@code FALSE}
     *         otherwise.
     */
    public final boolean rename(final String newName) {
        WriteFile file = new WriteFile();
        file.setFile(newName);
        if (!file.write(this.toString())) {
            System.err.print("Could not create file with new name\""
                             + newName + "\"" + this.nL);
            return false;
        }
        File file2delete = this.file;
        this.file = new File(newName);
        return file2delete.delete();
    }
    //--------------------------------------------------------------------------
    /**
     * Copy file.
     * @param newName
     *        - String object holding the path and filename of the file copy
     * @return {@code TRUE} if file copying is successful, {@code FALSE}
     *         otherwise.
     */
    public final boolean copy(final String newName) {
        WriteFile file = new WriteFile();
        file.setFile(newName);
        return file.write(this.toString());
    }
}

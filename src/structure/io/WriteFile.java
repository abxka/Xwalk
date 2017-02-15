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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;

import structure.constants.Constants;



/**
 * This class provides methods to write out String objects into files onto hard
 * drives.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class WriteFile {
    /**
     * BufferedWrite object to write out String objects.
     */
    private BufferedWriter fileWriter;
    /**
     * Path of the file to be created on the hard drive.
     */
    private String filePath;
    /**
     * System-dependent new line carrier.
     */
    private final String nL = Constants.LINE_SEPERATOR;

    //--------------------------------------------------------------------------

    /**
     * Sets the FileWriter object.
     * In order to create a file, this method is one of the first to be invoked.
     * @param fileName
     *        - String holding the path to the file to be created.
     * @return {@code TRUE} if creating file is successful, {@code FALSE} if
     *         IOException is thrown and caught.
     * @see #setFile(String, boolean)
     */
    public final boolean setFile(final String fileName) {
        return this.setFile(fileName, false);
    }
    //--------------------------------------------------------------------------

    /**
     * Sets the FileWriter object.
     * In order to create a file, this method is one of the first to be invoked.
     * @param fileName
     *        - String holding the path to the file to be created.
     * @param append
     *        - If {@code TRUE} the new content will be appended, if
     *          {@code FALSE} a new file will be created overwriting any
     *          existing file.
     * @return {@code TRUE} if creating/writing file is successful,
     *         {@code FALSE} if IOException is thrown and caught.
     * @see #setFile(String)
     */
    public final boolean setFile(final String fileName, final boolean append) {
        try {
            fileWriter = new BufferedWriter(new FileWriter(fileName, append));
        } catch (IOException e) {
            System.err.print("ERROR: Could not create file with file name \""
                           + fileName + "\"" + this.nL + e.getMessage()
                           + this.nL);
            return false;
        }
    return true;
    }
    //--------------------------------------------------------------------------

    /**
     * Writes some data into the file.
     * @param data
     *           - String holding the information to be written into the file.
     * @return {@code TRUE} if creating/writing file is successful,
     *         {@code FALSE} if IOException is thrown and caught.
     */
    public final boolean write(final String data) {
        try {
            fileWriter.write(data);
            fileWriter.close();
            fileWriter = null;
        } catch (IOException e) {
            System.err.print("ERROR: Could not write into file \""
                           + filePath + "\"" + this.nL + e.getMessage()
                           + this.nL);
            return false;
        }
    return true;
    }
    //--------------------------------------------------------------------------

    /**
     * Creates a directory.
     * @param dir
     *        - String holding the path to the directory to be created.
     * @return {@code TRUE} if creating directory is successful, {@code FALSE}
     *         otherwise.
     * @see #deleteDir(String)
     */
    public static boolean createDir(final String dir) {
        return new File(dir).mkdir();
    }
    //--------------------------------------------------------------------------

    /**
     * Deletes a directory with all its sub-directories.
     * @param dirPath
     *        - String holding the path of the directory to be deleted.
     * @return {@code TRUE} if all deletions are successful, {@code FALSE}
     *          otherwise.
     * @see #createDir(String)
     * @see #deleteFile(String)
     */
    public static boolean deleteDir(final String dirPath) {
        File dir = new File(dirPath);
        if (dir.isDirectory()) {
            String[] children = dir.list();
            for (String child : children) {
                boolean success = new File(dirPath
                                         + File.separator
                                         + child
                                          ).delete();
                if (!success) {
                    return false;
                }
            }
        }
        // The directory is now empty, thus can be deleted.
        return dir.delete();
    }
    //--------------------------------------------------------------------------

    /**
     * Deletes a file.
     * @param filePath
     *        - String holding the path of the file to be deleted.
     * @return {@code TRUE} if deletion is successful, {@code FALSE} otherwise.
     * @see #deleteDir(String)
     */
    public static boolean deleteFile(final String filePath) {
        File file = new File(filePath);
        if (file.isFile()) {
            return file.delete();
        } else {
            return false;
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Return the filename of the file to be written.
     * @return String object holding the name of the file that was set with
     *          the class's setFile() method.
     */
    public final String getFileName() {
        return filePath;
    }
}

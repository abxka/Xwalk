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

package external;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import structure.constants.Constants;

/**
 * Class for checking the status of the external applications.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
class StreamGobbler extends Thread {
    /**
     * InputStream of external application.
     */
    private InputStream is;
    /**
     * String representation InputStream.
     */
    private StringBuffer put = new StringBuffer();

    /**
     * Constructor.
     * @param inputStream
     *        InputStream of external application.
     */
    public StreamGobbler(final InputStream inputStream) {
        this.is = inputStream;
    }

    /**
     * Read InputStram.
     */
    public void run() {
        try {
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line = null;
            while ((line = br.readLine()) != null) {
                put.append(line + "\n");
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Outputs a string representation of the InputStream object.
     * @return String object holding the InputStream's string representation.
     */
    public String toString() {
        return put.toString();
    }
}
//------------------------------------------------------------------------------
/**
 * Class for running external programs on the commandline.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class ExternCommand {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected ExternCommand() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Executes a command on the commandline shell.
     * @param command
     *        String object holding the command to be executed.
     * @param verbose
     *        if {@code TRUE} than information on the executed command will be
     *        printed on STDERR, {@code FALSE} otherwise.
     * @return String object holding the STDOUT and STDERR outputs of the
     *         external command.
     */
    public static final String execute(final String command,
                                       final boolean verbose) {

        // Runtime environment for the process
        Runtime runtime = Runtime.getRuntime();

        String nL = Constants.LINE_SEPERATOR;
        try {
            if (verbose) {
                System.err.print("Executing: " + command + nL);
            }
            // External application/process to be executed.
            Process process = runtime.exec(command);

            // any error message?
            StreamGobbler error = new StreamGobbler(process.getErrorStream());
            // any output?
            StreamGobbler output = new StreamGobbler(process.getInputStream());
            // kick them off
            error.start();
            output.run();

            // any error???
            process.waitFor();
            return "STDOUT: " + Constants.LINE_SEPERATOR + output.toString()
                  + nL
                  + "STDERR:" + Constants.LINE_SEPERATOR + error.toString()
                  + nL;
        } catch (Throwable t) {
            t.printStackTrace();
            return "";
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether command execution causes any error. Can be used to check
     * for the existence of an application.
     * @param command
     *        String object holding the command to be executed.
     * @return {@code TRUE} if command executes on the system, {@code FALSE}
     *         otherwise.
     */
    public static final boolean exists(final String command) {
        Runtime runtime = Runtime.getRuntime();
        Process process = null;
        try {
            process = runtime.exec(command);
            process.destroy();
            process = null;
            runtime = null;
            return true;
        } catch (Exception e) {
            runtime = null;
            process.destroy();
            process = null;
            return false;
        }
    }
    //--------------------------------------------------------------------------
}

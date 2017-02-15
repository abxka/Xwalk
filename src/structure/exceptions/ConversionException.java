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

package structure.exceptions;

/**
 * Class that is thrown by the Mathematics class if inconsistencies occur during
 * transformation from Cartesian to spherical coordinates and vice versa.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class ConversionException extends Exception {

    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param errorMessage
     *        - String object holding the error message, ideally some text that
     *          indicates the origin of the error.
     */
    public ConversionException(final String errorMessage) {
        super(errorMessage);
    }
}

////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

package jing.rxnSys;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;

/**
 * The Logger class encapsulating logging functionality for RMG.
 *
 * During a normal RMG job, logging information will be printed both to the
 * console (stdout and stderr) and to a file. The amount of detail printed in
 * each can be set independently.
 *
 * This class contains only static methods (since we only one one instance of
 * the logger at runtime). To use, call the initialize
 */
public class Logger {

    // These constants map integers to various levels of detail
    // Larger integers correspond to more detail
    public static final int CRITICAL = 0;       // Critical (fatal) errors
    public static final int ERROR = 10;         // Regular (non-fatal) errors
    public static final int WARNING = 20;       // Warnings
    public static final int INFO = 30;          // Normal information
    public static final int VERBOSE = 40;       // Detailed information
    public static final int DEBUG = 50;         // Debug information

    /** The level of detail to use for log messages printed to stdout. */
    private static int consoleLevel = INFO;

    /** The level of detail to use for log messages printed to the file. */
    private static int fileLevel = VERBOSE;

    /** The object representing the log file. */
    private static BufferedWriter logFile = null;

    /**
     * Initialize the logger. The log file will be opened; if this is not
     * successful, the program will abort.
     */
    public static void initialize() {
        try {
            // Open the log file (throws IOException if unsuccessful)
            logFile = new BufferedWriter(new FileWriter("RMG.log"));
        }
        catch (IOException e) {
            // Log information is important, so we better stop if we're not
            // saving any!
            System.out.println("Unable to open file \"RMG.log\" for logging.");
            System.exit(0);
        }
    }

    /**
     * Finish the logger. The log file will be closed; if this is not
     * successful, the program will abort.
     */
    public static void finish() {
        try {
            // Close the log file (throws IOException if unsuccessful)
            logFile.close();
        }
        catch (IOException e) {

        }
    }

    /**
     * Set the level of detail to use for log messages printed to the console
     * (stdout and stderr). Generally you should try to use one of the
     * predefined levels if possible, but this is not required.
     * @param level The level of detail to use for log messages printed to the console
     */
    public static void setConsoleLevel(int level) {
        consoleLevel = level;
    }

    /**
     * Set the level of detail to use for log messages printed to the log file.
     * Generally you should try to use one of the predefined levels if
     * possible, but this is not required.
     * @param level The level of detail to use for log messages printed to the log file
     */
    public static void setFileLevel(int level) {
        fileLevel = level;
    }

    /**
     * Log a message with a specified level of detail. Wrapper methods have
     * been provided for cases when level is one of the predefined values; you
     * generally should use those instead.
     * @param level The level of detail to use for this message
     * @param message The message to log
     */
    public static void log(int level, String message) {
        // Do nothing if neither log target accepts messages at that level of detail
        if (consoleLevel < level && fileLevel < level)
            return;

        // Prefix the message if an error or warning
        if (level == CRITICAL)
            message = "CRITICAL: " + message;
        else if (level > CRITICAL && level <= ERROR)
            message = "ERROR: " + message;
        else if (level > ERROR && level <= WARNING)
            message = "Warning: " + message;

        // Print the log message to the console
        // Errors and warnings are printed to stderr, others to stdout
        if (consoleLevel >= level) {
            if (level > WARNING)
                System.out.println(message);
            else
                System.err.println(message);
        }
        
        // Print the log message to the file
        try {
            if (logFile != null && fileLevel >= level)
                logFile.write(message + "\n");
        }
        catch (IOException e) {
            // What should we do here?
        }

    }

    /**
     * Log a message describing a critical (fatal) error. Note that this will 
     * *not* cause the program to halt, but simply defines a unique level for 
     * fatal errors. The prefix "CRITICAL: " will automatically be added to
     * the message.
     * @param message The message to log
     */
    public static void critical(String message) {
        log(CRITICAL, message);
    }

    /**
     * Log a message describing a general error. The prefix "ERROR: " will
     * automatically be added to the message.
     * @param message The message to log
     */
    public static void error(String message) {
        log(ERROR, message);
    }

    /**
     * Log a message describing a warning. The prefix "Warning: " will
     * automatically be added to the message.
     * @param message The message to log
     */
    public static void warning(String message) {
        log(WARNING, message);
    }

    /**
     * Log an informational message. This is the method to use for "normal"
     * log messages.
     * @param message The message to log
     */
    public static void info(String message) {
        log(INFO, message);
    }

    /**
     * Log a detailed informational message.
     * @param message The message to log
     */
    public static void verbose(String message) {
        log(VERBOSE, message);
    }

    /**
     * Log a message containing debug information.
     * @param message The message to log
     */
    public static void debug(String message) {
        log(DEBUG, message);
    }

}

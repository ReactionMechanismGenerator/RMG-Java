// //////////////////////////////////////////////////////////////////////////////
//
// RMG - Reaction Mechanism Generator
//
// Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
// RMG Team (rmg_dev@mit.edu)
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// //////////////////////////////////////////////////////////////////////////////
package jing.rxnSys;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import jing.param.VersionInfo;

/**
 * The Logger class encapsulating logging functionality for RMG. During a normal RMG job, logging information will be
 * printed both to the console (stdout and stderr) and to a file. The amount of detail printed in each can be set
 * independently. This class contains only static methods (since we only one one instance of the logger at runtime). To
 * use, call the initialize
 */
public class Logger {
    // These constants map integers to various levels of detail
    // Larger integers correspond to more detail
    public static final int CRITICAL = 0; // Critical (fatal) errors
    public static final int ERROR = 10; // Regular (non-fatal) errors
    public static final int WARNING = 20; // Warnings
    public static final int INFO = 30; // Normal information
    public static final int VERBOSE = 40; // Detailed information
    public static final int DEBUG = 50; // Debug information
    /** The level of detail to use for log messages printed to stdout. */
    private static int consoleLevel = INFO;
    /** The level of detail to use for log messages printed to the file. */
    private static int fileLevel = VERBOSE;
    /** The object representing the log file. */
    private static BufferedWriter logFile = null;
    /** The newline character to use. */
    private static String newLine = System.getProperty("line.separator");

    /**
     * Initialize the logger. The log file will be opened; if this is not successful, the program will abort. If called
     * with no log file path, then "RMG.log" is assumed.
     */
    public static void initialize() {
        initialize("RMG.log");
    }

    /**
     * Initialize the logger. The specified log file will be opened; if this is not successful, the program will abort.
     */
    public static void initialize(String logFilePath) {
        try {
            // Open the log file (throws IOException if unsuccessful)
            logFile = new BufferedWriter(new FileWriter(logFilePath));
        } catch (IOException e) {
            // Log information is important, so we better stop if we're not
            // saving any!
            System.out.println(String.format(
                    "Unable to open file \"%s\" for logging.", logFilePath));
            System.exit(0);
        }
        // Set stderr to redirect to stdout
        // At the moment RMG's errors and warnings need to be placed in the
        // context of when they occur, so this is necessary
        System.setErr(System.out);
        // Print an initialization timestamp
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        info("RMG execution initiated at "
                + sdf.format(Calendar.getInstance().getTime()));
        info("");
    }

    /**
     * Finish the logger. The log file will be closed; if this is not successful, the program will abort.
     */
    public static void finish() {
        // Print a termination timestamp
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        info("");
        info("RMG execution terminated at "
                + sdf.format(Calendar.getInstance().getTime()));
        try {
            // Close the log file (throws IOException if unsuccessful)
            logFile.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static void flush() {
        // Flush the log file
        try {
            logFile.flush();
        } catch (IOException e) {
            System.err
                    .println("Couldn't flush RMG.log file. Did you initialize the Logger?");
            throw new RuntimeException(e);
        }
    }

    /**
     * Set the level of detail to use for log messages printed to the console (stdout and stderr). Generally you should
     * try to use one of the predefined levels if possible, but this is not required.
     * 
     * @param level
     *            The level of detail to use for log messages printed to the console
     */
    public static void setConsoleLevel(int level) {
        consoleLevel = level;
    }

    /**
     * Set the level of detail to use for log messages printed to the log file. Generally you should try to use one of
     * the predefined levels if possible, but this is not required.
     * 
     * @param level
     *            The level of detail to use for log messages printed to the log file
     */
    public static void setFileLevel(int level) {
        fileLevel = level;
    }

    /**
     * Log a message with a specified level of detail. Wrapper methods have been provided for cases when level is one of
     * the predefined values; you generally should use those instead.
     * 
     * @param level
     *            The level of detail to use for this message
     * @param message
     *            The message to log
     */
    public static void log(int level, String message) {
        // Do nothing if neither log target accepts messages at that level of detail
        if (consoleLevel < level && fileLevel < level)
            return;
        if (message == null)
            return;
        // Replace all "\n" with the proper platform end-of-line character
        message = message.replaceAll("\n", newLine);
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
                logFile.write(message + newLine);
        } catch (IOException e) {
            // What should we do here?
            throw new RuntimeException(e);
        }
        // If it was an error message, make sure the log file is up to date.
        if (level < WARNING && logFile != null) { // i.e. ERROR and CRITICAL
            flush();
        }
    }

    /**
     * Log a message describing a critical (fatal) error. Note that this will *not* cause the program to halt, but
     * simply defines a unique level for fatal errors. The prefix "CRITICAL: " will automatically be added to the
     * message.
     * 
     * @param message
     *            The message to log
     */
    public static void critical(String message) {
        log(CRITICAL, message);
    }

    /**
     * Log a message describing a general error. The prefix "ERROR: " will automatically be added to the message.
     * 
     * @param message
     *            The message to log
     */
    public static void error(String message) {
        log(ERROR, message);
    }

    /**
     * Log a message describing a warning. The prefix "Warning: " will automatically be added to the message.
     * 
     * @param message
     *            The message to log
     */
    public static void warning(String message) {
        log(WARNING, message);
    }

    /**
     * Log an informational message. This is the method to use for "normal" log messages.
     * 
     * @param message
     *            The message to log
     */
    public static void info(String message) {
        log(INFO, message);
    }

    /**
     * Log a detailed informational message.
     * 
     * @param message
     *            The message to log
     */
    public static void verbose(String message) {
        log(VERBOSE, message);
    }

    /**
     * Log a message containing debug information.
     * 
     * @param message
     *            The message to log
     */
    public static void debug(String message) {
        log(DEBUG, message);
    }

    /**
     * Log the current stack trace to both the console and the log file.
     */
    public static void logStackTrace(Throwable e) {
        String message = e.toString() + "\n";
        StackTraceElement[] stackTrace = e.getStackTrace();
        for (int index = 0; index < stackTrace.length; index++) {
            message += "    at " + stackTrace[index].toString() + "\n";
        }
        error(message);
    }

    /**
     * Logs a header containing information about RMG.
     */
    public static void logHeader() {
        String versionHash = VersionInfo.getVersionHash();
        info("######################################################################");
        info("#                                                                    #");
        info("#                 RMG - Reaction Mechanism Generator                 #");
        info("#                                                                    #");
        info("#                    http://rmg.sourceforge.net/                     #");
        info("#                                                                    #");
        info("#  This java code was compiled by ant at:                            #");
        info(String.format("#    %-60s    #", VersionInfo.getBuildDate()));
        info("#  The git repository was on the branch:                             #");
        info(String.format("#    %-60s    #", VersionInfo.getBranchName()));
        info("#  And at the commit with the hash:                                  #");
        info(String.format("#    %-60s    #", VersionInfo.getVersionHash()));
        info("#                                                                    #");
        info("#  For details visit:                                                #");
        if (versionHash.startsWith("*")) // error messages should start with a *
            info("#    http://github.com/GreenGroup/RMG-Java/                          #");
        else {
            info(String
                    .format("#    http://github.com/GreenGroup/RMG-Java/tree/%-17s    #",
                            versionHash.substring(0, 6)));
            info("#  To see changes since then visit:                                  #");
            info(String
                    .format("#    http://github.com/GreenGroup/RMG-Java/compare/%-6s...master   #",
                            versionHash.substring(0, 6)));
        }
        info("#                                                                    #");
        info("#  Copyright (c) 2002-2011                                           #");
        info("#  Prof. William H. Green and the RMG Team:                          #");
        info("#    Joshua W. Allen, Dr. Robert W. Ashcraft, Dr. Gregory J. Beran,  #");
        info("#    Caleb A. Class, Connie Gao, Dr. C. Franklin Goldsmith,          #");
        info("#    Michael R. Harper, Amrit Jalan, Gregory R. Magoon,              #");
        info("#    Dr. David M. Matheu, Shamel S. Merchant, Jeffrey D. Mo,         #");
        info("#    Sarah Petway, Sumathy Raman, Dr. Sandeep Sharma,                #");
        info("#    Prof. Kevin M. Van Geem, Dr. Jing Song, Dr. John Wen,           #");
        info("#    Dr. Richard H. West, Andrew Wong, Dr. Hsi-Wu Wong,              #");
        info("#    Dr. Paul E. Yelvington, Dr. Joanna Yu                           #");
        info("#                                                                    #");
        info("#  The RMGVE graphical user interface to the RMG database            #");
        info("#  was written by John Robotham.                                     #");
        info("#                                                                    #");
        info("#  This software package incorporates parts of the following         #");
        info("#  software packages:                                                #");
        info("#    DASSL    - Written by Prof. Linda Petzold et al                 #");
        info("#      http://www.cs.ucsb.edu/~cse/software.html                     #");
        info("#    CDK      - Written by Prof. Cristoph Steinbeck et al            #");
        info("#      http://cdk.sourceforge.net/                                   #");
        info("#    InChI    - Available from IUPAC                                 #");
        info("#      http://www.iupac.org/inchi/                                   #");
        info("#    cclib                                                           #");
        info("#      http://cclib.sourceforge.net                                  #");
        info("#                                                                    #");
        info("#  For more information, including how to properly cite this         #");
        info("#  program, see http://rmg.sourceforge.net/.                         #");
        info("#                                                                    #");
        info("######################################################################");
        info("");
    }
}

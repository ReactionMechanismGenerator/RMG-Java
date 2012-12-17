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
package jing.chem;

import java.util.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;
import jing.chem.ThermoData;
import jing.rxnSys.Logger;
import jing.chem.ChemGraph;
import jing.chemParser.ChemParser;

/**
 * This QMTPThermoWriter class will write the results of QMTP calculations to a on the fly generated thermo library to
 * avoid the expensive parsing and saving of all QMTP results.
 */
public class QMLibraryEditor {
    // data fields
    /** The object representing the dictionary. */
    private static File dictionaryFile = null;
    /** The object representing the library */
    private static File libraryFile = null;
    /** The newline character to use. */
    private static String newLine = System.getProperty("line.separator");

    // Constructor
    // methods
    /** Initializes the two files Dictionary.txt and Library.txt */
    public static void initialize() {
        dictionaryFile = new File(System.getProperty("RMG.qmLibraryDir"),
                "Dictionary.txt");
        if (!dictionaryFile.exists()) {
            Logger.info(String
                    .format("Creating new QM Thermo Dictionary file at %s because it does not yet exist.",
                            dictionaryFile.getPath()));
            try {
                dictionaryFile.createNewFile();
            } catch (IOException e) {
                // This is pretty essential to new QMTP regime. Stop if it fails
                Logger.critical("Could not create QM Thermo Dictionary");
                Logger.logStackTrace(e);
                System.exit(1);
            }
        }
        libraryFile = new File(System.getProperty("RMG.qmLibraryDir"),
                "Library.txt");
        if (!libraryFile.exists()) {
            Logger.info(String
                    .format("Creating new QM Thermo Library file at %s because it does not yet exist.",
                            libraryFile.getPath()));
            try {
                libraryFile.createNewFile();
            } catch (IOException e) {
                // This is pretty essential to new QMTP regime. Stop if it fails
                Logger.critical("Could not create QM Thermo Library file");
                Logger.logStackTrace(e);
                System.exit(1);
            }
        }
        if (!dictionaryFile.canWrite())
            throw new RuntimeException(
                    "Can't write to QM Thermo Dictionary file at "
                            + dictionaryFile.getAbsolutePath());
        if (!libraryFile.canWrite())
            throw new RuntimeException(
                    "Can't write to QM Thermo Dictionary file at "
                            + libraryFile.getAbsolutePath());
    }

    /**
     * Adds species to Dictionary.txt and Library.txt
     */
    // Will probably require p_graph also
    // Using InChI as name as per suggestion by gmagoon
    public static void addQMTPThermo(ChemGraph p_graph, String inChI,
            ThermoData qmResult, String qmMethod, String qmProgram) {
        // Create comment for logger. I need to double check if naming convention will interfere with model
        String comment = "Adding " + inChI + " to the QMTP thermo Library";
        Logger.info(comment);
        addLibraryEntry(inChI, qmResult, qmMethod, qmProgram);
        addDictionaryEntry(inChI, p_graph);
    }

    /**
     * Adds an entry to the Library of the QMTP thermo library using the ThermData calculated from mm4 or pm3.
     */
    public static void addLibraryEntry(String inChI, ThermoData qmResult,
            String qmMethod, String qmProgram) {
        // Assumes name in thermo library doesn't have to match the name in RMG
        // Not sure how true this is..figures crossed
        // line will be the entry to add to Library.txt using the format used in RMG thermo librarys
        String line = new String(inChI + " ");
        line = (line + Double.toString(qmResult.getH298()) + " ");
        line = (line + Double.toString(qmResult.getS298()) + " ");
        line = (line + Double.toString(qmResult.getCp300()) + " ");
        line = (line + Double.toString(qmResult.getCp400()) + " ");
        line = (line + Double.toString(qmResult.getCp500()) + " ");
        line = (line + Double.toString(qmResult.getCp600()) + " ");
        line = (line + Double.toString(qmResult.getCp800()) + " ");
        line = (line + Double.toString(qmResult.getCp1000()) + " ");
        line = (line + Double.toString(qmResult.getCp1500()) + " ");
        line = (line + Double.toString(qmResult.getDH()) + " ");
        line = (line + Double.toString(qmResult.getDS()) + " ");
        line = (line + Double.toString(qmResult.getDCp()) + " ");
        // include the source of the calculation
        line = (line + qmProgram + " " + qmMethod + " Calculation");
        // write to Library.txt
        // write to Dictionary.txt. Try up to five times before giving up
        int remainingTries = 5;
        boolean succeeded = false;
        while (remainingTries > 0 && !succeeded) {
            try {
                BufferedWriter libraryFileWriter = new BufferedWriter(
                        new FileWriter(libraryFile, true));
                libraryFileWriter.write(line + newLine);
                libraryFileWriter.close();
                succeeded = true;
            } catch (IOException e) {
                remainingTries -= 1;
                try {
                    Logger.warning(String
                            .format("Couldn't write to %s. Waiting 1s then trying again.",
                                    libraryFile.getPath()));
                    Thread.sleep(1000); // wait 1s in case another copy of RMG has the file open
                } catch (InterruptedException e1) {
                    Thread.currentThread().interrupt();
                }
            }
        }
        if (!succeeded) {
            Logger.error(String.format(
                    "Could not save species %s to QMLibrary library file %s",
                    inChI, libraryFile.getAbsolutePath()));
        }
    }

    /**
     * Adds an entry to the Dictionary.txt of the QMTP thermo library using
     */
    public static void addDictionaryEntry(String inChI, ChemGraph p_graph) {
        // Assumes name in thermo library doesn't have to match the name in RMG
        // Not sure how true this is..figures crossed
        // definition is a string of the name and adjacency list to be added to Dictionary.txt
        int i = 1;
        String definition = new String(inChI + newLine);
        // toString(i) is used as opposed to toString because that seems to give the adjacency list
        // without prefacing it with the chemical formula
        definition = (definition + p_graph.toString(i) + newLine + newLine);
        // write to Dictionary.txt. Try up to five times before giving up
        int remainingTries = 5;
        boolean succeeded = false;
        while (remainingTries > 0 && !succeeded) {
            try {
                BufferedWriter dictionaryFileWriter = new BufferedWriter(
                        new FileWriter(dictionaryFile, true));
                dictionaryFileWriter.write(definition);
                dictionaryFileWriter.close();
                succeeded = true;
            } catch (IOException e) {
                remainingTries -= 1;
                try {
                    Logger.warning(String
                            .format("Couldn't write to %s. Waiting 1s then trying again.",
                                    dictionaryFile.getPath()));
                    Thread.sleep(1000); // wait 1s in case another copy of RMG has the file open
                } catch (InterruptedException e1) {
                    Thread.currentThread().interrupt();
                }
            }
        }
        if (!succeeded) {
            Logger.error(String
                    .format("Could not save species %s to QMLibrary dictionary file %s",
                            inChI, dictionaryFile.getAbsolutePath()));
        }
    }

    /**
     * 
     */
    public static LinkedHashMap readLibrary(String p_qmThermoFileName)
            throws IOException {
        try {
            FileReader in = new FileReader(p_qmThermoFileName);
            BufferedReader data = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(data, true);
            LinkedHashMap<String, ThermoData> library = new LinkedHashMap<String, ThermoData>();
            while (line != null) {
                String[] result = line.split("\\s");
                // Construct information to create instance of ThermData
                String inChI = result[0];
                String thermo = result[1];
                for (int i = 2; i < 13; i++) {
                    thermo = thermo + " " + result[i];
                }
                int listLength = result.length;
                String comments = result[13];
                // check this upper bound might need to be listLength + 1 -nyee
                for (int i = 14; i < listLength; i++) {
                    comments = comments + " " + result[i];
                }
                // Parse out thermoData, and include name and comments
                ThermoData thermoData = ChemParser
                        .parseThermoFromLibrary(thermo);
                ThermoData newThermoData = new ThermoData(inChI, thermoData,
                        comments);
                library.put(inChI, newThermoData);
                line = ChemParser.readMeaningfulLine(data, true);
            }
            in.close();
            return library;
        } catch (FileNotFoundException e) {
            Logger.info("QMTP thermo library file could not be found at "
                    + p_qmThermoFileName);
            throw e;
        } catch (IOException e) {
            Logger.logStackTrace(e);
            throw new IOException("Can't read QMTP thermo library.");
        }
    }
}

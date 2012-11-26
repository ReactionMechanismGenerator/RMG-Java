////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
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

package jing.chem;

import java.util.*;
import java.io.BufferedWriter;
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
 * This QMTPThermoWriter class will write the results of QMTP calculations to a
 * on  the fly generated thermo library to avoid the expensive parsing and saving
 * of all QMTP results.
 */

public class QMLibraryEditor {

    //data fields
    /** The object representing the dictionary. */
    private static BufferedWriter dictionaryFile = null;

    /** The object representing the library     */
    private static BufferedWriter libraryFile = null;

    /** The newline character to use. */
    private static String newLine = System.getProperty("line.separator");

    /** Number of species added to Library used for naming */
    //may have to get rid of this if I find naming is important
    private static int counter = 0;

    //Constructor

    //methods
    /** Initializes the two files Dictionary.txt and Library.txt  */
    public static void initialize() {
        try {
            // Open the log file (throws IOException if unsuccessful)
            dictionaryFile = new BufferedWriter(new FileWriter("QMThermoLibrary/Dictionary.txt", true));
            Logger.info("Creating Dictionary for QM Thermo Library");
        }
        catch (IOException e) {
            // This is pretty essential to new QMTP regime. Stop if it fails
            Logger.critical("Could not create QM Thermo Dictionary");
            System.exit(0);
        }

        try {
            // Open the log file (throws IOException if unsuccessful)
            // Boolean Argument in FileWriter means that it will append to existing files.
            libraryFile = new BufferedWriter(new FileWriter("QMThermoLibrary/Library.txt", true));
            Logger.info("Creating Library for QM Thermo Library");
        }
        catch (IOException e) {
            // This is pretty essential to new QMTP regime. Stop if it fails
            Logger.critical("Could not create QM Thermo Library");
            System.exit(0);
        }

    }

    /**
     * Close the Dictionary and Library when finished. 
     */
    public static void finish() {
        try {
            // Close the log file (throws IOException if unsuccessful) 
            if (dictionaryFile != null) dictionaryFile.close();
            if (libraryFile != null) libraryFile.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Adds species to Dictionary.txt and Library.txt
     */
    //Will probably require p_graph also
    //Using InChI as name as per suggestion by gmagoon
    public static void addQMTPThermo(ChemGraph p_graph, String inChI, ThermoData qmResult, String qmMethod, String qmProgram){
        //Count the number of species added for naming
        counter ++;
        //Create comment for logger. I need to double check if naming convention will interfere with model
        String comment = "Adding " + inChI + " to the QMTP thermo Library";
        //    	String name = "SQM(";
        //    	name.concat(Integer.toString(counter));
        //    	name.concat(")");

        Logger.info(comment);
        //finish implementation of this method
        addLibraryEntry(inChI, qmResult, qmMethod, qmProgram);
        //finish implementation of this method
        addDictionaryEntry(inChI, p_graph);
    }
    /**
     * Adds an entry to the Library of the QMTP thermo library using 
     * the ThermData calculated from mm4 or pm3. 
     */
    public static void addLibraryEntry(String inChI, ThermoData qmResult, String qmMethod, String qmProgram) {
        //Assumes name in thermo library doesn't have to match the name in RMG
        //Not sure how true this is..figures crossed
        //line will be the entry to add to Library.txt using the format used in RMG thermo librarys
        String line = new String(inChI + " ");
        line = (line +Double.toString(qmResult.getH298()) + " ");
        line = (line +Double.toString(qmResult.getS298())+ " ");
        line = (line +Double.toString(qmResult.getCp300())+ " ");
        line = (line +Double.toString(qmResult.getCp400())+ " ");
        line = (line +Double.toString(qmResult.getCp500())+ " ");
        line = (line +Double.toString(qmResult.getCp600())+ " ");
        line = (line +Double.toString(qmResult.getCp800())+ " ");
        line = (line +Double.toString(qmResult.getCp1000())+ " ");
        line = (line +Double.toString(qmResult.getCp1500())+ " ");
        line = (line +Double.toString(qmResult.getDH())+ " ");
        line = (line +Double.toString(qmResult.getDS())+ " ");
        line = (line +Double.toString(qmResult.getDCp())+ " ");

        //include the source of the calculation 
        line = (line + qmProgram + " " + qmMethod + " Calculation");

        //write to Library.txt
        try {
            libraryFile.write(line + newLine);
            libraryFile.flush();
        }
        catch (IOException e) {
            // What should we do here?
            throw new RuntimeException(e);
        }

    }
    /**
     * Adds an entry to the Dictionary.txt of the QMTP thermo library using
     */
    public static void addDictionaryEntry(String inChI, ChemGraph p_graph) {
        //Assumes name in thermo library doesn't have to match the name in RMG
        //Not sure how true this is..figures crossed
        //definition is a string of the name and adjacency list to be added to Dictionary.txt
        int i = 1;
        String definition = new String(inChI + newLine);
        //toString(i) is used as opposed to toString because that seems to give the adjacency list 
        //without prefacing it with the chemical formula
        definition= (definition + p_graph.toString(i) + newLine + newLine);

        //write to Dictionary.txt
        try {
            dictionaryFile.write(definition);
            dictionaryFile.flush();
        }
        catch (IOException e) {
            // What should we do here?
            throw new RuntimeException(e);
        }
    }

    /**
     * Adds an entry to the Dictionary.txt of the QMTP thermo library using
     * @throws IOException 

     */
    public static HashMap readLibrary(String p_qmThermoFileName) throws IOException {
        try {
            FileReader in = new FileReader(p_qmThermoFileName);
            BufferedReader data = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(data, true);
            HashMap <String, ThermoData> library = new HashMap<String, ThermoData>();
            while (line!=null){
                String [] result = line.split("\\s");

                // Construct information to create instance of ThermData
                String inChI = result[0];
                String thermo = result[1];

                for(int i=2; i<13; i++){
                    thermo = thermo + " " + result[i];
                }

                int listLength = result.length;
                String comments = result[13];
                //check this upper bound might need to be listLength + 1 -nyee
                for(int i=14; i < listLength; i++){
                    comments = comments + " " + result[i];
                }

                //Parse out thermoData, and include name and comments
                ThermoData thermoData = ChemParser.parseThermoFromLibrary(thermo);
                ThermoData newThermoData = new ThermoData(inChI, thermoData, comments);

                library.put(inChI, newThermoData);
                line = ChemParser.readMeaningfulLine(data, true);
            }

            in.close();
            return library;

        }
        catch(FileNotFoundException e){
            Logger.info("QMTP thermo library file could not be found at "+p_qmThermoFileName);
            throw e;
        }
        catch(IOException e) {
            Logger.logStackTrace(e);
            throw new IOException("Can't read QMTP thermo library.");
        }
    }
}

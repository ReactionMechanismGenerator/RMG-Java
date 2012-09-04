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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;


import jing.chem.ThermoData;
import jing.rxnSys.Logger;
import jing.chem.ChemGraph;
import jing.chem.QMTP;
import jing.chemUtil.Graph;

/**
 * This QMTPThermoWriter class will write the results of QMTP calculations to a
 * on  the fly generated thermo library to avoid the expensive parsing and saving
 * of all QMTP results.
 */

public class QMTPThermoWriter {
	
	//data fields
    /** The object representing the dictionary. */
    private static BufferedWriter dictionary = null;
    
    /** The object representing the library     */
    private static BufferedWriter library = null;
	
	/** The newline character to use. */
	private static String newLine = System.getProperty("line.separator");
	
	/** The tab character to use. */
	// Is there one using System.getProperty("tab") like above?  -nyee
	private static String tab = ("\t");
	
	/** Number of species added to Library used for naming */
	//may have to get rid of this if I find naming is important
	private static int counter = 0;

	//Constructor
	
	//methods
	/** Initializes the two files Dictionary.txt and Library.txt  */
	public static void initialize() {
		try {
            // Open the log file (throws IOException if unsuccessful)
            dictionary = new BufferedWriter(new FileWriter("QMTPThermoLibrary/Dictionary.txt"));
            Logger.info("Creating Dictionary for QMTP Thermo Library");
        }
        catch (IOException e) {
            // This is pretty essential to new QMTP regime. Stop if it fails
        	Logger.critical("Could not create QMTP Thermo Dictionary");
            System.exit(0);
        }
		
		try {
            // Open the log file (throws IOException if unsuccessful)
            library = new BufferedWriter(new FileWriter("QMTPThermoLibrary/Library.txt"));
            Logger.info("Creating Library for QMTP Thermo Library");
        }
        catch (IOException e) {
        	// This is pretty essential to new QMTP regime. Stop if it fails
        	Logger.critical("Could not create QMTP Thermo Library");
            System.exit(0);
        }
		
	}
	
    /**
     * Close the Dictionary and Library when finished. If unsuccessful will abort 
     * ...I think
     */
    public static void finish() {
        try {
            // Close the log file (throws IOException if unsuccessful)
            dictionary.close();
            library.close();
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
    	String line = new String(inChI + tab);
    	line = (line +Double.toString(qmResult.getH298()) + tab);
    	line = (line +Double.toString(qmResult.getS298())+ tab);
    	line = (line +Double.toString(qmResult.getCp300())+ tab);
    	line = (line +Double.toString(qmResult.getCp400())+ tab);
    	line = (line +Double.toString(qmResult.getCp500())+ tab);
    	line = (line +Double.toString(qmResult.getCp600())+ tab);
    	line = (line +Double.toString(qmResult.getCp800())+ tab);
    	line = (line +Double.toString(qmResult.getCp1000())+ tab);
    	line = (line +Double.toString(qmResult.getCp1500())+ tab);
    	line = (line +Double.toString(qmResult.getDH())+ tab);
    	line = (line +Double.toString(qmResult.getDS())+ tab);
    	line = (line +Double.toString(qmResult.getDCp())+ tab);
    	
    	//include the source of the calculation 
    	line = (line +" Calculated using " + qmMethod + "  in " + qmProgram);
    	
    	//write to Library.txt
        try {
        	library.write(line + newLine);
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
        	dictionary.write(definition);
        }
        catch (IOException e) {
            // What should we do here?
			throw new RuntimeException(e);
        }
    }
}
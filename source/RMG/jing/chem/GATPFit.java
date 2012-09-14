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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import jing.mathTool.MathTool;

import java.io.*;
import jing.chem.*;
import java.util.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.mathTool.*;
import jing.chemUtil.*;
import jing.chemParser.*;
import jing.rxnSys.Logger;
/**
 * Contains methods used to interact with GATPFit.
 */
public class GATPFit {

    private static Process GATPFit; 
    private static BufferedReader errorStream, dataOutput; 
    private static PrintWriter commandInput; 

    static {
      try {
		    String workingDirectory = System.getProperty("RMG.workingDirectory");
        String[] command = {workingDirectory +  "/bin/GATPFit.exe"};
        File runningDir = new File("GATPFit");
        GATPFit = Runtime.getRuntime().exec(command, null, runningDir);
        errorStream = new BufferedReader(new InputStreamReader(GATPFit.getErrorStream()));
        commandInput = new PrintWriter(GATPFit.getOutputStream(), true);
        BufferedInputStream in = new BufferedInputStream(GATPFit.getInputStream());
        dataOutput = new BufferedReader(new InputStreamReader(in));

        Runtime.getRuntime().addShutdownHook(new Thread() {
          public void run() {
            GATPFit.destroy(); 
          }
        } ); 

        Thread Terr = new Thread(new Runnable(){
          public void run(){
            try {
              String errline = errorStream.readLine();
              if (errline!=null){
                String error_message="GATPFit Error: ";
                while (errline!=null){
                  error_message+=errline;
                  errline=errorStream.readLine();
                }
                throw new GATPFitException(error_message);
              }
            } catch (Exception e){
              Logger.logStackTrace(e); 
              throw new GATPFitException(e.toString()); 
            }
          }
        } );
        Terr.start(); 
      } catch (Exception e) {
        Logger.logStackTrace(e);
        String ls = System.getProperty("line.separator");
        String err = "Error running GATPFit" + ls;
        err += e.toString();
        throw new GATPFitException(err);
      }
    }


    private static NASAThermoData callGATPFit(Species species, String p_directory) {

        NASAThermoData nasaThermoData = null;

        if (p_directory == null) throw new NullPointerException("No direcotry specified to run GATPFit in.");

        // Construct GATPFit input
        // write species name
        String ls = System.getProperty("line.separator");
        StringBuilder result = new StringBuilder(1024);
		result.append( "SPEC " + species.getChemkinName() + ls );

        // write the element
        ChemGraph cg = species.getChemGraph();
        int Hn = cg.getHydrogenNumber();
        int Cn = cg.getCarbonNumber();
        int On = cg.getOxygenNumber();
        int Sin = cg.getSiliconNumber();
        int Sn = cg.getSulfurNumber();
        int Cln = cg.getChlorineNumber();

        int numUniqueElements = 0;
        if (Hn > 0) ++numUniqueElements;
        if (Cn > 0) ++numUniqueElements;
        if (On > 0) ++numUniqueElements;
        if (Sin > 0) ++numUniqueElements;
        if (Sn > 0) ++numUniqueElements;
        if (Cln > 0) ++numUniqueElements;

// GATPFit.exe requires at least two elements but no more than five
//        if (numUniqueElements > 4) {
//        	System.err.println("Species contains more than four unique elements.");
//        }

		result.append( "ELEM C " + MathTool.formatInteger(Cn,3,"L") + ls );
		result.append( "ELEM H " + MathTool.formatInteger(Hn,3,"L") + ls );
        if (On>0)  result.append( "ELEM O " + MathTool.formatInteger(On,3,"L") + ls );
        if (Sin>0) result.append( "ELEM Si " + MathTool.formatInteger(Sin,3,"L") + ls );
        if (Sn>0)  result.append( "ELEM S " + MathTool.formatInteger(Sn,3,"L") + ls );
        if (Cln>0) result.append( "ELEM Cl " + MathTool.formatInteger(Cln,3,"L") + ls );

        // write H and S at 298
        ThermoData td = species.getThermoData();
		result.append( "H298 " + Double.toString(td.getH298()) + "\n" );
		result.append( "S298 " + Double.toString(td.getS298()) + "\n" );
		result.append( "DLTH " + Double.toString(td.getH298()) + "\n" );
		result.append( "MWEI " + Double.toString(species.getMolecularWeight()) + "\n" );
		
		//result.append( "H298 " + MathTool.formatDouble(td.getH298(), 10, 2).trim() + ls );
        //result.append( "S298 " + MathTool.formatDouble(td.getS298(), 10, 2).trim() + ls );
        //result.append( "DLTH " + MathTool.formatDouble(td.getH298(), 10, 2).trim() + ls );
        // write MW, temperature, ouput format, etc
        //result.append( "MWEI " + MathTool.formatDouble(species.getMolecularWeight(), 6, 1).trim() + ls );
        result.append( "TEMP 1000.0" + ls );
		//note:the TMIN and TMAX values below only affect the valid range in the CHEMKIN output file;
		// the points used to deteremine the NASA-to-Wilhoit fit are determined independently inside the GATPFit fortran code.
		// Too large a range causes problems fitting the transport data when pre-processing in CHEMKIN.
		// (see https://github.com/GreenGroup/RMG-Java/issues/issue/18)
		result.append( "TMIN 250.0"+ls );
		result.append( "TMAX 5000.0" + ls );
        result.append( "CHEM" + ls );
        result.append( "TEM2 2000.0" + ls );
		if (species.getChemGraph().isLinear())  result.append( "LINEAR" + ls );
		else result.append( "NONLINEAR" + ls );
        result.append( String.valueOf(cg.getAtomNumber()) + ls );
        result.append( String.valueOf(species.getInternalRotor()) + ls );
		result.append( "TECP 300 " + Double.toString(td.Cp300) + "\n" );
		result.append( "TECP 400 " + Double.toString(td.Cp400) + "\n" );
		result.append( "TECP 500 " + Double.toString(td.Cp500) + "\n" );
		result.append( "TECP 600 " + Double.toString(td.Cp600) + "\n" );
		result.append( "TECP 800 " + Double.toString(td.Cp800) + "\n" );
		result.append( "TECP 1000 " + Double.toString(td.Cp1000) + "\n" );
		result.append( "TECP 1500 " + Double.toString(td.Cp1500) + "\n" );
		
		//result.append( "TECP 300 " + MathTool.formatDouble(td.Cp300,10,2).trim() + ls );
		//result.append( "TECP 400 " + MathTool.formatDouble(td.Cp400,10,2).trim() + ls );
		//result.append( "TECP 500 " + MathTool.formatDouble(td.Cp500,10,2).trim() + ls );
		//result.append( "TECP 600 " + MathTool.formatDouble(td.Cp600,10,2).trim() + ls );
		//result.append( "TECP 800 " + MathTool.formatDouble(td.Cp800,10,2).trim() + ls );
		//result.append( "TECP 1000 " + MathTool.formatDouble(td.Cp1000,10,2).trim() +ls );
		//result.append( "TECP 1500 " + MathTool.formatDouble(td.Cp1500,10,2).trim() + ls );
		result.append( "END" + ls );

        // finished writing text for input file, now save result to INPUT.txt
        String GATPFit_input_name = null;
        File GATPFit_input = null;

        // call GATPFit
        final String inputString = result.toString(); 
        boolean error = false;
        try {
        	commandInput.println(inputString);
            commandInput.flush();
            if (commandInput.checkError()) throw new GATPFitException("Error writing input to GATPFit buffer");
            
            String line = dataOutput.readLine();
            if (line==null) {
                throw new GATPFitException("no output from GATPFit");
            }
            line = dataOutput.readLine(); // skip first line (just says "The Chemkin polynomical coefficients calculated:")
            String nasaString = "";
            while ( (line != null) && !(line.contains("GATPFIT_HAS_FINISHED_ONE_INPUT"))) {
                nasaString += line + System.getProperty("line.separator");
                line = dataOutput.readLine();
            }
            //Logger.info(String.format("GATP string read: " + nasaString));
            nasaThermoData = new NASAThermoData(nasaString);
        }
        catch (Exception e) {
            Logger.logStackTrace(e);
            String err = "Error running GATPFit" + ls;
            err += e.toString();
            GATPFit_input_name = "GATPFit/INPUT.txt";
            err += ls + "To help diagnosis, writing GATPFit input to file "+GATPFit_input_name+ls;
            try {
                GATPFit_input = new File(GATPFit_input_name);
                FileWriter fw = new FileWriter(GATPFit_input);
                fw.write(inputString);
                fw.close();
            }
            catch (IOException e2) {
                err+= "Couldn't write to file "+ GATPFit_input_name + ls;
                err += e2.toString();
            }
            throw new GATPFitException(err);
        }

		/*
		// temporarily save all GATPFit files for debugging purposes
		GATPFit_input_name = "GATPFit/INPUT."+species.getChemkinName()+".txt";
		GATPFit_input = new File(GATPFit_input_name);
		try {
			FileWriter fw = new FileWriter(GATPFit_input);
			fw.write(inputString);
			fw.close();		
		}
		catch (IOException e2) {
			String err = "Couldn't write to file "+ GATPFit_input_name + ls;
			err += e2.toString();
			throw new GATPFitException(err);
		}
		 */

        return nasaThermoData;
    }
	
	
    public static NASAThermoData generateNASAThermoData(Species species) {
        // get working directory
        String dir = System.getProperty("RMG.workingDirectory");
        NASAThermoData nasaThermoData = null;
        try {
        	// prepare GATPFit input file and execute system call
            nasaThermoData = callGATPFit(species, dir);
        }
        catch (GATPFitException e) {
        	throw new NASAFittingException("Error in running GATPFit: " + e.toString());
        }
		
		return nasaThermoData;
    }
	
}

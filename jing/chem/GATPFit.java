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



package jing.chem;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import jing.mathTool.MathTool;

/**
 * Contains methods used to interact with GATPFit.
 * @author jwallen
 */
public class GATPFit {

		 //## operation callGATPFit(String)
    private static boolean callGATPFit(Species species, String p_directory) {
        //#[ operation callGATPFit(String)
        if (p_directory == null) throw new NullPointerException();

        // write GATPFit input file
		String workingDirectory = System.getProperty("RMG.workingDirectory");
        // write species name
        String ls = System.getProperty("line.separator");
        String result = "SPEC " + species.getChemkinName() + ls;

        // write the element
        ChemGraph cg = species.getChemGraph();
        int Hn = cg.getHydrogenNumber();
        int Cn = cg.getCarbonNumber();
        int On = cg.getOxygenNumber();
        int Sin = cg.getSiliconNumber();
        int Sn = cg.getSulfurNumber();
        
        int numUniqueElements = 0;
        if (Hn > 0) ++numUniqueElements;
        if (Cn > 0) ++numUniqueElements;
        if (On > 0) ++numUniqueElements;
        if (Sin > 0) ++numUniqueElements;
        if (Sn > 0) ++numUniqueElements;
        
        // GATPFit.exe requires at least two elements but no more than five
//        if (numUniqueElements > 4) {
//        	System.err.println("Species contains more than four unique elements.");
//        }
        
		result += "ELEM C " + MathTool.formatInteger(Cn,3,"L") + ls;
		result += "ELEM H " + MathTool.formatInteger(Hn,3,"L") + ls;
        if (On>0) result += "ELEM O " + MathTool.formatInteger(On,3,"L") + ls;
        if (Sin>0) result += "ELEM Si " + MathTool.formatInteger(Sin,3,"L") + ls;
        if (Sn>0) result += "ELEM S " + MathTool.formatInteger(Sn,3,"L") + ls;
        
        /*if (Cn>0) result += "ELEM C " + MathTool.formatInteger(Cn,3,"L") + ls;
        if (Hn>0) result += "ELEM H " + MathTool.formatInteger(Hn,3,"L") + ls;
        if (On>0) result += "ELEM O " + MathTool.formatInteger(On,3,"L") + ls;*/
		
//		result += "ELEM C " + MathTool.formatInteger(Cn,3,"L") + ls;
//        result += "ELEM H " + MathTool.formatInteger(Hn,3,"L") + ls;
//        result += "ELEM O " + MathTool.formatInteger(On,3,"L") + ls;
		
        // write H and S at 298
        ThermoData td = species.getThermoData();
		result += "H298 " + String.format("%4.2e \n",td.getH298());
		result += "S298 " + String.format("%4.2e \n",td.getS298());
		result += "DLTH " + String.format("%4.2e \n",td.getH298());
		result += "MWEI " + String.format("%6.1e \n",species.getMolecularWeight());
		
		//result += "H298 " + MathTool.formatDouble(td.getH298(), 10, 2).trim() + ls;
        //result += "S298 " + MathTool.formatDouble(td.getS298(), 10, 2).trim() + ls;
		
        //result += "DLTH " + MathTool.formatDouble(td.getH298(), 10, 2).trim() + ls;

        // write MW, temperature, ouput format, etc
        //result += "MWEI " + MathTool.formatDouble(species.getMolecularWeight(), 6, 1).trim() + ls;
        result += "TEMP 1000.0" + ls;
		result += "TMIN 300.0"+ls;
		result += "TMAX 5000.0" + ls;
        result += "CHEM" + ls;
        result += "TEM2 2000.0" + ls;
		if (species.getChemGraph().isLinear())  result += "LINEAR" + ls;
		else result += "NONLINEAR" + ls;
        result += String.valueOf(cg.getAtomNumber()) + ls;
        result += String.valueOf(species.getInternalRotor()) + ls;
		result += "TECP 300 " + String.format("%4.2e \n",td.Cp300);
		result += "TECP 400 " + String.format("%4.2e \n",td.Cp400);
		result += "TECP 500 " + String.format("%4.2e \n",td.Cp500);
		result += "TECP 600 " + String.format("%4.2e \n",td.Cp600);
		result += "TECP 800 " + String.format("%4.2e \n",td.Cp800);
		result += "TECP 1000 " + String.format("%4.2e \n",td.Cp1000);
		result += "TECP 1500 " + String.format("%4.2e \n",td.Cp1500);
		
		//result += "TECP 300 " + MathTool.formatDouble(td.Cp300,10,2).trim() + ls;
		//result += "TECP 400 " + MathTool.formatDouble(td.Cp400,10,2).trim() + ls;
		//result += "TECP 500 " + MathTool.formatDouble(td.Cp500,10,2).trim() + ls;
		//result += "TECP 600 " + MathTool.formatDouble(td.Cp600,10,2).trim() + ls;
		//result += "TECP 800 " + MathTool.formatDouble(td.Cp800,10,2).trim() + ls;
		//result += "TECP 1000 " + MathTool.formatDouble(td.Cp1000,10,2).trim() +ls;
		//result += "TECP 1500 " + MathTool.formatDouble(td.Cp1500,10,2).trim() + ls;
		result += "END" + ls;

        // finished writing text for input file, now save result to fort.1
        String GATPFit_input_name = null;
        File GATPFit_input = null;

        GATPFit_input_name = "GATPFit/INPUT.txt";
        try {
        	GATPFit_input = new File(GATPFit_input_name);
        	FileWriter fw = new FileWriter(GATPFit_input);
        	fw.write(result);
        	fw.close();
        	}
        catch (IOException e) {
        	String err = "GATPFit input file error: " + ls;
        	err += e.toString();
        	throw new GATPFitException(err);
        }

        // call GATPFit
        boolean error = false;
        try {
        	 // system call for therfit
        	String[] command = {workingDirectory +  "/software/GATPFit/GATPFit.exe"};
			File runningDir = new File("GATPFit");
        	Process GATPFit = Runtime.getRuntime().exec(command, null, runningDir);
        	int exitValue = GATPFit.waitFor();
        }
        catch (Exception e) {
        	String err = "Error in running GATPFit!" + ls;
        	err += e.toString();
        	throw new GATPFitException(err);
        }

        // return error = true, if there was a problem
        return error;


        //#]
    }
	
		  //## operation generateNASAThermoData()
    public static NASAThermoData generateNASAThermoData(Species species) {
        //#[ operation generateNASAThermoDatabyGATPFit()
        // get working directory
        String dir = System.getProperty("RMG.workingDirectory");
        
        try {
        	// prepare GATPFit input file and execute system call
        	boolean error = callGATPFit(species, dir);
        }
        catch (GATPFitException e) {
        	throw new NASAFittingException("Error in running GATPFit: " + e.toString());
        }

        // parse output from GATPFit, "output.txt" is the output file name
        String therfit_nasa_output = "GATPFit/OUTPUT.txt";

		NASAThermoData nasaThermoData = null;
		
        try {
        	FileReader in = new FileReader(therfit_nasa_output);
        	BufferedReader data = new BufferedReader(in);

        	String line = data.readLine();
        	line = data.readLine();
        	String nasaString = "";

        	while (line != null) {
        		nasaString += line + System.getProperty("line.separator");
        		line = data.readLine();
        	}
        	nasaThermoData = new NASAThermoData(nasaString);
        }
        catch (Exception e) {
        	throw new NASAFittingException("Error in reading in GATPFit output file: " + System.getProperty("line.separator") + e.toString());
        }
		
		return nasaThermoData;
        //#]
    }
	
}

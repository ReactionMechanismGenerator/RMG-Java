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
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import jing.chemParser.ChemParser;
import jing.mathTool.MathTool;

/**
 * Contains data members and methods for interacting with THERFIT, a module that
 * models the degrees of freedom of a molecule using three frequency-degeneracy
 * pairs developed by Chang, Bozzelli, and Dean.
 * This was the original method of modeling the degrees of freedom, but has
 * since been made obsolete in favor of a more general model contained in the
 * class SpectroscopicData.
 * <p>
 * For more information on THERFIT, see the following reference:
 * <p>
 * A. Y. Chang, J. W. Bozzelli, and A. M. Dean. "Kinetic Analysis of Complex
 * Chemical Activation and Unimolecular Dissociation Reactions using QRRK
 * Theory and the Modified Strong Collision Approximation." Z. Phys. Chem
 * 214 (11), p. 1533-1568 (2000).
 * 
 * @author jwallen
 */
public class Therfit {

		//## operation generateThreeFrequencyModel()
	public static ThreeFrequencyModel generateThreeFrequencyModel(Species species) {
        
		if (species.isTriatomicOrSmaller()) {
			System.out.println("Warning: Attempted to fit three frequency model to triatomic or smaller.");
			return null;
		}
		
        //String directory = System.getProperty("RMG.workingDirectory");
        String mode = "HOE"; // use Harmonic Oscillator Eqn to calc psuedo-freqs
		String dir = "therfit";
        // prepare therfit input file and execute system call
        boolean error = callTherfit(species, dir, mode);
        if (error) {
			System.out.println("Warning: Error detected in THERFIT execution. Three frequency model not fitted.");
			return null;
		}
		
        double [] degeneracy = new double [3];
        double [] frequency = new double [3];

        try {
        	String therfit_output = dir+"/fort.25";

        	FileReader in = new FileReader(therfit_output);
        	BufferedReader data = new BufferedReader(in);

        	String line = ChemParser.readMeaningfulLine(data);
        	if (!line.startsWith(species.getChemkinName())) {
        		System.out.println("Wrong output from therfit!");
        		System.exit(0);
        	}

        	int index = 0;
        	line = ChemParser.readMeaningfulLine(data);
        	while (line != null) {
        		line = line.trim();
        		if (line.startsWith("VIBRATION")) {
        			if (index > 2) break;
        			StringTokenizer token = new StringTokenizer(line);
        			String temp = token.nextToken();
        			temp = token.nextToken();
        			temp = token.nextToken();
        			temp = token.nextToken();
        			String deg = token.nextToken();
        			temp = token.nextToken();
        			temp = token.nextToken();
        	        String freq = token.nextToken();
        	        degeneracy[index] = Double.parseDouble(deg);
        	        frequency[index] = Double.parseDouble(freq);
        	        index++;
        	    }
        		line = ChemParser.readMeaningfulLine(data);
        	}

        	in.close();
        }
        catch (Exception e) {
        	System.out.println("Wrong output fort.25 from therfit!");
        	System.exit(0);
        }
        
		// create threeFrequencyModel for this species
        return new ThreeFrequencyModel(frequency, degeneracy);
    }
	
	private static boolean callTherfit(Species species, String p_directory, String p_mode) {
    	
        if (p_directory == null || p_mode == null) 
			throw new NullPointerException();
		
		String workingDirectory = System.getProperty("RMG.workingDirectory");
        
		// write therfit input file
        String result = p_mode.toUpperCase() + '\n';
        result += species.getChemkinName() + '\n';
        ThermoData td = species.getThermoData();

        result += MathTool.formatDouble(td.getH298(), 8, 2) + '\n';
        result += MathTool.formatDouble(td.getS298(), 8, 2) + '\n';
        result += MathTool.formatDouble(td.getCp300(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp400(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp500(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp600(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp800(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp1000(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp1500(), 7, 3) + '\n';

        ChemGraph cg = species.getChemGraph();
        result += "C\n";
        result += MathTool.formatInteger(cg.getCarbonNumber(),3,"R") + '\n';
        result += "H\n";
        result += MathTool.formatInteger(cg.getHydrogenNumber(),3,"R") + '\n';
        int oNum = cg.getOxygenNumber();
        if (oNum == 0) {
        	result += "0\n0\n";
        }
        else {
        	result += "O\n";
        	result += MathTool.formatInteger(oNum,3,"R") + '\n';
        }
        result += "0\n0\n";

        result += "G\n";
        result += Integer.toString(species.getInternalRotor()) + '\n';

        // finished writing text for input file, now save result to fort.1
        String therfit_input_name = null;
        File therfit_input = null;

		try {
        	// prepare therfit input file, "fort.1" is the input file name
        	therfit_input_name = p_directory+"/fort.1";//p_directory+"/software/therfit/fort.1";
        	therfit_input = new File(therfit_input_name);
        	FileWriter fw = new FileWriter(therfit_input);
        	fw.write(result);
        	fw.close();
        	}
        catch (IOException e) {
        	System.out.println("Wrong input file for therfit!");
        	System.exit(0);
        }

        // call therfit
        boolean error = false;
		try {
       	 // system call for therfit
			String[] command = {workingDirectory + "/bin/therfit.exe"};
			File runningDir = new File(p_directory );//+ "/therfit");// "/software/therfit");
			Process therfit = Runtime.getRuntime().exec(command, null, runningDir);
			InputStream is = therfit.getErrorStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				//System.out.println(line);
				line = line.trim();
				if (!line.contains("Therfit Job Complete")) {
					String speName = species.getName();
					System.out.println("therfit error for species: " + speName+"\n"+species.toString());
					File newfile = new File(therfit_input_name+"."+speName);
					therfit_input.renameTo(newfile);
					error = true;
				}
			}
			int exitValue = therfit.waitFor();
			br.close();
			isr.close();
		}
		catch (Exception e) {
			System.out.println("Error in run therfit!");
			System.exit(0);
		}

        // return error = true, if there was a problem
        return error;

        //#]
    }
	
		   //## operation generateNASAThermoData()
    public static NASAThermoData generateNASAThermoData(Species species) {
        ///#[ operation generateNASAThermoData()
       
		// get working directory
        //String dir = System.getProperty("RMG.workingDirectory");
        String mode = "WILHOI"; // use Wilhoit polynomial to extrapolate Cp
		String dir = "therfit";
        // prepare therfit input file and execute system call
        boolean error = callTherfit(species, dir, mode);

        // parse output from therfit, "fort.25" is the output file name
        if (error) {
        	System.out.println("Error in generating NASA thermodata!");
        	System.exit(0);
        }

		NASAThermoData nasaThermoData = null;
		
        try {
        	String therfit_nasa_output = dir+"/fort.54";

        	FileReader in = new FileReader(therfit_nasa_output);
        	BufferedReader data = new BufferedReader(in);

        	String line = data.readLine();
        	String nasaString = "";

        	while (line != null) {
        		nasaString += line + "\n";
        		line = data.readLine();
        	}
        	nasaThermoData = new NASAThermoData(nasaString);
        }
        catch (Exception e) {
        	System.out.println("Wrong output from therfit!");
        	System.exit(0);
        }

		return nasaThermoData;
        //#]
    }
}

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

package jing.rxn;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.StringReader;
import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

import Jama.Matrix;
import jing.chem.Species;
import jing.chem.SpectroscopicData;
import jing.mathTool.UncertainDouble;
import jing.param.GasConstant;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxnSys.CoreEdgeReactionModel;
import jing.rxnSys.ReactionModelGenerator;
import jing.rxnSys.ReactionSystem;

/**
 * Used to estimate pressure-dependent rate coefficients k(T, P) for a
 * PDepNetwork object. The estimation is done by a call to the Fortran module
 * FAME. There are two methods for estimating k(T, P): the modified strong
 * collision and the reservoir state methods.
 * @author jwallen
 */
public class FastMasterEqn implements PDepKineticsEstimator {

	/**
	 * The number of times the FAME module has been called for any network
	 * since the inception of this RMG execution. Used to be used to number the
	 * FAME input and output files, but now the networks have individual IDs
	 * that are used for this purpose.
	 */
	private static int runCount = 0;

	/**
	 * An enumeraction of pressure-dependent kinetics estimation methods:
	 * <ul>
	 * <li>NONE - The method could not be assessed.
	 * <li>STRONGCOLLISION - The steady state/modified strong collision method of Chang, Bozzelli, and Dean.
	 * <li>RESERVOIRSTATE - The steady state/reservoir state method of N. J. B. Green and Bhatti.
	 * </ul>
	 */
	public enum Mode { NONE, STRONGCOLLISION, RESERVOIRSTATE };

	/**
	 * The mode to use for estimating the pressure-dependent kinetics. The mode
	 * represents the set of approximations used to estimate the k(T, P) values.
	 */
	private Mode mode;

	/**
	 * A list of the temperatures at which the rate coefficient has been
	 * explicitly calculated..
	 */
	private static Temperature[] temperatures;

	/**
	 * A list of the pressures at which the rate coefficient has been
	 * explicitly calculated.
	 */
	private static Pressure[] pressures;
	
	private static int numTBasisFuncs = -1;
	private static int numPBasisFuncs = -1;
	
	/**
	 * The number of grains to use in a fame calculation (written in input file)
	 */
	private static int numGrains = 251;
	/**
	 * boolean, detailing whether the high-P-limit is greater than all of the
	 * 	fame-computed k(T,P)
	 */
	private static boolean pdepRatesOK = true;

	//==========================================================================
	//
	//	Constructors
	//

	/**
	 * Creates a new object with the desired mode. The mode represents the
	 * set of approximations used to estimate the k(T, P) values.
	 * @param m The mode to use for estimating the pressure-dependent kinetics.
	 */
	public FastMasterEqn(Mode m) {
		setMode(m);
	}

	//==========================================================================
	//
	//	Accessors
	//

	/**
	 * Returns the current pressure-dependent kinetics estimation mode.
	 * @return The current pressure-dependent kinetics estimation mode.
	 */
	public Mode getMode() {
		return mode;
	}

	/**
	 * Sets the pressure-dependent kinetics estimation mode.
	 * @param m The new pressure-dependent kinetics estimation mode.
	 */
	public void setMode(Mode m) {
		mode = m;
	}

	//==========================================================================
	//
	//	Other methods
	//

	/**
	 * Executes a pressure-dependent rate coefficient calculation.
	 * @param pdn The pressure-dependent reaction network of interest
	 * @param rxnSystem The reaction system of interest
	 * @param cerm The current core/edge reaction model
	 */
	public void runPDepCalculation(PDepNetwork pdn, ReactionSystem rxnSystem,
			CoreEdgeReactionModel cerm) {

		// No update needed if network is not altered
		if (pdn.getAltered() == false)
			return;

		// Don't do networks with isomers made up of only monatomic species
		boolean shouldContinue = true;
		for (ListIterator iter = pdn.getIsomers().listIterator(); iter.hasNext(); ) {
			PDepIsomer isomer = (PDepIsomer) iter.next();
			boolean allMonatomic = true;
			for (int i = 0; i < isomer.getNumSpecies(); i++) {
				if (!isomer.getSpecies(i).isMonatomic())
					allMonatomic = false;
			}
			if (allMonatomic)
				shouldContinue = false;
		}
		// No update needed if network is only two wells and one reaction (?)
		/*if (pdn.getUniIsomers().size() + pdn.getMultiIsomers().size() == 2 &&
				pdn.getPathReactions().size() == 1)
			shouldContinue = false;*/

		if (!shouldContinue)
		{
			/*LinkedList<PDepReaction> paths = pdn.getPathReactions();
			LinkedList<PDepReaction> net = pdn.getNetReactions();
			net.clear();
			for (int i = 0; i < paths.size(); i++)
				net.add(paths.get(i));*/
			return;
		}

		// Get working directory (to find FAME executable)
		String dir = System.getProperty("RMG.workingDirectory");

		// Determine wells and reactions; skip if no reactions in network
		LinkedList<Species> speciesList = pdn.getSpeciesList();
		LinkedList<PDepReaction> pathReactionList = pdn.getPathReactions();
		if (pathReactionList.size() == 0) {
			System.out.println("Warning: Empty pressure-dependent network detected. Skipping.");
			return;
		}

		// Sort isomers such that the order is:
		//	1. Explored unimolecular isomers
		//	2. Bimolecular reactant/product channels
		//	3. Unexplored unimolecular isomers (treated as product channels)
		LinkedList<PDepIsomer> isomerList = new LinkedList<PDepIsomer>();
		int nIsom = 0, nReac = 0, nProd = 0;
		for (int i = 0; i < pdn.getIsomers().size(); i++) {
			PDepIsomer isom = pdn.getIsomers().get(i);
			if (isom.isUnimolecular() && isom.getIncluded())
				isomerList.add(isom);
		}
		nIsom = isomerList.size();
		if (nIsom == 0) {
			// We need at least one unimolecular isomer in order to perform a
			// P-dep calculation
			for (int i = 0; i < pdn.getIsomers().size(); i++) {
				PDepIsomer isom = pdn.getIsomers().get(i);
				if (isom.isUnimolecular() && !isom.getIncluded())
					isomerList.add(isom);
			}
			nIsom = isomerList.size();
			for (int i = 0; i < pdn.getIsomers().size(); i++) {
				PDepIsomer isom = pdn.getIsomers().get(i);
				if (isom.isMultimolecular())
					isomerList.add(isom);
			}
			nReac = isomerList.size() - nIsom;
			nProd = 0;
		}
		else {
			for (int i = 0; i < pdn.getIsomers().size(); i++) {
				PDepIsomer isom = pdn.getIsomers().get(i);
				if (isom.isMultimolecular())
					isomerList.add(isom);
			}
			nReac = isomerList.size() - nIsom;
			for (int i = 0; i < pdn.getIsomers().size(); i++) {
				PDepIsomer isom = pdn.getIsomers().get(i);
				if (isom.isUnimolecular() && !isom.getIncluded())
					isomerList.add(isom);
			}
			nProd = isomerList.size() - nIsom - nReac;
		}

		// Make sure all species have spectroscopic data
		for (ListIterator<PDepIsomer> iter = isomerList.listIterator(); iter.hasNext(); ) {
			PDepIsomer isomer = iter.next();
			for (int i = 0; i < isomer.getNumSpecies(); i++) {
				Species species = isomer.getSpecies(i);
				if (!species.hasSpectroscopicData())
					species.generateSpectroscopicData();
			}
		}
		
		
		// Run garbage collection before we submit to FAME to prevent jobs from crashing with large P-Dep networks.
		// JDM January 22, 2010
		//if (runTime.freeMemory() < runTime.totalMemory()/3) {
		if (isomerList.size() >= 10) {
			// System.out.println("Number of isomers in PDepNetwork to follow: " + isomerList.size());
			// Find current time BEFORE running garbage collection.
			// JDM January 27, 2010
			// System.out.println("Running Time is: " + String.valueOf((System.currentTimeMillis())/1000/60) + " minutes before garbage collection.");
			Runtime runTime = Runtime.getRuntime();
			runTime.gc();
			// Find current time AFTER running garbage collection.
			// JDM January 27, 2010
			// System.out.println("Running Time is: " + String.valueOf((System.currentTimeMillis())/1000/60) + " minutes after garbage collection.");
		}
		//}

		// Create FAME input files
		String input = writeInputString(pdn, rxnSystem, speciesList, isomerList, pathReactionList, nIsom, nReac, nProd);
        String output = "";

		// DEBUG only: Write input file
		/*try {
			FileWriter fw = new FileWriter(new File("fame/" + Integer.toString(pdn.getID()) + "_input.txt"));
			fw.write(input);
			fw.close();
		} catch (IOException ex) {
			System.out.println("Unable to write FAME input file.");
			System.exit(0);
		}*/

		// Touch FAME output file
		//touchOutputFile(); //no longer needed with standard input/output

		// FAME system call
		try {
           	//System.out.println("Running FAME...");
			String[] command = {dir + "/bin/fame.exe"};
           	File runningDir = new File("fame/");
			Process fame = Runtime.getRuntime().exec(command, null, runningDir);

            BufferedReader stdout = new BufferedReader(new InputStreamReader(fame.getInputStream()));
            BufferedReader stderr = new BufferedReader(new InputStreamReader(fame.getErrorStream()));
			PrintStream stdin = new PrintStream(  new BufferedOutputStream( fame.getOutputStream(), 1024), true);

			stdin.print( input );
			stdin.flush();
			if (stdin.checkError()) { // Flush the stream and check its error state.
				System.out.println("ERROR sending input to fame.exe!!");
			}
			stdin.close();
			
            if (stdout == null)
				throw new PDepException("FAME output file is empty; FAME job was likely unsuccessful.");
			String line = stdout.readLine().trim();
			
			// advance to first line that's not for debug purposes
			while ( line.startsWith("#IN:") || line.contains("#DEBUG:") ) {
				output += line + "\n";
				line = stdout.readLine().trim();
			}
			
            if (line.startsWith("#####")) { // correct output begins with ######...
                output += line + "\n";
                while ( (line = stdout.readLine()) != null) {
                    output += line + "\n";
                }
            }
            else { // erroneous output does not begin with ######...
                // Print FAME stdout and error
				System.out.println("FAME Error::");
                System.out.println(line);
				output += line + "\n";
                while ( ((line = stdout.readLine()) != null) || ((line = stderr.readLine()) != null) ) {
					output += line + "\n";
                    System.out.println(line);
                }
				// clean up i/o streams (we may be trying again with ModifiedStrongCollision and carrying on)
				stdout.close();
				stderr.close();
                throw new Exception();
            }

            int exitValue = fame.waitFor();

			// Clean up i/o streams
			// This may be needed to release memory, which is especially 
			// important for FAME since it can easily be called tens of
			// thousands of times in a single job
			stdout.close();
			stderr.close();

			// Parse FAME output file and update accordingly
			if (parseOutputString(output, pdn, rxnSystem, cerm, isomerList)) {

				// Reset altered flag
				pdn.setAltered(false);

				// Clean up files
				int id = pdn.getID();
				/*String path = "fame/";
				if (id < 10)			path += "000";
				else if (id < 100)	path += "00";
				else if (id < 1000)	path += "0";
				path += Integer.toString(id);

				File input = new File("fame/fame_input.txt");
				File newInput = new File(path +  "_input.txt");
				input.renameTo(newInput);
				File output = new File("fame/fame_output.txt");
				File newOutput = new File(path +  "_output.txt");
				output.renameTo(newOutput);*/

				// Write finished indicator to console
				//System.out.println("FAME execution for network " + Integer.toString(id) + " complete.");
				String formula = pdn.getSpeciesType();
				System.out.println("PDepNetwork #" + Integer.toString(id) +
					" (" + formula + "): " +
					pdn.getNetReactions().size() + " included and " +
					pdn.getNonincludedReactions().size() + " nonincluded net reactions.");
				// Find time required to run FAME for large networks.
				// JDM January 27, 2010
				// if (isomerList.size() >= 10) {
				//	System.out.println("Running Time is: " + String.valueOf((System.currentTimeMillis())/1000/60) + " minutes after running fame.");
				// }


			}

        }
        catch (Exception e) {
			e.printStackTrace();
			
        	System.out.println(e.getMessage());
			// Save bad input to file
            try {
                FileWriter fw = new FileWriter(new File("fame/" + Integer.toString(pdn.getID()) + "_input.txt"));
                fw.write(input);
                fw.close();
				System.out.println("Troublesome FAME input saved to ./fame/" + Integer.toString(pdn.getID()) + "_input.txt");
                FileWriter fwo = new FileWriter(new File("fame/" + Integer.toString(pdn.getID()) + "_output.txt"));
                fwo.write(output);
                fwo.close();
				System.out.println("Troublesome FAME result saved to ./fame/" + Integer.toString(pdn.getID()) + "_output.txt");
            } catch (IOException ex) {
                System.out.println("Unable to save FAME input that caused the error.");
                System.exit(0);
            }
            // If using RS method, fall back to MSC
			if (mode == Mode.RESERVOIRSTATE) {  /// mode is not defined if running modified strong collision
				System.out.println("Falling back to modified strong collision mode for this network.");
				mode = Mode.STRONGCOLLISION;
				runPDepCalculation(pdn, rxnSystem, cerm);
				mode = Mode.RESERVOIRSTATE;
				return;
			}
			else {
				System.out.println("Error running FAME.");
				System.exit(0);				
			}
        }
        
        /*
         * MRH 26Feb2010:
         * Checking whether pdep rates exceed the high-P-limit
         * 
         * Although fame converges, the computed k(T,P) may exceed the high-P-limit,
         * 	due to the number of grains being too small.  If any of the pdep rates
         * 	exceed the high-P-limit by greater than a factor of 2, the "pdepRatesOK" boolean
         * 	is set to false and fame will be re-executed, using an increased number
         * 	of grains
         * After all pdep rates are below the high-P-limit (or the number of grains
         * 	exceeds 1000), we exit the while loop.  The number of grains is then
         * 	reset to 251 (which was the standard before)
         */
        while (!pdepRatesOK) {
        	numGrains = numGrains + 200;
        	runPDepCalculation(pdn, rxnSystem, cerm);
        }
        numGrains = 251;

		runCount++;

	}

	/**
	 * Creates the input file needed by FAME that represents a pressure-
	 * dependent reaction network.
	 * @param pdn The reaction network of interest
	 * @param rxnSystem The reaction system of interest
	 * @param speciesList The set of species in the network
	 * @param isomerList The set of isomers in the network
	 * @param pathReactionList The set of path reactions in the network
	 */
	public String writeInputString(PDepNetwork pdn, ReactionSystem rxnSystem,
			LinkedList<Species> speciesList,
			LinkedList<PDepIsomer> isomerList,
			LinkedList<PDepReaction> pathReactionList,
			int nIsom, int nReac, int nProd) {

		String input = "";

        Temperature stdTemp = new Temperature(298, "K");

		// Collect simulation parameters
        // MRH 28Feb: Commented temperature/pressure lines (variables are never read)
//		Temperature temperature = rxnSystem.getPresentTemperature();
//		Pressure pressure = rxnSystem.getPresentPressure();
		BathGas bathGas = new BathGas(rxnSystem);

		int numChebTempPolys, numChebPressPolys;
		if (numTBasisFuncs == -1) {
			numChebTempPolys = 4;
			numChebPressPolys = 4;
		}
		else {
			numChebTempPolys = numTBasisFuncs;
			numChebPressPolys = numPBasisFuncs;
		}

		// Determine reference energies
		double grainMinEnergy = getGrainMinEnergy(isomerList); // [=] kJ/mol
		double grainMaxEnergy = getGrainMaxEnergy(isomerList, 2100); // [=] kJ/mol
		double grainSize = getGrainSize(grainMinEnergy, grainMaxEnergy, pdn.getNumUniIsomers()); // [=] kJ/mol

		// Create the input string
		try {

			// Header
			input += "################################################################################\n";
			input += "#\n";
			input += "#	FAME input file\n";
			input += "#\n";
			input += "################################################################################\n";
			input += "\n";
			input += "# All syntax in this file is case-insensitive\n";
			input += "\n";

			// Method
			input += "# The method to use to extract the phenomenological rate coefficients k(T, P)\n";
			input += "# 	Options: ModifiedStrongCollision, ReservoirState\n";
			if (mode == Mode.STRONGCOLLISION)			input += "ModifiedStrongCollision\n";
			else if (mode == Mode.RESERVOIRSTATE)		input += "ReservoirState\n";
			else throw new Exception("Unable to determine method to use to estimate phenomenological rate coefficients.");
			input += "\n";

			// Temperatures
			input += "# The temperatures at which to estimate k(T, P)\n";
			input += "# 	First item is the number of temperatures \n";
			input += "#	Second item is the units; options are K, C, F, or R\n";
			input += "#	Remaining items are the temperature values in the specified units\n";
			input += Integer.toString(temperatures.length) + " K\n";
			for (int i = 0; i < temperatures.length; i++)
				input += Double.toString(temperatures[i].getK()) + "\n";
			input += "\n";

			// Pressures
			input += "# The pressures at which to estimate k(T, P)\n";
			input += "# 	First item is the number of pressures \n";
			input += "#	Second item is the units; options are bar, atm, Pa, or torr\n";
			input += "#	Remaining items are the temperature values in the specified units\n";
			input += Integer.toString(pressures.length) + " Pa\n";
			for (int i = 0; i < pressures.length; i++)
				input += Double.toString(pressures[i].getPa()) + "\n";
			input += "\n";

			// Interpolation model
			input += "# The interpolation model to use to fit k(T, P)\n";
			input += "#	Option 1: No interpolation\n";
			input += "#		Example: None\n";
			input += "#	Option 2: Chebyshev polynomials\n";
			input += "#		Option must be accompanied by two numbers, indicating the number of\n";
			input += "#		terms in the Chebyshev polynomials for temperature and pressure, \n";
			input += "#		respectively\n";
			input += "#		Example: Chebyshev 4 4\n";
			input += "#	Option 3: Pressure-dependent Arrhenius\n";
			input += "#		Example: PDepArrhenius\n";
			if (PDepRateConstant.getMode() == PDepRateConstant.Mode.CHEBYSHEV)
				input += "Chebyshev " + Integer.toString(numChebTempPolys) +
						" " + Integer.toString(numChebPressPolys) + "\n";
			else if (PDepRateConstant.getMode() == PDepRateConstant.Mode.PDEPARRHENIUS)
				input += "PDepArrhenius\n";
			else
				input += "None\n";
			input += "\n";

			// Number of energy grains to use (determines to an extent the accuracy and precision of the results)
			input += "# A method for determining the number of energy grains to use\n";
			input += "# 	Option 1: Specifying the number to use directly\n";
			input += "#		Example: NumGrains 201\n";
			input += "# 	Option 2: Specifying the grain size in J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n";
			input += "#		Example: GrainSize J/mol 4.184\n";
			input += "NumGrains " + numGrains;
			input += "\n\n";

			// Collisional transfer probability model to use
			input += "# Collisional transfer probability model\n";
			input += "# 	Option 1: Single exponential down\n";
			input += "#		Option must also be accompanied by unit and value of the parameter\n";
			input += "#		Allowed units are J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n";
			input += "#		Example: SingleExpDown kJ/mol 7.14\n";
			input += "SingleExpDown J/mol " + Double.toString(bathGas.getExpDownParam() * 1000);
			input += "\n\n";

			// Other parameters for bath gas
			input += "# Bath gas parameters\n";
			input += "# 	Molecular weight; allowed units are g/mol or u\n";
			input += "# 	Lennard-Jones sigma parameter; allowed units are m or A\n";
			input += "# 	Lennard-Jones epsilon parameter; allowed units are J or K\n";
			input += "u " + Double.toString(bathGas.getMolecularWeight()) + "\n";
			input += "m " + Double.toString(bathGas.getLJSigma()) + "\n";
			input += "J " + Double.toString(bathGas.getLJEpsilon()) + "\n";
			input += "\n";

			input += "# The number of species in the network (minimum of 2)\n";
			input += Integer.toString(speciesList.size()) + "\n";
			input += "\n";

			for (int i = 0; i < speciesList.size(); i++) {

				Species spec = speciesList.get(i);
				spec.calculateTransportParameters();

				input += "# Species identifier (128 characters or less, no spaces)\n";
				input += spec.getName() + "(" + Integer.toString(spec.getID()) + ")\n";

				input += "# Ground-state energy; allowed units are J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n";
				input += "J/mol " + Double.toString(spec.calculateH(stdTemp) * 4184) + "\n";

				input += "# Thermodynamics data:\n";
				input += "# 	Standard enthalpy of formation; allowed units are J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n";
				input += "# 	Standard entropy of formation; allowed units are permutations of energy (J, kJ, cal, or kcal) and temperature (K, C, F, or R)\n";
				input += "# 	Heat capacity at 300, 400, 500, 600, 800, 1000, and 1500 K\n";
				input += "J/mol " + Double.toString(spec.calculateH(stdTemp) * 4184) + "\n";
				input += "J/mol*K " + Double.toString(spec.calculateS(stdTemp) * 4.184) + "\n";
				input += "7 J/mol*K\n";
				input += Double.toString(spec.calculateCp(new Temperature(300, "K")) * 4.184) + "\n";
				input += Double.toString(spec.calculateCp(new Temperature(400, "K")) * 4.184) + "\n";
				input += Double.toString(spec.calculateCp(new Temperature(500, "K")) * 4.184) + "\n";
				input += Double.toString(spec.calculateCp(new Temperature(600, "K")) * 4.184) + "\n";
				input += Double.toString(spec.calculateCp(new Temperature(800, "K")) * 4.184) + "\n";
				input += Double.toString(spec.calculateCp(new Temperature(1000, "K")) * 4.184) + "\n";
				input += Double.toString(spec.calculateCp(new Temperature(1500, "K")) * 4.184) + "\n";
				
				input += "# Species gas parameters\n";
				input += "# 	Molecular weight; allowed units are g/mol or u\n";
				input += "# 	Lennard-Jones sigma parameter; allowed units are m or A\n";
				input += "# 	Lennard-Jones epsilon parameter; allowed units are J or K\n";
				input += "u " + Double.toString(spec.getMolecularWeight()) + "\n";
				input += "m " + Double.toString(spec.getChemkinTransportData().getSigma() * 1e-10) + "\n";
				input += "J " + Double.toString(spec.getChemkinTransportData().getEpsilon() * 1.380665e-23) + "\n";

				input += "# Harmonic oscillators; allowed units are Hz and cm^-1\n";
				SpectroscopicData data = spec.getSpectroscopicData();
				input += Integer.toString(data.getVibrationCount()) + " cm^-1\n";
				for (int j = 0; j < data.getVibrationCount(); j++)
					input += Double.toString(data.getVibration(j)) + "\n";
				
				input += "# Rigid rotors; allowed units are Hz and cm^-1\n";
				input += Integer.toString(data.getRotationCount()) + " cm^-1\n";
				for (int j = 0; j < data.getRotationCount(); j++)
					input += Double.toString(data.getRotation(j));
				
				input += "# Hindered rotor frequencies and barriers\n";
				input += Integer.toString(data.getHinderedCount()) + " cm^-1\n";
				for (int j = 0; j < data.getHinderedCount(); j++)
					input += Double.toString(data.getHinderedFrequency(j)) + "\n";
				input += Integer.toString(data.getHinderedCount()) + " cm^-1\n";
				for (int j = 0; j < data.getHinderedCount(); j++)
					input += Double.toString(data.getHinderedBarrier(j)) + "\n";
				
				input += "# Symmetry number\n";
				input += Integer.toString(spec.getChemGraph().calculateSymmetryNumber()) + "\n";

				input += "\n";
			}
			
			input += "# The number of isomers in the network (minimum of 2)\n";
			input += Integer.toString(nIsom) + "\n";
			input += "# The number of reactant channels in the network\n";
			input += Integer.toString(nReac) + "\n";
			input += "# The number of product channels in the network\n";
			input += Integer.toString(nProd) + "\n";
			input += "\n";

			for (int i = 0; i < isomerList.size(); i++) {
				PDepIsomer isomer = isomerList.get(i);

				input += "# The number and identifiers of each species in the isomer\n";
				input += Integer.toString(isomer.getNumSpecies());
				for (int j = 0; j < isomer.getNumSpecies(); j++) {
					Species spec = isomer.getSpecies(j);
					input += " " + spec.getName() + "(" + Integer.toString(spec.getID()) + ")";
				}
				input += "\n\n";
			}

			input += "# The number of reactions in the network (minimum of 2)\n";
			input += Integer.toString(pathReactionList.size()) + "\n";
			input += "\n";

			for (int i = 0; i < pathReactionList.size(); i++) {
				
				PDepReaction rxn = pathReactionList.get(i);
				
				double A = 0.0, Ea = 0.0, n = 0.0;
				if (rxn.isForward()) {
					Temperature stdtemp = new Temperature(298,"K");
					double Hrxn = rxn.calculateHrxn(stdtemp);
					Kinetics[] k_array = rxn.getKinetics();
					Kinetics kin = computeKUsingLeastSquares(k_array, Hrxn);
					A = kin.getAValue();
					n = kin.getNValue();
					Ea = kin.getEValue();//kin should be ArrheniusKinetics (rather than ArrheniusEPKinetics), so it should be correct to use getEValue here (similarly for other uses in this file)
				}
				else {
					Temperature stdtemp = new Temperature(298,"K");
					double Hrxn = rxn.calculateHrxn(stdtemp);
					((Reaction) rxn).generateReverseReaction();
					Kinetics[] k_array = ((Reaction)rxn).getFittedReverseKinetics();
					Kinetics kin = computeKUsingLeastSquares(k_array, -Hrxn);//gmagoon: I'm not sure, with forward/reverse reactions here whether it is correct to use Hrxn or -Hrxn, but in any case, getFittedReverseKinetics should return an ArrheniusKinetics (not ArrheniusEPKinetics) object, so it will not be used in computeKUsingLeastSquares anyway
//					Kinetics kin = ((Reaction) rxn).getFittedReverseKinetics();
					A = kin.getAValue();
					Ea = kin.getEValue();
					n = kin.getNValue();
				}

				input += "# The reaction equation, in the form A + B --> C + D\n";
				input += rxn.toString() + "\n";

				input += "# Indices of the reactant and product isomers, starting with 1\n";
				input += Integer.toString(isomerList.indexOf(rxn.getReactant()) + 1) + " ";
				input += Integer.toString(isomerList.indexOf(rxn.getProduct()) + 1) + "\n";

				input += "# Ground-state energy; allowed units are J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n";
				if (Ea < 0.0)
					input += "J/mol " + Double.toString((rxn.getReactant().calculateH(stdTemp)) * 4184) + "\n";
				else
					input += "J/mol " + Double.toString((Ea + rxn.getReactant().calculateH(stdTemp)) * 4184) + "\n";

				input += "# High-pressure-limit kinetics model k(T):\n";
				input += "#	Option 1: Arrhenius\n";
				input += "# 	Arrhenius preexponential factor; allowed units are combinations of volume {m^3, L, or cm^3} and time {s^-1}\n";
				input += "# 	Arrhenius activation energy; allowed units are J/mol, kJ/mol, cal/mol, or kcal/mol\n";
				input += "# 	Arrhenius temperature exponent\n";
				input += "Arrhenius\n";
				if (rxn.getReactant().isUnimolecular())
					input += "s^-1 " + Double.toString(A) + "\n";
				else if (rxn.getReactant().isMultimolecular())
					input += "cm^3/mol*s " + Double.toString(A) + "\n";
				input += "J/mol " + Double.toString(Ea * 4184) + "\n";
				input += Double.toString(n) + "\n";

				input += "\n";
			}

		}
		catch(IndexOutOfBoundsException e) {
			System.out.println("Error: IndexOutOfBoundsException thrown.");
			e.printStackTrace(System.out);
		}
		catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace(System.out);
		}

        return input;
	}

	/**
	 * Parses one meaningful line from a FAME output buffer.
	 */
	public String readMeaningfulLine(BufferedReader br) throws IOException {

		String str = "";
		boolean found = false;
		
		while (!found) {
			str = br.readLine().trim();
			found = !(str.length() == 0 || str.startsWith("#"));
		}

		return str;
	}

	/**
	 * Parses a FAME output file and updates the reaction network and system
	 * accordingly.
	 * @param pdn The pressure-dependent reaction network of interest
	 * @param rxnSystem The reaction system of interest
	 * @param cerm The current core/edge reaction model
	 * @param isomerList The set of isomers in the network
	 */
	public boolean parseOutputString(String output, PDepNetwork pdn,
            ReactionSystem rxnSystem, CoreEdgeReactionModel cerm,
			LinkedList<PDepIsomer> isomerList) throws PDepException {

	    String dir = System.getProperty("RMG.workingDirectory");
	    
		double Tmin = 0, Tmax = 0, Pmin = 0, Pmax = 0;

		// JWA January 28 2010
		// Clear the included reactions from the last successful fame execution
		LinkedList<PDepReaction> netReactionList = pdn.getNetReactions();
		netReactionList.clear();
		// Also need to clear the nonincluded reactions - we will regenerate them
		// based on the new output from fame
		pdn.getNonincludedReactions().clear();
		// Below we are going to use netReactionList to temporarily store all
		// reactions, both included and nonincluded, and will then sort them
		// into included and nonincluded

		boolean ignoredARate = false;
		/*
		 * MRH 26Feb2010:
		 * Checking whether pdep rates exceed the high-P-limit
		 * 
		 * This boolean reflects whether all pdep rates are below the high-P-limit.
		 * 	The variable is initialized to true for each call (in hope of avoiding
		 * 	any infinite loop scenarios).  If any of the pdep rates exceed the
		 * 	high-P-limit, the boolean is changed to false (see code below)
		 */
		pdepRatesOK = true;

		try {

			BufferedReader br = new BufferedReader(new StringReader(output));

			String str = "";
			StringTokenizer tkn;

			String method = readMeaningfulLine(br);
			
			str = readMeaningfulLine(br);
			tkn = new StringTokenizer(str);
			int numTemp = Integer.parseInt(tkn.nextToken());
			tkn.nextToken();
			Tmin = Double.parseDouble(tkn.nextToken());
			for (int i = 0; i < numTemp-2; i++) tkn.nextToken();
			Tmax = Double.parseDouble(tkn.nextToken());


			str = readMeaningfulLine(br);
			tkn = new StringTokenizer(str);
			int numPress = Integer.parseInt(tkn.nextToken());
			tkn.nextToken();
			Pmin = Double.parseDouble(tkn.nextToken());
			for (int i = 0; i < numPress-2; i++) tkn.nextToken();
			// Added by MRH on 10Feb2010 to handle the case where user only specifies one pressure
			if (numPress == 1) Pmax = Pmin;
			else Pmax = Double.parseDouble(tkn.nextToken());

			str = readMeaningfulLine(br);
			tkn = new StringTokenizer(str);
			String model = tkn.nextToken();
			int numChebT = 0, numChebP = 0;
			if (model.toLowerCase().contains("chebyshev")) {
				numChebT = Integer.parseInt(tkn.nextToken());
				numChebP = Integer.parseInt(tkn.nextToken());
			}

			str = readMeaningfulLine(br);
			int numSpecies = Integer.parseInt(str);

			str = readMeaningfulLine(br);
			int numIsomers = Integer.parseInt(str);
			str = readMeaningfulLine(br);
			int numReactants = Integer.parseInt(str);
			str = readMeaningfulLine(br);
			int numProducts = Integer.parseInt(str);

			str = readMeaningfulLine(br);
			int numReactions = Integer.parseInt(str);

			str = readMeaningfulLine(br);
			int numKinetics = Integer.parseInt(str);

			for (int i = 0; i < numKinetics; i++) {

				double[][] rates = new double[numTemp][numPress];
				double[][] chebyshev = null;
				PDepArrheniusKinetics pDepArrhenius = null;

				// Reactant and product isomers
				str = readMeaningfulLine(br);
				tkn = new StringTokenizer(str);
				int reac = Integer.parseInt(tkn.nextToken()) - 1;
				int prod = Integer.parseInt(tkn.nextToken()) - 1;
				
				// Table of phenomenological rate coefficients
				boolean valid = true;
				str = readMeaningfulLine(br);	// First row is list of pressures
				for (int t = 0; t < numTemp; t++) {
					str = readMeaningfulLine(br);
					tkn = new StringTokenizer(str);
					tkn.nextToken();
					for (int p = 0; p < numPress; p++) {
						rates[t][p] = Double.parseDouble(tkn.nextToken());
						if (rates[t][p] < 0 ||
								Double.isNaN(rates[t][p]) ||
								Double.isInfinite(rates[t][p]))
							valid = false;
					}
				}

				PDepRateConstant pDepRate = new PDepRateConstant(rates);

				// Chebyshev interpolation model
				if (model.toLowerCase().contains("chebyshev")) {
					chebyshev = new double[numChebT][numChebP];
					for (int t = 0; t < numChebT; t++) {
						str = readMeaningfulLine(br);
						tkn = new StringTokenizer(str);
						for (int p = 0; p < numChebP; p++) {
							chebyshev[t][p] = Double.parseDouble(tkn.nextToken());
						}
					}
					pDepRate.setChebyshev(chebyshev);
				}

				// PDepArrhenius interpolation model
				else if (model.toLowerCase().contains("pdeparrhenius")) {
					pDepArrhenius = new PDepArrheniusKinetics(numPress);
					for (int p = 0; p < numPress; p++) {
						str = readMeaningfulLine(br);
						tkn = new StringTokenizer(str);
						tkn.nextToken();
						UncertainDouble A = new UncertainDouble(Double.parseDouble(tkn.nextToken()), 0.0, "A");
						UncertainDouble Ea = new UncertainDouble(Double.parseDouble(tkn.nextToken()) / 4184.0, 0.0, "A");
						UncertainDouble n = new UncertainDouble(Double.parseDouble(tkn.nextToken()), 0.0, "A");
						String Trange = Double.toString(Tmin) + "-" + Double.toString(Tmax) + " K";
						String Prange = Double.toString(Pmin / 1.0e5) + "-" + Double.toString(Pmax / 1.0e5) + " bar";
						ArrheniusKinetics kinetics = new ArrheniusKinetics(A, n, Ea, Trange, 0, "Result of FAME calculation", "Prange = " + Prange);
						pDepArrhenius.setKinetics(p, pressures[p], kinetics);
					}
					pDepRate.setPDepArrheniusKinetics(pDepArrhenius);
				}

				// If the fitted rate coefficient is not valid, then don't add the net reaction
				if (!valid) {
					ignoredARate = true;
					continue;
				}

				// Initialize net reaction
				PDepIsomer reactant = isomerList.get(reac);
				PDepIsomer product = isomerList.get(prod);
				PDepReaction rxn = new PDepReaction(reactant, product, pDepRate);
				
				/*
				 * MRH 26Feb2010:
				 * Checking whether pdep rates exceed the high-P-limit
				 * 
				 * We grab all of the reactions in the pathReactionList.  These are the
				 * 	RMG-generated reactions (not chemically-activated reactions) and thus
				 * 	have "natural" high-P-limit kinetics.
				 * If the current PDepReaction "rxn"'s structure matches one of the structures
				 * 	of the PDepReactions located in the pathReactionList, we compare the high-P-limit
				 * 	kinetics of the pathReactionList (these values either come from the RMG database
				 * 	or are "fitted" parameters, based on the reverse kinetics + equilibrium constant)
				 * 	with the k(T,P_max) for each T in the "temperatures" array.  NOTE: Not every "rxn"
				 * 	will have a match in the pathReactionList.
				 * If the pdep rate is greater than 2x the high-P-limit, we consider this to be different.
				 * 	If the number of grains is less than 1000, we set the pdepRatesOK boolean to false,
				 * 		so that another fame calculation will ensue
				 * 	If the number of grains exceeds 1000, we continue on with the simulation, but alert the
				 * 		user of the discrepancy.
				 * 
				 * The value of 2 was somewhat randomly chosen by MRH.
				 * For a toy case of tBuOH pyrolysis (with 1e-6 reaction time and 0.9 error tolerance):
				 * 	Before code addition:	iC4H8/H2O final concentration = 7.830282E-10, run time < 1 min
				 * 	2x : 					iC4H8/H2O final concentration = 6.838976E-10, run time ~ 2 min
				 *  1.5x :					iC4H8/H2O final concentration = 6.555548E-10, run time ~ 4 min
				 *  1.1x :                  iC4H8/H2O final concentration = 6.555548E-10, run time ~ 10 min
				 *  
				 * The value of 2 (for the toy model) seems to be a good balance between speed and accuracy
				 * 
				 * P.S. Want to keep the name fame (fast approximate me) instead of having to change to
				 * 	smame (slow more accurate me).  ;)
				 *  
				 */
				LinkedList pathReactionList = pdn.getPathReactions();
				boolean foundHighPLimitRxn = false;
				for (int HighPRxNum = 0; HighPRxNum < pathReactionList.size(); HighPRxNum++) {
					PDepReaction rxnWHighPLimit = (PDepReaction)pathReactionList.get(HighPRxNum);
					Temperature stdtemp = new Temperature(298,"K");
					double Hrxn = rxnWHighPLimit.calculateHrxn(stdtemp);
					if (rxn.getStructure().equals(rxnWHighPLimit.getStructure())) {
						foundHighPLimitRxn = true;
						double A = 0.0, Ea = 0.0, n = 0.0;
						if (rxnWHighPLimit.isForward()) {
							Kinetics[] k_array = rxnWHighPLimit.getKinetics();
							Kinetics kin = computeKUsingLeastSquares(k_array, Hrxn);
							A = kin.getAValue();
							Ea = kin.getEValue();
							n = kin.getNValue();
							// While I'm here, and know which reaction was the High-P limit, set the comment in the P-dep reaction 
							rxn.setComments("NetReaction from PDepNetwork #" + Integer.toString(pdn.getID()) + " (" + pdn.getSpeciesType() + ")" + 
									" High-P Limit: " + kin.getSource().toString() + " " + kin.getComment().toString() );
						}
						else {
							Kinetics[] k_array = rxnWHighPLimit.getFittedReverseKinetics();
							Kinetics kin = computeKUsingLeastSquares(k_array, -Hrxn);//gmagoon: I'm not sure, with forward/reverse reactions here whether it is correct to use Hrxn or -Hrxn, but in any case, getFittedReverseKinetics should return an ArrheniusKinetics (not ArrheniusEPKinetics) object, so it will not be used in computeKUsingLeastSquares anyway
							A = kin.getAValue();
							Ea = kin.getEValue();
							n = kin.getNValue();
							//  While I'm here, and know which reaction was the High-P limit, set the comment in the P-dep reaction 
							Kinetics[] fwd_kin = rxnWHighPLimit.getKinetics();
							String commentsForForwardKinetics = "";
							if (fwd_kin.length > 1) commentsForForwardKinetics += "Summation of kinetics:\n!";
							for (int numKs=0; numKs<fwd_kin.length; ++numKs) {
								commentsForForwardKinetics += "High-P Limit Reverse: " + fwd_kin[numKs].getSource().toString() +fwd_kin[numKs].getComment().toString();
								if (numKs != fwd_kin.length-1) commentsForForwardKinetics += "\n!";
							}
							rxn.setComments("NetReaction from PDepNetwork #" + Integer.toString(pdn.getID()) + " (" + pdn.getSpeciesType() + ")" +
									" " + commentsForForwardKinetics);
						}
						if (ReactionModelGenerator.rerunFameWithAdditionalGrains()) {
							double[][] all_ks = rxn.getPDepRate().getRateConstants();
							for (int numTemps=0; numTemps<temperatures.length; numTemps++) {
								double T = temperatures[numTemps].getK();
								double k_highPlimit = A * Math.pow(T,n) * Math.exp(-Ea/GasConstant.getKcalMolK()/T);
								if (all_ks[numTemps][pressures.length-1] > 2*k_highPlimit) {
									if (numGrains > 1000) {
										System.out.println("Pressure-dependent rate exceeds high-P-limit rate: " +
											"Number of grains already exceeds 1000.  Continuing RMG simulation " +
											"with results from current fame run.");
										break;	// No need to continue checking the rates for this one reaction
									}
									else {
										System.out.println("Pressure-dependent rate exceeds high-P-limit rate: " +
											"Re-running fame with additional number of grains");
										pdepRatesOK = false;
										return false;
									}
								}
							}
						//break;	// Why did MRH put this break here?
						}
					}
				}
				// If not found, we have a "nonIncluded" (pressure-dependent) reaction
				if (!foundHighPLimitRxn) {
					rxn.setComments("NetReaction from PDepNetwork #" + Integer.toString(pdn.getID()) + " (" + pdn.getSpeciesType() + ")");
				}
				
				// Add net reaction to list
				netReactionList.add(rxn);


			}

			// Close file when finished
			br.close();

			// Set reverse reactions
			int i = 0;
			while (i < netReactionList.size()) {
				PDepReaction rxn1 = netReactionList.get(i);
				if (rxn1.getReverseReaction() == null) {
					boolean found = false;
					for (int j = 0; j < netReactionList.size(); j++) {
						PDepReaction rxn2 = netReactionList.get(j);
						if (rxn1.getReactant().equals(rxn2.getProduct()) &&
							rxn2.getReactant().equals(rxn1.getProduct())) {
							rxn1.setReverseReaction(rxn2);
							rxn2.setReverseReaction(rxn1);
							netReactionList.remove(rxn2);
							found = true;
							break;
						}
					}
					if (!found) {
						rxn1.setReverseReaction(null);
						i++;
					}
				}
				else
					i++;
			}

			// Update reaction lists (sort into included and nonincluded)
			pdn.updateReactionLists(cerm);

		}
		catch (PDepException e) {
			System.out.println(pdn);
			System.out.println(e.getMessage());
			e.printStackTrace();
			throw new PDepException("Unable to parse FAME output file.");
		}
		catch(IOException e) {
			System.out.println("Error: Unable to read from file \"fame_output.txt\".");
			System.exit(0);
		}
		
		return true;
    }

	/**
	 * Determines the maximum energy grain to use in the calculation. The
	 * maximum energy grain is chosen to be 25 * R * T above the highest
	 * potential energy on the potential energy surface, which hopefully will
	 * capture the majority of the equilibrium distributions even at the
	 * high temperature end of the calculation.
	 *
	 * @param uniIsomers The set of unimolecular isomers in the network
	 * @param multiIsomers The set of multimolecular isomers in the network
	 * @param T The temperature of the calculation in K
	 * @return The maximum grain energy in kJ/mol
	 */
	private double getGrainMaxEnergy(LinkedList<PDepIsomer> isomerList, double T) {

		double Emax = 0.0;
		Temperature stdTemp = new Temperature(298, "K");

		for (int i = 0; i < isomerList.size(); i++) {
			PDepIsomer isomer = isomerList.get(i);
			double E = isomer.calculateH(stdTemp);
			if (E > Emax)
				Emax = E;
		}

		Emax *= 4.184;
		//Emax += 100 * 0.008314 * T;
		Emax += 50 * 0.008314 * T;

		// Round up to nearest ten kJ/mol
		return Math.ceil(Emax / 10) * 10.0;	// [=] kJ/mol
	}

	/**
	 * Determines the minimum energy grain to use in the calculation. The
	 * maximum energy grain is chosen to be at or just below the energy of the
	 * lowest energy isomer in the system.
	 *
	 * @param uniIsomers The set of unimolecular isomers in the network
	 * @param multiIsomers The set of multimolecular isomers in the network
	 * @return The minimum grain energy in kJ/mol
	 */
	private double getGrainMinEnergy(LinkedList<PDepIsomer> isomerList) {

		double Emin = 1000000.0;
		Temperature stdTemp = new Temperature(298, "K");

		for (int i = 0; i < isomerList.size(); i++) {
			PDepIsomer isomer = isomerList.get(i);
			double E = isomer.calculateH(stdTemp);
			if (E < Emin)
				Emin = E;
		}

		Emin *= 4.184;

		// Round down to nearest ten kJ/mol
		return Math.floor(Emin / 10) * 10.0;	// [=] kJ/mol
	}

	/**
	 * Determines a suitable grain size for the calculation. Too small a grain
	 * size will cause the simulation to take a very long time to conduct; too
	 * large a grain size will cause the simulation results to not be very
	 * precise. The choice is hopefully made such that the grains get larger
	 * as the system does, so as to not get completely bogged down in the
	 * larger networks.
	 *
	 * @param grainMinEnergy The minimum energy grain in kJ/mol
	 * @param grainMaxEnergy The maximum energy grain in kJ/mol
	 * @param numUniWells The number of unimolecular isomers in the network.
	 * @return
	 */
	private double getGrainSize(double grainMinEnergy, double grainMaxEnergy, int numUniWells) {
		if (numUniWells < 5)
			return (grainMaxEnergy - grainMinEnergy) / 200;
		else if (numUniWells < 10)
			return (grainMaxEnergy - grainMinEnergy) / 100;
		else if (numUniWells < 20)
			return (grainMaxEnergy - grainMinEnergy) / 50;
		else
			return (grainMaxEnergy - grainMinEnergy) / 20;
	}

	/**
	 * Creates an empty file on the hard disk for the FORTRAN executable to
	 * use to write its output. If this file is not created, the FORTRAN
	 * executable will write to the file 'fort.2', which is checked by the
	 * readOutputFile() function but still discouraged.
	 */
	public void touchOutputFile() {
		try {
			File output = new File("fame/fame_output.txt");
			if (output.exists())
				output.delete();
			output.createNewFile();
		}
		catch(IOException e) {
			System.out.println("Error: Unable to touch file \"fame_output.txt\".");
			System.out.println(e.getMessage());
			System.exit(0);
		}
		catch(Exception e) {
			System.out.println(e.getMessage());
		}
	}

	/**
	 * Returns the mode as a string suitable for writing the FAME input file.
	 * @return The mode as a string suitable for writing the FAME input file.
	 */
	public String getModeString() {
		if (mode == Mode.STRONGCOLLISION)
			return "ModifiedStrongCollision";
		else if (mode == Mode.RESERVOIRSTATE)
			return "ReservoirState";
		else
			return "";

	}

	public static Temperature[] getTemperatures() {
		return temperatures;
	}

	public static void setTemperatures(Temperature[] temps) {
		temperatures = temps;
	}

	public static Pressure[] getPressures() {
		return pressures;
	}

	public static void setPressures(Pressure[] press) {
		pressures = press;
	}
	
	public static void setNumTBasisFuncs(int n) {
		numTBasisFuncs = n;
	}
	
	public static void setNumPBasisFuncs(int m) {
		numPBasisFuncs = m;
	}
	
	public static Kinetics computeKUsingLeastSquares(Kinetics[] k_array, double Hrxn) {
        /*
         * MRH 24MAR2010:
         * 	If the reaction has more than one set of Arrhenius kinetics,
         * 	sum all kinetics and re-fit for a single set of Arrhenius kinetics
         */
		Kinetics k = null;
        if (k_array.length > 1) {
        	double[] T = new double[5];
        	T[0] = 300; T[1] = 600; T[2] = 900; T[3] = 1200; T[4] = 1500;
        	double[] k_total = new double[5];
        	for (int numKinetics=0; numKinetics<k_array.length; ++numKinetics) {
        		for (int numTemperatures=0; numTemperatures<T.length; ++numTemperatures) {
				double Ea = 0.0;
				if (k_array[numKinetics] instanceof ArrheniusEPKinetics){
				    Ea = ((ArrheniusEPKinetics)k_array[numKinetics]).getEaValue(Hrxn);
				}
				else{
				    Ea = k_array[numKinetics].getEValue();
				}
        			k_total[numTemperatures] += k_array[numKinetics].getAValue() * 
        			Math.pow(T[numTemperatures], k_array[numKinetics].getNValue()) * 
        			Math.exp(-Ea/GasConstant.getKcalMolK()/T[numTemperatures]);
        		}
        	}
        	// Construct matrix X and vector y
        	double[][] y = new double[k_total.length][1];
        	double[][] X = new double[k_total.length][3];
        	for (int i=0; i<5; i++) {
        		y[i][0] = Math.log(k_total[i]);
        		X[i][0] = 1;
        		X[i][1] = Math.log(T[i]);
        		X[i][2] = 1/T[i];
        	}
        	// Solve least-squares problem using inv(XT*X)*(XT*y)
        	Matrix X_matrix = new Matrix(X);
        	Matrix y_matrix = new Matrix(y);
        	Matrix b_matrix = X_matrix.solve(y_matrix);
        	UncertainDouble uA = new UncertainDouble(Math.exp(b_matrix.get(0,0)),0.0,"Adding");
        	UncertainDouble un = new UncertainDouble(b_matrix.get(1,0),0.0,"Adding");
        	UncertainDouble uE = new UncertainDouble(-GasConstant.getKcalMolK()*b_matrix.get(2,0),0.0,"Adding");
			String commentsForFittedKinetics = "";
			for (int numKs=0; numKs<k_array.length; ++numKs) {
				commentsForFittedKinetics += "!" + k_array[numKs].getSource().toString();
				if (k_array[numKs].getComment() != null) commentsForFittedKinetics += k_array[numKs].getComment().toString();
				if (numKs != k_array.length-1) commentsForFittedKinetics += "\n";
			}
        	k = new ArrheniusKinetics(uA,un,uE,"300-1500K",5,"Summation of kinetics:",commentsForFittedKinetics);
        } else {
        	k = k_array[0];
        }
        return k;
	}
	
	public static int getNumTBasisFuncs() {
		return numTBasisFuncs;
	}
	
	public static int getNumPBasisFuncs() {
		return numPBasisFuncs;
	}

}

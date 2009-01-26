////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

package jing.rxn;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.StringTokenizer;
import jing.chem.Species;
import jing.chem.SpectroscopicData;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxnSys.CoreEdgeReactionModel;
import jing.rxnSys.ReactionSystem;

/**
 *
 * @author jwallen
 */
public class FastMasterEqn implements PDepKineticsEstimator {

	/**
	 * The number of times the FAME module has been called for any network
	 * since the inception of this RMG execution.
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
		
		// No update needed if network is only two wells and one reaction (?)
		/*if (pdn.getUniIsomers().size() + pdn.getMultiIsomers().size() == 2 &&
				pdn.getPathReactions().size() == 1) {
			System.out.println("No pDep calculation needed for network with only two isomers and one reaction.");
			LinkedList<PDepReaction> paths = pdn.getPathReactions();
			LinkedList<PDepReaction> net = pdn.getNetReactions();
			net.clear();
			net.add(paths.get(0));
			return;
		}*/
		
		// Get working directory (to find FAME executable)
		String dir = System.getProperty("RMG.workingDirectory");
        
		// Determine wells and reactions; skip if no reactions in network
		LinkedList<PDepIsomer> uniIsomers = pdn.getUniIsomers();
		LinkedList<PDepIsomer> multiIsomers = pdn.getMultiIsomers();
		LinkedList<PDepReaction> pathReactions = pdn.getPathReactions();
		if (pathReactions.size() == 0) {
			System.out.println("Warning: Empty pressure-dependent network detected. Skipping.");
			return;
		}
		/*for (int i = 0; i < pathReactions.size(); i++) {
			Reaction reaction = (Reaction) pathReactions.get(i).getTheReaction();
			System.out.println(reaction.toString() + 
				"    A = " + reaction.getKinetics().getAValue() + " s^-1 " +
				"    Ea = " + reaction.getKinetics().getEValue() + " kcal/mol " +
				"    n = " + reaction.getKinetics().getNValue());
			reaction = (Reaction) pathReactions.get(i);
			System.out.println(reaction.toString() + 
				"    A = " + reaction.getKinetics().getAValue() + " s^-1 " +
				"    Ea = " + reaction.getKinetics().getEValue() + " kcal/mol " +
				"    n = " + reaction.getKinetics().getNValue());
			
		}*/
		
		
		// Create FAME input files
		writeInputFile(pdn, rxnSystem, uniIsomers, multiIsomers, pathReactions);
		
		// Touch FAME output file
		touchOutputFile();
		
		// FAME system call
		try {
           	System.out.println("Running FAME...");
			String[] command = {dir + "/software/fame/fame.exe"};
           	File runningDir = new File("fame/");
			Process fame = Runtime.getRuntime().exec(command, null, runningDir);                     
            /*InputStream ips = fame.getInputStream();
            InputStreamReader is = new InputStreamReader(ips);
            BufferedReader br = new BufferedReader(is);
            // Print FAME stdout
			String line = null;
            while ( (line = br.readLine()) != null) {
            	System.out.println(line);
            }*/
            int exitValue = fame.waitFor();
			
        }
        catch (Exception e) {
        	System.out.println("Error while executing FAME!");
        	System.exit(0);
        }
        
		// Parse FAME output file and update accordingly
        readOutputFile(pdn, rxnSystem, cerm, uniIsomers, multiIsomers);
        
		// Reset altered flag
        pdn.setAltered(false);
		
		// Clean up files
		String path = "fame/";
		runCount++;
		if (runCount < 10)			path += "000";
		else if (runCount < 100)	path += "00";
		else if (runCount < 1000)	path += "0";
		path += Integer.toString(runCount);
		
		File input = new File("fame/fame_input.txt");
		File newInput = new File(path +  "_input.txt");
		input.renameTo(newInput);
		File output = new File("fame/fame_output.txt");
		File newOutput = new File(path +  "_output.txt");
		output.renameTo(newOutput);
		
		// Write finished indicator to console
		System.out.println("FAME execution for network " + Integer.toString(runCount) + " complete.");
		
	}
	
	/**
	 * Creates the input file needed by FAME that represents a pressure-
	 * dependent reaction network.
	 * @param pdn The reaction network of interest
	 * @param rxnSystem The reaction system of interest
	 * @param uniIsomers The set of unimolecular isomers in the network
	 * @param multiIsomers The set of multimolecular isomers in the network
	 * @param pathReactions The set of path reactions in the network
	 */
	public void writeInputFile(PDepNetwork pdn,	ReactionSystem rxnSystem,
			LinkedList<PDepIsomer> uniIsomers,
			LinkedList<PDepIsomer> multiIsomers, 
			LinkedList<PDepReaction> pathReactions) {
		
		Temperature stdTemp = new Temperature(298, "K");
                
		// Collect simulation parameters
		Temperature temperature = rxnSystem.getPresentTemperature();
		Pressure pressure = rxnSystem.getPresentPressure();
		BathGas bathGas = new BathGas(rxnSystem);
		
		int numChebTempPolys = 4, numChebPressPolys = 4;
		
		// Determine reference energies
		double Eref = getReferenceEnergy(uniIsomers, multiIsomers, stdTemp);		// [=] kcal/mol
		double grainMaxEnergy = getGrainMaxEnergy(uniIsomers, multiIsomers, 2100, Eref); // [=] kJ/mol
		double grainSize = getGrainSize(grainMaxEnergy, uniIsomers.size()); // [=] kJ/mol
		
		// Create the simulation parameters file fame/simData.txt
		try {
        	
			File simData = new File("fame/fame_input.txt");
			FileWriter fw = new FileWriter(simData);

			fw.write("# FAME input for RMG-generated pressure dependent network #" + runCount + "\n");
			fw.write("Mode                                " + getModeString() + "\n");
			fw.write("Temperatures                        7\n");
			fw.write("300 K\n600 K\n900 K\n1200 K\n1500 K\n1800 K\n2100 K\n");
			fw.write("Pressures                           5\n");
			fw.write("0.01 bar\n0.1 bar\n1 bar\n10 bar\n100 bar\n");
            fw.write("Grain size                          " + grainSize + " kJ/mol\n");
			fw.write("Maximum grain energy                " + grainMaxEnergy + " kJ/mol\n");
			fw.write("Number of unimolecular wells        " + uniIsomers.size() + "\n");
			fw.write("Number of multimolecular wells      " + multiIsomers.size() + "\n");
			fw.write("Number of reactions                 " + pathReactions.size() + "\n");
			fw.write("Exponential down parameter          " + bathGas.getExpDownParam() + " kJ/mol\n");
			fw.write("Bath gas LJ sigma parameter         " + bathGas.getLJSigma() + " m\n");
			fw.write("Bath gas LJ epsilon parameter       " + bathGas.getLJEpsilon() + " J\n");
			fw.write("Bath gas molecular weight           " + bathGas.getMolecularWeight() + " g/mol\n");
			fw.write("Number of Chebyshev temperatures    " + numChebTempPolys + "\n");
			fw.write("Number of Chebyshev pressures       " + numChebPressPolys + "\n");
			fw.write("\n");

			for (int i = 0; i < uniIsomers.size(); i++) {
				
				Species species = uniIsomers.get(i).getSpecies(0);
				species.calculateLJParameters();
				
				fw.write("# Unimolecular well " + Integer.toString(i+1) + ": " + species.getName() + "(" + Integer.toString(species.getID()) + ")" + "\n");
				fw.write("Ground-state energy                 " + (species.calculateH(stdTemp) - Eref) * 4.184 + " kJ/mol\n");
				fw.write("Enthalpy of formation               " + (species.calculateH(stdTemp) - Eref) * 4.184 + " kJ/mol\n");
				fw.write("Free energy of formation            " + (species.calculateG(stdTemp) - Eref) * 4.184 + " kJ/mol\n");
				fw.write("LJ sigma parameter                  " + (species.getLJ().getSigma() * 1e-10) + " m\n");
				fw.write("LJ epsilon parameter                " + (species.getLJ().getEpsilon() * 1.381e-23) + " J\n");
				fw.write("Molecular weight                    " + species.getMolecularWeight() + " g/mol\n");
				
				SpectroscopicData data = species.getSpectroscopicData();
				if (data.getVibrationCount() > 0) {
					fw.write("Harmonic oscillators                " + data.getVibrationCount() + "\n");
					for (int j = 0; j < data.getVibrationCount(); j++)
						fw.write(data.getVibration(j) + " cm^-1\n");
				}
				if (data.getRotationCount() > 0) {
					fw.write("Rigid rotors                        " + data.getRotationCount() + "\n");
					for (int j = 0; j < data.getRotationCount(); j++)
						fw.write(data.getRotation(j) + " cm^-1\n");
				}
				if (data.getHinderedCount() > 0) {
					fw.write("Hindered rotors                     " + data.getHinderedCount() + "\n");
					for (int j = 0; j < data.getHinderedCount(); j++)
						fw.write(data.getHinderedFrequency(j) + " cm^-1\n");
					for (int j = 0; j < data.getHinderedCount(); j++)
						fw.write(data.getHinderedBarrier(j) + " cm^-1\n");
				}
				fw.write("Symmetry number                         " + data.getSymmetryNumber() + "\n");
					
				fw.write("\n");
			}
			
			for (int i = 0; i < multiIsomers.size(); i++) {
				
				PDepIsomer isomer = multiIsomers.get(i);

				fw.write("# Multimolecular well " + Integer.toString(i+1) + ": ");
				fw.write(isomer.getSpecies(0).getName() + "(" + Integer.toString(isomer.getSpecies(0).getID()) + ")");
				for (int j = 1; j < isomer.getNumSpecies(); j++)
					fw.write(" + " + isomer.getSpecies(j).getName() + "(" + Integer.toString(isomer.getSpecies(j).getID()) + ")");
				fw.write("\n");
				fw.write("Number of species                   " + Integer.toString(isomer.getNumSpecies()) + "\n");
				fw.write("Ground-state energy                 ");
				for (int j = 0; j < isomer.getNumSpecies(); j++)
					fw.write((isomer.calculateH(j, stdTemp) - Eref) * 4.184 + " kJ/mol    ");
				fw.write("\n");
				fw.write("Enthalpy of formation               ");
				for (int j = 0; j < isomer.getNumSpecies(); j++)
					fw.write((isomer.calculateH(j, stdTemp) - Eref) * 4.184 + " kJ/mol    ");
				fw.write("\n");
				fw.write("Free energy of formation            ");
				for (int j = 0; j < isomer.getNumSpecies(); j++)
					fw.write((isomer.calculateG(j, stdTemp) - Eref) * 4.184 + " kJ/mol    ");
				fw.write("\n");
				
				fw.write("Harmonic oscillators                ");
				for (int j = 0; j < isomer.getNumSpecies(); j++)
					fw.write(isomer.getSpecies(j).getSpectroscopicData().getVibrationCount() + " ");
				fw.write("\n");
				for (int j = 0; j < isomer.getNumSpecies(); j++) {
					SpectroscopicData data = isomer.getSpecies(j).getSpectroscopicData();
					for (int k = 0; k < data.getVibrationCount(); k++)
						fw.write(data.getVibration(k) + " cm^-1\n");
				}
				
				fw.write("Rigid rotors                        ");
				for (int j = 0; j < isomer.getNumSpecies(); j++)
					fw.write(isomer.getSpecies(j).getSpectroscopicData().getRotationCount() + " ");
				fw.write("\n");
				for (int j = 0; j < isomer.getNumSpecies(); j++) {
					SpectroscopicData data = isomer.getSpecies(j).getSpectroscopicData();
					for (int k = 0; k < data.getRotationCount(); k++)
						fw.write(data.getRotation(k) + " cm^-1\n");
				}
				
				fw.write("Hindered rotors                     ");
				for (int j = 0; j < isomer.getNumSpecies(); j++)
					fw.write(isomer.getSpecies(j).getSpectroscopicData().getHinderedCount() + " ");
				fw.write("\n");
				for (int j = 0; j < isomer.getNumSpecies(); j++) {
					SpectroscopicData data = isomer.getSpecies(j).getSpectroscopicData();
					for (int k = 0; k < data.getHinderedCount(); k++)
						fw.write(data.getHinderedFrequency(k) + " cm^-1\n");
					for (int k = 0; k < data.getHinderedCount(); k++)
						fw.write(data.getHinderedBarrier(k) + " cm^-1\n");
				}
				fw.write("Symmetry number                     ");
				for (int j = 0; j < isomer.getNumSpecies(); j++)
					fw.write(isomer.getSpecies(j).getSpectroscopicData().getSymmetryNumber() + " ");
				fw.write("\n");
				
				fw.write("\n");
			}
			
			for (int i = 0; i < pathReactions.size(); i++) {
				
				PDepReaction rxn = pathReactions.get(i);
				
				// Determine isomers associated with reactant and product
				PDepIsomer reactant = rxn.getReactant();
				PDepIsomer product = rxn.getProduct();
				int isomer1 = 0, isomer2 = 0;
				if (reactant == null || product == null) {
					continue;
				}
				if (reactant.isUnimolecular())
					isomer1 = uniIsomers.indexOf(reactant) + 1;
				else
					isomer1 = multiIsomers.indexOf(reactant) + uniIsomers.size() + 1;
				if (product.isUnimolecular())
					isomer2 = uniIsomers.indexOf(product) + 1;
				else
					isomer2 = multiIsomers.indexOf(product) + uniIsomers.size() + 1;
					
				// Arrhenius parameters
				double A = rxn.getHighPKinetics().getAValue();
				double Ea = rxn.getHighPKinetics().getEValue();
				double n = rxn.getHighPKinetics().getNValue();
				if (Ea < 0) {
					System.out.println("Warning: Adjusted activation energy of reaction " +
							rxn.toString() + " from " + Double.toString(Ea) +
							" kcal/mol to 0 kcal/mol for FAME calculation.");
					Ea = 0;
				}
				
				// Calculate transition state energy = ground-state energy of reactant (isomer 1) + activation energy
				double E0 = Ea + reactant.calculateH(stdTemp);
						
				fw.write("# Reaction " + rxn.toString() + ":\n");
				fw.write("Isomer 1                            " + isomer1 + "\n");
				fw.write("Isomer 2                            " + isomer2 + "\n");
				fw.write("Ground-state energy                 " + (E0 - Eref) * 4.184 + " kJ/mol\n");
				fw.write("Arrhenius preexponential            " + A + " s^-1\n");
				fw.write("Arrhenius activation energy         " + Ea * 4.184 + " kJ/mol\n");
				fw.write("Arrhenius temperature exponent      " + n + "\n");
				fw.write("\n");
			}
			
			fw.write("\n");
            fw.close();
		}
		catch(IOException e) {
			System.out.println("Error: Unable to create file \"fame_input.txt\".");
			System.exit(0);
		}
		catch(IndexOutOfBoundsException e) {
			System.out.println("Error: IndexOutOfBoundsException thrown.");
			e.printStackTrace(System.out);
		}
		catch(Exception e) {
			System.out.println(e.getMessage());
		}
		
	}
	
	/**
	 * Parses a FAME output file and updates the reaction network and system
	 * accordingly.
	 * @param pdn The pressure-dependent reaction network of interest
	 * @param rxnSystem The reaction system of interest
	 * @param cerm The current core/edge reaction model
	 * @param uniIsomers The set of unimolecular isomers in the network
	 * @param multiIsomers The set of multimolecular isomers in the network
	 */
	public void readOutputFile(PDepNetwork pdn, ReactionSystem rxnSystem,
			CoreEdgeReactionModel cerm,
			LinkedList<PDepIsomer> uniIsomers, 
			LinkedList<PDepIsomer> multiIsomers) {
        
	    String dir = System.getProperty("RMG.workingDirectory");
	    int numUniWells = 0;
		int numMultiWells = 0;
		int numReactions = 0;
		int numTemperatures = 0;
		int numPressures = 0;
		int entryWell = 0;
		
		double Tmin = 0, Tmax = 0, Pmin = 0, Pmax = 0;
		
		try {
        	
			File output = new File("fame/fame_output.txt");
			if (!output.exists()) {
				output = new File("fame/fort.2");
				if (!output.exists()) 
					throw new Exception("FAME output file not found!");	
			}
			
			if (output.length() == 0)
				throw new IOException();
			
			BufferedReader br = new BufferedReader(new FileReader(output));
			
			String str = "";
			
			// Read output file header
			boolean found = false;
			while (!found) {
				
				str = br.readLine().trim();
				
				if (str.length() == 0)
					found = true;
				else if (str.charAt(0) == '#')
					continue;
				else if (str.substring(0, 28).equals("Number of unimolecular wells"))
					numUniWells = Integer.parseInt(str.substring(29).trim());
				else if (str.substring(0, 30).equals("Number of multimolecular wells"))
					numMultiWells = Integer.parseInt(str.substring(31).trim());
				else if (str.substring(0, 32).equals("Number of Chebyshev temperatures"))
					numTemperatures = Integer.parseInt(str.substring(33).trim());
				else if (str.substring(0, 29).equals("Number of Chebyshev pressures"))
					numPressures = Integer.parseInt(str.substring(30).trim());
				else if (str.substring(0, 24).equals("Temperature range of fit")) {
					StringTokenizer strtok = new StringTokenizer(str.substring(25).trim());
					Tmin = Double.parseDouble(strtok.nextToken());
					strtok.nextToken();
					Tmax = Double.parseDouble(strtok.nextToken());
				}
				else if (str.substring(0, 21).equals("Pressure range of fit")) {
					StringTokenizer strtok = new StringTokenizer(str.substring(22).trim());
					Pmin = Double.parseDouble(strtok.nextToken());
					strtok.nextToken();
					Pmax = Double.parseDouble(strtok.nextToken());
				}
			}
			
			// Initialize temperature and pressure variables based on read values
			Temperature tLow = new Temperature(Tmin, "K");
			Temperature tHigh = new Temperature(Tmax, "K");
			Pressure pLow = new Pressure(Pmin, "bar");
			Pressure pHigh = new Pressure(Pmax, "bar");
		
			LinkedList<PDepReaction> netReactionList = pdn.getNetReactions();
			netReactionList.clear();
        	
			// Read Chebyshev coefficients for each reaction
			boolean ignoredARate = false;
			for (int i = 0; i < numUniWells + numMultiWells; i++) {
				for (int j = 0; j < numUniWells + numMultiWells; j++) {
					if (i != j && br.ready()) {
					
						double[][] alpha = new double[numTemperatures][numPressures];

						// Comment line at start of reaction
						str = br.readLine().trim();
						
						// Read Chebyshev coefficients from file
						boolean valid = true;
						for (int t = 0; t < numTemperatures; t++) {
							str = br.readLine().trim();
							StringTokenizer tkn = new StringTokenizer(str);
							for (int p = 0; p < numPressures; p++) {
								alpha[t][p] = Double.parseDouble(tkn.nextToken());
								if (alpha[t][p] == 0 ||
										Double.isNaN(alpha[t][p]) ||
										Double.isInfinite(alpha[t][p])) 
									valid = false;
							}
						}
						
						// Skip blank line between records
						str = br.readLine().trim();
							
						// If the fitted rate coefficient is not valid, then don't add the net reaction
						if (!valid) {
							ignoredARate = true;
							continue;
						}
						
						// Create Chebyshev polynomial object for current rate coefficient
						ChebyshevPolynomials cp = new ChebyshevPolynomials(
								numTemperatures, tLow, tHigh, numPressures, pLow, pHigh,
								alpha);

						// Determine reactants
						PDepIsomer reactant = null;
						if (i < numUniWells) 
							reactant = uniIsomers.get(i);
						else
							reactant = multiIsomers.get(i - numUniWells);
						
						// Determine products
						PDepIsomer product = null;
						if (j < numUniWells) 
							product = uniIsomers.get(j);
						else
							product = multiIsomers.get(j - numUniWells);
						
						// Initialize net reaction
						PDepReaction rxn = new PDepReaction(reactant, product, cp);

						// Add net reaction to list
						netReactionList.add(rxn);

					}
				}
			}
			
			// Close file when finished
			br.close();
			
			// Set reverse reactions
			for (int i = 0; i < netReactionList.size(); i++) {
				PDepReaction rxn1 = netReactionList.get(i);
				for (int j = 0; j < netReactionList.size(); j++) {
					PDepReaction rxn2 = netReactionList.get(j);
					if (rxn1.getReactant().equals(rxn2.getProduct()) &&
						rxn2.getReactant().equals(rxn1.getProduct())) {
						rxn1.setReverseReaction(rxn2);
						rxn2.setReverseReaction(rxn1);
						netReactionList.remove(rxn2);
					}
				}
			}
			
			if (ignoredARate)
				System.out.println("Warning: One or more rate coefficients in FAME output was ignored due to invalid values.");

			// Update reaction lists (sort into included and nonincluded)
			pdn.updateReactionLists(cerm);
			
		}
		catch(IOException e) {
			System.out.println("Error: Unable to read from file \"fame_output.txt\".");
			System.exit(0);
		}
		catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
        
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
	 * @param Eref The reference energy (the energy used to shift all 
	 * energies of the system such that the lowest one is zero) in kJ/mol
	 * @return The maximum grain energy in kJ/mol
	 */
	private double getGrainMaxEnergy(LinkedList<PDepIsomer> uniIsomers, 
			LinkedList<PDepIsomer> multiIsomers, double T, double Eref) {

		double Emax = 0.0;
		Temperature stdTemp = new Temperature(298, "K");
        
		for (int i = 0; i < uniIsomers.size(); i++) {
			PDepIsomer isomer = uniIsomers.get(i);
			double E = isomer.calculateH(stdTemp);
			if (E > Emax)
				Emax = E;
		}
		for (int n = 0; n < multiIsomers.size(); n++) {
			PDepIsomer isomer = multiIsomers.get(n);
			double E = isomer.calculateH(stdTemp);
			if (E > Emax)
				Emax = E;
		}
		
		Emax -= Eref;
		Emax *= 4.184;
		Emax += 25 * 0.008314 * T;
		
		// Round up to nearest ten kJ/mol
		return Math.ceil(Emax / 10) * 10.0;	// [=] kJ/mol
	}
	
	/**
	 * Determines the reference energy: the energy used to shift all 
	 * energies of the system such that the lowest one is zero. This reference
	 * energy comes from the lowest enthalpy or free energy of formation in 
	 * the network.
	 * 
	 * @param uniIsomers The set of unimolecular isomers in the network
	 * @param multiIsomers The set of multimolecular isomers in the network
	 * @param stdTemp The reference temperature in K
	 * @return The reference energy in kJ/mol
	 */
	private double getReferenceEnergy(LinkedList<PDepIsomer> uniIsomers,
			LinkedList<PDepIsomer> multiIsomers, Temperature stdTemp) {
		
		double Eref = 1000000.0;
		for (ListIterator<PDepIsomer> iter = uniIsomers.listIterator(); iter.hasNext(); ) {
			PDepIsomer isomer = iter.next();
			for (int i = 0; i < isomer.getNumSpecies(); i++) {
				double E = isomer.calculateH(i, stdTemp);
				if (E < Eref) Eref = E; 
			}
		
		}
		for (ListIterator<PDepIsomer> iter = multiIsomers.listIterator(); iter.hasNext(); ) {
			PDepIsomer isomer = iter.next();
			for (int i = 0; i < isomer.getNumSpecies(); i++) {
				double E = isomer.calculateH(i, stdTemp);
				if (E < Eref) Eref = E; 
			}
		
		}
		return Eref;
	}

	/**
	 * Determines a suitable grain size for the calculation. Too small a grain
	 * size will cause the simulation to take a very long time to conduct; too 
	 * large a grain size will cause the simulation results to not be very
	 * precise. The choice is hopefully made such that the grains get larger
	 * as the system does, so as to not get completely bogged down in the
	 * larger networks.
	 * 
	 * @param grainMaxEnergy The maximum energy grain in kJ/mol
	 * @param numUniWells The number of unimolecular isomers in the network.
	 * @return
	 */
	private double getGrainSize(double grainMaxEnergy, int numUniWells) {
		if (numUniWells < 5)
			return grainMaxEnergy / 200;
		else if (numUniWells < 10)
			return grainMaxEnergy / 100;
		else if (numUniWells < 20)
			return grainMaxEnergy / 50;
		else
			return grainMaxEnergy / 20;
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

}

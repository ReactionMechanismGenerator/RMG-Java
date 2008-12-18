//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator                                            
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!     
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!  
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED 
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR  
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;  
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//! 
//!******************************************************************************

package jing.rxn;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.StringTokenizer;
import jing.chem.LennardJones;
import jing.chem.Species;
import jing.chem.SpectroscopicData;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxnSys.ReactionSystem;

/**
 * Contains data members and methods for interacting with FAME, a pressure-
 * dependent rate coefficient estimator developed by Allen. The method is based
 * on the steady-state/reservoir-state method of Green and Bhatti, and replaces
 * the Chemdis module as the primary pressure-dependent kinetics estimator.
 * <p>
 * For more information on FAME, see the following references:
 * <p>
 * A. Y. Chang, J. W. Bozzelli, and A. M. Dean. "Kinetic Analysis of Complex
 * Chemical Activation and Unimolecular Dissociation Reactions using QRRK
 * Theory and the Modified Strong Collision Approximation." Z. Phys. Chem
 * 214 (11), p. 1533-1568 (2000).
 * 
 * @author jwallen
 */
public class FastMasterEqn implements PDepKineticsEstimator {

	private static int runCount = 0;
	
	/**
	 * Runs a pressure-dependent calculation by preparing the input file,
	 * calling the FAME executable, parsing the output file, and updating the
	 * network/system accordingly.
	 * @param pdn The pressure-dependent reaction network of interest
	 * @param rxnSystem The reaction system of interest
	 */
	public void runPDepCalculation(PDepNetwork pdn, ReactionSystem rxnSystem) {
		
		// Get working directory (to find FAME executable)
		String dir = System.getProperty("RMG.workingDirectory");
		
		// Skip if network does not contain any unimolecular isomer wells
		if (!pdn.isActive() && pdn.getIsChemAct()) {
        	Temperature t = rxnSystem.getPresentTemperature();
  		    pdn.setKLeak(rxnSystem.getIndex(), pdn.getEntryReaction().calculateTotalRate(t));
        	return;
        }
        
		// Determine wells and reactions
		LinkedList uniWells = new LinkedList();
		LinkedList multiWells = new LinkedList();
		LinkedList reactions = new LinkedList();
		
        // Collect well and reaction parameters
		getUniWells(pdn, uniWells);
		getMultiWells(pdn, multiWells);
		getReactions(pdn, reactions, uniWells, multiWells);

		if (reactions.size() == 0) {
			System.out.println("Warning: Empty pressure-dependent network detected. Skipping.");
			return;
		}
		
		// Create FAME input files
		writeInputFile(pdn, rxnSystem, uniWells, multiWells, reactions);
		
		// Touch FAME output file
		touchOutputFile();
		
		// FAME system call
		try {
           	System.out.println("Running FAME...");
			String[] command = {dir + "/software/fame/fame.exe"};
           	File runningDir = new File("fame/");
            Process fame = Runtime.getRuntime().exec(command, null, runningDir);                     
            InputStream ips = fame.getInputStream();
            InputStreamReader is = new InputStreamReader(ips);
            BufferedReader br = new BufferedReader(is);
            String line = null;
            while ( (line = br.readLine()) != null) {
            	System.out.println(line);
            }
            int exitValue = fame.waitFor();
        }
        catch (Exception e) {
        	System.out.println("Error while executing FAME!");
        	System.exit(0);
        }
        
		// Parse FAME output file and update accordingly
        readOutputFile(pdn, rxnSystem, uniWells, multiWells);
        pdn.updateKLeak();
        
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
	 * Creates the input file needed by Chemdis that represents a pressure-
	 * dependent reaction network.
	 * @param pdn The reaction network of interest
	 * @param rxnSystem The reaction system of interest
	 */
	public void writeInputFile(PDepNetwork pdn, ReactionSystem rxnSystem,
			LinkedList uniWells, LinkedList multiWells, LinkedList reactions) {
		
		Temperature stdTemp = new Temperature(298, "K");
                
		// Collect simulation parameters
		Temperature temperature = rxnSystem.getPresentTemperature();
		Pressure pressure = rxnSystem.getPresentPressure();
		BathGas bathGas = new BathGas(rxnSystem);
		
		int numChebTempPolys = 4, numChebPressPolys = 4;
		
		// Determine reference energies (units are kcal/mol)
		double Eref = 1000000.0, Href = 1000000.0, Gref = 1000000.0;
		for (int i = 0; i < uniWells.size(); i++) {
			Species isomer = ((PDepWell) uniWells.get(i)).getIsomer();
			double E = isomer.calculateH(stdTemp);
			double H = isomer.calculateH(stdTemp);
			double G = isomer.calculateG(stdTemp);
			if (E < Eref) {
				Eref = E; Href = H; Gref = G;
			}
		}
		for (int n = 0; n < multiWells.size(); n++) {
			MultiWell isomer = ((MultiWell) multiWells.get(n));
			for (int i = 0; i < isomer.getNumberOfSpecies(); i++) {
				double E = isomer.calculateH(i, stdTemp);
				double H = isomer.calculateH(i, stdTemp);
				double G = isomer.calculateG(i, stdTemp);
				if (E < Eref) {
					Eref = E; Href = H; Gref = G;
				}
			}
		}
		
		double grainMaxEnergy = getGrainMaxEnergy(uniWells, multiWells, 2100, Eref); // [=] kJ/mol
		double grainSize = getGrainSize(grainMaxEnergy, uniWells.size()); // [=] kJ/mol
		
		// Create the simulation parameters file fame/simData.txt
		try {
        	
			File simData = new File("fame/fame_input.txt");
			FileWriter fw = new FileWriter(simData);

			fw.write("# FAME input for RMG-generated partial network " + String.valueOf(PDepNetwork.getID()) + "\n");
			fw.write("Temperatures                        7\n");
			fw.write("300 K\n600 K\n900 K\n1200 K\n1500 K\n1800 K\n2100 K\n");
			fw.write("Pressures                           5\n");
			fw.write("0.01 bar\n0.1 bar\n1 bar\n10 bar\n100 bar\n");
            fw.write("Grain size                          " + grainSize + " kJ/mol\n");
			fw.write("Maximum grain energy                " + grainMaxEnergy + " kJ/mol\n");
			fw.write("Number of unimolecular wells        " + uniWells.size() + "\n");
			fw.write("Number of multimolecular wells      " + multiWells.size() + "\n");
			fw.write("Number of reactions                 " + reactions.size() + "\n");
			fw.write("Exponential down parameter          " + bathGas.getExpDownParam() + " kJ/mol\n");
			fw.write("Bath gas LJ sigma parameter         " + bathGas.getLJSigma() + " m\n");
			fw.write("Bath gas LJ epsilon parameter       " + bathGas.getLJEpsilon() + " J\n");
			fw.write("Bath gas molecular weight           " + bathGas.getMolecularWeight() + " g/mol\n");
			fw.write("Number of Chebyshev temperatures    " + numChebTempPolys + "\n");
			fw.write("Number of Chebyshev pressures       " + numChebPressPolys + "\n");
			fw.write("\n");

			for (int i = 0; i < uniWells.size(); i++) {
				
				Species isomer = ((PDepWell) uniWells.get(i)).getIsomer();
				isomer.calculateLJParameters();
				
				fw.write("# Unimolecular well " + Integer.toString(i+1) + ": " + isomer.getName() + "(" + Integer.toString(isomer.getID()) + ")" + "\n");
				fw.write("Ground-state energy                 " + (isomer.calculateH(stdTemp) - Eref) * 4.184 + " kJ/mol\n");
				fw.write("Enthalpy of formation               " + (isomer.calculateH(stdTemp) - Href) * 4.184 + " kJ/mol\n");
				fw.write("Free energy of formation            " + (isomer.calculateG(stdTemp) - Gref) * 4.184 + " kJ/mol\n");
				fw.write("LJ sigma parameter                  " + (isomer.getLJ().getSigma() * 1e-10) + " m\n");
				fw.write("LJ epsilon parameter                " + (isomer.getLJ().getEpsilon() * 1.381e-23) + " J\n");
				fw.write("Molecular weight                    " + isomer.getMolecularWeight() + " g/mol\n");
				
				SpectroscopicData data = isomer.getSpectroscopicData();
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
			
			for (int i = 0; i < multiWells.size(); i++) {
				
				MultiWell isomer = (MultiWell) multiWells.get(i);

				fw.write("# Multimolecular well " + Integer.toString(i+1) + ": ");
				fw.write(isomer.getSpecies(0).getName() + "(" + Integer.toString(isomer.getSpecies(0).getID()) + ")");
				for (int j = 1; j < isomer.getNumberOfSpecies(); j++)
					fw.write(" + " + isomer.getSpecies(j).getName() + "(" + Integer.toString(isomer.getSpecies(j).getID()) + ")");
				fw.write("\n");
				fw.write("Number of species                   " + Integer.toString(isomer.getNumberOfSpecies()) + "\n");
				fw.write("Ground-state energy                 ");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++)
					fw.write((isomer.calculateH(j, stdTemp) - Href) * 4.184 + " kJ/mol    ");
				fw.write("\n");
				fw.write("Enthalpy of formation               ");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++)
					fw.write((isomer.calculateH(j, stdTemp) - Href) * 4.184 + " kJ/mol    ");
				fw.write("\n");
				fw.write("Free energy of formation            ");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++)
					fw.write((isomer.calculateG(j, stdTemp) - Gref) * 4.184 + " kJ/mol    ");
				fw.write("\n");
				
				fw.write("Harmonic oscillators                ");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++)
					fw.write(isomer.getSpecies(j).getSpectroscopicData().getVibrationCount() + " ");
				fw.write("\n");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++) {
					SpectroscopicData data = isomer.getSpecies(j).getSpectroscopicData();
					for (int k = 0; k < data.getVibrationCount(); k++)
						fw.write(data.getVibration(k) + " cm^-1\n");
				}
				
				fw.write("Rigid rotors                        ");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++)
					fw.write(isomer.getSpecies(j).getSpectroscopicData().getRotationCount() + " ");
				fw.write("\n");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++) {
					SpectroscopicData data = isomer.getSpecies(j).getSpectroscopicData();
					for (int k = 0; k < data.getRotationCount(); k++)
						fw.write(data.getRotation(k) + " cm^-1\n");
				}
				
				fw.write("Hindered rotors                     ");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++)
					fw.write(isomer.getSpecies(j).getSpectroscopicData().getHinderedCount() + " ");
				fw.write("\n");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++) {
					SpectroscopicData data = isomer.getSpecies(j).getSpectroscopicData();
					for (int k = 0; k < data.getHinderedCount(); k++)
						fw.write(data.getHinderedFrequency(k) + " cm^-1\n");
					for (int k = 0; k < data.getHinderedCount(); k++)
						fw.write(data.getHinderedBarrier(k) + " cm^-1\n");
				}
				fw.write("Symmetry number                     ");
				for (int j = 0; j < isomer.getNumberOfSpecies(); j++)
					fw.write(isomer.getSpecies(j).getSpectroscopicData().getSymmetryNumber() + " ");
				fw.write("\n");
				
				fw.write("\n");
			}
			
			for (int i = 0; i < reactions.size(); i++) {
				
				PDepPathReaction pdpr = (PDepPathReaction) reactions.get(i);
				
				// Determine isomer(s) associated with reactant
				int isomer1 = getAssociatedIsomer(pdpr.getReactantList(), uniWells, multiWells); 
				int isomer2 = getAssociatedIsomer(pdpr.getProductList(), uniWells, multiWells); 
				if (isomer1 < 1 || isomer2 < 1) {
					continue;
				}
				
				String comment = "# Reaction " + Integer.toString(i+1) + ": ";
				if (isomer1 > uniWells.size())
					comment += ((MultiWell) multiWells.get(isomer1-uniWells.size()-1)).getSpeciesNames();
				else
					comment += ((PDepWell) uniWells.get(isomer1-1)).getIsomer().getName();
				comment += " <---> ";
				if (isomer2 > uniWells.size())
					comment += ((MultiWell) multiWells.get(isomer2-uniWells.size()-1)).getSpeciesNames();
				else
					comment += ((PDepWell) uniWells.get(isomer2-1)).getIsomer().getName();
				comment += "\n";
				
				// Calculate transition state energy = ground-state energy of reactant (isomer 1) + activation energy
				double E0 = 0.0;
				if (isomer1 > uniWells.size())
					E0 = ((MultiWell) multiWells.get(isomer1-uniWells.size()-1)).calculateH(stdTemp);
				else
					E0 = ((PDepWell) uniWells.get(isomer1-1)).getIsomer().calculateH(stdTemp);
				E0 += pdpr.getKinetics().getEValue();
						
				fw.write(comment);
				fw.write("Isomer 1                            " + isomer1 + "\n");
				fw.write("Isomer 2                            " + isomer2 + "\n");
				fw.write("Ground-state energy                 " + (E0 - Eref) * 4.184 + " kJ/mol\n");
				fw.write("Arrhenius preexponential            " + pdpr.getKinetics().getAValue() + " s^-1\n");
				fw.write("Arrhenius activation energy         " + pdpr.getKinetics().getEValue() * 4.184 + " kJ/mol\n");
				fw.write("Arrhenius temperature exponent      " + pdpr.getKinetics().getNValue() + "\n");
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
	 * Parses a Chemdis output file and updates the reaction network and system
	 * accordingly.
	 * @param pdn The pressure-dependent reaction network of interest
	 * @param rxnSystem The reaction system of interest
	 */
	public static void readOutputFile(PDepNetwork pdn, ReactionSystem rxnSystem,
			LinkedList uniWells, LinkedList multiWells) {
        
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
			
			if (numUniWells > 5) {
				numUniWells--;
				numUniWells++;
			}
			
			// Initialize temperature and pressure variables based on read values
			Temperature tLow = new Temperature(Tmin, "K");
			Temperature tHigh = new Temperature(Tmax, "K");
			Pressure pLow = new Pressure(Pmin, "bar");
			Pressure pHigh = new Pressure(Pmax, "bar");
							
			// Determine entrant reaction
			LinkedList entryReactant = pdn.getReactant();
			if (entryReactant.size() == 1) {
				for (int i = 0; i < numUniWells; i++) {
					if (((PDepWell) uniWells.get(i)).getIsomer().equals((Species) entryReactant.get(0)))
						entryWell = i + 1;
				}
			}
			else {
				for (int n = 0; n < numMultiWells; n++) {
					if (((MultiWell) multiWells.get(n)).getSpeciesList().equals(entryReactant))
						entryWell = numUniWells + n + 1;
				}
			}
		
			// Reset reaction lists
			LinkedList pDepNetReactionList = pdn.getPDepNetReactionList();
			pDepNetReactionList.clear();
        	LinkedList pDepNonincludedReactionList = pdn.getPDepNonincludedReactionList();
			pDepNonincludedReactionList.clear();
        
			int netReactionCount = 0;
			int nonincludedReactionCount = 0;
			
			// Read Chebyshev coefficients for each reaction
			for (int i = 0; i < numUniWells + numMultiWells; i++) {
				for (int j = 0; j < numUniWells + numMultiWells; j++) {
					if (i != j) {
					
						double[][] alpha = new double[numTemperatures][numPressures];

						// Comment line at start of reaction
						str = br.readLine().trim();

						// Read Chebyshev coefficients from file
						for (int t = 0; t < numTemperatures; t++) {
							str = br.readLine().trim();
							StringTokenizer tkn = new StringTokenizer(str);
							for (int p = 0; p < numPressures; p++)
								alpha[t][p] = Double.parseDouble(tkn.nextToken());

						}

						// Create Chebyshev polynomial object for current rate coefficient
						ChebyshevPolynomials cp = new ChebyshevPolynomials(
								numTemperatures, tLow, tHigh, numPressures, pLow, pHigh,
								alpha);

						// Determine reactants
						LinkedList reactant = new LinkedList();
						if (i < numUniWells) {
							reactant.add( ((PDepWell) uniWells.get(i)).getIsomer() );
						}
						else {
							MultiWell mw = (MultiWell) multiWells.get(i - numUniWells);
							for (int l = 0; l < mw.getNumberOfSpecies(); l++) {
								reactant.add(mw.getSpecies(l));
							}
						}

						// Determine products
						LinkedList product = new LinkedList();
						if (j < numUniWells) {
							product.add( ((PDepWell) uniWells.get(j)).getIsomer() );
						}
						else {
							MultiWell mw = (MultiWell) multiWells.get(j - numUniWells);
							for (int l = 0; l < mw.getNumberOfSpecies(); l++) {
								product.add(mw.getSpecies(l));
							}
						}

						// Initialize net reaction
						PDepNetReaction pdnr = new PDepNetReaction(reactant, product, cp);

						// Add net reaction to appropriate linked list
						boolean includeReactants = false;
						boolean includeProducts = false;
						
						if (i < numUniWells) {
							PDepWell well1 = (PDepWell) uniWells.get(i);
							includeReactants = pdn.includeAsIsomer(well1.getIsomer());
						}
						else {
							MultiWell well1 = (MultiWell) multiWells.get(i - numUniWells);
							includeReactants = well1.includeAsIsomer(pdn);
						}
						
						if (j < numUniWells) {
							PDepWell well2 = (PDepWell) uniWells.get(j);
							includeProducts = pdn.includeAsIsomer(well2.getIsomer());
						}
						else {
							MultiWell well2 = (MultiWell) multiWells.get(j - numUniWells);
							includeProducts = well2.includeAsIsomer(pdn);
						}
						
						// Only add forward reactions
						//if (i > j) {
						if (j + 1 == entryWell) {
							if (includeReactants && includeProducts) {
								pDepNetReactionList.add(pdnr);
								netReactionCount++;
							}
							else {
								pDepNonincludedReactionList.add(pdnr);
								nonincludedReactionCount++;
							}
						}
						
						// Skip blank line between records
						str = br.readLine().trim();

					}
				}
			}
			
			// Close file when finished
			br.close();
			
			System.out.println(Integer.toString(netReactionCount) + " net reactions and " + 
					Integer.toString(nonincludedReactionCount) + " nonincluded reactions in network");
		}
		catch(IOException e) {
			System.out.println("Error: Unable to read from file \"fame_output.txt\".");
			System.exit(0);
		}
		catch(Exception e) {
			System.out.println(e.getMessage());
		}
        
    }

	private int getAssociatedIsomer(LinkedList speciesList, LinkedList uniWells, 
			LinkedList multiWells) {
		if (speciesList.size() == 1) {
			for (int j = 0; j < uniWells.size(); j++) {
				PDepWell pdw = (PDepWell) uniWells.get(j);
				if (pdw.getIsomer().equals(speciesList.get(0)))
					return (j + 1);
			}
		}
		else if (speciesList.size() > 1) {
			for (int j = 0; j < multiWells.size(); j++) {
				MultiWell mw = (MultiWell) multiWells.get(j);
				if (mw.getSpeciesList().equals(speciesList))
					return (uniWells.size() + j + 1);
			}
		}
		return 0;
	}

	private double getGrainMaxEnergy(LinkedList uniWells, LinkedList multiWells, double T, double Eref) {

		double Emax = 0.0;
		Temperature stdTemp = new Temperature(298, "K");
        
		for (int i = 0; i < uniWells.size(); i++) {
			Species isomer = ((PDepWell) uniWells.get(i)).getIsomer();
			double E = isomer.calculateH(stdTemp);
			if (E > Emax)
				Emax = E;
		}
		for (int n = 0; n < multiWells.size(); n++) {
			MultiWell isomer = ((MultiWell) multiWells.get(n));
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
	 * Returns the unimolecular wells in the pressure-dependent network.
	 * @param pdn The pressure-dependent network of interest
	 * @param uniWells A list in which to place the unimolecular wells
	 */
	private void getUniWells(PDepNetwork pdn, LinkedList uniWells) {
		uniWells.clear();
		for (ListIterator iter = pdn.getPDepWellListIterator(); iter.hasNext(); ) {
			PDepWell pdw = (PDepWell) iter.next();
			uniWells.add(pdw);
		}
		for (Iterator iter = pdn.getNonIncludedSpecies(); iter.hasNext(); ) {
			Species spe = (Species) iter.next();
			PDepWell pdw = new PDepWell(spe);
			uniWells.add(pdw);
		}
		
	}
	
	/**
	 * Returns the multimolecular wells in the pressure-dependent network.
	 * @param pdn The pressure-dependent network of interest
	 * @param multiWells A list in which to place the multimolecular wells
	 */
	private void getMultiWells(PDepNetwork pdn, LinkedList multiWells) {
		multiWells.clear();
		for (ListIterator iter = pdn.getPDepWellListIterator(); iter.hasNext(); ) {
			PDepWell pdw = (PDepWell) iter.next();
			for (Iterator iter2 = pdw.getPathsIterator(); iter2.hasNext(); ) {
				PDepPathReaction pdpr = (PDepPathReaction) iter2.next();
				if (pdpr.getProductList().size() > 1)
					multiWells.add(new MultiWell(pdpr.getProductList()));
				else if (pdpr.getReactantList().size() > 1)
					multiWells.add(new MultiWell(pdpr.getProductList()));
			}
		}
	}

	/**
	 * Determines the unique reactions in the pressure-dependent network.
	 * @param pdn The pressure-dependent network of interest
	 * @param reactions A list in which to place the reactions
	  */
	private void getReactions(PDepNetwork pdn, LinkedList reactions,
			LinkedList uniWells, LinkedList multiWells) {
		reactions.clear();
		for (ListIterator iter = pdn.getPDepWellListIterator(); iter.hasNext(); ) {
			PDepWell pdw = (PDepWell) iter.next();
			for (Iterator iter2 = pdw.getPathsIterator(); iter2.hasNext(); ) {
				PDepPathReaction pdpr = (PDepPathReaction) iter2.next();
				// Only add if reactants and products are both included in isomer lists
				int isomer1 = getAssociatedIsomer(pdpr.getReactantList(), uniWells, multiWells); 
				int isomer2 = getAssociatedIsomer(pdpr.getProductList(), uniWells, multiWells); 
				if (isomer1 > 0 && isomer2 > 0) 
					reactions.add(pdpr);
			}
		}
	}

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
	
	//==========================================================================
	//
	//	MultiWell class
	//
	
	/**
	 * Represents a multimolecular pressure-dependent well.
	 * @author jwallen
	 */
	private class MultiWell {
		
		private LinkedList speciesList;
		
		public MultiWell(LinkedList spe) {
			speciesList = spe;
		}
		
		public LinkedList getSpeciesList() {
			return speciesList;
		}
		
		public Species getSpecies(int i) {
			return (Species) speciesList.get(i);
		}
		
		public int getNumberOfSpecies() {
			return speciesList.size();
		}
		
		@Override
		public String toString() {
			String str = ((Species) speciesList.get(0)).toString();
			for (int i = 1; i < speciesList.size(); i++) {
				str += " + " + ((Species) speciesList.get(i)).toString();
			}
			return str;
		}

		public double calculateG(Temperature temperature) {
			double G = 0;
			for (int i = 0; i < speciesList.size(); i++) {
				G += ((Species) speciesList.get(i)).calculateG(temperature);
			}
			return G;
		}

		public double calculateG(int species, Temperature temperature) {
			return ((Species) speciesList.get(species)).calculateG(temperature);
		}
		
		public double calculateH(Temperature temperature) {
			double H = 0;
			for (int i = 0; i < speciesList.size(); i++) {
				H += ((Species) speciesList.get(i)).calculateH(temperature);
			}
			return H;
		}
		
		public double calculateH(int species, Temperature temperature) {
			return ((Species) speciesList.get(species)).calculateH(temperature);
		}

		public String getSpeciesNames() {
			String str = ((Species) speciesList.get(0)).getName();
			for (int i = 1; i < speciesList.size(); i++) {
				str += " + " + ((Species) speciesList.get(i)).getName();
			}
			return str;
		}
		
		public boolean includeAsIsomer(PDepNetwork pdn) {
			for (ListIterator iter = pdn.getPDepWellListIterator(); iter.hasNext(); ) {
				PDepWell pdw = (PDepWell) iter.next();
				for (Iterator iter2 = pdw.getPathsIterator(); iter2.hasNext(); ) {
					PDepPathReaction pdpr = (PDepPathReaction) iter2.next();
					if (pdpr.getProductList().equals(speciesList))
						return (pdpr.isProduct() || pdpr.isReactant());
				}
			}
			return false;
		}
		
		
	}
	
	//==========================================================================
	//
	//	BathGas class
	//
	
	/**
	 * Represents the bath gas. Calculates parameters of the bath gas as a
	 * weighted average of the components in the bath gas, where the weights
	 * are the mole fractions.
	 * <p>
	 * This class is contained within FastMasterEqn because it is set up to 
	 * primarily interact with it. In the future this class could become
	 * a standalone class.
	 * @author jwallen
	 */
	private class BathGas {
		private double expDownParam = 0.0;
		private double ljSigma = 0.0;
		private double ljEpsilon = 0.0;
		private double molWt = 0.0;
		private HashMap colliders;
		
		public BathGas() {
			colliders = null;
		}
		
		public BathGas(ReactionSystem rxnSystem) {
			colliders = rxnSystem.identifyColliders();
			update();
		}
		
		public double getExpDownParam() {
			return expDownParam;
		}
		
		public double getLJSigma() {
			return ljSigma;
		}
		
		public double getLJEpsilon() {
			return ljEpsilon;
		}
		
		public double getMolecularWeight() {
			return molWt;
		}
		
		public HashMap getColliders() {
			return colliders;
		}
		
		public void setColliders(ReactionSystem rxnSystem) {
			colliders = rxnSystem.identifyColliders();
			update();
		}
		
		/**
		 * Updates the bath gas parameters via a weighted average.
		 */
		private void update() {
			
			// Check for null pointer (i.e. colliders not set)
			if (colliders == null)
				throw new NullPointerException();
			
			// Clear parameters
			expDownParam = 0.0;
			ljSigma = 0.0;
			ljEpsilon = 0.0;
			molWt = 0.0;
			
			// Determine bath gas concentration (i.e. total concentration of colliders)
			double totalConc = 0;
			for (Iterator iter = colliders.values().iterator(); iter.hasNext();) {
				totalConc += ((Double) iter.next()).doubleValue();
			}

			// Calculate values as weighted average
			for (Iterator iter = colliders.keySet().iterator(); iter.hasNext();) {

				Object key = iter.next();
				if (key instanceof Species) {
					Species spe = (Species) key;
					double conc = ((Double) colliders.get(spe)).doubleValue();
					double mf = conc/totalConc;
					
					molWt += mf * spe.getMolecularWeight();
					expDownParam += mf * spe.getDeltaEDown();
					ljSigma += mf * spe.getLJ().getSigma();
					ljEpsilon += mf * spe.getLJ().getEpsilon();
				}
				else if (key instanceof String) {
					String name = (String) key;
					LennardJones lj = new LennardJones();
					ljSigma = lj.getSigma();
					ljEpsilon = lj.getSigma();
					
					if (name.equals("Ar") || name.equals("AR")) {
						expDownParam = 374.0;
						molWt = 39.95;
					}
					else if (name.equals("N2")) {
						expDownParam = 461.0;
						molWt = 28.01;
					}
					else if (name.equals("He") || name.equals("HE")) {
						expDownParam = 291.0;
						molWt = 4.00;
					}
					else {
						System.out.println("unknown colliders: " + name);
						System.exit(0);
					}
				}
				else {
					System.out.println("unknown colliders: " + key.toString());
					System.exit(0);
				}
			}
			
			// Convert to units used by FAME
			expDownParam *= 2.9979e10 * 6.626e-34 * 6.022e23 / 1000; // cm^-1 --> kJ/mol
			ljSigma *= 1e-10; // A --> m
			ljEpsilon *= 1.381e-23; // K --> J
			
		}
		
	}
	
}

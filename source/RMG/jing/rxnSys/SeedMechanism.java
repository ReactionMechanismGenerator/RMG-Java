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



package jing.rxnSys;

import java.io.*;

import jing.mathTool.UncertainDouble;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxn.*;
import jing.chem.*;
import java.util.*;

import jing.chemUtil.*;
import jing.chemParser.*;

/**
 * This is a new class called SeedMechanism.  SeedMechanism is the same class
 * as the old PrimaryKineticLibrary, just with a new (and more appropriate)
 * name.  RMG will automatically include every species and reaction contained
 * in a Seed Mechanism.  Furthermore, the user has the option to pass multiple
 * Seed Mechanisms to RMG.  In the event of a duplicate species/reaction, RMG
 * will use the first instance it finds (i.e. the order of the Seed Mechanisms
 * listed in the condition.txt file is important).
 * 
 * MRH 9-Jun-2009
 */

/*
 * Comments from old PrimaryKineticLibrary:
 * 
 * This is the primary reaction set that any reaction system has to include 
 * into its model.  For example, in combustion system, we build a primary small 
 * molecule reaction set, and every combustion/oxidation system should include 
 * such a primary reaction library.  The reaction / rates are basically from Leeds 
 * methane oxidation mechanism.
 */

public class SeedMechanism {
    
    protected String name; 
    
    protected LinkedHashSet reactionSet = new LinkedHashSet(); 
    
    protected HashMap speciesSet = new HashMap(); 
    
    private boolean generateReactions = false;
    
    public LinkedList allPdepNetworks = new LinkedList();

    // Constructors
    
    public  SeedMechanism(String p_mechName, String p_directoryPath, boolean p_generateReactions, boolean p_fromRestart) throws IOException {
    	name = p_mechName;
		generateReactions = p_generateReactions;
        if ( p_directoryPath == null) throw new NullPointerException("RMG does not recognize Seed Mechanism directory path: Value is null");
        try {
        	read(p_directoryPath,p_fromRestart,p_mechName);
        }
        catch (IOException e) {
        	throw new IOException("Error in reading Seed Mechanism: " + p_mechName + '\n' + e.getMessage());
        }
    }
    
    public SeedMechanism() {
	}

	public void appendSeedMechanism(String new_mechName, String new_directoryPath, boolean p_generateReactions, boolean p_fromRestart) throws IOException {
     	if (p_generateReactions)
			setGenerateReactions(p_generateReactions);
		setName(name + "/" + new_mechName);
    	try {
    		read(new_directoryPath,p_fromRestart,new_mechName);	
    	}
        catch (IOException e) {
        	throw new IOException("Error in reading Seed Mechanism: " + new_mechName + '\n' + e.getMessage());
        }
    }
    
    public LinkedHashSet getSpeciesSet() {
        return new LinkedHashSet(speciesSet.values());
    }
    
    public void read(String p_directoryName, boolean p_fromRestart, String seedMechName) throws IOException {
        System.out.println("Reading seed mechanism from directory " + p_directoryName);
    	HashMap localSpecies = null;
    	LinkedHashSet localReactions = null;
		try {
        	if (!p_directoryName.endsWith("/")) p_directoryName = p_directoryName + "/";
        	
        	if (!p_fromRestart) {
	            String speciesFile = p_directoryName + "species.txt";
	            String reactionFile = p_directoryName + "reactions.txt";
	            String pdepreactionFile = p_directoryName + "pdepreactions.txt";

	        	speciesSet.putAll(readSpecies(speciesFile,seedMechName,"Seed Mechanism: "));
	        	reactionSet.addAll(readReactions(reactionFile,seedMechName,speciesSet,"Seed Mechanism: ",false));
	        	reactionSet.addAll(readPdepReactions(pdepreactionFile,seedMechName,speciesSet,"Seed Mechanism: ",false));
        	}
        	else {
	            String speciesFile = p_directoryName + "coreSpecies.txt";
	            String pdepreactionFile = p_directoryName + "pdepreactions.txt";
	        	
	        	speciesSet.putAll(readSpecies(speciesFile,seedMechName,"Seed Mechanism: "));
	        	reactionSet.addAll(readPdepReactions(pdepreactionFile,seedMechName,speciesSet,"Seed Mechanism: ",false));
        	}
        	return;
        }
        catch (Exception e) {
        	throw new IOException("RMG cannot read entire Seed Mechanism: " 
        			+ p_directoryName + "\n" + e.getMessage());
        }
    }
    

    public LinkedHashSet readReactions(String p_reactionFileName, String p_name, HashMap allSpecies, String source, boolean pkl) throws IOException {
    	LinkedHashSet localReactions = new LinkedHashSet();
        try {
        	FileReader in = new FileReader(p_reactionFileName);
        	BufferedReader data = new BufferedReader(in);
        	
        	double[] multipliers = parseReactionRateUnits(data);
        	double A_multiplier = multipliers[0];
        	double E_multiplier = multipliers[1];
            
        	String line = ChemParser.readMeaningfulLine(data);
        	read: while (line != null) {
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(allSpecies, line, A_multiplier, E_multiplier);
					r.setKineticsSource(source+ p_name,0);
					r.setKineticsComments(" ",0);
					if (pkl) {
						r.setIsFromPrimaryKineticLibrary(true);
						(r.getKinetics())[0].setFromPrimaryKineticLibrary(true);
					}
				}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(line + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(line);
        		
        		localReactions = updateReactionList(r,localReactions,true);
        		
        		line = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return localReactions;
        }
        catch (Exception e) {
        	System.out.println("RMG did not read the following " + source + p_name + " file: " 
        			+ p_reactionFileName + " because " + e.getMessage() );
			e.printStackTrace();
        	return null;
        }
    }
    
    public HashMap readSpecies(String p_speciesFileName, String p_name, String source) throws IOException {
    	HashMap localSpecies = new HashMap();
        try {
        	FileReader in = new FileReader(p_speciesFileName);
        	BufferedReader data = new BufferedReader(in);
        	
        	// step 1: read in structure
        	String line = ChemParser.readMeaningfulLine(data);
        	read: while (line != null) {
        		// GJB allow unreactive species
        		StringTokenizer st = new StringTokenizer(line);
				String name = st.nextToken().trim();
				boolean IsReactive = true;
				if (st.hasMoreTokens()) {
					String reactive = st.nextToken().trim();
					if ( reactive.equalsIgnoreCase("unreactive") ) IsReactive = false;
				}
            	Graph graph;
        		try {
        			graph = ChemParser.readChemGraph(data);
					if (graph == null) throw new IOException("Graph was null");
        		}
        		catch (IOException e) {
        			throw new InvalidChemGraphException("Cannot read species '" + name + "': " + e.getMessage());
        		}
				ChemGraph cg = ChemGraph.make(graph);	
        		Species spe = Species.make(name, cg);
        		// GJB: Turn off reactivity if necessary, but don't let code turn it on
        		// again if was already set as unreactive from input file
        		if(IsReactive==false) spe.setReactivity(IsReactive);
        		localSpecies.put(name, spe);
        		line = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return localSpecies;
        }
        catch (Exception e) {
			throw new IOException("RMG cannot read the \"species.txt\" file in the " + source + p_name + "\n" + e.getMessage());
        }
    }
    
    public LinkedHashSet readPdepReactions(String pdepFileName, String p_name, HashMap allSpecies, String source, boolean pkl) throws IOException {
    	LinkedHashSet localReactions = new LinkedHashSet();
		LinkedList pdepNetworks = getPDepNetworks();
        try {
        	FileReader in = new FileReader(pdepFileName);
        	BufferedReader data = new BufferedReader(in);

        	double[] multipliers = parseReactionRateUnits(data);
        	double A_multiplier = multipliers[0];
        	double E_multiplier = multipliers[1];
            
        	String nextLine = ChemParser.readMeaningfulLine(data);
        	read: while (nextLine != null) {	
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(allSpecies, nextLine, A_multiplier, E_multiplier);
        		}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(nextLine + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(nextLine);
                
        		/*
        		 * Read the next line and determine what to do based on the
        		 * 	presence/absence of keywords
        		 */
        		nextLine = ChemParser.readMeaningfulLine(data);
        		boolean continueToReadRxn = true;
        		
        		// Initialize all of the possible pdep variables
        		HashMap thirdBodyList = new HashMap();
        		UncertainDouble uA = new UncertainDouble(0.0, 0.0, "Adder");
        		UncertainDouble un = new UncertainDouble(0.0, 0.0, "Adder");
        		UncertainDouble uE = new UncertainDouble(0.0, 0.0, "Adder");
        		ArrheniusKinetics low = new ArrheniusKinetics(uA, un, uE, "", 0, "", "");
        		double a = 0.0;
        		double T3star = 0.0;
        		double Tstar = 0.0;
        		double T2star = 0.0;
        		boolean troe7 = false;
        		int numPLOGs = 0;
        		PDepArrheniusKinetics pdepkineticsPLOG = new PDepArrheniusKinetics(numPLOGs);
        		
        		
        		/*
        		 * When reading in the auxillary information for the pdep reactions,
        		 * 	let's not assume the order is fixed (i.e. third-bodies and
        		 * 	their efficiencies, then lindemann, then troe)
        		 * The order of the if statement is important as the "troe" and
        		 * 	"low" lines will also contain a "/"; thus, the elseif contains
        		 * 	"/" needs to be last.
        		 */
        		/*
        		 * Additions by MRH on 26JUL2010:
        		 * 	Allowing RMG to read-in PLOG and CHEB formatting
        		 */
        		while (continueToReadRxn) {
        			if (nextLine == null) {
        				continueToReadRxn = false;
        			} else if (nextLine.toLowerCase().contains("troe")) {
	        			// read in troe parameters
	            		StringTokenizer st = new StringTokenizer(nextLine, "/");
	            		String temp = st.nextToken().trim();	// TROE
	            		String troeString = st.nextToken().trim();	// List of troe parameters
	            		st = new StringTokenizer(troeString);
	                    int n = st.countTokens();
	                    if (n != 3 && n != 4) throw new InvalidKineticsFormatException("Troe parameter number = "+n + " for reaction: " + r.toString());
	            
	              		a = Double.parseDouble(st.nextToken().trim());
	            		T3star = Double.parseDouble(st.nextToken().trim());
	            		Tstar = Double.parseDouble(st.nextToken().trim());

	            		if (st.hasMoreTokens()) {
	            			troe7 = true;
	            			T2star = Double.parseDouble(st.nextToken().trim());
	               		}
	            		nextLine = ChemParser.readMeaningfulLine(data);
	        		} else if (nextLine.toLowerCase().contains("low")) {
	        			// read in lindemann parameters
	            		StringTokenizer st = new StringTokenizer(nextLine, "/");
	            		String temp = st.nextToken().trim();	// LOW
	            		String lowString = st.nextToken().trim();	// Modified Arrhenius parameters
	            		/*
	            		 * MRH 17Feb2010:
	            		 * 	The units of the k_zero (LOW) Arrhenius parameters are different from the units of
	            		 * 	k_inf Arrhenius parameters by a factor of cm3/mol, hence the getReactantNumber()+1
	            		 */
	            		low = ChemParser.parseSimpleArrheniusKinetics(lowString, A_multiplier, E_multiplier, r.getReactantNumber()+1);
	            		nextLine = ChemParser.readMeaningfulLine(data);
	        		} else if (nextLine.contains("CHEB")) {
        				// Read in the Tmin/Tmax and Pmin/Pmax information
        				StringTokenizer st_cheb = new StringTokenizer(nextLine,"/");
        				String nextToken = st_cheb.nextToken(); // Should be TCHEB or PCHEB
        				StringTokenizer st_minmax = new StringTokenizer(st_cheb.nextToken());
        				double Tmin = 0.0; double Tmax = 0.0; double Pmin = 0.0; double Pmax = 0.0;
        				if (nextToken.trim().equals("TCHEB")) {
        					Tmin = Double.parseDouble(st_minmax.nextToken());
        					Tmax = Double.parseDouble(st_minmax.nextToken());
        				} else {
        					Pmin = Double.parseDouble(st_minmax.nextToken());
        					Pmax = Double.parseDouble(st_minmax.nextToken());
        				}
        				nextToken = st_cheb.nextToken(); // Should be PCHEB or TCHEB
        				st_minmax = new StringTokenizer(st_cheb.nextToken());
        				if (nextToken.trim().equals("TCHEB")) {
        					Tmin = Double.parseDouble(st_minmax.nextToken());
        					Tmax = Double.parseDouble(st_minmax.nextToken());
        				} else {
        					Pmin = Double.parseDouble(st_minmax.nextToken());
        					Pmax = Double.parseDouble(st_minmax.nextToken());
        				}
        				// Read in the N/M values (number of polynomials in the Temp and Press domain)
        				nextLine = ChemParser.readMeaningfulLine(data);
        				st_cheb = new StringTokenizer(nextLine,"/");
        				nextToken = st_cheb.nextToken(); // Should be CHEB
        				st_minmax = new StringTokenizer(st_cheb.nextToken());
        				int numN = Integer.parseInt(st_minmax.nextToken());
        				int numM = Integer.parseInt(st_minmax.nextToken());
        				// Read in the coefficients
        				nextLine = ChemParser.readMeaningfulLine(data);
        				double[] unorderedChebyCoeffs = new double[numN*numM];
        				int chebyCoeffCounter = 0;
        				while (nextLine != null && nextLine.contains("CHEB")) {
        					st_cheb = new StringTokenizer(nextLine,"/");
        					nextToken = st_cheb.nextToken(); // Should be CHEB
        					st_minmax = new StringTokenizer(st_cheb.nextToken());
        					while (st_minmax.hasMoreTokens()) {
        						unorderedChebyCoeffs[chebyCoeffCounter] = Double.parseDouble(st_minmax.nextToken());
        						++chebyCoeffCounter;
        					}
        					nextLine = ChemParser.readMeaningfulLine(data);
        				}
        				// Order the chebyshev coefficients
        				double[][] chebyCoeffs = new double[numN][numM];
        				for (int numRows=0; numRows<numN; numRows++) {
        					for (int numCols=0; numCols<numM; numCols++) {
        						chebyCoeffs[numRows][numCols] = unorderedChebyCoeffs[numM*numRows+numCols];
        					}
        				}
        				// Make the ChebyshevPolynomials, PDepReaction, and add to list of PDepReactions
        				ChebyshevPolynomials coeffs = new ChebyshevPolynomials(
        						numN,new Temperature(Tmin,"K"),new Temperature(Tmax,"K"),
        						numM,new Pressure(Pmin,"bar"),new Pressure(Pmax,"bar"),
        						chebyCoeffs);
	        			PDepIsomer reactants = new PDepIsomer(r.getStructure().getReactantList());
	        			PDepIsomer products = new PDepIsomer(r.getStructure().getProductList());
	        			PDepRateConstant pdepRC = new PDepRateConstant(coeffs);
	        			PDepReaction pdeprxn = new PDepReaction(reactants,products,pdepRC);
	        			allPdepNetworks.add(pdeprxn);
	        			continueToReadRxn = false;
	        		} else if (nextLine.contains("PLOG")) {
	        			while (nextLine != null && nextLine.contains("PLOG")) {
		        			// Increase the PLOG counter
		        			++numPLOGs;
		        			// Store the previous PLOG information in temporary arrays
		        			Pressure[] previousPressures = new Pressure[numPLOGs];
		        			ArrheniusKinetics[] previousKinetics = new ArrheniusKinetics[numPLOGs];
		        			for (int previousNumPLOG=0; previousNumPLOG<numPLOGs-1; ++previousNumPLOG) {
		        				previousPressures[previousNumPLOG] = pdepkineticsPLOG.getPressures(previousNumPLOG);
		        				previousKinetics[previousNumPLOG] = pdepkineticsPLOG.getKinetics(previousNumPLOG);
		        			}
		        			// Read in the new PLOG information, and add this to the temporary array
		        			PDepArrheniusKinetics newpdepkinetics = parsePLOGline(nextLine);
		        			previousPressures[numPLOGs-1] = newpdepkinetics.getPressures(0);
		        			previousKinetics[numPLOGs-1] = newpdepkinetics.getKinetics(0);
		        			// Re-initialize pdepkinetics and populate with stored information
		        			pdepkineticsPLOG = new PDepArrheniusKinetics(numPLOGs);
		        			pdepkineticsPLOG.setPressures(previousPressures);
		        			pdepkineticsPLOG.setRateCoefficients(previousKinetics);
		        			// Read the next line
		        			nextLine = ChemParser.readMeaningfulLine(data);
	        			}
	        			// Make the PDepReaction
	        			PDepIsomer reactants = new PDepIsomer(r.getStructure().getReactantList());
	        			PDepIsomer products = new PDepIsomer(r.getStructure().getProductList());
	        			PDepRateConstant pdepRC = new PDepRateConstant(pdepkineticsPLOG);
	        			PDepReaction pdeprxn = new PDepReaction(reactants,products,pdepRC);
	        			// Add to the list of PDepReactions
	        			allPdepNetworks.add(pdeprxn);
	        			// Re-initialize the pdepkinetics variable
	            		numPLOGs = 0;
	            		pdepkineticsPLOG = new PDepArrheniusKinetics(numPLOGs);
	            		continueToReadRxn = false;
	        		} else if (nextLine.contains("/")) {
	        			// read in third body colliders + efficiencies
	            		thirdBodyList = ChemParser.parseThirdBodyList(nextLine,allSpecies);
	            		nextLine = ChemParser.readMeaningfulLine(data);
	        		} else {
	        			// the nextLine is a "new" reaction, hence we need to exit the while loop
	        			continueToReadRxn = false;
	        		}
        		}
        		   
        		/*
        		 * Make the pdep reaction, according to which parameters
        		 * 	are present
        		 */
        		
        		if ((a==0.0) && (T3star==0.0) && (Tstar==0.0)) {
        			// Not a troe reaction
        			if (low.getAValue() == 0.0) {
        				// thirdbody reaction
                		ThirdBodyReaction tbr = ThirdBodyReaction.make(r,thirdBodyList);
        				tbr.setKineticsSource(source+ p_name,0);
        				tbr.setKineticsComments(" ",0);
    					if (pkl) {
    						tbr.setIsFromPrimaryKineticLibrary(true);
    						(tbr.getKinetics())[0].setFromPrimaryKineticLibrary(true);
    					}
        				localReactions.add(tbr);
                		Reaction reverse = tbr.getReverseReaction();
        				if (reverse != null) localReactions.add(reverse);
        			} else {
        				// lindemann reaction
                		LindemannReaction tbr = LindemannReaction.make(r,thirdBodyList,low);
        				tbr.setKineticsSource(source+ p_name,0);
        				tbr.setKineticsComments(" ",0);
    					if (pkl) {
    						tbr.setIsFromPrimaryKineticLibrary(true);
    						(tbr.getKinetics())[0].setFromPrimaryKineticLibrary(true);
    					}
        				localReactions.add(tbr);
                		Reaction reverse = tbr.getReverseReaction();
        				if (reverse != null) localReactions.add(reverse);
        			}
        		} else {
        			// troe reaction
            		TROEReaction tbr = TROEReaction.make(r,thirdBodyList, low, a, T3star, Tstar, troe7, T2star);
    				tbr.setKineticsSource(source+ p_name,0);
    				tbr.setKineticsComments(" ",0);
					if (pkl) {
						tbr.setIsFromPrimaryKineticLibrary(true);
						(tbr.getKinetics())[0].setFromPrimaryKineticLibrary(true);
					}
    				localReactions.add(tbr);
            		Reaction reverse = tbr.getReverseReaction();
    				if (reverse != null) localReactions.add(reverse);
        		}
        	}
        	   
            in.close();
        	return localReactions;
        }
        catch (Exception e) {

			/*
			 * 25Jun2009-MRH: When reading the Primary Kinetic Library, we should not require the user to supply
			 * 		troe reactions.  In the instance that no "troeReactions.txt" file exists, inform
			 * 		user of this but continue simulation.
			 */
        	System.out.println("RMG did not find/read pressure-dependent reactions (pdepreactions.txt) " +
        			"in the " + source + p_name + "\n" + e.getMessage());
        	return new LinkedHashSet();
        }
    }
	
    public int size() {
        return speciesSet.size();
    }
    
    public String getName() {
        return name;
    }
    
    public void setName(String p_name) {
        name = p_name;
    }
    
    public LinkedHashSet getReactionSet() {
        return reactionSet;
    }

	public boolean shouldGenerateReactions() {
		return generateReactions;
	}

	public void setGenerateReactions(boolean generateReactions) {
		this.generateReactions = generateReactions;
	}
	
	public double[] parseReactionRateUnits(BufferedReader data) {
		double[] multipliers = new double[2];
    	String line = ChemParser.readMeaningfulLine(data);
    	if (line.startsWith("Unit")) {
    		line = ChemParser.readMeaningfulLine(data);
    		unit: while(!(line.startsWith("Reaction"))) {
    			if (line.startsWith("A")) {
    				StringTokenizer st = new StringTokenizer(line);
    				String temp = st.nextToken();
    				String unit = st.nextToken().trim();
    				if (unit.compareToIgnoreCase("mol/cm3/s") == 0) {
    					multipliers[0] = 1;
    				}
    				else if (unit.compareToIgnoreCase("mol/liter/s") == 0) {
       					multipliers[0] = 1e-3;
    				}
    				else if (unit.compareToIgnoreCase("molecule/cm3/s") == 0) {
    					multipliers[0] = 6.022e23;
    				}
    			}
    			else if (line.startsWith("E")) {
    				StringTokenizer st = new StringTokenizer(line);
    				String temp = st.nextToken();
    				String unit = st.nextToken().trim();
    				if (unit.compareToIgnoreCase("kcal/mol") == 0) {
    					multipliers[1] = 1;
    				}
    				else if (unit.compareToIgnoreCase("cal/mol") == 0) {
       					multipliers[1] = 1e-3;
    				}
    				else if (unit.compareToIgnoreCase("kJ/mol") == 0) {
       					multipliers[1] = 1/4.186;
    				}
    				else if (unit.compareToIgnoreCase("J/mol") == 0) {
       					multipliers[1] = 1/4186;
    				}
    				else if (unit.compareToIgnoreCase("Kelvin") == 0) {
    					multipliers[1] = 1.987e-3;
    				}
    			}
    			line = ChemParser.readMeaningfulLine(data);
    		}
    	}
    	return multipliers;
	}
	
	public LinkedHashSet updateReactionList(Reaction r, LinkedHashSet listOfRxns, boolean generateReverse) {
		Iterator allRxnsIter = listOfRxns.iterator();
		boolean foundRxn = false;
		while (allRxnsIter.hasNext()) {
			Reaction old = (Reaction)allRxnsIter.next();
			if (old.equals(r)) {
				old.addAdditionalKinetics(r.getKinetics()[0],1);
				foundRxn = true;
				break;
			}
		}
		if (!foundRxn) {
			listOfRxns.add(r);
			if (generateReverse) {
	    		Reaction reverse = r.getReverseReaction();
	    		if (reverse != null) {
					listOfRxns.add(reverse);
	    		}
			}
		}
		return listOfRxns;
	}
	
	public PDepArrheniusKinetics parsePLOGline(String line) {
		PDepRateConstant pdepk = new PDepRateConstant();
		
		// Delimit the PLOG line by "/"
		StringTokenizer st_slash = new StringTokenizer(line,"/");
		String dummyString = st_slash.nextToken();
		// Delimit the data of the PLOG line by whitespace
		StringTokenizer st = new StringTokenizer(st_slash.nextToken());
		
		// If reading a chemkin plog line, pressure MUST be in atmospheres
		Pressure p = new Pressure(Double.parseDouble(st.nextToken().trim()),"atm");
		UncertainDouble dA = new UncertainDouble(Double.parseDouble(st.nextToken().trim()),0.0,"A");
		UncertainDouble dn = new UncertainDouble(Double.parseDouble(st.nextToken().trim()),0.0,"A");
		double Ea = Double.parseDouble(st.nextToken().trim());
		// If reading a chemkin plog line, Ea MUST be in cal/mol
		Ea = Ea / 1000;
		UncertainDouble dE = new UncertainDouble(Ea,0.0,"A");
		ArrheniusKinetics k = new ArrheniusKinetics(dA, dn, dE, "", 1, "", "");
		PDepArrheniusKinetics pdepAK = new PDepArrheniusKinetics(1);
		pdepAK.setKinetics(0, p, k);
		
		return pdepAK;
	}
	
	public LinkedList getPDepNetworks() {
		return allPdepNetworks;
	}
    
}


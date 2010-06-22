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
import jing.rxn.*;
import jing.chem.*;
import java.util.*;

import jing.chemUtil.*;
import jing.chemParser.*;



/**
This is the primary kinetics set that any reaction system has to include into its model.  For example, in combustion system, we build a primary small molecule reaction set, and every combustion/oxydation system should include such a primary kinetic library.  The reaction / rates are basically from Leeds methane oxidation mechanism.
*/
public class PrimaryKineticLibrary {
    
    protected String name;		//## attribute name 
    
    protected static LinkedHashSet reactionSet = new LinkedHashSet();		//## attribute reactionSet 
    
    protected HashMap speciesSet = new HashMap();		//## attribute speciesSet 
    
    
    // Constructors
    
    public  PrimaryKineticLibrary(String p_libraryName, String p_directoryPath) throws IOException {
        name = p_libraryName;

        if ( p_directoryPath == null) throw new NullPointerException("PrimaryKineticLibrary directory path");
        try {
        	read(p_directoryPath,p_libraryName);
        }
        catch (IOException e) {
        	throw new IOException("Error reading Primary Kinetic Library: " + name + '\n' + e.getMessage());
        }
    }
    public  PrimaryKineticLibrary() {
    }
    
    public void appendPrimaryKineticLibrary(String new_p_libraryName, String p_directoryPath) throws IOException {
    	// Appends the current PKLib with an additional one, allowing the user
    	// to combine separate PKLibs easily.  GJB 10/05.
     	setName(name+"/"+new_p_libraryName);
    	try {
    		read(p_directoryPath,new_p_libraryName);	
    	}
        catch (IOException e) {
        	throw new IOException("Error reading Primary Kinetic Library: " + new_p_libraryName + '\n' + e.getMessage());
        }
    }
    
    
    
    public LinkedHashSet getSpeciesSet() {
        
        return new LinkedHashSet(speciesSet.values());
        
    }
    
    public void read(String p_directoryName, String p_name) throws IOException {
        try {
        	if (!p_directoryName.endsWith("/")) p_directoryName = p_directoryName + "/";
			System.out.println("Reading Primary Kinetic Library from: "+p_directoryName);
        	
            String speciesFile = p_directoryName + "species.txt";
            String reactionFile = p_directoryName + "reactions.txt";
            String pdepreactionFile = p_directoryName + "pdepreactions.txt";
        	readSpecies(speciesFile);
        	readReactions(reactionFile,p_name);
        	readPdepReactions(pdepreactionFile,p_name);
        	return;
        }
        catch (Exception e) {
        	throw new IOException("Can't read primary kinetic library.\n" + e.getMessage());
        }
    }
    
    public void readReactions(String p_reactionFileName, String p_name) throws IOException {
        try {
        	FileReader in = new FileReader(p_reactionFileName);
        	BufferedReader data = new BufferedReader(in);
        	
        	double A_multiplier = 1;
        	double E_multiplier = 1;
        	
        	String line = ChemParser.readMeaningfulLine(data);
        	if (line.startsWith("Unit")) {
        		line = ChemParser.readMeaningfulLine(data);
        		unit: while(!(line.startsWith("Reaction"))) {
        			if (line.startsWith("A")) {
        				StringTokenizer st = new StringTokenizer(line);
        				String temp = st.nextToken();
        				String unit = st.nextToken().trim();
        				if (unit.compareToIgnoreCase("mol/cm3/s") == 0) {
        					A_multiplier = 1;
        				}
        				else if (unit.compareToIgnoreCase("mol/liter/s") == 0) {
           					A_multiplier = 1e-3;
        				}
        				else if (unit.compareToIgnoreCase("molecule/cm3/s") == 0) {
        					A_multiplier = 6.022e23;
        				}
        			}
        			else if (line.startsWith("E")) {
        				StringTokenizer st = new StringTokenizer(line);
        				String temp = st.nextToken();
        				String unit = st.nextToken().trim();
        				if (unit.compareToIgnoreCase("kcal/mol") == 0) {
        					E_multiplier = 1;
        				}
        				else if (unit.compareToIgnoreCase("cal/mol") == 0) {
           					E_multiplier = 1e-3;
        				}
        				else if (unit.compareToIgnoreCase("kJ/mol") == 0) {
           					E_multiplier = 1/4.186;
        				}
        				else if (unit.compareToIgnoreCase("J/mol") == 0) {
           					E_multiplier = 1/4186;
        				}
        				else if (unit.compareToIgnoreCase("Kelvin") == 0) {
        					E_multiplier = 1.987e-3;
        				}
        			}
        			line = ChemParser.readMeaningfulLine(data);
        		}
        	}
            
        	LinkedHashSet currentPKLReactions = new LinkedHashSet();
        	line = ChemParser.readMeaningfulLine(data);
        	read: while (line != null) {
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(speciesSet, line, A_multiplier, E_multiplier);
        			r.setIsFromPrimaryKineticLibrary(true);
        			(r.getKinetics())[0].setFromPrimaryKineticLibrary(true);
        			// Changed source from "Seed Mechanism" to "Primary Kinetic Library"
					r.setKineticsSource("Primary Kinetic Library: "+ p_name,0);
					r.setKineticsComments(" ",0);
				}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(line + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(line);
        		
        		Iterator prlRxnIter = currentPKLReactions.iterator();
        		boolean foundRxn = false;
        		while (prlRxnIter.hasNext()) {
        			Reaction old = (Reaction)prlRxnIter.next();
        			if (old.equals(r)) {
        				old.addAdditionalKinetics(r.getKinetics()[0],1);
        				foundRxn = true;
        				break;
        			}
        		}
        		if (!foundRxn) {
        			currentPKLReactions.add(r);
	        		Reaction reverse = r.getReverseReaction();
					
	        		if (reverse != null) {
						//reverse.getKinetics().setSource("Seed Mechanism: " + name);
						currentPKLReactions.add(reverse);
	        		}
        		}
        		line = ChemParser.readMeaningfulLine(data);
        	}
        	   
        	reactionSet.addAll(currentPKLReactions);
            in.close();
        	return;
        }
        catch (Exception e) {
//        	throw new IOException("Can't read reaction in primary kinetic library.\n" + e.getMessage());
			/*
			 * 25Jun2009-MRH: When reading the Primary Kinetic Library, we should not require the user to supply
			 * 		non pressure-dependent reactions.  In the instance that no "reactions.txt" file exists, inform
			 * 		user of this but continue simulation.
			 */
        	System.out.println("RMG did not find/read non pressure-dependent reactions (reaction.txt) " +
        			"in the Primary Kinetic Library: " + name + "\n" + e.getMessage());
        }
    }

    public void readSpecies(String p_speciesFileName) throws IOException {
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
        		}
        		catch (InvalidGraphFormatException e) {
        			throw new InvalidChemGraphException(name + ": " + e.getMessage());
        		}
        		if (graph == null) throw new InvalidChemGraphException(name);
        		ChemGraph cg = ChemGraph.make(graph);	
        		Species spe = Species.make(name, cg);
        		// GJB: Turn off reactivity if necessary, but don't let code turn it on
        		// again if was already set as unreactive from input file
        		if(IsReactive==false) spe.setReactivity(IsReactive);
        		speciesSet.put(name, spe);
        		line = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {

        	/*
        	 * 25Jun2009-MRH: When reading the Primary Kinetic Library, it IS NECESSARY 
        	 * 		to have a species file.  In the instance that no "species.txt" file exists,
        	 * 		inform user of this and terminate simulation.
        	 */
        	throw new IOException("RMG did not find/read species (species.txt) " +
        			"in the Primary Kinetic Library: " + name + "\n" + e.getMessage());
        }
    }

    public void readPdepReactions(String pdepFileName, String p_name) throws IOException { 
        try {
        	FileReader in = new FileReader(pdepFileName);
        	BufferedReader data = new BufferedReader(in);
        	
        	double A_multiplier = 1;
        	double E_multiplier = 1;
        	
        	String line = ChemParser.readMeaningfulLine(data);
        	if (line.startsWith("Unit")) {
        		line = ChemParser.readMeaningfulLine(data);
        		unit: while(!(line.startsWith("Reaction"))) {
        			if (line.startsWith("A")) {
        				StringTokenizer st = new StringTokenizer(line);
        				String temp = st.nextToken();
        				String unit = st.nextToken().trim();
        				if (unit.compareToIgnoreCase("mol/cm3/s") == 0) {
        					A_multiplier = 1;
        				}
        				else if (unit.compareToIgnoreCase("mol/liter/s") == 0) {
           					A_multiplier = 1e-3;
        				}
        				else if (unit.compareToIgnoreCase("molecule/cm3/s") == 0) {
        					A_multiplier = 6.022e23;
        				}
        			}
        			else if (line.startsWith("E")) {
        				StringTokenizer st = new StringTokenizer(line);
        				String temp = st.nextToken();
        				String unit = st.nextToken().trim();
        				if (unit.compareToIgnoreCase("kcal/mol") == 0) {
        					E_multiplier = 1;
        				}
        				else if (unit.compareToIgnoreCase("cal/mol") == 0) {
           					E_multiplier = 1e-3;
        				}
        				else if (unit.compareToIgnoreCase("kJ/mol") == 0) {
           					E_multiplier = 1/4.186;
        				}
        				else if (unit.compareToIgnoreCase("J/mol") == 0) {
           					E_multiplier = 1/4186;
        				}
        				else if (unit.compareToIgnoreCase("Kelvin") == 0) {
        					E_multiplier = 1.987e-3;
        				}
        			}
        			line = ChemParser.readMeaningfulLine(data);
        		}
        	}
            
        	String nextLine = ChemParser.readMeaningfulLine(data);
        	read: while (nextLine != null) {	
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(speciesSet, nextLine, A_multiplier, E_multiplier);
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
        		
        		/*
        		 * When reading in the auxillary information for the pdep reactions,
        		 * 	let's not assume the order is fixed (i.e. third-bodies and
        		 * 	their efficiencies, then lindemann, then troe)
        		 * The order of the if statement is important as the "troe" and
        		 * 	"low" lines will also contain a "/"; thus, the elseif contains
        		 * 	"/" needs to be last.
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
	        		} else if (nextLine.contains("/")) {
	        			// read in third body colliders + efficiencies
	            		thirdBodyList = ChemParser.parseThirdBodyList(nextLine);
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
        				tbr.setKineticsSource("Primary Kinetic Library: "+ p_name,0);
        				tbr.setKineticsComments(" ",0);
        				reactionSet.add(tbr);
                		Reaction reverse = tbr.getReverseReaction();
        				if (reverse != null) reactionSet.add(reverse);
        			} else {
        				// lindemann reaction
                		LindemannReaction tbr = LindemannReaction.make(r,thirdBodyList,low);
        				tbr.setKineticsSource("Primary Kinetic Library: "+ p_name,0);
        				tbr.setKineticsComments(" ",0);
        				reactionSet.add(tbr);
                		Reaction reverse = tbr.getReverseReaction();
        				if (reverse != null) reactionSet.add(reverse);
        			}
        		} else {
        			// troe reaction
            		TROEReaction tbr = TROEReaction.make(r,thirdBodyList, low, a, T3star, Tstar, troe7, T2star);
    				tbr.setKineticsSource("Primary Kinetic Library: "+ p_name,0);
    				tbr.setKineticsComments(" ",0);
    				reactionSet.add(tbr);
            		Reaction reverse = tbr.getReverseReaction();
    				if (reverse != null) reactionSet.add(reverse);
        		}
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {
//        	throw new IOException("Can't read reaction in primary kinetic library: troe reaction list.\n" + e.getMessage());
			/*
			 * 25Jun2009-MRH: When reading the Primary Kinetic Library, we should not require the user to supply
			 * 		troe reactions.  In the instance that no "troeReactions.txt" file exists, inform
			 * 		user of this but continue simulation.
			 */
        	System.out.println("RMG did not find/read pressure-dependent reactions (pdepreactions.txt) " +
        			"in the Primary Kinetic Library: " + name + "\n" + e.getMessage());
        }
    }
    
    public static int size() {
        return reactionSet.size();
    }
    
    public String getName() {
        return name;
    }
    
    public void setName(String p_name) {
        name = p_name;
    }
    
    public static LinkedHashSet getReactionSet() {
        return reactionSet;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\PrimaryKineticLibrary.java
*********************************************************************/


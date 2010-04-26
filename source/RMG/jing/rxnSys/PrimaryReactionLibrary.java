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
import jing.rxn.*;
import jing.chem.*;
import java.util.*;
import jing.chemUtil.*;
import jing.chemParser.*;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\PrimaryReactionLibrary.java                                                                  
//----------------------------------------------------------------------------

/**
This is the primary reaction set that any reaction system has to include into its model.  For example, in combustion system, we build a primary small molecule reaction set, and every combustion/oxydation system should include such a primary reaction library.  The reaction / rates are basically from Leeds methane oxidation mechanism.
*/
//## class PrimaryReactionLibrary 
public class PrimaryReactionLibrary {
    
    protected String name;		//## attribute name 
    
    protected static LinkedHashSet reactionSet = new LinkedHashSet();		//## attribute reactionSet 
    
    protected HashMap speciesSet = new HashMap();		//## attribute speciesSet 
    
    
    // Constructors
    
    public  PrimaryReactionLibrary(String p_libraryName, String p_directoryPath) throws IOException {
        name = p_libraryName;

        if ( p_directoryPath == null) throw new NullPointerException("PrimaryReactionLibrary directory path");
        try {
        	read(p_directoryPath);
        }
        catch (IOException e) {
        	throw new IOException("Error reading Primary Reaction Library: " + name + '\n' + e.getMessage());
        }
    }
    public  PrimaryReactionLibrary() {
    }
    
    public void appendPrimaryReactionLibrary(String new_p_libraryName, String p_directoryPath) throws IOException {
    	// Appends the current PRLib with an additional one, allowing the user
    	// to combine separate PRLibs easily.  GJB 10/05.
     	setName(name+"/"+new_p_libraryName);
    	try {
    		read(p_directoryPath);	
    	}
        catch (IOException e) {
        	throw new IOException("Error reading Primary Reaction Library: " + new_p_libraryName + '\n' + e.getMessage());
        }
    }
    
    
    //## operation getSpeciesSet() 
    public LinkedHashSet getSpeciesSet() {
        //#[ operation getSpeciesSet() 
        return new LinkedHashSet(speciesSet.values());
        //#]
    }
    
    public void read(String p_directoryName) throws IOException {
        try {
        	if (!p_directoryName.endsWith("/")) p_directoryName = p_directoryName + "/";
			System.out.println("Reading Primary Reaction Library from: "+p_directoryName);
        	
            String speciesFile = p_directoryName + "species.txt";
            String reactionFile = p_directoryName + "reactions.txt";
            String pdepreactionFile = p_directoryName + "pdepreactions.txt";
//            String thirdBodyReactionFile = p_directoryName + "3rdBodyReactions.txt";
//            String troeReactionsFile = p_directoryName + "troeReactions.txt"; 
//            String lindemannReactionsFile = p_directoryName + "lindemannReactions.txt";
        	readSpecies(speciesFile);
        	readReactions(reactionFile);
        	readPdepReactions(pdepreactionFile);
//            readThirdBodyReactions(thirdBodyReactionFile);
//            readTroeReactions(troeReactionsFile);
//            readLindemannReactions(lindemannReactionsFile);
        	return;
        }
        catch (Exception e) {
        	throw new IOException("Can't read primary reaction library.\n" + e.getMessage());
        }
    }
    
    public void readReactions(String p_reactionFileName) throws IOException {
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
            
        	LinkedHashSet currentPRLReactions = new LinkedHashSet();
        	line = ChemParser.readMeaningfulLine(data);
        	read: while (line != null) {
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(speciesSet, line, A_multiplier, E_multiplier);
        			r.setIsFromPrimaryReactionLibrary(true);
        			(r.getKinetics())[0].setFromPrimaryReactionLibrary(true);
        			// Changed source from "Seed Mechanism" to "Primary Reaction Library"
					r.setKineticsSource("Primary Reaction Library: "+ name,0);
					r.setKineticsComments(" ",0);
				}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(line + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(line);
        		
        		Iterator prlRxnIter = currentPRLReactions.iterator();
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
        			currentPRLReactions.add(r);
	        		Reaction reverse = r.getReverseReaction();
					
	        		if (reverse != null) {
						//reverse.getKinetics().setSource("Seed Mechanism: " + name);
						currentPRLReactions.add(reverse);
	        		}
        		}
        		line = ChemParser.readMeaningfulLine(data);
        	}
        	   
        	reactionSet.addAll(currentPRLReactions);
            in.close();
        	return;
        }
        catch (Exception e) {
//        	throw new IOException("Can't read reaction in primary reaction library.\n" + e.getMessage());
			/*
			 * 25Jun2009-MRH: When reading the Primary Reaction Library, we should not require the user to supply
			 * 		non pressure-dependent reactions.  In the instance that no "reactions.txt" file exists, inform
			 * 		user of this but continue simulation.
			 */
        	System.out.println("RMG did not find/read non pressure-dependent reactions (reaction.txt) " +
        			"in the Primary Reaction Library: " + name + "\n" + e.getMessage());
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
//        	throw new IOException("Can't read the species in primary reaction library: " + '\n' + e.getMessage());
        	/*
        	 * 25Jun2009-MRH: When reading the Primary Reaction Library, it IS NECESSARY 
        	 * 		to have a species file.  In the instance that no "species.txt" file exists,
        	 * 		inform user of this and terminate simulation.
        	 */
        	throw new IOException("RMG did not find/read species (species.txt) " +
        			"in the Primary Reaction Library: " + name + "\n" + e.getMessage());
        }
    }

    public void readThirdBodyReactions(String p_thirdBodyReactionFileName) throws IOException {
        try {
        	FileReader in = new FileReader(p_thirdBodyReactionFileName);
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
        			}
        			line = ChemParser.readMeaningfulLine(data);
        		}
        	}
            
        	String reactionLine = ChemParser.readMeaningfulLine(data);
        	read: while (reactionLine != null) {	
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(speciesSet, reactionLine, A_multiplier, E_multiplier);
        		}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(reactionLine + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(reactionLine);
        
        		String thirdBodyLine = ChemParser.readMeaningfulLine(data);
        		HashMap thirdBodyList = ChemParser.parseThirdBodyList(thirdBodyLine);
        		
        		ThirdBodyReaction tbr = ThirdBodyReaction.make(r,thirdBodyList);
        		tbr.setIsFromPrimaryReactionLibrary(true);
        		tbr.getKinetics()[0].setFromPrimaryReactionLibrary(true);
        		// Changed source from "Seed Mechanism" to "Primary Reaction Library"
				tbr.setKineticsSource("Primary Reaction Library: "+ name,0);
				tbr.setKineticsComments(" ",0);
        		reactionSet.add(tbr);
        		
        		Reaction reverse = tbr.getReverseReaction();
				
        		if (reverse != null) {
					//reverse.getKinetics().setSource("Seed Mechanism: "+ name);
					reactionSet.add(reverse);
        		}
        		
                /*
                 * MRH 23APR2010:
                 * Allowing RMG to handle general third-body reactions,
                 * 	e.g. A+M=B+C+M (where no specific third-body colliders are given)
                 * 
                 * The thirdBodyLine should either contain a list of third-bodies and their
                 * 	collision efficiencies (e.g. H2/2.0/ H2O/15.4/ CO/0.75/ etc.)
                 * OR
                 * the next reaction line (e.g. O+O+M=O2+M)
                 * 
                 * If the former, the thirdBodyList will not be empty and the next
                 * reactionLine string needs to be read from the BufferedReader data;
                 * if the latter, the thirdBodyList will be empty and the next
                 * reactionLine string has already been read (and is currently
                 * stored in the thirdBodyLine).
                 */
        		if (thirdBodyList.isEmpty()) reactionLine = thirdBodyLine;
        		else reactionLine = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {
//        	throw new IOException("Can't read reaction in primary reaction library.\n" + e.getMessage());
			/*
			 * 25Jun2009-MRH: When reading the Primary Reaction Library, we should not require the user to supply
			 * 		third body reactions.  In the instance that no "3rdBodyReactions.txt" file exists, inform
			 * 		user of this but continue simulation.
			 */
        	System.out.println("RMG did not find/read third body reactions (3rdBodyReactions.txt) " +
        			"in the Primary Reaction Library: " + name + "\n" + e.getMessage());
        }
    }
    
	   //## operation readTroeReactions(String) 
    public void readTroeReactions(String p_troeReactionFileName) throws IOException {
        //#[ operation readTroeReactions(String) 
        try {
        	FileReader in = new FileReader(p_troeReactionFileName);
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
        			}
        			line = ChemParser.readMeaningfulLine(data);
        		}
        	}
            
        	String reactionLine = ChemParser.readMeaningfulLine(data);
        	read: while (reactionLine != null) {	
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(speciesSet, reactionLine, A_multiplier, E_multiplier);
        		}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(reactionLine + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(reactionLine);
                
        		String thirdBodyLine = ChemParser.readMeaningfulLine(data);
        		HashMap thirdBodyList = ChemParser.parseThirdBodyList(thirdBodyLine);
        		
        		// parse the K at low limit
        		String lowLine = ChemParser.readMeaningfulLine(data);
        		StringTokenizer st = new StringTokenizer(lowLine, "/");
        		String temp = st.nextToken().trim();
        		String lowString = st.nextToken().trim();
        		/*
        		 * MRH 17Feb2010:
        		 * 	The units of the k_zero (LOW) Arrhenius parameters are different from the units of
        		 * 	k_inf Arrhenius parameters by a factor of cm3/mol, hence the getReactantNumber()+1
        		 */
        		ArrheniusKinetics low = ChemParser.parseSimpleArrheniusKinetics(lowString, A_multiplier, E_multiplier, r.getReactantNumber()+1);
        		
        		// parse Troe parameters
        		String troeLine = ChemParser.readMeaningfulLine(data);
        		st = new StringTokenizer(troeLine, "/");
        		temp = st.nextToken().trim();
        		String troeString = st.nextToken().trim();
        		st = new StringTokenizer(troeString);
                int n = st.countTokens();
                if (n != 3 && n != 4) throw new InvalidKineticsFormatException("Troe parameter number = "+n);
        
          		double a = Double.parseDouble(st.nextToken().trim());
        		double T3star = Double.parseDouble(st.nextToken().trim());
        		double Tstar = Double.parseDouble(st.nextToken().trim());
        		boolean troe7 = false;
        		double T2star = 0;
        		if (st.hasMoreTokens()) {
        			troe7 = true;
        			T2star = Double.parseDouble(st.nextToken().trim());
           		}
        		
        		TROEReaction tbr = TROEReaction.make(r,thirdBodyList, low, a, T3star, Tstar, troe7, T2star);
        		tbr.setIsFromPrimaryReactionLibrary(true);
        		tbr.getKinetics()[0].setFromPrimaryReactionLibrary(true);
        		// Changed source from "Seed Mechanism" to "Primary Reaction Library"
				tbr.setKineticsSource("Primary Reaction Library: "+ name,0);
				tbr.setKineticsComments(" ",0);
				
        		reactionSet.add(tbr);
        		Reaction reverse = tbr.getReverseReaction();
				
        		if (reverse != null) {
					//reverse.getKinetics().setSource("Seed Mechanism: "+ name);
					reactionSet.add(reverse);
        		}
        		
        		reactionLine = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {
//        	throw new IOException("Can't read reaction in primary reaction library: troe reaction list.\n" + e.getMessage());
			/*
			 * 25Jun2009-MRH: When reading the Primary Reaction Library, we should not require the user to supply
			 * 		troe reactions.  In the instance that no "troeReactions.txt" file exists, inform
			 * 		user of this but continue simulation.
			 */
        	System.out.println("RMG did not find/read troe reactions (troeReactions.txt) " +
        			"in the Primary Reaction Library: " + name + "\n" + e.getMessage());
        }
    }
    
    public void readLindemannReactions(String p_lindemannReactionFileName) throws IOException {
        try {
        	FileReader in = new FileReader(p_lindemannReactionFileName);
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
        			}
        			line = ChemParser.readMeaningfulLine(data);
        		}
        	}
            
        	String reactionLine = ChemParser.readMeaningfulLine(data);
        	read: while (reactionLine != null) {	
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(speciesSet, reactionLine, A_multiplier, E_multiplier);
        		}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(reactionLine + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(reactionLine);
                
        		String thirdBodyLine = ChemParser.readMeaningfulLine(data);
        		HashMap thirdBodyList = ChemParser.parseThirdBodyList(thirdBodyLine);
        		
        		// parse the K at low limit
        		String lowLine = ChemParser.readMeaningfulLine(data);
        		StringTokenizer st = new StringTokenizer(lowLine, "/");
        		String temp = st.nextToken().trim();
        		String lowString = st.nextToken().trim();
        		/*
        		 * MRH 17Feb2010:
        		 * 	The units of the k_zero (LOW) Arrhenius parameters are different from the units of
        		 * 	k_inf Arrhenius parameters by a factor of cm3/mol, hence the getReactantNumber()+1
        		 */
        		ArrheniusKinetics low = ChemParser.parseSimpleArrheniusKinetics(lowString, A_multiplier, E_multiplier, r.getReactantNumber()+1);        		
        		LindemannReaction lr = LindemannReaction.make(r,thirdBodyList, low);
				lr.setKineticsSource("Primary Reaction Library: "+ name,0);
				lr.setKineticsComments(" ",0);
				
        		reactionSet.add(lr);
        		Reaction reverse = lr.getReverseReaction();
				
        		if (reverse != null) {
					//reverse.getKinetics().setSource("Seed Mechanism: "+ name);
					reactionSet.add(reverse);
        		}
        		
        		reactionLine = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {
        	System.out.println("RMG did not read the following Primary Reaction Library file:"
        			+ p_lindemannReactionFileName);
        }
    }
    
    public void readPdepReactions(String pdepFileName) throws IOException { 
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
            
        	String reactionLine = ChemParser.readMeaningfulLine(data);
        	read: while (reactionLine != null) {	
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(speciesSet, reactionLine, A_multiplier, E_multiplier);
        		}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(reactionLine + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(reactionLine);
                
        		/*
        		 * Read the next line and determine what to do based on the
        		 * 	presence/absence of keywords
        		 */
        		String nextLine = ChemParser.readMeaningfulLine(data);
        		boolean continueToReadRxn = true;
        		
        		// Initialize all of the possible pdep variables
        		HashMap thirdBodyList = new HashMap();
        		ArrheniusKinetics low = new ArrheniusKinetics();
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
	            		String lowLine = ChemParser.readMeaningfulLine(data);
	            		StringTokenizer st = new StringTokenizer(lowLine, "/");
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
	            		String thirdBodyLine = ChemParser.readMeaningfulLine(data);
	            		thirdBodyList = ChemParser.parseThirdBodyList(thirdBodyLine);
	            		nextLine = ChemParser.readMeaningfulLine(data);
	        		} else {
	        			// the nextLine is a "new" reaction, hence we need to exit the while loop
	        			continueToReadRxn = false;
	        		}
        		}
        		   
        		// Make the "troe" reaction
        		TROEReaction tbr = TROEReaction.make(r,thirdBodyList, low, a, T3star, Tstar, troe7, T2star);
        		tbr.setIsFromPrimaryReactionLibrary(true);
        		tbr.getKinetics()[0].setFromPrimaryReactionLibrary(true);
        		// Changed source from "Seed Mechanism" to "Primary Reaction Library"
				tbr.setKineticsSource("Primary Reaction Library: "+ name,0);
				tbr.setKineticsComments(" ",0);
				
        		reactionSet.add(tbr);
        		Reaction reverse = tbr.getReverseReaction();
				
        		if (reverse != null) {
					//reverse.getKinetics().setSource("Seed Mechanism: "+ name);
					reactionSet.add(reverse);
        		}
        		
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {
//        	throw new IOException("Can't read reaction in primary reaction library: troe reaction list.\n" + e.getMessage());
			/*
			 * 25Jun2009-MRH: When reading the Primary Reaction Library, we should not require the user to supply
			 * 		troe reactions.  In the instance that no "troeReactions.txt" file exists, inform
			 * 		user of this but continue simulation.
			 */
        	System.out.println("RMG did not find/read pressure-dependent reactions (pdepReactions.txt) " +
        			"in the Primary Reaction Library: " + name + "\n" + e.getMessage());
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
	File Path	: RMG\RMG\jing\rxnSys\PrimaryReactionLibrary.java
*********************************************************************/


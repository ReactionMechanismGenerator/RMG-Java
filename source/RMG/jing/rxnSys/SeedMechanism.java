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

/**
 * This is a new class called SeedMechanism.  SeedMechanism is the same class
 * as the old PrimaryReactionLibrary, just with a new (and more appropriate)
 * name.  RMG will automatically include every species and reaction contained
 * in a Seed Mechanism.  Furthermore, the user has the option to pass multiple
 * Seed Mechanisms to RMG.  In the event of a duplicate species/reaction, RMG
 * will use the first instance it finds (i.e. the order of the Seed Mechanisms
 * listed in the condition.txt file is important).
 * 
 * MRH 9-Jun-2009
 */

/*
 * Comments from old PrimaryReactionLibrary:
 * 
 * This is the primary reaction set that any reaction system has to include 
 * into its model.  For example, in combustion system, we build a primary small 
 * molecule reaction set, and every combustion/oxidation system should include 
 * such a primary reaction library.  The reaction / rates are basically from Leeds 
 * methane oxidation mechanism.
 */

public class SeedMechanism {
    
    protected static String name;		//## attribute name 
    
    protected static LinkedHashSet reactionSet = new LinkedHashSet();		//## attribute reactionSet 
    
    protected static HashMap speciesSet = new HashMap();		//## attribute speciesSet 
    
    
    private boolean generateReactions = false;

    // Constructors
    
    public  SeedMechanism(String p_mechName, String p_directoryName, boolean p_generateReactions) throws IOException {
    	name = p_mechName;
		generateReactions = p_generateReactions;
        String dir = System.getProperty("RMG.workingDirectory");
        if (dir==null || p_directoryName == null) throw new NullPointerException("RMG " +
        		"does not recognize Seed Mechanism file name");
        try {
        	read(dir+"/databases/"+p_directoryName);
        }
        catch (IOException e) {
        	throw new IOException("Error in reading Seed Mechanism: " + name + '\n' + e.getMessage());
        }
    }
    
    public void appendSeedMechanism(String new_mechName, String new_directoryName, boolean p_generateReactions) throws IOException {
    	String dir = System.getProperty("RMG.workingDirectory");
     	if (p_generateReactions)
			setGenerateReactions(p_generateReactions);
		setName(name + "/" + new_mechName);
    	try {
    		read(dir + "/databases/" + new_directoryName);	
    	}
        catch (IOException e) {
        	throw new IOException("Error in reading Seed Mechanism: " + new_mechName + '\n' + e.getMessage());
        }
    }
    
    public LinkedHashSet getSpeciesSet() {
        return new LinkedHashSet(speciesSet.values());
    }
    
    public void read(String p_directoryName) throws IOException {
        System.out.println("Reading seed mechanism from directory " + p_directoryName);
		try {
        	if (!p_directoryName.endsWith("/")) p_directoryName = p_directoryName + "/";
        	
            String speciesFile = p_directoryName + "species.txt";
            String reactionFile = p_directoryName + "reactions.txt";
            String thirdBodyReactionFile = p_directoryName + "3rdBodyReactions.txt";
            String troeReactionsFile = p_directoryName + "troeReactions.txt";
            String lindemannReactionsFile = p_directoryName + "lindemannReactions.txt";
        	readSpecies(speciesFile);
        	readReactions(reactionFile);
            readThirdBodyReactions(thirdBodyReactionFile);
            readTroeReactions(troeReactionsFile);
            readLindemannReactions(lindemannReactionsFile);
        	return;
        }
        catch (Exception e) {
        	throw new IOException("RMG cannot read entire Seed Mechanism: " 
        			+ p_directoryName + "\n" + e.getMessage());
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
            
        	line = ChemParser.readMeaningfulLine(data);
        	read: while (line != null) {
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(speciesSet, line, A_multiplier, E_multiplier);
					r.setKineticsSource("Seed Mechanism: "+ name);
					r.setKineticsComments(" ");
				}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(line + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(line);
        		
        		Iterator prlRxnIter = reactionSet.iterator();
        		boolean foundRxn = false;
        		while (prlRxnIter.hasNext()) {
        			Reaction old = (Reaction)prlRxnIter.next();
        			if (old.equals(r)) {
        				old.addAdditionalKinetics(r.getKinetics(),1);
        				foundRxn = true;
        				break;
        			}
        		}
        		if (!foundRxn) {
        			reactionSet.add(r);
	        		Reaction reverse = r.getReverseReaction();
					
	        		if (reverse != null) {
						//reverse.getKinetics().setSource("Seed Mechanism: " + name);
						reactionSet.add(reverse);
	        		}
        		}
        		
        		line = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {
        	System.out.println("RMG did not read the following Seed Mechanism file: " 
        			+ p_reactionFileName);
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
        		speciesSet.put(name, spe);
        		line = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {
			throw new IOException("RMG cannot read the \"species.txt\" file in the Seed Mechanism\n" + e.getMessage());
        }
    }
    
    public static void readThirdBodyReactions(String p_thirdBodyReactionFileName) throws IOException {
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
				tbr.setKineticsSource("Seed Mechanism: "+ name);
				tbr.setKineticsComments(" ");
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
        	System.out.println("RMG did not read the following Seed Mechanism file:" 
        			+ p_thirdBodyReactionFileName);
        }
    }
     
    public static void readTroeReactions(String p_troeReactionFileName) throws IOException {
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
				tbr.setKineticsSource("Seed Mechanism: "+ name);
				tbr.setKineticsComments(" ");
				
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
        	System.out.println("RMG did not read the following Seed Mechanism file:" 
        			+ p_troeReactionFileName);
        }
    }
    
    public static void readLindemannReactions(String p_lindemannReactionFileName) throws IOException {
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
				lr.setKineticsSource("Seed Mechanism: "+ name);
				lr.setKineticsComments(" ");
				
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
        	System.out.println("RMG did not read the following Seed Mechanism file:" 
        			+ p_lindemannReactionFileName);
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
    
}


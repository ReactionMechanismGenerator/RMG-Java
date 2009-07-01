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
    
    //## operation PrimaryReactionLibrary(String,String) 
    public  PrimaryReactionLibrary(String p_libraryName, String p_directoryName) throws IOException {
        //#[ operation PrimaryReactionLibrary(String,String) 
        name = p_libraryName;
        String dir = System.getProperty("RMG.workingDirectory");
        if (dir==null || p_directoryName == null) throw new NullPointerException("PrimaryReactionLibrary file name");
        try {
        	read(dir+"/databases/"+p_directoryName);
        }
        catch (IOException e) {
        	throw new IOException("error in read primary library: " + name + '\n' + e.getMessage());
        }
        //#]
    }
    public  PrimaryReactionLibrary() {
    }
    
    //## operation appendPrimaryReactionLibrary(String, String) 
    public void appendPrimaryReactionLibrary(String new_p_libraryName, String new_p_directoryName) throws IOException {
    	//#[ operation appendPrimaryReactionLibrary(String, String) 

    	// Appends the current PRLib with an additional one, allowing the user
    	// to combine separate PRLibs easily.  GJB 10/05.
    	String dir = System.getProperty("RMG.workingDirectory");
     	setName(name+"/"+new_p_libraryName);
    	try {
    		read(dir+"/databases/"+new_p_directoryName);	
    	}
        catch (IOException e) {
        	throw new IOException("error in read primary library: " + new_p_libraryName + '\n' + e.getMessage());
        }
        //#]
    }
    
    
    //## operation getSpeciesSet() 
    public LinkedHashSet getSpeciesSet() {
        //#[ operation getSpeciesSet() 
        return new LinkedHashSet(speciesSet.values());
        //#]
    }
    
    //## operation read(String) 
    public void read(String p_directoryName) throws IOException {
        //#[ operation read(String) 
        try {
        	if (!p_directoryName.endsWith("/")) p_directoryName = p_directoryName + "/";
        	
            String speciesFile = p_directoryName + "species.txt";
            String reactionFile = p_directoryName + "reactions.txt";
            String thirdBodyReactionFile = p_directoryName + "3rdBodyReactions.txt";
            String troeReactions = p_directoryName + "troeReactions.txt"; 
        	readSpecies(speciesFile);
        	readReactions(reactionFile);
            readThirdBodyReactions(thirdBodyReactionFile);
            readTroeReactions(troeReactions);
        	return;
        }
        catch (Exception e) {
        	throw new IOException("Can't read primary reaction library.\n" + e.getMessage());
        }
        
        
        
        
        //#]
    }
    
    //## operation readReactions(String) 
    public void readReactions(String p_reactionFileName) throws IOException {
        //#[ operation readReactions(String) 
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
            
        	LinkedHashSet currentPRLReactions = new LinkedHashSet();
        	line = ChemParser.readMeaningfulLine(data);
        	read: while (line != null) {
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(speciesSet, line, A_multiplier, E_multiplier);
        			r.setIsFromPrimaryReactionLibrary(true);
        			r.getKinetics().setFromPrimaryReactionLibrary(true);
        			// Changed source from "Seed Mechanism" to "Primary Reaction Library"
					r.setKineticsSource("Primary Reaction Library: "+ name);
					r.setKineticsComments(" ");
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
        				old.addAdditionalKinetics(r.getKinetics(),1);
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
        
        
        
        
        //#]
    }
    
    //## operation readSpecies(String) 
    public void readSpecies(String p_speciesFileName) throws IOException {
        //#[ operation readSpecies(String) 
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
        
        
        
        
        //#]
    }
    
    //## operation readThirdBodyReactions(String) 
    public void readThirdBodyReactions(String p_thirdBodyReactionFileName) throws IOException {
        //#[ operation readThirdBodyReactions(String) 
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
        		tbr.getKinetics().setFromPrimaryReactionLibrary(true);
        		// Changed source from "Seed Mechanism" to "Primary Reaction Library"
				tbr.setKineticsSource("Primary Reaction Library: "+ name);
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
//        	throw new IOException("Can't read reaction in primary reaction library.\n" + e.getMessage());
			/*
			 * 25Jun2009-MRH: When reading the Primary Reaction Library, we should not require the user to supply
			 * 		third body reactions.  In the instance that no "3rdBodyReactions.txt" file exists, inform
			 * 		user of this but continue simulation.
			 */
        	System.out.println("RMG did not find/read third body reactions (3rdBodyReactions.txt) " +
        			"in the Primary Reaction Library: " + name + "\n" + e.getMessage());
        }
        
        
        
        
        //#]
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
        		ArrheniusKinetics low = ChemParser.parseSimpleArrheniusKinetics(lowString, A_multiplier, E_multiplier);
        		
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
        		tbr.getKinetics().setFromPrimaryReactionLibrary(true);
        		// Changed source from "Seed Mechanism" to "Primary Reaction Library"
				tbr.setKineticsSource("Primary Reaction Library: "+ name);
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
//        	throw new IOException("Can't read reaction in primary reaction library: troe reaction list.\n" + e.getMessage());
			/*
			 * 25Jun2009-MRH: When reading the Primary Reaction Library, we should not require the user to supply
			 * 		troe reactions.  In the instance that no "troeReactions.txt" file exists, inform
			 * 		user of this but continue simulation.
			 */
        	System.out.println("RMG did not find/read troe reactions (troeReactions.txt) " +
        			"in the Primary Reaction Library: " + name + "\n" + e.getMessage());
        }
        
        
        
        
        //#]
    }
	
    //## operation size() 
    public static int size() {
        //#[ operation size() 
        return reactionSet.size();
        //#]
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


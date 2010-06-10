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


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import jing.chem.ChemGraph;
import jing.chem.ForbiddenStructureException;
import jing.chem.InvalidChemGraphException;
import jing.chem.Species;
import jing.chemParser.ChemParser;
import jing.chemParser.InvalidReactionFormatException;
import jing.chemUtil.Graph;
import jing.mathTool.UncertainDouble;
import jing.rxnSys.PrimaryReactionLibrary;
import jing.rxnSys.SeedMechanism;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ReactionLibrary.java                                                                  
//----------------------------------------------------------------------------

/**
This is a super user-defined reaction library.  Here, "super" means that all the library reactions whose reactant(s)/product(s) appears in the reaction system will be included in the final reaction system model.  Also, the kinetics in this library has the higher priority than the kinetics in our default reaction/kinetics template.  So, be careful when make this ReactionLibrary.  we will include things in this library without any checking.  But, if there is an alternative kinetics we find in our kinetics template, we will output it to warn user.
*/
 
public class ReactionLibrary {
    
    private static ReactionLibrary INSTANCE = new ReactionLibrary();		//attribute INSTANCE 
    
    protected HashMap dictionary = new HashMap();
    
    protected HashSet library = new HashSet();
    
    protected String name;		//attribute name
    
    // Constructors
    
    /**
    Requires:
    Effects: this is the only constructor in this singleton reaction library.  it is protected, which means no user can construct it.  the construction should go through instance() method to check if there is only one instance of this class.
    Modifies: itsLibraryReaction
    */
    
    public  ReactionLibrary(String p_libraryName, String p_directoryPath) throws IOException {
    	// Called by readANDMakeReaction Library when it reads the first folder mentioned in condition file in Reaction Library 
    	// section
    	
    	// Picks the library name in condition file which needs to be attached to the library reactions
        name = p_libraryName;
        
        // Check if directory path is given/exists
        if ( p_directoryPath == null) throw new NullPointerException("New ReactionLibrary directory path");
        try {
        	read(p_directoryPath);
        }
        catch (IOException e) {
        	// Throws error if cant read in Reaction Library Path
        	throw new IOException("Error reading New Reaction Library: " + name + '\n' + e.getMessage());
        }
        
        
    }
    
     public  ReactionLibrary() {
      // This might be not called anywhere is this Redundant? CHECK 
     }
     
        
     public void appendReactionLibrary(String new_p_libraryName, String p_directoryPath) throws IOException {
     	// Appends the current Reaction Library with an additional one, allowing the user
     	// to combine separate Reaction Library easily. 
    	 
    	// Adds additional new name to Reaction Library  
      	setName(name+"/"+new_p_libraryName);
      	// Read in the directory path
     	try {
     		read(p_directoryPath);	
     	}
         catch (IOException e) {
         	throw new IOException("Error reading New Reaction Library: " + new_p_libraryName + '\n' + e.getMessage());
         }
     }
     
     public void read(String p_directoryName) throws IOException {
    	 // Reads in the directory and creates the Dictionary and Library
    	 // The Library is composed of Pdep and Non Pdep Reactions
    	 // shamel 6/10/2010 can handle only Lindemann, TROE and Third Body form
    	 
         try {
         	if (!p_directoryName.endsWith("/")) p_directoryName = p_directoryName + "/";
 			System.out.println("Reading New Reaction Library from: "+p_directoryName);
         	
 			 
             String dictionaryFile = p_directoryName + "species.txt";
             String libraryFile = p_directoryName + "reaction.txt";
             String pdeplibraryFile = p_directoryName + "pdepreaction.txt";
            // Read in Dictionary File, species list
             readDictionary(dictionaryFile);
            // Read in Non Pdep Reactions  
         	readLibrary(libraryFile);
         	// Read in Pdep Reactions
         	readPdepReactions(pdeplibraryFile);
         	return;
         }
         catch (Exception e) {
         	throw new IOException("Can't read New reaction library.\n" + e.getMessage());
         }
         
         
     }

          
   
    
    public void readPdepReactions(String pdepFileName) throws IOException {
    	// shamel As of 6/10/2010 This function is EXACTLY like the seed mechanism function
    	// Created by M.R.Harper and I am keeping  some of his comments below which will match with 
    	// seed mechanism function
    	
    	// To read in Pdep Reactions
    	try {
    		// Read in the file name
        	FileReader in = new FileReader(pdepFileName);
        	// Buffer the file
        	BufferedReader data = new BufferedReader(in);
        	
        	double A_multiplier = 1;
        	double E_multiplier = 1;
        	
        	// Here we want the Reaction Library to be in certain form -- as standard library form in RMG. The
        	// first line must begin with units followed by Reaction
        	
        	String line = ChemParser.readMeaningfulLine(data);
        	if (line.startsWith("Unit")) {
        		line = ChemParser.readMeaningfulLine(data);
        		unit: while(!(line.startsWith("Reaction"))) {
        			if (line.startsWith("A")) {
        				// To get units of Arrehnius Constant and get correct conversion multiplier
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
        				// To grab E and get correct conversion multiplier
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
        			// To grab the reaction from the pdeplibrary file
        			r = ChemParser.parseArrheniusReaction(dictionary, nextLine, A_multiplier, E_multiplier);
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
        				tbr.setKineticsSource("ReactionLibrary: " + name,0);
        				tbr.setKineticsComments(" ",0);
        				library.add(tbr);
                		Reaction reverse = tbr.getReverseReaction();
        				if (reverse != null) library.add(reverse);
        			} else {
        				// lindemann reaction
                		LindemannReaction tbr = LindemannReaction.make(r,thirdBodyList,low);
        				tbr.setKineticsSource("ReactionLibrary: " + name,0);
        				tbr.setKineticsComments(" ",0);
        				library.add(tbr);
                		Reaction reverse = tbr.getReverseReaction();
        				if (reverse != null) library.add(reverse);
        			}
        		} else {
        			// troe reaction
            		TROEReaction tbr = TROEReaction.make(r,thirdBodyList, low, a, T3star, Tstar, troe7, T2star);
    				tbr.setKineticsSource("ReactionLibrary: " + name,0);
    				tbr.setKineticsComments(" ",0);
    				library.add(tbr);
            		Reaction reverse = tbr.getReverseReaction();
    				if (reverse != null) library.add(reverse);
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
        	System.out.println("RMG did not find/read pressure-dependent reactions (PdepLibrary.txt) " +
        			"in the Reaction Library: "  + "\n" + e.getMessage());
        }

    	
    }
    
    
    
    
    public void readLibrary(String p_reactionFileName) throws IOException {
    	// shamel As of 6/10/2010 This function is EXACTLY like the seed mechanism function
    	// Created by M.R.Harper and I am keeping  some of his comments below which will match with 
    	// seed mechanism function
    	
    	// To read in Non Pdep Reactions
    	
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
        			r = ChemParser.parseArrheniusReaction(dictionary, line, A_multiplier, E_multiplier);
					r.setKineticsSource("ReactionLibrary: "+ name,0);
					r.setKineticsComments(" ",0);
				}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(line + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(line);
        		
        		Iterator prlRxnIter = library.iterator();
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
        			library.add(r);
	        		Reaction reverse = r.getReverseReaction();
					
	        		if (reverse != null) {
						//reverse.getKinetics().setSource("Seed Mechanism: " + name);
						library.add(reverse);
	        		}
        		}
        		
        		line = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {
        	System.out.println("RMG did not read the following New Reaction Library file: " 
        			+ p_reactionFileName);
        }
    }

    
    
    
    private HashMap readDictionary(String p_fileName) throws FileNotFoundException, IOException{
    	  //Function to read in dictionary file
    	  try{
    	    FileReader in = new FileReader(p_fileName);
    	    BufferedReader data = new BufferedReader(in);

    	    String line = ChemParser.readMeaningfulLine(data);

    	    read: while(line != null){
    	      StringTokenizer st = new StringTokenizer(line);
    	      String name = st.nextToken();

    	      data.mark(10000);
    	      line = ChemParser.readMeaningfulLine(data);
    	      if (line == null) break read;
    	      line = line.trim();
    	      data.reset();
    	      Graph graph = null;

    	        graph = ChemParser.readChemGraph(data);
    	        if (graph == null) throw new InvalidChemGraphException(name);
        		ChemGraph cg ;
        		Species spe  ;
				try {
					cg = ChemGraph.make(graph);
					spe = Species.make(name, cg);
					
					Object old = dictionary.get(name);
		    	      if (old == null){
		    	        dictionary.put(name, spe);
		    	      }
		    	      else{
		    	        Species oldGraph = (Species)old;
		    	        if (!oldGraph.getChemGraph().getGraph().equals(graph)) {
		    	          System.out.println("Can't replace graph in Reaction Library!");
		    	          System.exit(0);
		    	        }
		    	    }
					
				} catch (InvalidChemGraphException e) {
					System.err.println("The graph \n" + graph.toString() + " is not valid");
					e.printStackTrace();
				} catch (ForbiddenStructureException e) {
					System.err.println("The graph \n" + graph.toString() + " is a forbidden structure");
					e.printStackTrace();
				}	
        		
    	      
    	      line = ChemParser.readMeaningfulLine(data);
    	    }
    	    in.close();
    	    return dictionary;
    	  }
    	  catch (FileNotFoundException e){
    	      throw new FileNotFoundException(p_fileName);
    	  }
    	  catch (IOException e){
    	    throw new IOException(p_fileName + ":" + e.getMessage());
    	  }
    	  
    	}
    
 
    // clearLibraryReaction() 
    public void clearLibraryReaction() {
     
        library.clear();
         }
    
        // getLibraryReaction() 
    public Iterator getLibraryReaction() {
        
        Iterator iter=library.iterator();
        return iter;
        
    }
    
        /**
    Requires:
    Effects: check if all the library reaction in this reaction library are valid
    Modifies:
    */
    
    public boolean repOk() {
    
        Iterator iter = getLibraryReaction();
        Reaction reaction;
        while (iter.hasNext()) {
        	reaction = (Reaction)iter.next();
        	if (!reaction.repOk()) return false;
        }
        
        return true;
        //#]
    }

    /**
    Requires:
    Effects: check if this reaction library is empty.  if it is, return true; otherwise, return false;
    Modifies:
    */
    
    public boolean isEmpty() {
        
        return (size() == 0);
        
        
        
        
    }
    
    
    public void removeLibraryReaction(LibraryReaction p_LibraryReaction) {
        
        library.remove(p_LibraryReaction);
        
    }
    
 
    /**
    Requires:
    Effects: return the size of itsLibraryReaction.
    Modifies:
    */
    
    public int size() {
    
        return library.size();
    
    }
    
    // shamel Asof 6/10/2010 this method is not needed
    //public static ReactionLibrary getINSTANCE() {
    //    return INSTANCE;
   // }
    
   
    public String getName() {
        return name;
    }
    
    public void setName(String p_name) {
        name = p_name;
    }
    
    
}    

/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ReactionLibrary.java
*********************************************************************/


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

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ReactionLibrary.java                                                                  
//----------------------------------------------------------------------------

/**
This is a super user-defined reaction library.  Here, "super" means that all the library reactions whose reactant(s)/product(s) appears in the reaction system will be included in the final reaction system model.  Also, the kinetics in this library has the higher priority than the kinetics in our default reaction/kinetics template.  So, be careful when make this ReactionLibrary.  we will include things in this library without any checking.  But, if there is an alternative kinetics we find in our kinetics template, we will output it to warn user.
*/
//## class ReactionLibrary 
public class ReactionLibrary {
    
    private static ReactionLibrary INSTANCE = new ReactionLibrary();		//## attribute INSTANCE 
    
    protected HashMap dictionary;
    
    protected HashSet library;
    
    // Constructors
    
    /**
    Requires:
    Effects: this is the only constructor in this singleton reaction library.  it is protected, which means no user can construct it.  the construction should go through instance() method to check if there is only one instance of this class.
    Modifies: itsLibraryReaction
    */
    //## operation ReactionLibrary() 
    private  ReactionLibrary() {
        {
            library=new HashSet();
            dictionary = new HashMap();
        }
        String directory = System.getProperty("jing.rxn.ReactionLibrary.pathName");
        String dictionaryFile = directory + "/Dictionary.txt";
        String libraryFile = directory + "/Library.txt";
			
			System.out.println("\nReading reaction library from "+directory);
        try{
            read(dictionaryFile, libraryFile);
          }
          catch (IOException e){
            System.out.println("Can't read reaction library files!");
        }
    }
    
    
    private void read(String p_dictionary, String p_library) throws IOException, FileNotFoundException {
    	  //#[ operation read(String)
    	    dictionary = readDictionary(p_dictionary);
    	    readLibrary(p_library);
    	    //#]
    }
    
    private void readLibrary(String p_reactionFileName) throws IOException {
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
            
        	line = ChemParser.readMeaningfulLine(data);
        	read: while (line != null) {
        		Reaction r;
        		try {
        			r = ChemParser.parseArrheniusReaction(dictionary, line, A_multiplier, E_multiplier);
					r.kinetics[0].setComments(" ");
					r.kinetics[0].setSource("ReactionLibrary");
				}
        		catch (InvalidReactionFormatException e) {
        			throw new InvalidReactionFormatException(line + ": " + e.getMessage());
        		}
        		if (r == null) throw new InvalidReactionFormatException(line);
        		
        		HashSet reactants = new HashSet();
        		reactants.addAll(r.getReactantList());
        		LibraryReaction fLR = LibraryReaction.makeLibraryReaction(r);
        		
        		library.add(fLR);
        		
        		Reaction reverse = r.getReverseReaction();
				LibraryReaction rLR = LibraryReaction.makeLibraryReaction(reverse);
				
        		if (rLR != null) {
        			reactants.addAll(reverse.getReactantList());
					library.add(rLR);
					fLR.setReverseReaction(rLR);
					rLR.setReverseReaction(fLR);						
			    }
        		
        		line = ChemParser.readMeaningfulLine(data);
        	}
        	   
            in.close();
        	return;
        }
        catch (Exception e) {
        	throw new IOException("Can't read reaction in primary reaction library.\n" + e.getMessage());
        }
        
    }
    
    
    
    private HashMap readDictionary(String p_fileName) throws FileNotFoundException, IOException{
    	  //#[ operation readDictionary(String)
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
		    	        if (!oldGraph.equals(graph)) {
		    	          System.out.println("Can't replace graph in Reaction Library!");
		    	          System.exit(0);
		    	        }
		    	    }
					
				} catch (InvalidChemGraphException e) {
					System.err.println("The graph \n" + graph.toString() + " is not valid");
					e.printStackTrace();
				} catch (ForbiddenStructureException e) {
					System.err.println("The graph \n" + graph.toString() + " is not a forbidden structure");
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
    	  //#]
    	}
    
 
    //## operation clearLibraryReaction() 
    public void clearLibraryReaction() {
        //#[ operation clearLibraryReaction() 
        library.clear();
        //#]
    }
    
        //## operation getLibraryReaction() 
    public Iterator getLibraryReaction() {
        //#[ operation getLibraryReaction() 
        Iterator iter=library.iterator();
        return iter;
        //#]
    }
    
        /**
    Requires:
    Effects: check if all the library reaction in this reaction library are valid
    Modifies:
    */
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
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
    //## operation isEmpty() 
    public boolean isEmpty() {
        //#[ operation isEmpty() 
        return (size() == 0);
        
        
        
        //#]
    }
    
    //## operation removeLibraryReaction(LibraryReaction) 
    public void removeLibraryReaction(LibraryReaction p_LibraryReaction) {
        //#[ operation removeLibraryReaction(LibraryReaction) 
        library.remove(p_LibraryReaction);
        //#]
    }
    
 
    /**
    Requires:
    Effects: return the size of itsLibraryReaction.
    Modifies:
    */
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return library.size();
        //#]
    }
    
    protected static ReactionLibrary getINSTANCE() {
        return INSTANCE;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ReactionLibrary.java
*********************************************************************/


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
import jing.rxnSys.PrimaryKineticLibrary;
import jing.rxnSys.SeedMechanism;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ReactionLibrary.java                                                                  
//----------------------------------------------------------------------------

/**
This is a super user-defined reaction library.  Here, "super" means that all the library reactions whose reactant(s)/product(s) appears in the reaction system will be included in the final reaction system model.  Also, the kinetics in this library has the higher priority than the kinetics in our default reaction/kinetics template.  So, be careful when make this ReactionLibrary.  we will include things in this library without any checking.  But, if there is an alternative kinetics we find in our kinetics template, we will output it to warn user.
*/
 
public class ReactionLibrary {
    
    private static ReactionLibrary INSTANCE = new ReactionLibrary(); 
    
    protected HashMap dictionary = new HashMap();
    
    protected LinkedHashSet library = new LinkedHashSet();
    
    protected String name;
    
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
        if ( p_directoryPath == null) throw new NullPointerException("RMG cannot locate Reaction Library: directory path is null");
        try {
        	read(p_directoryPath,p_libraryName);
        }
        catch (IOException e) {
        	// Throws error if cant read in Reaction Library Path
        	throw new IOException("Error reading Reaction Library: " + name + '\n' + e.getMessage());
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
     		read(p_directoryPath,new_p_libraryName);	
     	}
         catch (IOException e) {
         	throw new IOException("Error reading Reaction Library: " + new_p_libraryName + '\n' + e.getMessage());
         }
     }
     
     public void read(String p_directoryName, String p_name) throws IOException {
    	 // Reads in the directory and creates the Dictionary and Library
    	 // The Library is composed of Pdep and Non Pdep Reactions
    	 // shamel 6/10/2010 can handle only Lindemann, TROE and Third Body form
    	 
    	 // MRH making some changes to remove duplicate code in 
    	 //	PrimaryKineticsLibrary / ReactionLibrary 
    	 
         try {
         	if (!p_directoryName.endsWith("/")) p_directoryName = p_directoryName + "/";
 			System.out.println("Reading Reaction Library from: "+p_directoryName);
         	
 			 
             String dictionaryFile = p_directoryName + "species.txt";
             String libraryFile = p_directoryName + "reactions.txt";
             String pdeplibraryFile = p_directoryName + "pdepreactions.txt";
             
             SeedMechanism sm = new SeedMechanism();
             dictionary.putAll(sm.readSpecies(dictionaryFile,p_name,"ReactionLibrary: "));
             library.addAll(sm.readReactions(libraryFile,p_name,dictionary,"ReactionLibrary: ",false));
             library.addAll(sm.readPdepReactions(pdeplibraryFile,p_name,dictionary,"ReactionLibrary: ",false));
            
             return;
         }
         catch (Exception e) {
         	throw new IOException("Can't read Reaction library.\n" + e.getMessage());
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


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


import jing.chem.*;

import java.util.*;

import jing.chem.Species;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\LibraryReactionGenerator.java                                                                  
//----------------------------------------------------------------------------

/**
Generate reaction system from a list of species according to the reaction library.  i.e., find out matched reactions in reaction library. 
This ReactionLibraryGenerator should be a singleton since if there are two ReactionLibraryGenerators, we can't maintain the consistency.
*/
//## class LibraryReactionGenerator 
public class LibraryReactionGenerator implements ReactionGenerator {
    
    protected ReactionLibrary reactionLibrary;
    
    
    public LibraryReactionGenerator(ReactionLibrary p_reactionlibrary) {
    	reactionLibrary = p_reactionlibrary;
    }
    
    /**
    This is the overridden operation of the reactable interface.   Here, we will define the meaning of "match" the pass-in species to the library reactions.  
    some possible meanings:
    (1) if any of the reactant of a Library reaction appears in the species list, it is considered matched;
    (2) if all the reactants of a Library reaction appear in the species list, it is considered matched.
    Now, we will use the second definition.
    
    Why not the First ? As the other reactants might neither be  edge nor core species, hence that reaction is not considered for time being
    **/
    
    
    private boolean match_reactants(LinkedHashSet p_speciesSeed, Reaction p_libraryReaction) {
    	/*
    	 * This function will run through the Species Set and find "true" for those reactions in which all the reactants or all the  products 
    	 * are present in the Species Set  
    	 */
    
    	// Iterator to iterate over all reactants of the current library reaction
    	Iterator iter_reactants = p_libraryReaction.getReactants();
        
       int count_reactants =0; // Counter for No of Reactants in Reaction
       int count_reactantcontainedinseed =0; // Counter for No of Reactants in Reaction Present in Species Set
       
        while(iter_reactants.hasNext()){
        	
        	// Cast the reactants into RMG species
        	Species spe=(Species)iter_reactants.next();
        	
        	count_reactants =count_reactants +1;
        	
        	// Check if species is present in Species Set 
        	if(p_speciesSeed.contains(spe)){
        		// shamel: 6/11/2010 Line for me to debug
        		//System.out.println("Found a Species Match"+spe);
        		count_reactantcontainedinseed = count_reactantcontainedinseed +1;
        	}
        	
        }
        if(count_reactants == count_reactantcontainedinseed){
        	// If all the reactants in the reaction are present in species set return true
    		return true;
    	}
     // If any of the reactants in the reaction are not present in species set return false
        return false;
    }
   
     
    
    /*  
     // This portion of code was written to match chemgraphs of individual species instead of using the .contains(spe),might be useful 
     // later on. More rigorous than using .contains() 
   	    int i =0 ;
        while (iter_reactants.hasNext()) {
        	Species spe = (Species)iter_reactants.next();
        	
        	Iterator iter_speciesSeed = p_speciesSeed.iterator();
        	
        	while(iter_speciesSeed.hasNext()){
        		Species coreSpecies_indiv = (Species)iter_speciesSeed.next();
        		
        		if (coreSpecies_indiv.equals(spe)){
            		System.out.println("Found a Species Match"+spe);
            		return true;
            		
            	}	
        	}
        }
        
        return false;    		
    }
    
     */
        
    /**
    Search all the reaction in the itsReactionLibrary.  If a library reaction has a species in the pass-in SpeciesList, add this reaction into present ReactionSystem.  Return the final ReactionSystem, when the search for checking all the library reactions are done.
    */ 
    /**
    the species involved in reaction system.
    */

    public LinkedHashSet react(LinkedHashSet p_speciesSeed) {
        
    	LinkedHashSet reaction_set = new LinkedHashSet();
        
    	// Check whether Reaction Library is OK or empty or the Species Set is Empty
    	
        if (reactionLibrary == null || !reactionLibrary.repOk() || reactionLibrary.isEmpty() || p_speciesSeed.size()==0) {
        	return reaction_set;
        }
        
        // Algorithm for search reaction library to find out proper library reactions                                                  
        Reaction current_reaction; 
        
        Iterator iter = reactionLibrary.getLibraryReaction();
        
        // Run through all the reactions in current instance of Reaction Library
        while (iter.hasNext()) {
        	
        	// Cast it into a  Reaction ( i.e pick the reaction )
        	// As Reaction type encompasses all types i.e Library Reaction, Template Reaction 
        	current_reaction = (Reaction)iter.next();
        	
        	// check if the reaction in reaction library is a library reaction, if it is not, return an empty reaction system
        	
        	/*
        	 *  We only need to search through the reactants as the ReactionLibrary returns both forward and backward reactions 
        	 *  for reversible reactions as two separate reactions in the Set.   
        	 */
        	
        	if (match_reactants(p_speciesSeed, current_reaction))
        		reaction_set.add(current_reaction);	
        }
        
        // add in reaction set into reaction system
        return reaction_set;
        
    }
    
     
    public LinkedHashSet react(LinkedHashSet p_speciesSet, Species p_species, String p_rxnFamilyName) {
               if (!p_speciesSet.contains(p_species)) p_speciesSet.add(p_species);
        LinkedHashSet species = (LinkedHashSet)p_speciesSet.clone();
        
        return react(species);
    }
    
    public LinkedHashSet generatePdepReactions(Species p_species){
    	LinkedHashSet speciesSet = new LinkedHashSet();
    	speciesSet.add(p_species);
    	LinkedHashSet reactionSet = react(speciesSet); 
    	// Iterate through reaction set and for a lindemann / TROE / Third body reaction print out
    	// Reaction and a warning saying only using the high pressure limit rate 
    	Iterator iter_reaction = reactionSet.iterator();
    	while (iter_reaction.hasNext()){
    		Reaction rxn =(Reaction)iter_reaction.next();
    		if (rxn instanceof ThirdBodyReaction || rxn instanceof TROEReaction || rxn instanceof LindemannReaction)
    			System.out.println("RMG is only utilizing the high-pressure limit parameters for Reaction Library, reaction: " + rxn.toString());
    	}
        return reactionSet;
    }
    
    public ReactionLibrary getReactionLibrary() {
    	return reactionLibrary;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\LibraryReactionGenerator.java
*********************************************************************/


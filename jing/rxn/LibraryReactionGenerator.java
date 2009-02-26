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
    
    // Constructors
    
    //## operation LibraryReactionGenerator() 
    public  LibraryReactionGenerator() {
        //#[ operation LibraryReactionGenerator() 
        reactionLibrary = ReactionLibrary.getINSTANCE();
        //#]
    }
    
    /**
    This is the overridden operation of the reactable interface.   Here, we will define the meaning of "match" the pass-in species to the library reactions.  
    some possible meanings:
    (1) if any of the reactant of a Library reaction appears in the species list, it is considered matched;
    (2) if both reactants of a Library reaction appear in the species list, it is considered matched;
    (3) if all the reactants and products of a Library reaction appear in the species list, it is considered matched.
    Now, we use the first definition.
    
    Requires:
    Effects:
    Modifies:
    */
    // Argument HashSetp_speciesSeed : 
    /**
    pass-in species as key to search the proper reaction.
    */
    // Argument LibraryReactionp_libraryReaction : 
    /**
    the pass-in structure of a library reaction.
    */
    //## operation match(HashSet,LibraryReaction) 
    private boolean match(LinkedHashSet p_speciesSeed, LibraryReaction p_libraryReaction) {
        //#[ operation match(HashSet,LibraryReaction) 
        Iterator iter = p_libraryReaction.getReactants();
        
        while (iter.hasNext()) {
        	Species spe = (Species)iter.next();
        	if (p_speciesSeed.contains(spe)) return true;
        }
        
        return false;
        
        		
        
        
        
        //#]
    }
    
    /**
    Search all the reaction in the itsReactionLibrary.  If a library reaction has a species in the pass-in SpeciesList, add this reaction into present ReactionSystem.  Return the final ReactionSystem, when the search for checking all the library reactions are done.
    */
    // Argument HashSetp_speciesSeed : 
    /**
    the species involved in reaction system.
    */
    //## operation react(HashSet) 
    public LinkedHashSet react(LinkedHashSet p_speciesSeed) {
        //#[ operation react(HashSet) 
    	LinkedHashSet reaction_set = new LinkedHashSet();
        
        if (!reactionLibrary.repOk() || reactionLibrary.isEmpty() || p_speciesSeed.size()==0) {
        	return reaction_set;
        }
        
        // add algorithm for search reaction library to find out proper library reactions                                                  
        LibraryReaction current_reaction;
        Iterator iter = ReactionLibrary.getINSTANCE().getLibraryReaction();
        while (iter.hasNext()) {
        	current_reaction = (LibraryReaction)iter.next();
        	// check if the reaction in reaction library is a library reaction, if it is not, return an empty reaction system
        	if (match(p_speciesSeed, current_reaction)) reaction_set.add(current_reaction);
        }
        // add in reaction set into reaction system
        return reaction_set;
        
        
        //#]
    }
    
    
    //## operation react(HashSet,Species) 
    public LinkedHashSet react(LinkedHashSet p_speciesSet, Species p_species) {
        //#[ operation react(HashSet,Species) 
        if (!p_speciesSet.contains(p_species)) p_speciesSet.add(p_species);
        LinkedHashSet species = (LinkedHashSet)p_speciesSet.clone();
        
        return react(species);
        //#]
    }
    
    public LinkedHashSet generatePdepReactions(Species p_species){
    	LinkedHashSet speciesSet = new LinkedHashSet();
    	speciesSet.add(p_species);
    	LinkedHashSet reactionSet = react(speciesSet);
    	
        return reactionSet;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\LibraryReactionGenerator.java
*********************************************************************/


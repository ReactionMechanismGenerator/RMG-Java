//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator                                            
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!     
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!  
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED 
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR  
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;  
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//! 
//!******************************************************************************



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
    
    public void generatePdepReactions(Species p_species){
    	LinkedHashSet speciesSet = new LinkedHashSet();
    	speciesSet.add(p_species);
    	LinkedHashSet reactionSet = react(speciesSet);
    	
        p_species.addPdepPaths(reactionSet);
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\LibraryReactionGenerator.java
*********************************************************************/


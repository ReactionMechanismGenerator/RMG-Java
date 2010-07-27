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
    			System.out.println("RMG is only using the high-pressure limit parameters from Reaction Library, reaction: " + rxn.toString());
    	}
        return reactionSet;
    }
    
    public ReactionLibrary getReactionLibrary() {
    	return reactionLibrary;
    }

    public LinkedHashSet RemoveDuplicateReac(LinkedHashSet reaction_set){
    	
   	 // Get the reactants and products of a reaction and check with other reaction if both reactants and products
   	 // match - delete duplicate entry, give preference to Seed Mechanism > Reaction Library >  Reaction Template 
   	 // this information might be available from the comments 
   	
   	LinkedHashSet newreaction_set = new LinkedHashSet();
   	
   	Iterator iter_reaction =reaction_set.iterator();
   	
   	Reaction current_reaction;
   	
   	while(iter_reaction.hasNext()){
   		// Cast it into a  Reaction ( i.e pick the reaction )
       	current_reaction = (Reaction)iter_reaction.next();
       	
       	// To remove current reaction from reaction_set
       	reaction_set.remove(current_reaction);
       	
       	// Match Current Reaction with the reaction set and if a duplicate reaction is found remove that reaction 
              LinkedHashSet dupreaction_set = dupreaction(reaction_set,current_reaction);
           // Remove the duplicate reaction from reaction set
              reaction_set.removeAll(dupreaction_set);
           
           // If duplicate reaction set was not empty 
              if(!dupreaction_set.isEmpty()){
 
           // Add current reaction to duplicate set and from among this choose reaction according to
           // following hierarchy Seed > Reaction Library > Template. Add that reaction to the newreaction_set
              
           // Add current_reaction to duplicate set 
              dupreaction_set.add(current_reaction);
           
           // Get Reaction according to hierarchy
              LinkedHashSet reaction_toadd = reaction_add(dupreaction_set);
              
           // Add all the Reactions to be kept to new_reaction set     
              newreaction_set.addAll(reaction_toadd);
              }
              else{
           	   // If no duplicate reaction was found add the current reaction to the newreaction set
           	   newreaction_set.add(current_reaction);
              }
              
              
           // Need to change iterate over counter here 
              iter_reaction =reaction_set.iterator();
       	}
   	return newreaction_set;
   }
  
   
    public LinkedHashSet dupreaction(LinkedHashSet reaction_set, Reaction test_reaction){
    	// Iterate over the reaction set and find if a duplicate reaction exist for the the test reaction 
 
    	LinkedHashSet dupreaction_set = new LinkedHashSet();	
    	
    	Iterator iter_reaction =reaction_set.iterator();
    	
    	Reaction current_reaction;
    	
    	// we will test if reaction are equal by structure test here, structure dosent require kinetics
    	
    	// Get Structure of test reaction
    	Structure test_reactionstructure = test_reaction.getStructure();
    	
    	// Get reverse structure of test reaction
    	Structure test_reactionstructure_rev = test_reactionstructure.generateReverseStructure();
    	
    	   	    	
    	while(iter_reaction.hasNext()){
    		// Cast it into a  Reaction ( i.e pick the reaction )
        	current_reaction = (Reaction)iter_reaction.next();
        	
        	// Get Structure of current reaction to be tested for equality to test reaction
        	Structure current_reactionstructure = current_reaction.getStructure();
        	
        	// Check if Current Reaction Structure matches the Fwd Structure of Test Reaction
        	if(current_reactionstructure.equals(test_reactionstructure)){
        		dupreaction_set.add(current_reaction);
        	}
        	
        	// Check if Current Reaction Structure matches the Reverse Structure of Test Reaction
        	if(current_reactionstructure.equals(test_reactionstructure_rev)){
        		dupreaction_set.add(current_reaction);
        	}
        	
        	        		
    	}
    	
    	// Print out the dupreaction set if not empty
    	if(!dupreaction_set.isEmpty()){
    		// shamel 07-23-2010 closed annoying printing of duplicate reactions found
    	//System.out.println("dupreaction_set (Generated By RMG and also Present in Reaction Library, keeping the one from RL)" + dupreaction_set);
    	}
    	// Return the duplicate reaction set
    	return dupreaction_set;
    }

    public LinkedHashSet reaction_add(LinkedHashSet reaction_set){
   	
   	Reaction current_reaction;
   	
   	Iterator iter_reaction = reaction_set.iterator();
   	
   	LinkedHashSet reaction_seedset = new LinkedHashSet();
   	
   	LinkedHashSet reaction_rlset = new LinkedHashSet();
   	
   	LinkedHashSet reaction_trset = new LinkedHashSet();
   	
   	
   	while(iter_reaction.hasNext()){
   		// Cast it into a  Reaction ( i.e pick the reaction )
       	current_reaction = (Reaction)iter_reaction.next();
       	
       	// As I cant call the instance test as I have casted my reaction as a Reaction 
       	// I will use the source (comments) to find whether a reaction is from Seed Mechanism
       	// Reaction Library or Template Reaction
       	
       	String source = current_reaction.getKineticsSource(0);
       	//System.out.println("Source"+source);
       	
       	if (source == null){
       		// If source is null I am assuming that its not a Reaction from Reaction Library or Seed Mechanism
       		source = "TemplateReaction:";
       	}
       	
       	// To grab the First word from the source of the comment
       	// As we have Reaction_Type:, we will use : as our string tokenizer
       	StringTokenizer st = new StringTokenizer(source,":");
       	String reaction_type = st.nextToken();
       	
       	// shamel: Cant think of more elegant way for now
       	// Splitting the set into Reactions from Seed Mechanism/Reaction Library and otherwise Template Reaction
       	if (reaction_type.equals( "SeedMechanism")){
       		// Add to seed mechanism set
       		reaction_seedset.add(current_reaction);
       	}        	
       	else if (reaction_type.equals("ReactionLibrary") ){
       		// Add to reaction library set
       		reaction_rlset.add(current_reaction);
       	}
       	else{
       		// Add to template reaction set
       		reaction_trset.add(current_reaction);
       	}
       		
       	
       	
   	}
   	 if(!reaction_seedset.isEmpty()){
   		 // shamel: 6/10/2010 Debug lines
   		 //System.out.println("Reaction Set Being Returned"+reaction_seedset);
   		 return reaction_seedset;
   	 }
   	 else if(!reaction_rlset.isEmpty()){
   		 //System.out.println("Reaction Set Being Returned in RatebasedRME"+reaction_rlset);
   		 return reaction_rlset;
   	 }
   	 else{
   		 //System.out.println("Reaction Set Being Returned"+reaction_trset); 
   		return reaction_trset; 
   	 }
   	
   }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\LibraryReactionGenerator.java
*********************************************************************/


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
import jing.param.Global;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\TemplateReactionGenerator.java                                                                  
//----------------------------------------------------------------------------

//## class TemplateReactionGenerator 
public class TemplateReactionGenerator implements ReactionGenerator {
    
    protected ReactionTemplateLibrary reactionTemplateLibrary;
    
    // Constructors
    
    //## operation TemplateReactionGenerator() 
    public  TemplateReactionGenerator() {
        //#[ operation TemplateReactionGenerator() 
        reactionTemplateLibrary = ReactionTemplateLibrary.getINSTANCE();
        
        //#]
    }
    
    /**
    Requires:
    Effects: Check all the pass-in species to see if there is a match between the species and any reaction template's reactant's pattern, if there is a match, generate the new reaction, find out the kinetics from first kinetics library, then kineticstemplatelibrary, then add the new reaction 
    into the reaction set.  keep doing this until there is no more new match.  return the final reaction set.
    Modifies: 
    */
    // Argument HashSetp_speciesSeed : 
    /**
    the set of species in the reaction system.
    */
    //## operation react(HashSet) 
    public HashSet react(HashSet p_speciesSeed) {
        //#[ operation react(HashSet) 
        if (p_speciesSeed.size() == 0) {
        	return null;
        }
        
        HashSet reaction_set = new HashSet();
        HashSet species_set = new HashSet();
        
        // add here the algorithm to generate reaction
        // loop over all the reaction template to find any possible match between the species seed set and the reaction template library
        Iterator template_iter = reactionTemplateLibrary.getReactionTemplate();
        while (template_iter.hasNext()) {
        	ReactionTemplate current_template = (ReactionTemplate)template_iter.next();
        	// the reaction template has only one reactant, we only need to loop over the whole species seed set to find a match
        	if (current_template.hasOneReactant()) {
        		Iterator species_iter1 = p_speciesSeed.iterator();
        		while (species_iter1.hasNext()) {
        			Species first_reactant = (Species)species_iter1.next();
        			HashSet current_reactions = current_template.reactOneReactant(first_reactant);
        			reaction_set.addAll(current_reactions);
        		}                                
        	}
        	// the reaction template has two reactants, we need to check all the possible combination of two species
        	else if (current_template.hasTwoReactants()) {
				
        		Iterator species_iter1 = p_speciesSeed.iterator();
        		while (species_iter1.hasNext()) {
        			Species first_reactant = (Species)species_iter1.next();
        			Iterator species_iter2 = p_speciesSeed.iterator();
        			while (species_iter2.hasNext()) {
        				Species second_reactant = (Species)species_iter2.next();
        				HashSet current_reactions = current_template.reactTwoReactants(first_reactant,second_reactant);
        				reaction_set.addAll(current_reactions);
        			}
        		}
        	}
        }
        
        return reaction_set;
        
        
        
        //#]
    }
    
    //## operation react(HashSet,Species) 
    public HashSet react(HashSet p_speciesSet, Species p_species) {
        //#[ operation react(HashSet,Species) 
        HashSet reaction_set = new HashSet();
        
        if (p_speciesSet.size() == 0 && p_species == null) {
        	return reaction_set;
        }
        
		double singleReaction = 0, doubleReaction = 0;
		double longestTime = 0;
		String longestTemplate = "";
		StringBuffer HAbs = new StringBuffer();//"H_Abstraction");
		
        // add here the algorithm to generate reaction
        // loop over all the reaction template to find any possible match between the species seed set and the reaction template library
        Iterator template_iter = reactionTemplateLibrary.getReactionTemplate();
        while (template_iter.hasNext()) {
        	ReactionTemplate current_template = (ReactionTemplate)template_iter.next();
        	// the reaction template has only one reactant, we only need to loop over the whole species seed set to find a match
        	double startTime = System.currentTimeMillis();
			if (current_template.hasOneReactant()) {
        		HashSet current_reactions = current_template.reactOneReactant(p_species);
        		reaction_set.addAll(current_reactions);   
				singleReaction = singleReaction + ((System.currentTimeMillis()-startTime)/1000/60);
        	}
			
        	// the reaction template has two reactants, we need to check all the possible combination of two species
			
			
        	else if (current_template.hasTwoReactants()) {
				double tempStart = System.currentTimeMillis();
				double species_StartTime = 0;
        		Iterator species_iter = p_speciesSet.iterator();
        		while (species_iter.hasNext()) {
					if (current_template.name.equals("H_Abstraction"))
						species_StartTime = System.currentTimeMillis();
					
        			Species first_reactant = (Species)species_iter.next();
        			HashSet current_reactions = current_template.reactTwoReactants(first_reactant,p_species);
        			reaction_set.addAll(current_reactions);
        			current_reactions = current_template.reactTwoReactants(p_species,first_reactant);
        			reaction_set.addAll(current_reactions);
					if (current_template.name.equals("H_Abstraction")) {
						double speciesTime = ((System.currentTimeMillis()-species_StartTime)/1000/60);
						if (speciesTime >= 1)
							HAbs.append(first_reactant.getChemkinName() + "\t" + speciesTime + "\t" + current_reactions.size() + "\t");
					}
					
        		}
				doubleReaction += ((System.currentTimeMillis()-startTime)/1000/60);
				double thisDoubleReaction = ((System.currentTimeMillis()-startTime)/1000/60);
				if (thisDoubleReaction >= longestTime){
					longestTime = thisDoubleReaction;
					longestTemplate = current_template.name;
				}
        		if (!p_speciesSet.contains(p_species)) {
        			HashSet current_reactions = current_template.reactTwoReactants(p_species, p_species);
        			reaction_set.addAll(current_reactions);
        		}
        	}
        }
        
		Global.enlargerInfo.append(p_species.getChemkinName() + "\t" + singleReaction + "\t" + doubleReaction + "\t" + longestTime +"\t" + longestTemplate + "\t" + HAbs.toString() + "\n");
		
        return reaction_set;
        
        
        
        //#]
    }
    
    public ReactionTemplateLibrary getReactionTemplateLibrary() {
        return reactionTemplateLibrary;
    }
    
    public void setReactionTemplateLibrary(ReactionTemplateLibrary p_ReactionTemplateLibrary) {
        reactionTemplateLibrary = p_ReactionTemplateLibrary;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\TemplateReactionGenerator.java
*********************************************************************/


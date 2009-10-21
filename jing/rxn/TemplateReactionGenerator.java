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
import jing.param.Global;
import jing.param.Temperature;

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
    public LinkedHashSet react(LinkedHashSet p_speciesSeed) {
        //#[ operation react(HashSet) 
        if (p_speciesSeed.size() == 0) {
        	return null;
        }
        
        LinkedHashSet reaction_set = new LinkedHashSet();
        //LinkedHashSet species_set = new LinkedHashSet();
        
        // add here the algorithm to generate reaction
        // loop over all the reaction template to find any possible match between the species seed set and the reaction template library
        Iterator template_iter = reactionTemplateLibrary.getReactionTemplate();
        while (template_iter.hasNext()) {
        	ReactionTemplate current_template = (ReactionTemplate)template_iter.next();
        	/*
        	 * Added by MRH on 15-Jun-2009
        	 * Display more information to the user:
        	 * 	This println command informs the user which rxn family template
        	 *		the new seed species is reacting against
        	 */
        	System.out.println("Reacting the species given in the input file: " + current_template.name);
        	// the reaction template has only one reactant, we only need to loop over the whole species seed set to find a match
//        	System.out.println(current_template.name);
        	if (current_template.hasOneReactant()) {
        		Iterator species_iter1 = p_speciesSeed.iterator();
        		while (species_iter1.hasNext()) {
        			Species first_reactant = (Species)species_iter1.next();
        			if (first_reactant.isReactive()) {
        				LinkedHashSet current_reactions = current_template.reactOneReactant(first_reactant);
        				first_reactant.addPdepPaths(current_reactions);
        				reaction_set.addAll(current_reactions);
        			}
    				Runtime runTime = Runtime.getRuntime();
    				if (runTime.freeMemory() < runTime.totalMemory()/3) {
    					runTime.gc();
    				}
    				//System.out.println("Free memory: " + runTime.freeMemory());        			
        		}                                
        	}
        	// the reaction template has two reactants, we need to check all the possible combination of two species
        	else if (current_template.hasTwoReactants()) {
				
        		Iterator species_iter1 = p_speciesSeed.iterator();
        		while (species_iter1.hasNext()) {
        			Species first_reactant = (Species)species_iter1.next();
        			if (first_reactant.isReactive()) {
        				Iterator species_iter2 = p_speciesSeed.iterator();
        				while (species_iter2.hasNext()) {
        					Species second_reactant = (Species)species_iter2.next();
        					if (second_reactant.isReactive()) {
 		       					LinkedHashSet current_reactions = current_template.reactTwoReactants(first_reactant,second_reactant);
        						reaction_set.addAll(current_reactions);
        					}
        				}
        				Runtime runTime = Runtime.getRuntime();
        				if (runTime.freeMemory() < runTime.totalMemory()/3) {
        					runTime.gc();
        				}
        				//System.out.println("Free memory: " + runTime.freeMemory());
        			}
        		}
        	}
        }
        
        //PDepNetwork.completeNetwork(p_speciesSeed);
		Runtime runTime = Runtime.getRuntime();
		System.out.print("Memory used: ");
		System.out.println(runTime.totalMemory());
		System.out.print("Free memory: ");
		System.out.println(runTime.freeMemory());
		if (runTime.freeMemory() < runTime.totalMemory()/3) {
			runTime.gc();
			System.out.println("After garbage collection:");
			System.out.print("Memory used: ");
			System.out.println(runTime.totalMemory());
			System.out.print("Free memory: ");
			System.out.println(runTime.freeMemory());
		}

        return reaction_set;
        
        
        
        //#]
    }
    
    public LinkedHashSet generatePdepReactions(Species p_species){
    	Iterator template_iter = reactionTemplateLibrary.getReactionTemplate();
    	LinkedHashSet pdepReactionSet = new LinkedHashSet();
        while (template_iter.hasNext()) {
        	ReactionTemplate current_template = (ReactionTemplate)template_iter.next();
        	// the reaction template has only one reactant, we only need to loop over the whole species seed set to find a match
        	double startTime = System.currentTimeMillis();
			if (current_template.hasOneReactant()) {
	        	/*
	        	 * Added by MRH on 15-Jun-2009
	        	 * Display more information to the user:
	        	 * 	This println command informs the user which rxn family template
	        	 *		the species is reacting against
	        	 */
	        	System.out.println("Generating pressure dependent network for " + p_species.getChemkinName() + ": " + current_template.name);
        		
				LinkedHashSet current_reactions = current_template.reactOneReactant(p_species);
				 
				pdepReactionSet.addAll(current_reactions);
				
        	}
        }
		Runtime runTime = Runtime.getRuntime();
		System.out.print("Memory used: ");
		System.out.println(runTime.totalMemory());
		System.out.print("Free memory: ");
		System.out.println(runTime.freeMemory());
		if (runTime.freeMemory() < runTime.totalMemory()/3) {
			runTime.gc();
			System.out.println("After garbage collection:");
			System.out.print("Memory used: ");
			System.out.println(runTime.totalMemory());
			System.out.print("Free memory: ");
			System.out.println(runTime.freeMemory());
		}
		
        return pdepReactionSet;
    }
    
//    //## operation react(HashSet,Species) 
//    public LinkedHashSet react(LinkedHashSet p_speciesSet, Species p_species) {
//        //#[ operation react(HashSet,Species) 
//    	LinkedHashSet reaction_set = new LinkedHashSet();
//    	if (!p_species.isReactive()){
//    		return reaction_set;
//    	}
//        
//        if (p_speciesSet.size() == 0 && p_species == null) {
//        	return reaction_set;
//        }
//        
//		double singleReaction = 0, doubleReaction = 0;
//		double longestTime = 0;
//		String longestTemplate = "";
//		StringBuffer HAbs = new StringBuffer();//"H_Abstraction");
//		
//		HashSet pdepReactionSet = new HashSet();
//		
//        // add here the algorithm to generate reaction
//        // loop over all the reaction template to find any possible match between the species seed set and the reaction template library
//        Iterator template_iter = reactionTemplateLibrary.getReactionTemplate();
//        while (template_iter.hasNext()) {
//        	ReactionTemplate current_template = (ReactionTemplate)template_iter.next();
//        	/*
//        	 * Added by MRH on 12-Jun-2009
//        	 * Display more information to the user:
//        	 * 	This println command informs the user which rxn family template
//        	 *		the new core species is reacting against
//        	 */
//        	System.out.println("Reacting " + p_species.getChemkinName() + " with the core: " + current_template.name);
//        	
//        	// the reaction template has only one reactant, we only need to loop over the whole species seed set to find a match
//        	double startTime = System.currentTimeMillis();
//			if (current_template.hasOneReactant()) {
//        		
//				LinkedHashSet current_reactions = current_template.reactOneReactant(p_species);
//				reaction_set.addAll(current_reactions);   
//				pdepReactionSet.addAll(current_reactions);
//				singleReaction = singleReaction + ((System.currentTimeMillis()-startTime)/1000/60);
//        		
//        	}
//			
//        	// the reaction template has two reactants, we need to check all the possible combination of two species
//			
//			
//        	else if (current_template.hasTwoReactants()) {
//				double tempStart = System.currentTimeMillis();
//				double species_StartTime = 0;
//        		Iterator species_iter = p_speciesSet.iterator();
//        		while (species_iter.hasNext()) {
//					if (current_template.name.equals("H_Abstraction"))
//						species_StartTime = System.currentTimeMillis();
//					
//        			Species first_reactant = (Species)species_iter.next();
//        			if (first_reactant.isReactive()){
//        				LinkedHashSet current_reactions = current_template.reactTwoReactants(first_reactant,p_species);
//            			reaction_set.addAll(current_reactions);
//            			current_reactions = current_template.reactTwoReactants(p_species,first_reactant);
//            			reaction_set.addAll(current_reactions);
//    					if (current_template.name.equals("H_Abstraction")) {
//    						double speciesTime = ((System.currentTimeMillis()-species_StartTime)/1000/60);
//    						if (speciesTime >= 1)
//    							HAbs.append(first_reactant.getChemkinName() + "\t" + speciesTime + "\t" + current_reactions.size() + "\t");
//    					}
//        			}
//        		}
//				doubleReaction += ((System.currentTimeMillis()-startTime)/1000/60);
//				double thisDoubleReaction = ((System.currentTimeMillis()-startTime)/1000/60);
//				if (thisDoubleReaction >= longestTime){
//					longestTime = thisDoubleReaction;
//					longestTemplate = current_template.name;
//				}
//        		if (!p_speciesSet.contains(p_species)) {
//        			LinkedHashSet current_reactions = current_template.reactTwoReactants(p_species, p_species);
//        			reaction_set.addAll(current_reactions);
//        		}
//        	}
//        }
//        
//        p_species.addPdepPaths(pdepReactionSet);
//        
//		Global.enlargerInfo.append(p_species.getChemkinName() + "\t" + singleReaction + "\t" + doubleReaction + "\t" + longestTime +"\t" + longestTemplate + "\t" + HAbs.toString() + "\n");
//		
//		//PDepNetwork.completeNetwork(p_species);
//		
//		Runtime runTime = Runtime.getRuntime();
//		System.out.print("Memory used: ");
//		System.out.println(runTime.totalMemory());
//		System.out.print("Free memory: ");
//		System.out.println(runTime.freeMemory());
//		runTime.gc();
//		System.out.println("After garbage collection:");
//		System.out.print("Memory used: ");
//		System.out.println(runTime.totalMemory());
//		System.out.print("Free memory: ");
//		System.out.println(runTime.freeMemory());
//		
//        return reaction_set;
//        
//        
//        
//        //#]
//    }
    
            //## operation react(HashSet,Species) 
    public LinkedHashSet react(LinkedHashSet p_speciesSet, Species newCoreSpecies) {
//    	double pT = System.currentTimeMillis();
        //#[ operation react(HashSet,Species) 
    	LinkedHashSet reaction_set = new LinkedHashSet();
    	if (!newCoreSpecies.isReactive()){
    		return reaction_set;
    	}
        
        if (p_speciesSet.size() == 0 && newCoreSpecies == null) {
        	return reaction_set;
        }
        
		double singleReaction = 0, doubleReaction = 0;
		double longestTime = 0;
		String longestTemplate = "";
		StringBuffer HAbs = new StringBuffer();//"H_Abstraction");
		
		HashSet pdepReactionSet = new HashSet();
		
        // add here the algorithm to generate reaction
        // loop over all the reaction template to find any possible match between the species seed set and the reaction template library
        Iterator template_iter = reactionTemplateLibrary.getReactionTemplate();
        while (template_iter.hasNext()) {
        	ReactionTemplate current_template = (ReactionTemplate)template_iter.next();
        	/*
        	 * Added by MRH on 12-Jun-2009
        	 * Display more information to the user:
        	 * 	This println command informs the user which rxn family template
        	 *		the new core species is reacting against
        	 */
        	System.out.println("Reacting " + newCoreSpecies.getChemkinName() + " with the core: " + current_template.name);
        	
        	// the reaction template has only one reactant, we only need to loop over the whole species seed set to find a match
        	double startTime = System.currentTimeMillis();
			if (current_template.hasOneReactant()) {
        		
				LinkedHashSet current_reactions = current_template.reactOneReactant(newCoreSpecies);
				reaction_set.addAll(current_reactions);   
				pdepReactionSet.addAll(current_reactions);
				singleReaction = singleReaction + ((System.currentTimeMillis()-startTime)/1000/60);
        		
        	}
			
        	// the reaction template has two reactants, we need to check all the possible combination of two species
			
			
        	else if (current_template.hasTwoReactants()) {
//        		LinkedHashSet current_reactions = new LinkedHashSet();
        		StructureTemplate structTemp = current_template.structureTemplate;
        		
        		LinkedHashSet site1_reactiveSites_sp1 = new LinkedHashSet();
        		LinkedHashSet site2_reactiveSites_sp1 = new LinkedHashSet();
        		LinkedHashSet site1_reactiveSites_sp2 = new LinkedHashSet();
        		LinkedHashSet site2_reactiveSites_sp2 = new LinkedHashSet();
        		
        		if (!newCoreSpecies.hasResonanceIsomers()) {
        			ChemGraph newCoreCG = newCoreSpecies.getChemGraph();
        			site1_reactiveSites_sp1 = structTemp.identifyReactedSites(newCoreCG,1);
        			site2_reactiveSites_sp1 = structTemp.identifyReactedSites(newCoreCG,2);
        			Iterator coreSpeciesIter = p_speciesSet.iterator();
        			while (coreSpeciesIter.hasNext()) {
        				Species coreSpecies = (Species)coreSpeciesIter.next();
        				if (coreSpecies.isReactive()) {
        					if (!coreSpecies.hasResonanceIsomers()) {
        						ChemGraph tempCG = coreSpecies.getChemGraph();
        						ChemGraph oldCoreCG = generateCGcopyIfNecessary(tempCG,newCoreCG);
        						if (!site1_reactiveSites_sp1.isEmpty()) site2_reactiveSites_sp2 = structTemp.identifyReactedSites(oldCoreCG,2);
        						if (!site2_reactiveSites_sp1.isEmpty()) site1_reactiveSites_sp2 = structTemp.identifyReactedSites(oldCoreCG,1);
        	        			//	React A + B
        	        			LinkedHashSet current_reactions = current_template.reactTwoReactants(newCoreCG,site1_reactiveSites_sp1,oldCoreCG,site2_reactiveSites_sp2);
        	        			reaction_set.addAll(current_reactions);
        	        			//	React B + A
        	        			current_reactions = current_template.reactTwoReactants(oldCoreCG,site1_reactiveSites_sp2,newCoreCG,site2_reactiveSites_sp1);
        	        			reaction_set.addAll(current_reactions);
        					} else {
        						Iterator coreSpeciesCGIter = coreSpecies.getResonanceIsomers();
        						while (coreSpeciesCGIter.hasNext()) {
            						ChemGraph tempCG = (ChemGraph)coreSpeciesCGIter.next();
            						ChemGraph oldCoreCG = generateCGcopyIfNecessary(tempCG,newCoreCG);
            						if (!site1_reactiveSites_sp1.isEmpty()) site2_reactiveSites_sp2 = structTemp.identifyReactedSites(oldCoreCG,2);
            						if (!site2_reactiveSites_sp1.isEmpty()) site1_reactiveSites_sp2 = structTemp.identifyReactedSites(oldCoreCG,1);
            	        			//	React A + B
            	        			LinkedHashSet current_reactions = current_template.reactTwoReactants(newCoreCG,site1_reactiveSites_sp1,oldCoreCG,site2_reactiveSites_sp2);
            	        			reaction_set.addAll(current_reactions);
            	        			//	React B + A
            	        			current_reactions = current_template.reactTwoReactants(oldCoreCG,site1_reactiveSites_sp2,newCoreCG,site2_reactiveSites_sp1);
            	        			reaction_set.addAll(current_reactions);
	        					}
        					}
        				}
        			}
        		} else {
        			Iterator newSpeciesCGIter = newCoreSpecies.getResonanceIsomers();
        			while (newSpeciesCGIter.hasNext()) {
        				ChemGraph newCoreCG = (ChemGraph)newSpeciesCGIter.next();
	        			site1_reactiveSites_sp1 = structTemp.identifyReactedSites(newCoreCG,1);
	        			site2_reactiveSites_sp1 = structTemp.identifyReactedSites(newCoreCG,2);
	        			Iterator coreSpeciesIter = p_speciesSet.iterator();
	        			while (coreSpeciesIter.hasNext()) {
	        				Species coreSpecies = (Species)coreSpeciesIter.next();
	        				if (coreSpecies.isReactive()) {
	        					if (!coreSpecies.hasResonanceIsomers()) {
	        						ChemGraph tempCG = coreSpecies.getChemGraph();
	        						ChemGraph oldCoreCG = generateCGcopyIfNecessary(tempCG,newCoreCG);
	        						if (!site1_reactiveSites_sp1.isEmpty()) site2_reactiveSites_sp2 = structTemp.identifyReactedSites(oldCoreCG,2);
	        						if (!site2_reactiveSites_sp1.isEmpty()) site1_reactiveSites_sp2 = structTemp.identifyReactedSites(oldCoreCG,1);
	        	        			//	React A + B
	        	        			LinkedHashSet current_reactions = current_template.reactTwoReactants(newCoreCG,site1_reactiveSites_sp1,oldCoreCG,site2_reactiveSites_sp2);
	        	        			reaction_set.addAll(current_reactions);
	        	        			//	React B + A
	        	        			current_reactions = current_template.reactTwoReactants(oldCoreCG,site1_reactiveSites_sp2,newCoreCG,site2_reactiveSites_sp1);
	        	        			reaction_set.addAll(current_reactions);
	        					} else {
	        						Iterator coreSpeciesCGIter = coreSpecies.getResonanceIsomers();
	        						while (coreSpeciesCGIter.hasNext()) {
	            						ChemGraph tempCG = (ChemGraph)coreSpeciesCGIter.next();
	            						ChemGraph oldCoreCG = generateCGcopyIfNecessary(tempCG,newCoreCG);
	            						if (!site1_reactiveSites_sp1.isEmpty()) site2_reactiveSites_sp2 = structTemp.identifyReactedSites(oldCoreCG,2);
	            						if (!site2_reactiveSites_sp1.isEmpty()) site1_reactiveSites_sp2 = structTemp.identifyReactedSites(oldCoreCG,1);
	            	        			//	React A + B
	            	        			LinkedHashSet current_reactions = current_template.reactTwoReactants(newCoreCG,site1_reactiveSites_sp1,oldCoreCG,site2_reactiveSites_sp2);
	            	        			reaction_set.addAll(current_reactions);
	            	        			//	React B + A
	            	        			current_reactions = current_template.reactTwoReactants(oldCoreCG,site1_reactiveSites_sp2,newCoreCG,site2_reactiveSites_sp1);
	            	        			reaction_set.addAll(current_reactions);
		        					}
	        					}
	        				}
	        			}
        			}
        		}
        	}
			doubleReaction += ((System.currentTimeMillis()-startTime)/1000/60);
			double thisDoubleReaction = ((System.currentTimeMillis()-startTime)/1000/60);
			if (thisDoubleReaction >= longestTime){
				longestTime = thisDoubleReaction;
				longestTemplate = current_template.name;
			}
        }
		
        newCoreSpecies.addPdepPaths(pdepReactionSet);
        
		Global.enlargerInfo.append(newCoreSpecies.getChemkinName() + "\t" + singleReaction + "\t" + doubleReaction + "\t" + longestTime +"\t" + longestTemplate + "\t" + HAbs.toString() + "\n");
		
		//PDepNetwork.completeNetwork(p_species);
		
		Runtime runTime = Runtime.getRuntime();
		System.out.print("Memory used: ");
		System.out.println(runTime.totalMemory());
		System.out.print("Free memory: ");
		System.out.println(runTime.freeMemory());
		if (runTime.freeMemory() < runTime.totalMemory()/3) {
			runTime.gc();
			System.out.println("After garbage collection:");
			System.out.print("Memory used: ");
			System.out.println(runTime.totalMemory());
			System.out.print("Free memory: ");
			System.out.println(runTime.freeMemory());
		}
		
//		double t = (System.currentTimeMillis()-pT)/1000/60;
//        Global.RT_reactTwoReactants += t;
		
        return reaction_set;
        
        
        
        //#]
    }
    

    
    public ChemGraph generateCGcopyIfNecessary(ChemGraph cg1, ChemGraph cg2) {
    	ChemGraph cg_copy = null;
		if (cg1 == cg2) {
			try {
				cg_copy = ChemGraph.copy(cg2);
			} catch (ForbiddenStructureException e) {
				System.out.println("Forbidden Structure encountered in react(): " + e.toString());
			}
		} else return cg1;
		return cg_copy;
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


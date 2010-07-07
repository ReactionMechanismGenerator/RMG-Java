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


import jing.rxn.*;
import jing.chem.*;
import java.util.*;

import jing.chem.Species;
import jing.rxn.Reaction;
import jing.param.Global;
import jing.param.Pressure;
import jing.param.Temperature;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\CoreEdgeReactionModel.java                                                                  
//----------------------------------------------------------------------------

//## class CoreEdgeReactionModel 
public class CoreEdgeReactionModel implements ReactionModel {
    
    protected Core core;
    protected Edge edge;
	
    
    // Constructors
    
    //## operation CoreEdgeReactionModel() 
    public CoreEdgeReactionModel() {
        initRelations();
        //#[ operation CoreEdgeReactionModel() 
        //#]
    }

	
    // Argument HashSetp_reactedSpeciesSet : 
    /**
    reacted species set
    */
    // Argument HashSetp_reactionSet : 
    /**
    All the reaction generated from the reacted species set
    */
    //## operation CoreEdgeReactionModel(HashSet,HashSet) 
    public  CoreEdgeReactionModel(LinkedHashSet p_reactedSpeciesSet, LinkedHashSet p_reactionSet) {
        initRelations();
        
        
        core.setSpeciesSet(p_reactedSpeciesSet);
        addReactionSet(p_reactionSet);

        //#]
    }

	/**
	 * This initializer just takes in the species and the reactions are added from the Pdep-networks
	 * 
	 */
	public CoreEdgeReactionModel(LinkedHashSet p_reactedSpeciesSet){
		initRelations();
		core.setSpeciesSet(p_reactedSpeciesSet);
		Iterator iter = PDepNetwork.getNetworks().iterator();
		while (iter.hasNext()){
			PDepNetwork pdn = (PDepNetwork)iter.next();
			ListIterator reactionIter = pdn.getNetReactions().listIterator();
			while (reactionIter.hasNext()) {
				addReaction(((Reaction)reactionIter.next()));
			}
			
		}
	}
    
    
    public void addPrimaryKineticSet(LinkedHashSet p_reactionSet) {
        
        for (Iterator iter = p_reactionSet.iterator(); iter.hasNext(); ) {
        	Reaction rxn = (Reaction)iter.next();
        	int rxnType = categorizeReaction(rxn);
        	// here the same reaction has been generated from template, we need to remove that one, and then add the one from PKL
        	if (rxnType == 1) {
        		if (containsAsReactedReaction(rxn)) {
        			// remove the present one
        			for (Iterator rIter = getReactedReactionSet().iterator(); rIter.hasNext(); ) {
        				Reaction r = (Reaction)rIter.next();
        				if (r.equals(rxn)) rIter.remove();
        			}
        		}
        		addReactedReaction(rxn);
        	}
        	else if (rxnType == -1) {
        		if (containsAsUnreactedReaction(rxn)) {
        			// remove the present one
        			for (Iterator rIter = getUnreactedReactionSet().iterator(); rIter.hasNext(); ) {
        				Reaction r = (Reaction)rIter.next();
        				if (r.equals(rxn)) rIter.remove();
        			}
        
        		}
        		addUnreactedReaction(rxn);
        	}
        	
        	if (rxn.hasReverseReaction()){
        		Reaction reverse = rxn.getReverseReaction();
            	rxnType = categorizeReaction(reverse);
            	// here the same reaction has been generated from template, we need to remove that one, and then add the one from PRL
            	if (rxnType == 1) {
            		if (containsAsReactedReaction(reverse)) {
            			// remove the present one
            			for (Iterator rIter = getReactedReactionSet().iterator(); rIter.hasNext(); ) {
            				Reaction r = (Reaction)rIter.next();
            				if (r.equals(reverse)) rIter.remove();
            			}
            		}
            		addReactedReaction(reverse);
            	}
            	else if (rxnType == -1) {
            		if (containsAsUnreactedReaction(reverse)) {
            			// remove the present one
            			for (Iterator rIter = getUnreactedReactionSet().iterator(); rIter.hasNext(); ) {
            				Reaction r = (Reaction)rIter.next();
            				if (r.equals(reverse)) rIter.remove();
            			}
            
            		}
            		addUnreactedReaction(reverse);
            	}
        	}
        
        }
        
        return;
        
        
        //#]
    }
    
    //## operation addReactedReaction(Reaction) 
    public void addReactedReaction(Reaction p_reaction) throws InvalidReactedReactionException {
        //#[ operation addReactedReaction(Reaction) 
        if (isReactedReaction(p_reaction)) {
			getReactedReactionSet().add(p_reaction);
			//if (p_reaction.isForward()) System.out.println(p_reaction.toRestartString());
        }
        else throw new InvalidReactedReactionException(p_reaction.toString());
        //#]
    }
    
    //## operation addReactedReactionSet(HashSet) 
    public void addReactedReactionSet(LinkedHashSet p_reactedReactionSet) throws InvalidReactedReactionException {
        //#[ operation addReactedReactionSet(HashSet) 
        try {
        	for (Iterator iter = p_reactedReactionSet.iterator(); iter.hasNext(); ) {
        		Reaction r = (Reaction)iter.next();
        		addReactedReaction(r);
        	}
        }
        catch (InvalidReactedReactionException e) {
        	throw new InvalidReactedReactionException(e.getMessage());
        }   
        
        
        //#]
    }
    
    //## operation addReactedSpecies(Species) 
    public void addReactedSpecies(Species p_species) {
        //#[ operation addReactedSpecies(Species) 
        if (containsAsUnreactedSpecies(p_species)) {
        	moveFromUnreactedToReactedSpecies(p_species);
        }
        else {
        	getReactedSpeciesSet().add(p_species);
        }
        
        moveFromUnreactedToReactedReaction();
        //#]
    }
    
    //## operation addReactedSpeciesSet(HashSet) 
    public void addReactedSpeciesSet(LinkedHashSet p_reactedSpeciesSet) {
        //#[ operation addReactedSpeciesSet(HashSet) 
        boolean added = false;
        LinkedHashSet rs = getReactedSpeciesSet();
        for (Iterator iter = p_reactedSpeciesSet.iterator(); iter.hasNext(); ) {
        	Species spe = (Species)iter.next();
        	if (containsAsUnreactedSpecies(spe)) {
        		moveFromUnreactedToReactedSpecies(spe);
        	}
        	else {
        		rs.add(spe);
				if (!added){
					Iterator originalSpeciesIterator = rs.iterator();
					while (originalSpeciesIterator.hasNext()){
						Species temp = (Species)originalSpeciesIterator.next();
						if(spe.equals(temp)){
							temp.setName(temp.getName());
						}
					}
				}
        	}
        }
        
        moveFromUnreactedToReactedReaction();
        
        
        //#]
    }
    
    //## operation addReactionSet(HashSet) 
    public void addReactionSet(LinkedHashSet p_reactionSet) {
        //#[ operation addReactionSet(HashSet) 
        for (Iterator iter = p_reactionSet.iterator(); iter.hasNext(); ) {
        	Reaction rxn = (Reaction)iter.next();
        	int rxnType = categorizeReaction(rxn);
        	if (rxnType == 1) {
        		addReactedReaction(rxn);
        	}
        	else if (rxnType == -1) {
        		addUnreactedReaction(rxn);
        	}
        	
        	//also add the reverse reaction 
        	if (rxn.hasReverseReaction()){
        		Reaction reverse = (Reaction)rxn.getReverseReaction();
        		rxnType = categorizeReaction(reverse);
        		if (rxnType == 1){
        			addReactedReaction(reverse);
        		}
        		else if (rxnType == -1){
        			addUnreactedReaction(reverse);
        		}
        	}
        	
        }
        
        return;
        //#]
    }
    
    /**
     * Add a given reaction. First Categorize the
     * reaction as either core or Edge and then add it to the list
     * 
     */
    public void addReaction(Reaction p_reaction){
    	int rxnType = categorizeReaction(p_reaction);
    	if (rxnType == 1){
    		addReactedReaction(p_reaction);
    	}
    	else
    		addUnreactedReaction(p_reaction);
    }
    
    //## operation addUnreactedReaction(Reaction) 
    public void addUnreactedReaction(Reaction p_reaction) throws InvalidUnreactedReactionException {
        //#[ operation addUnreactedReaction(Reaction) 
    	if (p_reaction instanceof PDepReaction)
    		System.out.println(p_reaction);
    	if (p_reaction.hasReverseReaction()) {
	        if (isUnreactedReaction(p_reaction)) getUnreactedReactionSet().add(p_reaction);
	        else if (isUnreactedReaction(p_reaction.getReverseReaction())) getUnreactedReactionSet().add(p_reaction.getReverseReaction());
	        else if (isUnreactedReversiblePathReaction(p_reaction))
			{ 
				// We need to run isUnreactedReversiblePathReaction() because it
				// does some work with the reactants and products of the path
				// reaction
				// However, we don't want to add the path reaction to the edge,
				// so we do nothing here
			}
	        else throw new InvalidUnreactedReactionException(p_reaction.toString());
    	} else {
    		if (isUnreactedReaction(p_reaction)) getUnreactedReactionSet().add(p_reaction);
	        else if (isUnreactedIrreversiblePathReaction(p_reaction))
			{
				// Again, don't add a pressure-dependent path reaction to the
				// edge here
			}
	        else throw new InvalidUnreactedReactionException(p_reaction.toString());
    	}
        //#]
    }
    
    public void addUnreactedReactionSet(LinkedHashSet reactions) {
		Iterator rxnIter = reactions.iterator();
		while (rxnIter.hasNext()){
			Reaction rxn = (Reaction)rxnIter.next();
			System.out.println(rxn.getStructure().toString());
			addUnreactedReaction(rxn);
		}
    }
    
    //## operation addUnreactedSpecies(Species) 
    public void addUnreactedSpecies(Species p_species) {
        //#[ operation addUnreactedSpecies(Species) 
        if (containsAsReactedSpecies(p_species)) {
        	// this is not a unreacted species
        	System.out.println("This is a reacted species " + p_species.getName());
        	System.out.println("Can't add it into unreacted species set!");
        }
        else {
        	getUnreactedSpeciesSet().add(p_species);
        }
        //#]
    }
	
//	## operation addUnreactedSpecies(Species) 
    public void addUnreactedSpeciesSet(LinkedHashSet p_species) {
        //#[ operation addUnreactedSpecies(Species)
		Iterator speciesIter = p_species.iterator();
		while (speciesIter.hasNext()){
			Species species = (Species)speciesIter.next();
			if (containsAsReactedSpecies(species)) {
	        	// this is not a unreacted species
	        	System.out.println("This is a reacted species " + species.getName());
	        	System.out.println("Can't add it into unreacted species set!");
	        }
	        else {
	        	getUnreactedSpeciesSet().add(species);
	        }
		}
        
        //#]
    }
    
    /**
    Requires: the reactionSpecies set and unreactedSpecies set has been defined properly.
    Effects: according to the reaction's reactants and products, categorize the pass-in reaction as reacted reaction (return 1), or unreacted reaction(return -1), or reaction not in the model (return 0).
    Modifies:
    */
    protected int categorizeReaction(Reaction p_reaction) {
		return categorizeReaction(p_reaction.getStructure());
    }
	
	/**
    Requires: the reactionSpecies set and unreactedSpecies set has been defined properly.
    Effects: according to the reaction's reactants and products, categorize the pass-in reaction as reacted reaction (return 1), or unreacted reaction(return -1), or reaction not in the model (return 0).
    Modifies:
    */
    //## operation categorizeReaction(Reaction) 
    public int categorizeReaction(Structure p_structure) {
        //#[ operation categorizeReaction(Reaction) 
	if (p_structure == null) throw new NullPointerException();
	
	if (!reactantsInCoreQ(p_structure)){
	    return 0;
	}
        
        int type = 1;
        Iterator iter = p_structure.getProducts();
        while (iter.hasNext()) {
			Species spe = (Species)iter.next();
        	if (!contains(spe)) {
        		// new unreacted species 
        		type = -1;
        		addUnreactedSpecies(spe);
        	}
        	else if (containsAsUnreactedSpecies(spe)) {
        		type = -1;
        	}
        }	                                          
        
        return type;
        //#]
    }

    public boolean reactantsInCoreQ(Structure p_structure){
	Iterator iter = p_structure.getReactants();
        while (iter.hasNext()) {
			Species spe = (Species)iter.next();
        	if (!containsAsReactedSpecies(spe))
        		return false;
        }
	return true;
    }
    
    //## operation contains(Reaction) 
    public boolean contains(Reaction p_reaction) {
        //#[ operation contains(Reaction) 
        return (containsAsReactedReaction(p_reaction) || containsAsUnreactedReaction(p_reaction));
        //#]
    }
    
    //## operation contains(Species) 
    public boolean contains(Species p_species) {
        //#[ operation contains(Species) 
        return (containsAsReactedSpecies(p_species) || containsAsUnreactedSpecies(p_species));
        //#]
    }
    
    //## operation containsAsReactedReaction(Reaction) 
    public boolean containsAsReactedReaction(Reaction p_reaction) {
        //#[ operation containsAsReactedReaction(Reaction) 
        return getReactedReactionSet().contains(p_reaction);
        //#]
    }
    
    //## operation containsAsReactedSpecies(Species) 
    public boolean containsAsReactedSpecies(Species p_species) {
        //#[ operation containsAsReactedSpecies(Species) 
        return getReactedSpeciesSet().contains(p_species);
        //#]
    }
    
    //## operation containsAsUnreactedReaction(Reaction) 
    public boolean containsAsUnreactedReaction(Reaction p_reaction) {
        //#[ operation containsAsUnreactedReaction(Reaction) 
        return getUnreactedReactionSet().contains(p_reaction);
        //#]
    }
    
    //## operation containsAsUnreactedSpecies(Species) 
    public boolean containsAsUnreactedSpecies(Species p_species) {
        //#[ operation containsAsUnreactedSpecies(Species) 
        return getUnreactedSpeciesSet().contains(p_species);
        //#]
    }
    
    //## operation getReactedReactionSet() 
    public LinkedHashSet getReactedReactionSet() {
        //#[ operation getReactedReactionSet() 
        return getCore().getReactionSet();
        //#]
    }
    
    //## operation getReactedSpeciesSet() 
    public LinkedHashSet getReactedSpeciesSet() {
        //#[ operation getReactedSpeciesSet() 
        return getCore().getSpeciesSet();
        //#]
    }
    
    //## operation getReaction() 
    public Iterator getReaction() {
        //#[ operation getReaction() 
        return core.getReaction();
        //#]
    }
  
	 public LinkedList generatePDepReactionSet() {
	        //#[ operation generatePDepReactionSet() 
	        LinkedList nonPDepList = new LinkedList();
	        LinkedList pDepList = new LinkedList();
	        LinkedList duplicates = new LinkedList();
	        
	        
	        
	        LinkedHashSet rSet = getReactionSet();
	        for (Iterator iter = rSet.iterator(); iter.hasNext(); ) {
	        	Reaction r = (Reaction)iter.next();
	        	if ((r instanceof ThirdBodyReaction || r instanceof TROEReaction || r instanceof LindemannReaction) && r.isForward()) {
	        		pDepList.add(r);
	        		
	        	}
	        }
	        
	        for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
	        	PDepNetwork pdn = (PDepNetwork)iter.next();
	        	for (ListIterator pdniter = pdn.getNetReactions().listIterator(); pdniter.hasNext();) {
	        		PDepReaction rxn = (PDepReaction) pdniter.next();
	        		Structure s = rxn.getStructure();
	        		if (!rxn.reactantEqualsProduct() && isReactedReaction(rxn)) {
	        			// if 2 reactants = 2 products then compare it with reaction in the core.
	        			// if the same reaction is present then the reactions should go to the duplicate list.
	        			if (rxn.getStructure().getReactantNumber() == 2 && rxn.getStructure().getProductNumber() ==2) {
	        				if (containsAsReactedReaction(rxn))
	        					duplicates.add(rxn);
	        				else
	        					pDepList.add(rxn);
	        			}
	        			else
	        				pDepList.add(rxn);
	        			
	        		}
	        	}
	        }
	        
	        Iterator iter = getReactionSet().iterator();
	        while (iter.hasNext()) {
	        	Reaction r = (Reaction)iter.next();
	        	if (!r.reactantEqualsProduct() && !(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction) && !(r instanceof LindemannReaction) && r.isForward()) {
	        		
	        		if (!pDepList.contains(r)  ) {
	        			nonPDepList.add(r);
	        		}
	        		else if (r.getStructure().getReactantNumber() == 2 && r.getStructure().getProductNumber() == 2 && !r.hasMultipleKinetics())
	        			duplicates.add(r);
	        		else if (r.getStructure().getReactantNumber() == 2 && r.getStructure().getProductNumber() == 2 && r.hasMultipleKinetics())
	        			nonPDepList.add(r);
	        	}
	        }
	        
	        LinkedList all = new LinkedList();
	        all.add(0,nonPDepList);
	        all.add(1,pDepList);
	        all.add(2,duplicates);
	        return all;
	        //#]
	    }
	 
    //## operation getReactionNumber() 
    public int getReactionNumber() {
        //#[ operation getReactionNumber() 
        return core.getReactionNumber();
        //#]
    }
    
	public int getUnreactedReactionSetIncludingReverseSize() {
		Iterator iter = edge.getReaction();
		int count = 0;
		while (iter.hasNext()){
			Reaction r = (Reaction) iter.next();
			if (r.isBackward()) {
				count = count + 2;
				
			}
			else 
				count++;
		}
		return count;
	}
	
    //## operation getReactionSet() 
    public LinkedHashSet getReactionSet() {
        //#[ operation getReactionSet() 
        return getReactedReactionSet();
        //#]
    }
    
    //## operation getSpecies() 
    public Iterator getSpecies() {
        //#[ operation getSpecies() 
        return core.getSpecies();
        //#]
    }
    
    //## operation getSpeciesNumber() 
    public int getSpeciesNumber() {
        //#[ operation getSpeciesNumber() 
        return core.getSpeciesNumber();
        //#]
    }
    
    //## operation getSpeciesSet() 
    public LinkedHashSet getSpeciesSet() {
        //#[ operation getSpeciesSet() 
        return getReactedSpeciesSet();
        //#]
    }
    
    //## operation getUnreactedReactionSet() 
    public LinkedHashSet getUnreactedReactionSet() {
        //#[ operation getUnreactedReactionSet() 
        return getEdge().getReactionSet();
        
        
        //#]
    }

    public void removeFromUnreactedReactionSet(Reaction rxn) {
        boolean success = edge.reaction.remove(rxn);
	//if (!success){
	//    System.out.println("Pruning debugging line: "+ rxn);
	//}
        return;
    }
    
    //## operation getUnreactedSpeciesSet() 
    public LinkedHashSet getUnreactedSpeciesSet() {
        //#[ operation getUnreactedSpeciesSet() 
        return getEdge().getSpeciesSet();
        
        
        //#]
    }
    
    //## operation isCoreEdgeConsistent() 
    public boolean isCoreEdgeConsistent() {
        //#[ operation isCoreEdgeConsistent() 
        return AbstractReactionModel.isDisjoint(core,edge);
        
        
        
        //#]
    }
    
    //## operation isEmpty() 
    public boolean isEmpty() {
        //#[ operation isEmpty() 
        return core.isEmpty();
        //#]
    }
    
    public boolean isEmpty(FinishController fc) {
        //#[ operation isEmpty() 
        return core.isEmpty(fc);
        //#]
    }
    
    //## operation isReactedReaction(Reaction) 
    public boolean isReactedReaction(Reaction p_reaction) {
        //#[ operation isReactedReaction(Reaction) 
        return (categorizeReaction(p_reaction) == 1) ;
        //#]
    }
    
    //## operation isUnreactedReaction(Reaction) 
    public boolean isUnreactedReaction(Reaction p_reaction) {
        //#[ operation isUnreactedReaction(Reaction) 
        return (categorizeReaction(p_reaction) == -1) ;
        //#]
    }
    
    /*
     * isUnreactedReversiblePathReaction() method:
     * 
     * In a pdepnetwork, suppose we have A+B=C (where A and B are core
     * 	species and C is an edge species).  If the leak flux from this
     * 	network is the largest "edge" flux in the reaction model, the
     * 	reactions of C will be explored.  Imagine C=D and C=E+F are two
     * 	such reactions formed from C.  If we try to add C=D (or C=E+F)
     * 	to the unreacted reaction set, RMG will throw an error.
     * 
     * The reason for this is that both hypothetical reactions are "edge"
     * 	species going to "edge" species, which would never occur in the
     * 	normaly RMG algorithm.
     * 
     * The categorizeReaction() method will return 0
     */
    public boolean isUnreactedReversiblePathReaction(Reaction p_reaction) {
    	int categorizeForwardRxn = categorizeReaction(p_reaction);
    	int categorizeReverseRxn = categorizeReaction(p_reaction.getReverseReaction());
    	if (categorizeForwardRxn == 0 && categorizeReverseRxn == 0) {
    		/*
    		 *  Loop over the reactants and label all non-"Reacted" and 
    		 *  	non-"Unreacted" species as "Unreacted" 
    		 */    		
            Iterator iter = p_reaction.getStructure().getReactants();
            while (iter.hasNext()) {
            	Species spe = (Species)iter.next();
            	if (!contains(spe)) {
            		// new unreacted species
            		addUnreactedSpecies(spe);
            	}
            }
            /*
             *	Loop over the products 
             */
            iter = p_reaction.getStructure().getProducts();
            while (iter.hasNext()) {
            	Species spe = (Species)iter.next();
            	if (!contains(spe)) {
            		// new unreacted species
            		addUnreactedSpecies(spe);
            	}
            }
    		return true;
    	}
    	else return false;
    }
    
    public boolean isUnreactedIrreversiblePathReaction(Reaction p_reaction) {
    	int categorizeForwardRxn = categorizeReaction(p_reaction);
    	if (categorizeForwardRxn == 0) {
    		/*
    		 *  Loop over the reactants and label all non-"Reacted" and 
    		 *  	non-"Unreacted" species as "Unreacted" 
    		 */    		
            Iterator iter = p_reaction.getStructure().getReactants();
            while (iter.hasNext()) {
            	Species spe = (Species)iter.next();
            	if (!contains(spe)) {
            		// new unreacted species
            		addUnreactedSpecies(spe);
            	}
            }
            /*
             *	Loop over the products 
             */
            iter = p_reaction.getStructure().getProducts();
            while (iter.hasNext()) {
            	Species spe = (Species)iter.next();
            	if (!contains(spe)) {
            		// new unreacted species
            		addUnreactedSpecies(spe);
            	}
            }
    		return true;
    	}
    	else return false;
    }
    
    //## operation mergeBasedOnReactedSpecies(CoreEdgeReactionModel,CoreEdgeReactionModel) 
    public static CoreEdgeReactionModel mergeBasedOnReactedSpecies(CoreEdgeReactionModel p_cerm1, CoreEdgeReactionModel p_cerm2) {
        //#[ operation mergeBasedOnReactedSpecies(CoreEdgeReactionModel,CoreEdgeReactionModel) 
        // notice!!! 
        // here p_cerm1 is the primary one, which means every reacted thing in it will be guaranteed showing in the merge result
        // while p_cerm2 is the secondary one, which means if the reacted thing in it is already in p_cerm1, it wont be in the merged one
    	LinkedHashSet rs1 = p_cerm1.getReactedSpeciesSet();
    	LinkedHashSet us1 = p_cerm1.getUnreactedSpeciesSet();    
    	LinkedHashSet rr1 = p_cerm1.getReactedReactionSet();
    	LinkedHashSet ur1 = p_cerm1.getUnreactedReactionSet(); 
        
    	LinkedHashSet rs2 = p_cerm2.getReactedSpeciesSet();
    	LinkedHashSet us2 = p_cerm2.getUnreactedSpeciesSet();    
    	LinkedHashSet rr2 = p_cerm2.getReactedReactionSet();
    	LinkedHashSet ur2 = p_cerm2.getUnreactedReactionSet(); 
        
    	LinkedHashSet rs = new LinkedHashSet(rs1);
        rs.addAll(rs2);
        LinkedHashSet us = new LinkedHashSet(us1);
        us.addAll(us2);
        
        // remove the one possibly in both reacted and unreacted species set from unreacted set
        for (Iterator iter = us.iterator(); iter.hasNext(); ) {
        	Species spe = (Species)iter.next();
        	if (rs.contains(spe)) {
        		iter.remove();
        	}
        }
        
        CoreEdgeReactionModel result = new CoreEdgeReactionModel();
        result.getCore().setSpeciesSet(rs);
        result.getEdge().setSpeciesSet(us);
        
        result.addReactionSet(rr1);
        result.addReactionSet(ur1);
        result.addReactionSet(rr2);
        result.addReactionSet(ur2);
        
        return result;
        
        
        //#]
    }
    
    //## operation moveFromUnreactedToReactedReaction() 
    public void moveFromUnreactedToReactedReaction() {
        //#[ operation moveFromUnreactedToReactedReaction() 
    	LinkedHashSet ur = getUnreactedReactionSet();
        Iterator iter = ur.iterator();
        while (iter.hasNext()) {
        	Reaction r = (Reaction)iter.next();
        	if (categorizeReaction(r) == 1) {
        		addReactedReaction(r);
                        //also add the reverse reaction (added by gmagoon 7/28/09: this should address potential issues where forward and reverse kinetics do not agree due to a non-deterministic identification of a different group for fused cyclic species and reaction families like birad recombination or intraRaddEndocyclic, for which the size of the ring is not always identified the same (e.g. one time, it might be 5-membered ring and another it might be 6-membered ring)); testing of this change suggests that there is a small effect on the numbers for the "official" hexadiene test case, but it appears that this is only due to the kinetics being written to the ODEsolver input file in a different order
                        if (r.hasReverseReaction()){
                            Reaction reverse = (Reaction)r.getReverseReaction();
                            if (categorizeReaction(reverse) == 1) addReactedReaction(reverse);
                        }
        		iter.remove();
        	}
        }
        
        return;
        //#]
    }
    
    //## operation moveFromUnreactedToReactedSpecies(Species) 
    public void moveFromUnreactedToReactedSpecies(Species p_species) {
        //#[ operation moveFromUnreactedToReactedSpecies(Species) 
        boolean rs = containsAsReactedSpecies(p_species);
        boolean us = containsAsUnreactedSpecies(p_species);
        
        if (rs && !us) return;
        else if (!rs && us) {
        	getUnreactedSpeciesSet().remove(p_species);
        	getReactedSpeciesSet().add(p_species);
        	return;
        }
        else throw new InvalidCoreEdgeRelationException(p_species.getName());
        //#]
    }
    
    //## operation printPDepModel(Temperature) 
    //10/25/07 gmagoon; changed to also take pressure
    public void printPDepModel(Temperature p_temperature,Pressure p_pressure) {
        //#[ operation printPDepModel(Temperature) 
        String modelInformation ="";
	System.out.print("This model include totally " + String.valueOf(getSpeciesNumber()) + " Species and ");
        System.out.println(String.valueOf(getReactionNumber()) + " Reactions.");
       
	//System.out.println("Species Set:");
        //System.out.println("Totally " + String.valueOf(getSpeciesNumber()) + " Species:");
        /*LinkedList sortedSpeList = new LinkedList();
        for (Iterator iter = getSpecies(); iter.hasNext(); ) {
       	Species spe = (Species)iter.next();
       	int id = spe.getID();
       	boolean added = false;
       	if (sortedSpeList.isEmpty()) sortedSpeList.add(spe);
       	else {
     		for (int i = 0; i<sortedSpeList.size(); i++) {
        			Species thisSpe = (Species)sortedSpeList.get(i);
        			if (thisSpe.getID()>id) {
        				sortedSpeList.add(i, spe);
        				added = true;
        				break;
        			}
        		}
        		if (!added) sortedSpeList.add(spe);
        	}
        }
        
        for (int i=0; i<sortedSpeList.size(); i++) {
        	Species spe = (Species)sortedSpeList.get(i);
        	System.out.println(spe.toStringWithoutH());
        }
        
        System.out.println("Thermo Properties:");
        System.out.println("SpeciesID\tNamw\tH298\tS298\tCp300\tCp400\tCp500\tCp600\tCp800\tCp1000\tCp1500");
        for (int i=0; i<sortedSpeList.size(); i++) {
        	Species spe = (Species)sortedSpeList.get(i);
        	System.out.println(String.valueOf(spe.getID()) + '\t' + spe.getName() + '\t' + spe.getThermoData().toString());
        }*/
        
        LinkedList nonPDepList = new LinkedList();
        LinkedList pDepList = new LinkedList();
       
        HashSet pDepStructureSet = new HashSet();
        for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
        	PDepNetwork pdn = (PDepNetwork)iter.next();
        	for (Iterator pdniter = pdn.getNetReactions().listIterator(); pdniter.hasNext();) {
        		PDepReaction rxn = (PDepReaction) pdniter.next();
        		if (isReactedReaction(rxn)) {
        			pDepList.add(rxn);
        			pDepStructureSet.add(rxn.getStructure());
        		}
        	}
        }
        
        for (Iterator iter = getReactionSet().iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	Structure s = r.getStructure();
        	if (!pDepStructureSet.contains(s)) {
        		nonPDepList.add(r);
        	} 
        }
        
        System.out.println("//non p_dep reactions:");
        for (Iterator iter = nonPDepList.iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	double rate = r.calculateTotalRate(p_temperature);
            if (r instanceof TemplateReaction) rate = ((TemplateReaction)r).calculateTotalPDepRate(p_temperature, p_pressure);//10/25/07 gmagoon: added pressure
            //System.out.println(r.toString()+"\t rate = \t"+ String.valueOf(rate));
			System.out.println(r.toChemkinString(p_temperature));
		//System.out.println(r.toChemkinString(Global.temperature));//10/25/07 gmagoon eliminating use of Global.temperature
        }
        
        System.out.println("//p_dep reactions:");
        for (Iterator iter = pDepList.iterator(); iter.hasNext(); ) {
        	PDepReaction r = (PDepReaction) iter.next();
            //System.out.println(r.getStructure().toString() + "\t rate = \t" + Double.toString(r.getRate()));
                System.out.println(r.toChemkinString(p_temperature));
		//System.out.println(r.toChemkinString(Global.temperature));//10/25/07 gmagoon eliminating use of Global.temperature
        }
        System.out.println("/////////////////////////////");
        return;
        //#]
    }

    
//	## operation printPDepModel(Temperature) 
    public String returnPDepModel(SystemSnapshot p_ss) {
		Temperature p_temperature = p_ss.temperature;
        //#[ operation printPDepModel(Temperature) 
        String modelInformation ="";
		modelInformation = modelInformation + "This model include totally " + String.valueOf(getSpeciesNumber()) + " Species and ";
		modelInformation = modelInformation + String.valueOf(getReactionNumber()) + " Reactions.\n";
        
		//System.out.println("Species Set:");
        //System.out.println("Totally " + String.valueOf(getSpeciesNumber()) + " Species:");
        /*LinkedList sortedSpeList = new LinkedList();
        for (Iterator iter = getSpecies(); iter.hasNext(); ) {
        	Species spe = (Species)iter.next();
        	int id = spe.getID();
        	boolean added = false;
        	if (sortedSpeList.isEmpty()) sortedSpeList.add(spe);
        	else {
        		for (int i = 0; i<sortedSpeList.size(); i++) {
        			Species thisSpe = (Species)sortedSpeList.get(i);
        			if (thisSpe.getID()>id) {
        				sortedSpeList.add(i, spe);
        				added = true;
        				break;
        			}
        		}
        		if (!added) sortedSpeList.add(spe);
        	}
        }
        
        for (int i=0; i<sortedSpeList.size(); i++) {
        	Species spe = (Species)sortedSpeList.get(i);
        	System.out.println(spe.toStringWithoutH());
        }
        
        System.out.println("Thermo Properties:");
        System.out.println("SpeciesID\tNamw\tH298\tS298\tCp300\tCp400\tCp500\tCp600\tCp800\tCp1000\tCp1500");
        for (int i=0; i<sortedSpeList.size(); i++) {
        	Species spe = (Species)sortedSpeList.get(i);
        	System.out.println(String.valueOf(spe.getID()) + '\t' + spe.getName() + '\t' + spe.getThermoData().toString());
        }*/
        
        LinkedList nonPDepList = new LinkedList();
        LinkedList pDepList = new LinkedList();
        
        HashSet pDepStructureSet = new HashSet();
        for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
        	PDepNetwork pdn = (PDepNetwork)iter.next();
        	for (Iterator pdniter = pdn.getNetReactions().listIterator(); pdniter.hasNext();) {
        		PDepReaction rxn = (PDepReaction) pdniter.next();
        		if (isReactedReaction(rxn)) {
        			pDepList.add(rxn);
        			pDepStructureSet.add(rxn.getStructure());
        		}
        	}
        }
        
        for (Iterator iter = getReactionSet().iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	Structure s = r.getStructure();
        	if (!pDepStructureSet.contains(s)) {
        		nonPDepList.add(r);
        	} 
        }
        
		modelInformation = modelInformation + "//non p_dep reactions:\n";
        for (Iterator iter = nonPDepList.iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	double rate = r.calculateTotalRate(p_temperature);
            if (r instanceof TemplateReaction) rate = ((TemplateReaction)r).calculateTotalPDepRate(p_temperature, p_ss.pressure);//10/25/07 gmagoon: added pressure
			if (r instanceof ThirdBodyReaction) rate = ((ThirdBodyReaction)r).calculateRate(p_ss);
            //System.out.println(r.toString()+"\t rate = \t"+ String.valueOf(rate));
		modelInformation = modelInformation + r.toChemkinString(p_temperature)+"\t"+rate+"\n";
		//modelInformation = modelInformation + r.toChemkinString(Global.temperature)+"\t"+rate+"\n";//10/25/07 gmagoon: eliminating use of Global.temperature
        }
        
		modelInformation = modelInformation + "//p_dep reactions:\n";
        for (Iterator iter = pDepList.iterator(); iter.hasNext(); ) {
        	PDepReaction r = (PDepReaction)iter.next();
            //System.out.println(r.getStructure().toString() + "\t rate = \t" + Double.toString(r.getRate()));
                        modelInformation = modelInformation + r.toChemkinString(p_temperature)+"\n";
			//modelInformation = modelInformation + r.toChemkinString(Global.temperature)+"\n";//10/25/07 gmagoon: eliminating use of Global.temperature
        }
		modelInformation = modelInformation + "/////////////////////////////";
        return modelInformation;
        //#]
    }
    
	
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return true;
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String s = "Model Core:\n";
        s = s + core.toString();
        s = s + "Model Edge:\n";
        s = s + edge.toString();
        
        return s;
        
        
        //#]
    }
    
    public Core getCore() {
        return core;
    }
    
    public Core newCore() {
        core = new Core();
        return core;
    }
    
    public void deleteCore() {
        core=null;
    }
    
    public Edge getEdge() {
        return edge;
    }
    
    public Edge newEdge() {
        edge = new Edge();
        return edge;
    }
    
    public void deleteEdge() {
        edge=null;
    }
    
    protected void initRelations() {
        core = newCore();
        edge = newEdge();
    }
	
	public int getMaxSpeciesID() {
		int maxID = 0;
		for (Iterator iter = core.getSpecies(); iter.hasNext(); ) {
			Species species = (Species) iter.next();
			if (species.getID() > maxID)
				maxID = species.getID();
		}
		for (Iterator iter = edge.getSpecies(); iter.hasNext(); ) {
			Species species = (Species) iter.next();
			if (species.getID() > maxID)
				maxID = species.getID();
		}
		return maxID;
	}
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\CoreEdgeReactionModel.java
*********************************************************************/


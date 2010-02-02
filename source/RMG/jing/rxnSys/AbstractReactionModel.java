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


import jing.param.Global;
import jing.rxn.*;
import jing.chem.*;
import java.util.*;
import jing.mathTool.*;
import jing.chem.Species;
import jing.rxn.Reaction;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\AbstractReactionModel.java                                                                  
//----------------------------------------------------------------------------

//## class AbstractReactionModel 
public class AbstractReactionModel implements ReactionModel {
    
    protected LinkedHashSet reaction;
    protected ReactionModelEnlarger reactionModelEnlarger;
    protected LinkedHashSet species;
    
    // Constructors
    
    //## operation AbstractReactionModel(ReactionModelEnlarger) 
    public  AbstractReactionModel(ReactionModelEnlarger p_reactionModelEnlarger) {
        {
            species=new LinkedHashSet();
        }
        {
            reaction=new LinkedHashSet();
        }
        //#[ operation AbstractReactionModel(ReactionModelEnlarger) 
        reactionModelEnlarger = p_reactionModelEnlarger;
        
        
        //#]
    }
    //## operation AbstractReactionModel(HashSet) 
    public  AbstractReactionModel(LinkedHashSet p_reactionSet) {
        {
            species=new LinkedHashSet();
        }
        {
            reaction=new LinkedHashSet();
        }
        //#[ operation AbstractReactionModel(HashSet) 
        reaction = p_reactionSet;
        
        Iterator iter = p_reactionSet.iterator();
        while (iter.hasNext()) {
        	Reaction rxn = (Reaction)iter.next();
        	addSpeciesFromReaction(rxn);
        }
        
        
        
        //#]
    }
    public  AbstractReactionModel() {
        {
            reaction=new LinkedHashSet();
        }
        {
            species=new LinkedHashSet();
        }
    }
    
    //## operation addReaction(Reaction) 
    public void addReaction(Reaction p_reaction) {
        //#[ operation addReaction(Reaction) 
        addReactionAndSpecies(p_reaction);
        //#]
    }
    
    /**
    Requires:
    Effects: add pass-in reaction into reaction set of this reaction model, meanwhile, if any of the reactants and products is not in speices set of this reaction model, add the species in.
    Modifies: this.reaction, this.species
    */
    //## operation addReactionAndSpecies(Reaction) 
    protected void addReactionAndSpecies(Reaction p_reaction) {
        //#[ operation addReactionAndSpecies(Reaction) 
        reaction.add(p_reaction);
        addSpeciesFromReaction(p_reaction);
        
        
        //#]
    }
    
    //## operation addSpecies(Species) 
    protected void addSpecies(Species p_species) {
        //#[ operation addSpecies(Species) 
        species.add(p_species);
        //#]
    }
    
    /**
    Requires:
    Effects: add any of the reactants and products in the pass-in reaction, if it is not in the species set, into the species set of this reaction model
    Modifies: this.species
    */
    //## operation addSpeciesFromReaction(Reaction) 
    protected void addSpeciesFromReaction(Reaction p_reaction) {
        //#[ operation addSpeciesFromReaction(Reaction) 
        Iterator r_iter = p_reaction.getReactants();
        while (r_iter.hasNext()) {
        	ChemGraph cg = (ChemGraph)r_iter.next();
        	Species spe = cg.getSpecies();
        	addSpecies(spe);
        }
        r_iter = p_reaction.getProducts();
        while (r_iter.hasNext()) {
        	ChemGraph cg = (ChemGraph)r_iter.next();
        	Species spe = cg.getSpecies();
        	addSpecies(spe);
        }
        
        
        //#]
    }
    
    //## operation clearSpecies() 
    protected void clearSpecies() {
        //#[ operation clearSpecies() 
        species.clear();
        //#]
    }
    
    //## operation contains(Species) 
    public boolean contains(Species p_species) {
        //#[ operation contains(Species) 
        return (species.contains(p_species));
        
        
        //#]
    }
    
    //## operation contains(Reaction) 
    public boolean contains(Reaction p_reaction) {
        //#[ operation contains(Reaction) 
        return (reaction.contains(p_reaction));
        //#]
    }
    
    //## operation getReactionNumber() 
    public int getReactionNumber() {
        //#[ operation getReactionNumber() 
        return reaction.size();
        //#]
    }
    
    //## operation getReactionSet() 
    public LinkedHashSet getReactionSet() {
        //#[ operation getReactionSet() 
        return reaction;
        //#]
    }
    
    //## operation getSpeciesNumber() 
    public int getSpeciesNumber() {
        //#[ operation getSpeciesNumber() 
        return species.size();
        //#]
    }
    
    //## operation getSpeciesSet() 
    public LinkedHashSet getSpeciesSet() {
        //#[ operation getSpeciesSet() 
        return species;
        //#]
    }
    
    //## operation isDisjoint(AbstractReactionModel,AbstractReactionModel) 
    public static boolean isDisjoint(AbstractReactionModel p_model1, AbstractReactionModel p_model2) {
        //#[ operation isDisjoint(AbstractReactionModel,AbstractReactionModel) 
    	LinkedHashSet s1 = p_model1.getSpeciesSet();
    	LinkedHashSet s2 = p_model2.getSpeciesSet();    
    	LinkedHashSet r1 = p_model1.getReactionSet();
    	LinkedHashSet r2 = p_model2.getReactionSet(); 
        
        return (MathTool.isCollectionDisjoint(s1,s2) && MathTool.isCollectionDisjoint(s1,s2));
        
        
        
        //#]
    }
    
    //## operation isEmpty() 
    public boolean isEmpty() {
        //#[ operation isEmpty() 
        return reaction.size()==0;
        //#]
    }
    
    public boolean isEmpty(FinishController fc){
    	boolean impSpeciesPresent = false;
    	if (fc.terminationTester instanceof ConversionTT){
    		Species sp = ((SpeciesConversion)((ConversionTT)fc.terminationTester).speciesGoalConversionSet.get(0)).species;
    		for (Iterator iter = reaction.iterator(); iter.hasNext();){
    			Reaction r = (Reaction)iter.next();
                        if (r.isForward() && r.calculateKeq(fc.getReactionSystem().getPresentTemperature()) <= 0.1 && r.getStructure().containsAsProducts(sp)){//10/26/07 gmagoon: changed to avoid use of Global.temperature; I am assuming finish controller has access to an up-to-date copy of reaction system
    	//		if (r.isForward() && r.calculateKeq(Global.temperature) <= 0.1 && r.getStructure().containsAsProducts(sp)){
    				impSpeciesPresent = true;
    				break;
    			}
    		}
    		return (reaction.size() == 0 || !impSpeciesPresent);
    	}
    	return reaction.size() == 0;
    }
    
    //## operation outputModelByTemplates() 
    public void outputModelByTemplates() {
        //#[ operation outputModelByTemplates() 
        String s = "Species Set:\n";
        s = s + "Total of " + String.valueOf(getSpeciesNumber()) + " Species:\n";
        LinkedList sortedSpeList = new LinkedList();
        for (Iterator iter = species.iterator(); iter.hasNext(); ) {
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
        	s = s + spe.toStringWithoutH() + '\n';
        }
         
        s = s + "\nReaction Set:\n";
        s = s + "Total of " + String.valueOf(getReactionNumber()) + " Reactions:\n";
        for (Iterator iter = reaction.iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	s = s + r.toString() + '\n';
        }
        return;
        //#]
    }
    
    //## operation printModel() 
    public void printModel() {
        //#[ operation printModel() 
        System.out.print("This model include a total of " + String.valueOf(getSpeciesNumber()) + " Species and ");
        System.out.println(String.valueOf(getReactionNumber()) + " Reactions.");
        System.out.println("Species Set:");
        System.out.println("Total of " + String.valueOf(getSpeciesNumber()) + " Species:");
        LinkedList sortedSpeList = new LinkedList();
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
         
        System.out.println("\nReaction Set:");
        System.out.println("Total of " + String.valueOf(getReactionNumber()) + " Reactions:");
        for (Iterator iter = getReaction(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	System.out.println(r);
        }
        return;	
        
        
        //#]
    }
    
    //## operation removeReaction(Reaction) 
    public void removeReaction(Reaction p_reaction) {
        //#[ operation removeReaction(Reaction) 
        reaction.remove(p_reaction);
        //#]
    }
    
    //## operation removeSpecies(Species) 
    protected void removeSpecies(Species p_Species) {
        //#[ operation removeSpecies(Species) 
        species.remove(p_Species);
        //#]
    }
    
    //## operation setReactionSet(HashSet) 
    public void setReactionSet(LinkedHashSet p_reactionSet) {
        //#[ operation setReactionSet(HashSet) 
        reaction = p_reactionSet;
        //#]
    }
    
    //## operation setSpeciesSet(HashSet) 
    public void setSpeciesSet(LinkedHashSet p_speciesSet) {
        //#[ operation setSpeciesSet(HashSet) 
        species = p_speciesSet;
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String s = "This model include a total of " + String.valueOf(getSpeciesNumber()) + " Species and ";
        s = s + String.valueOf(getReactionNumber()) + " Reactions.\n";
        s = s + "Species Set:\n";
        s = s + "Total of " + String.valueOf(getSpeciesNumber()) + " Species:\n";
        LinkedList sortedSpeList = new LinkedList();
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
        	s = s + spe.toStringWithoutH() + '\n';
        }
         
        s = s + "\nReaction Set:\n";
        s = s + "Total of " + String.valueOf(getReactionNumber()) + " Reactions:\n";
        for (Iterator iter = getReaction(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	s = s + r.toString() + '\n';
        }
        return s;	
        
        
        //#]
    }
    
    public Iterator getReaction() {
        Iterator iter=reaction.iterator();
        return iter;
    }
    
    public void clearReaction() {
        reaction.clear();
    }
    
    public ReactionModelEnlarger getReactionModelEnlarger() {
        return reactionModelEnlarger;
    }
    
    public void setReactionModelEnlarger(ReactionModelEnlarger p_ReactionModelEnlarger) {
        reactionModelEnlarger = p_ReactionModelEnlarger;
    }
    
    public Iterator getSpecies() {
        Iterator iter=species.iterator();
        return iter;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\AbstractReactionModel.java
*********************************************************************/


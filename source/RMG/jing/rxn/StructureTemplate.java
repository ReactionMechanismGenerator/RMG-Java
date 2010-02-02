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
import jing.chemUtil.*;
import jing.chem.Matchable;
import jing.chem.ChemGraph;
import jing.chemUtil.HierarchyTree;
import jing.param.Global;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\StructureTemplate.java                                                                  
//----------------------------------------------------------------------------

/**
This is the basic description of a flexible reaction structure.  
A (+ B) -> C (+ D)
in which, A, B, C, and D are not fixed chemical structures; instead, they are some functionalgroup, for example: R-H, R., R-COOH, etc.  (Here, R site can be changed, also, although not shown in those examples, the bond could vary, too.)
So, structure template actually defines a group of reactions.
Each reaction template could have a bunch of structure templates.  for each structure template, there will be one or more kinetics corresponding to it.


*/
//## class StructureTemplate 
public class StructureTemplate {
    
    protected LinkedList reactants = new LinkedList();		//## attribute reactants 
    
    protected ArrayList reactantTree;
    
    // Constructors
    
    //## operation StructureTemplate(Matchable,Matchable) 
    public  StructureTemplate(Matchable p_reactant1, Matchable p_reactant2) {
        {
            reactantTree=new ArrayList();
        }
        //#[ operation StructureTemplate(Matchable,Matchable) 
        reactants.add(p_reactant1);
        if (p_reactant2 != null) reactants.add(p_reactant2);
        //#]
    }
    //## operation StructureTemplate(Matchable) 
    public  StructureTemplate(Matchable p_reactant) {
        {
            reactantTree=new ArrayList();
        }
        //#[ operation StructureTemplate(Matchable) 
        reactants.add(p_reactant);
        //#]
    }
    //## operation StructureTemplate(LinkedList) 
    public  StructureTemplate(LinkedList p_reactants) {
        {
            reactantTree=new ArrayList();
        }
        //#[ operation StructureTemplate(LinkedList) 
        reactants = p_reactants;
        //#]
    }
    public  StructureTemplate() {
        {
            reactantTree=new ArrayList();
        }
    }
    
    //## operation generateProduct(ReactionAdjList) 
    public void generateProduct(ReactionAdjList p_reactionAdjList) {
        //#[ operation generateProduct(ReactionAdjList) 
        //#]
    }
    
    /**
    Requires:
    Effects: reverse this structure template by exchanging the reactants structure and products structure.  return the reversed structure template.  
    Modifies:
    */
    //## operation generateReverse(ReactionAdjList) 
    public StructureTemplate generateReverse(ReactionAdjList p_reactionAdjList) {
        //#[ operation generateReverse(ReactionAdjList) 
        int r_num = getReactantNumber();
        Matchable r1 = getAllowedFunctionalGroupAt(1);
        LinkedList reactant = new LinkedList();
        reactant.add(r1);
        if (r_num == 2) {
        	Matchable r2 = getAllowedFunctionalGroupAt(2);
        	reactant.add(r2);
        }
        
        LinkedList reverse = p_reactionAdjList.reactFunctionalGroup(reactant);
        
        StructureTemplate reverseRT = new StructureTemplate(reverse);
        reverseRT.reactantTree=this.reactantTree;
        
        return reverseRT;
        
        
        
        //#]
    }
    
    //## operation getAllowedFunctionalGroupAt(int) 
    public Matchable getAllowedFunctionalGroupAt(int p_index) {
        //#[ operation getAllowedFunctionalGroupAt(int) 
        Matchable t = (Matchable)(reactants.get(p_index-1));
        if (t == null) return null;
        return t;
        //#]
    }
    
    /**
    Requires: the central sites of each of the pass-in reactants have been set correctly; otherwise, no matching can be returned.  This is used when we search a kinetics.
    Effects: compare the central sites of the reactants and the reactants tree to find out the match leaf, and return the matched leaf as a collection.
    Modifies:
    */
    //## operation getMatchedFunctionalGroup(LinkedList) 
    public LinkedList getMatchedFunctionalGroup(LinkedList p_reactants) {
        //#[ operation getMatchedFunctionalGroup(LinkedList) 
        LinkedList fgCollection = new LinkedList();
        boolean found = false;
        Iterator r_iter = p_reactants.iterator();
        while (r_iter.hasNext()) {
        	found = false;
			Object o = r_iter.next();
			ChemGraph cg = null;
			if (o instanceof Species)
				cg = ((Species)r_iter.next()).getChemGraph();
			else
				cg = (ChemGraph)o;
        	 
        	Iterator t_iter = reactantTree.iterator();
        	while (t_iter.hasNext()) {
        		HierarchyTree t = (HierarchyTree)t_iter.next();
        		Stack s = t.findMatchedPath(cg); 
        	    if (s != null && !s.isEmpty() ) {
        	    	found = true;
        	    	fgCollection.add(s);
        	    }
        	}
        	if (!found) {
        		//System.out.println("can't find matched path: " + cg.toString());
        		//System.exit(0);
				return null;
        	}
        }
        
        return fgCollection;
        //#]
    }
    
    //## operation getReactantNumber() 
    public int getReactantNumber() {
        //#[ operation getReactantNumber() 
        return reactants.size();
        //#]
    }
    
    //## operation getReactantTreeNumber() 
    public int getReactantTreeNumber() {
        //#[ operation getReactantTreeNumber() 
        return reactantTree.size();
        //#]
    }
    
    //## operation getReactants() 
    public Iterator getReactants() {
        //#[ operation getReactants() 
        Iterator iter = reactants.iterator();
        return iter;
        //#]
    }
    
    /**
    Requires: 
    Effects: find out all the possible reacted sites of the pass-in reactants according to this structure template.
    Modifies: 
    */
    //## operation identifyReactedSites(ChemGraph,int) 
    public LinkedHashSet identifyReactedSites(ChemGraph p_reactant, int p_position) {
//    	Global.identifyReactiveSitesCount++;
		double pT = System.currentTimeMillis();
        //#[ operation identifyReactedSites(ChemGraph,int) 
        Matchable allowed = getAllowedFunctionalGroupAt(p_position);
        if (allowed == null) return null;
        
        LinkedHashSet hS =  p_reactant.identifyReactionMatchedSite(allowed);
		double t = (System.currentTimeMillis()-pT)/1000/60;
        Global.RT_identifyReactedSites += t;
		return hS;
        
        
        //#]
    }
    
    //## operation isReverse(StructureTemplate,ReactionAdjList) 
    public boolean isReverse(StructureTemplate p_structureTemplate, ReactionAdjList p_reactionAdjList) {
        //#[ operation isReverse(StructureTemplate,ReactionAdjList) 
        return (p_structureTemplate.equals(generateReverse(p_reactionAdjList)));
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return true;
        //#]
    }
    
    //## operation setReactant(int,HashSet) 
    public void setReactant(int p_position, LinkedHashSet p_reactantFGSet) {
        //#[ operation setReactant(int,HashSet) 
        reactants.add(p_position, p_reactantFGSet);
        //#]
    }
    
    //## operation setReactantTree(HashSet) 
    public void setReactantTree(LinkedHashSet p_treeSet) {
        //#[ operation setReactantTree(HashSet) 
        if (p_treeSet == null) throw new InvalidReactantTreeException();
        int size = p_treeSet.size();
        if (size == 0) throw new InvalidReactantTreeException();
        
        Iterator iter = p_treeSet.iterator();
        while (iter.hasNext()) {
        	HierarchyTree tree = (HierarchyTree)iter.next();
        	addReactantTree(tree);
        }
        
        return;
        //#]
    }
    
    public void setReactants(LinkedList p_reactants) {
        reactants = p_reactants;
    }
    
    public ListIterator getReactantTree() {
        ListIterator iter=reactantTree.listIterator();
        return iter;
    }
    
    public void addReactantTree(HierarchyTree p_HierarchyTree) {
        reactantTree.add(p_HierarchyTree);
    }
    
    public void removeReactantTree(HierarchyTree p_HierarchyTree) {
        reactantTree.remove(p_HierarchyTree);
    }
    
    public void clearReactantTree() {
        reactantTree.clear();
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\StructureTemplate.java
*********************************************************************/


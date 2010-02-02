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
import jing.chemUtil.Arc;
import jing.chemUtil.Graph;
import jing.param.Global;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ReactionAdjList.java                                                                  
//----------------------------------------------------------------------------

//## class ReactionAdjList 
public class ReactionAdjList {
    
    protected LinkedList actions = new LinkedList();		//## attribute actions 
    
    protected int productNumber;		//## attribute productNumber 
    
    protected int reactantNumber;		//## attribute reactantNumber 
    
    
    // Constructors
    
    //## operation ReactionAdjList(int,int) 
    public  ReactionAdjList(int p_reactantNumber, int p_productNumber) {
        //#[ operation ReactionAdjList(int,int) 
        reactantNumber = p_reactantNumber;
        productNumber = p_productNumber;
        
        
        
        //#]
    }
    public  ReactionAdjList() {
    }
    
    //## operation addAction(Action) 
    public void addAction(Action p_action) {
        //#[ operation addAction(Action) 
        if (actions == null) actions = new LinkedList();
        actions.add(p_action);
        //#]
    }
    
    //## operation convertBenzene(Arc) 
    public void convertBenzene(Arc p_arc) {
        //#[ operation convertBenzene(Arc) 
        //#]
    }
    
    //## operation generateReverse() 
    public ReactionAdjList generateReverse() {
        //#[ operation generateReverse() 
        ReactionAdjList reverse = new ReactionAdjList(productNumber,reactantNumber);
        Iterator iter = actions.iterator();
        while (iter.hasNext()) {
        	Action action = (Action)iter.next();
        	Action newAction = action.generateReverse();
        	reverse.addAction(newAction);
        }
        	
        return reverse;
        //#]
    }
    
    //## operation getActions() 
    public Iterator getActions() {
        //#[ operation getActions() 
        Iterator iter = actions.iterator();
        return iter;
        //#]
    }
    
    //## operation mutate(Graph) 
    public void mutate(Graph p_graph) throws InvalidActionException {
        //#[ operation mutate(Graph) 
        try {
        	Iterator act_iter = getActions();
        	LinkedHashSet changedAtom = new LinkedHashSet();
        	while (act_iter.hasNext()) {
        		Action act = (Action)act_iter.next();
        		switch (act.type) {
        			case Action.CHANGE_BOND:
        				{
        				// locate the atoms linked by changed bond
        				Iterator iter = act.getSite();
        		        Integer a1 = (Integer)iter.next();
        		        Node n1 = p_graph.getCentralNodeAt(a1.intValue());
        		        Integer a2 = (Integer)iter.next();
        		        Node n2 = p_graph.getCentralNodeAt(a2.intValue());
        		        // locate the changed bond
        		        Integer changedOrder = (Integer)act.getElement();
        		        Arc arc = p_graph.getArcBetween(n1,n2);
        		        Object oldBond = arc.getElement();
        		        if (oldBond instanceof Collection) {
        		        	LinkedHashSet b = new LinkedHashSet();
        		        	Iterator bond_iter = ((Collection)oldBond).iterator();
        		        	while (bond_iter.hasNext()) {
        		        		Bond thisBond = (Bond)bond_iter.next();
        		        		b.add(thisBond.changeBond(changedOrder.intValue()));
        		        	}
        		        	arc.setElement(b);
        		        }
        		        else if (oldBond instanceof Bond) {
        		        	Bond b = ((Bond)oldBond).changeBond(changedOrder.intValue());
        		        	arc.setElement(b);
        		        }
        		        else {
        		        	throw new InvalidBondException("change bond");
        		        }
        		        
        		        Node leftNode = null;
        		        if (n1.includeFgElementInChemNodeElement(FGElement.make("Cdd"))) {
        		        	Iterator a_iter = n1.getNeighbor();
        		        	while (a_iter.hasNext()) {
        		        		Arc a = (Arc)a_iter.next();
        		        		if (a != arc) leftNode = n1.getOtherNode(a);
        		        	}
        		        }
        		        n1.changeChemNodeElement(changedOrder.intValue(), leftNode);
        		        
        		        if (n2.includeFgElementInChemNodeElement(FGElement.make("Cdd"))) {
        		        	Iterator a_iter = n2.getNeighbor();
        		        	while (a_iter.hasNext()) {
        		        		Arc a = (Arc)a_iter.next();
        		        		if (a != arc) leftNode = n2.getOtherNode(a);
        		        	}
        		        }
        		        n2.changeChemNodeElement(changedOrder.intValue(),leftNode);
        		        
        		        changedAtom.add(n1);
        		        changedAtom.add(n2);
        		        
        		        break;
        		    	}
        			case Action.BREAK_BOND:
        				{
        				// locate the atoms linked by changed bond
        				Iterator iter = act.getSite();
        				Integer a1 = (Integer)iter.next();
        				Node n1 = p_graph.getCentralNodeAt(a1.intValue());
        				Integer a2 = (Integer)iter.next();
        				Node n2 = p_graph.getCentralNodeAt(a2.intValue());
        				// locate the broken bond
        				Bond b = (Bond)act.getElement();
        				Arc arc = p_graph.getArcBetween(n1,n2);
        				try {
        					if ((Bond)(arc.getElement()) != b) {
        						throw new InvalidActionException("break bond");
        					}
        				}
        				catch (ClassCastException e) {
        					throw new InvalidBondException("break bond");
        				}
        				
        				// break the bond element
        				p_graph.removeArc(arc);
        
        		        changedAtom.add(n1);
        		        changedAtom.add(n2);
        
        				break;
        			    }
        			case Action.FORM_BOND:
        				{
        				// locate the atoms linked by changed bond
        				Iterator iter = act.getSite();
        				Integer a1 = (Integer)iter.next();
        				Node n1 = p_graph.getCentralNodeAt(a1.intValue());
        				Integer a2 = (Integer)iter.next();
        				Node n2 = p_graph.getCentralNodeAt(a2.intValue());
        
        				// form the new bond and add it to graph
        				Bond b = (Bond)act.getElement();
        				p_graph.addArcBetween(n1,b,n2);
        
        		        changedAtom.add(n1);
        		        changedAtom.add(n2);
        
        				break;
        			    }
        			case Action.GAIN_RADICAL:
        				{
        				// locate the atom
        				Iterator iter = act.getSite();
        				Integer a = (Integer)iter.next();
        				Node n = p_graph.getCentralNodeAt(a.intValue());
        
        				// get the number of radical gained
        				Integer r = (Integer)act.getElement();
        				String spin = null;
        				if (iter.hasNext()) spin = (String)iter.next();
                        
        				// new a radical or set the new radical	order
        				Object oldAtom = n.getElement();
        				if (oldAtom instanceof Collection) {
        		        	LinkedHashSet atom = new LinkedHashSet();
        		        	Iterator atom_iter = ((Collection)oldAtom).iterator();
        		        	while (atom_iter.hasNext()) {
        		        		ChemNodeElement thisAtom = (ChemNodeElement)atom_iter.next();
        		        		atom.add(thisAtom.changeRadical(r.intValue(),spin));
        		        	}
        		        	n.setElement(atom);
        		        }
        		        else if (oldAtom instanceof ChemNodeElement) {
        		        	ChemNodeElement atom = ((ChemNodeElement)oldAtom).changeRadical(r.intValue(),spin);
        		        	n.setElement(atom);
        		        }
        		        else {
        		        	throw new InvalidActionException();
        		        }
        		        
                        changedAtom.add(n);
        
        				break;
        				}
        	        case Action.LOSE_RADICAL:
        	        	{
        	        	// locate the atom
        	        	Iterator iter = act.getSite();
        	        	Integer a = (Integer)iter.next();
        	        	Node n = p_graph.getCentralNodeAt(a.intValue());
        
        	        	// get the number of radical lost
        	        	Integer r = (Integer)act.getElement();
        	        	String spin = null;
        	        	if (iter.hasNext()) spin = (String)iter.next();
        
        	        	// new a radical or set the new radical	order
        				Object oldAtom = n.getElement();
        				if (oldAtom instanceof Collection) {
        		        	LinkedHashSet atom = new LinkedHashSet();
        		        	Iterator atom_iter = ((Collection)oldAtom).iterator();
        		        	while (atom_iter.hasNext()) {
        		        		ChemNodeElement thisAtom = (ChemNodeElement)atom_iter.next();
        		        		atom.add(thisAtom.changeRadical(-r.intValue(),spin));
        		        	}
        		        	n.setElement(atom);
        		        }
        		        else if (oldAtom instanceof ChemNodeElement) {
        		        	ChemNodeElement atom = ((ChemNodeElement)oldAtom).changeRadical(-r.intValue(),spin);
        		        	n.setElement(atom);
        		        }
        		        else {
        		        	throw new InvalidActionException();
        		        }
        		        
        	        	changedAtom.add(n);
        
        	        	break;
        	        	}
        	    	default:
        	    		throw new InvalidActionException("unknown action");
        		}
        	}
        	Iterator iter = changedAtom.iterator();
        	while (iter.hasNext()) {
        		Node node = (Node)iter.next();
        		node.updateFeElement();
        		node.updateFgElement();
        	}
        	return;
        }
        catch (UnknownSymbolException e) {
        	System.err.println("unknown symbols: " + e.getMessage());
        	System.exit(0);
        }
        //#]
    }
    
    /**
    Requires: p_reactants is an ordered linked list, where at index1 there is reactant1's graph, and at indext2, there is reactant2's graph.  the reactedSite ArrayList is the combined ordered reaction sites.
    Effects: according to the reaction adjlist, react at reaction site to form product's graph.
    Modifies:
    */
    //## operation react(LinkedList) 
    public LinkedList react(LinkedList p_reactants) throws InvalidReactantException, InvalidProductException, InvalidReactantNumberException, InvalidProductNumberException {
        //#[ operation react(LinkedList) 
        if (p_reactants.size() != reactantNumber) {
           	throw new InvalidReactantNumberException();
        }
           
        Graph graph = null;
                
        if (reactantNumber == 1) {
           	Graph g1 = (Graph)p_reactants.get(0);
           	if (g1 == null) {
           		throw new InvalidReactantException();
           	}
           	graph = Graph.copy(g1);
        }
        else if (reactantNumber == 2) {
           	Graph g1 = (Graph)p_reactants.get(0);
           	Graph g2 = (Graph)p_reactants.get(1);
           	if (g1 == null || g2 == null) {
           		throw new InvalidReactantException();
           	}
        	graph = Graph.combine(g1, g2);
        }
        else {
        	throw new InvalidReactantNumberException();
        }
                
        mutate(graph);
              
        if (graph == null) {
           	throw new InvalidProductException();
        }
                
        LinkedList productGraph = graph.partition();
        if (productGraph.size() != productNumber) {
           	throw new InvalidProductNumberException();
        }  
               
        return productGraph;
        
        
        //#]
    }
    
    //## operation reactChemGraph(LinkedList) 
    public LinkedList reactChemGraph(LinkedList p_reactants) throws InvalidChemGraphException, ForbiddenStructureException {
		double pT = System.currentTimeMillis();
        //#[ operation reactChemGraph(LinkedList) 
        LinkedList reactants = new LinkedList();
        LinkedList products = new LinkedList();
        
        for (Iterator iter = p_reactants.iterator(); iter.hasNext(); ) {
			Object o = iter.next();
			Graph g = null;
			if (o instanceof ChemGraph)
				 g = ((ChemGraph)o).getGraph();
			else
				g = ((Species)o).getChemGraph().getGraph();
			reactants.add(g);
        }
        
        LinkedList productGraph = null;
        try {
        	productGraph = react(reactants);
        }
        catch (InvalidProductNumberException e) {
        	throw new InvalidProductNumberException(e.getMessage());
        }
        
        for (Iterator iter = productGraph.iterator(); iter.hasNext(); ) {
           	Graph pg = (Graph)iter.next();
           	/*if (ChemGraph.isForbiddenStructure(pg)) 
           		throw new ForbiddenStructureException(pg.toString());*/
           	String name = null;
		
			ChemGraph pcg = ChemGraph.make(pg,true);
			//Species ps = Species.make(name, pcg);
			products.add(pcg);
		
			
        }
		Global.RT_reactChemGraph += (System.currentTimeMillis()-pT)/1000/60;
        return products;
        
        
        //#]
    }
    
    //## operation reactFunctionalGroup(LinkedList) 
    public LinkedList reactFunctionalGroup(LinkedList p_reactants) {
        //#[ operation reactFunctionalGroup(LinkedList) 
        if (p_reactants.size() != reactantNumber) {
           	throw new InvalidReactantNumberException();
        }
           
        LinkedList productGraph = null;
        LinkedList allProduct = new LinkedList();
        
        if (reactantNumber == 1) {
           	Matchable g1 = (Matchable)p_reactants.get(0);
           	if (g1 == null) {
           		throw new InvalidReactantException();
           	}
        
        	LinkedList reactant = new LinkedList();
           	if (g1 instanceof FunctionalGroupCollection) {
           		Iterator iter = ((FunctionalGroupCollection)g1).getFunctionalGroups();
           		while (iter.hasNext()) {
           			Graph thisGraph = ((FunctionalGroup)iter.next()).getGraph();
           			reactant.add(thisGraph);
           			productGraph = react(reactant);
           			allProduct.addAll(productGraph);
           			reactant.clear();
           		}
           	}
           	else if (g1 instanceof FunctionalGroup) {
           		Graph thisGraph = ((FunctionalGroup)g1).getGraph();
           		reactant.add(thisGraph);
           		productGraph = react(reactant);
           		allProduct.addAll(productGraph);
           	}
           	else {
           		throw new InvalidReactantException();
           	}
        }
        else if (reactantNumber == 2) {
           	Matchable g1 = (Matchable)p_reactants.get(0);
           	Matchable g2 = (Matchable)p_reactants.get(1);
           	if (g1 == null || g2 == null) {
           		throw new InvalidReactantException();
           	}
           	boolean r1fg = (g1 instanceof FunctionalGroup);
           	boolean r2fg = (g2 instanceof FunctionalGroup); 
           	boolean r1fgc = (g1 instanceof FunctionalGroupCollection);
           	boolean r2fgc = (g2 instanceof FunctionalGroupCollection);
           	
           	LinkedList reactant = new LinkedList();
           	if (r1fg && r2fg) {
           		Graph thisGraph1 = ((FunctionalGroup)g1).getGraph();
           		Graph thisGraph2 = ((FunctionalGroup)g2).getGraph();
           		reactant.add(thisGraph1);
           		reactant.add(thisGraph2);
           		allProduct = react(reactant);
           	}
           	else if (r1fg && r2fgc) {
           		Graph thisGraph1 = ((FunctionalGroup)g1).getGraph();
           		Iterator iter = ((FunctionalGroupCollection)g2).getFunctionalGroups();
        		reactant.add(thisGraph1);
           		while (iter.hasNext()) {
           			Graph thisGraph2 = ((FunctionalGroup)iter.next()).getGraph();
           			reactant.add(thisGraph2);
           			productGraph = react(reactant);
           			allProduct.addAll(productGraph);
           			reactant.remove(thisGraph2);
           		}
           	}
           	else if (r1fgc && r2fg) {
           		Graph thisGraph2 = ((FunctionalGroup)g2).getGraph();
           		Iterator iter = ((FunctionalGroupCollection)g1).getFunctionalGroups();
        		reactant.add(thisGraph2);
           		while (iter.hasNext()) {
           			Graph thisGraph1 = ((FunctionalGroup)iter.next()).getGraph();
           			reactant.add(thisGraph1);
           			productGraph = react(reactant);
           			allProduct.addAll(productGraph);
           			reactant.remove(thisGraph1);
           		}
           	}
           	else if (r1fgc && r2fgc) {
           		Iterator iter1 = ((FunctionalGroupCollection)g1).getFunctionalGroups();
           		while (iter1.hasNext()) {
        	   		Graph thisGraph1 = ((FunctionalGroup)iter1.next()).getGraph();
        			reactant.add(thisGraph1);
        	   		Iterator iter2 = ((FunctionalGroupCollection)g2).getFunctionalGroups();
        	   		while (iter2.hasNext()) {
        	   			Graph thisGraph2 = ((FunctionalGroup)iter2.next()).getGraph();
        	   			reactant.add(thisGraph2);
        	   			productGraph = react(reactant);
        	   			allProduct.addAll(productGraph);
        	   			reactant.remove(thisGraph2);
        	   		}
        	   		reactant.remove(thisGraph1);
        	   	}
           	}
        }
        else {
        	throw new InvalidReactantNumberException();
        }
        
        LinkedList p_collection = new LinkedList();
        
        if (productNumber == 1) {
        	FunctionalGroupCollection p1 = new FunctionalGroupCollection();
        	           
        	Iterator p_iter = allProduct.iterator();
        	while (p_iter.hasNext()) {
        	   	Graph pg = (Graph)p_iter.next();
        	   	String name = "";
        	   	FunctionalGroup fg = FunctionalGroup.make(name, pg);
        		p1.addFunctionalGroups(fg);
        	}
        	p_collection.add(p1);
        }
        else if (productNumber == 2) {
        	FunctionalGroupCollection p1 = new FunctionalGroupCollection();
        	FunctionalGroupCollection p2 = new FunctionalGroupCollection();
        	/*int highestCentralID1 = -1;
        	int highestCentralID2 = -1;*/
        	int lowestCentralID1 = 10000;
        	int lowestCentralID2 = 10000;
        	Iterator p_iter = allProduct.iterator(); 
        	while (p_iter.hasNext()) {
        	   	Graph pg = (Graph)p_iter.next();
        	   	String name = "";
        	   	FunctionalGroup fg = FunctionalGroup.make(name, pg);
        	   	/*int present_HCID = pg.getHighestCentralID();
        	   	if (present_HCID <= 0) throw new InvalidCentralIDException();
        	   	if (present_HCID == highestCentralID1) {
        	   		p1.addFunctionalGroups(fg);
        	   	}
        	   	else if (present_HCID == highestCentralID2) {
        	   		p2.addFunctionalGroups(fg);
        	   	}
        	   	else {
        	   		if (highestCentralID1 == -1) {
        	   			p1.addFunctionalGroups(fg);
        	   			highestCentralID1 = present_HCID;
        	   		}
        	   		else if (highestCentralID2 == -1) {
        	   			p2.addFunctionalGroups(fg);
        	   			highestCentralID2 = present_HCID;
        	   		}
        	   		else {
        	   			throw new InvalidCentralIDException();
        	   		}
        	   	}*/
        	   	int present_LCID = pg.getLowestCentralID();
        	   	if (present_LCID <= 0) throw new InvalidCentralIDException();
        	   	if (present_LCID == lowestCentralID1) {
        	   		p1.addFunctionalGroups(fg);
        	   	}
        	   	else if (present_LCID == lowestCentralID2) {
        	   		p2.addFunctionalGroups(fg);
        	   	}
        	   	else {
        	   		if (lowestCentralID1 == 10000) {
        	   			p1.addFunctionalGroups(fg);
        	   			lowestCentralID1 = present_LCID;
        	   		}
        	   		else if (lowestCentralID2 == 10000) {
        	   			p2.addFunctionalGroups(fg);
        	   			lowestCentralID2 = present_LCID;
        	   		}
        	   		else {
        	   			throw new InvalidCentralIDException();
        	   		}
        	   	}
        
        	}                   
        	p_collection.add(p1);
        	p_collection.add(p2);
        }
        
        return p_collection;
        
        
        //#]
    }
    
    public void setActions(LinkedList p_actions) {
        actions = p_actions;
    }
    
    public int getProductNumber() {
        return productNumber;
    }
    
    public void setProductNumber(int p_productNumber) {
        productNumber = p_productNumber;
    }
    
    public int getReactantNumber() {
        return reactantNumber;
    }
    
    public void setReactantNumber(int p_reactantNumber) {
        reactantNumber = p_reactantNumber;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ReactionAdjList.java
*********************************************************************/


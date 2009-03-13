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



package jing.chem;


import java.util.*;
import jing.chemUtil.*;
import jing.chemUtil.Graph;
import jing.chemUtil.GraphComponent;

//## package jing::chem 

//----------------------------------------------------------------------------
// jing\chem\FunctionalGroup.java                                                                  
//----------------------------------------------------------------------------

//## class FunctionalGroup 
public class FunctionalGroup implements Matchable {
    
    /**
    This is the name user defines for a functional group.
    */
    protected String name;		//## attribute name 
    
    protected Graph graph;
    
    protected int C_count = 0;
    protected int O_count = 0;
    protected int rad_count = 0;
    
    // Constructors
    
    //## operation FunctionalGroup() 
    private  FunctionalGroup() {
        //#[ operation FunctionalGroup() 
        //#]
    }
    //## operation FunctionalGroup(String,Graph) 
    private  FunctionalGroup(String p_name, Graph p_graph) {
        //#[ operation FunctionalGroup(String,Graph) 
        name = p_name;
        graph = p_graph;
        //#]
    }
    
    private FunctionalGroup(String p_name, Graph p_graph, int p_rad_count, int p_O_count, int p_C_count) {
        name = p_name;
        graph = p_graph;
        rad_count = p_rad_count;
        O_count = p_O_count;
        C_count = p_C_count;
    }
    
    /**
    Requires:
    Effects: if all the elements in p_atomList are instance of atom, add that at the p_position of this functional group
    Modifies: this.nodeList
    */
    //## operation addAtomAt(int,Collection) 
    public void addAtomAt(int p_position, Collection p_atomList) throws AtomOcuppiedException {
        //#[ operation addAtomAt(int,Collection) 
        Iterator iter = p_atomList.iterator();
        while (iter.hasNext()) {
        	if (!(iter.next() instanceof Atom)) {
        		throw new InvalidAtomListException();
        	}
        }
        graph.addNodeAt(p_position, p_atomList);
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: add a new node at p_position storing p_atom
    Modifies: this.graph.nodeList
    */
    //## operation addAtomAt(int,Atom) 
    public void addAtomAt(int p_position, Atom p_atom) {
        //#[ operation addAtomAt(int,Atom) 
        graph.addNodeAt(p_position, p_atom);
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: add a new bond between p_position1 and p_position2 to store a collection of bonds: p_bondList
    Modifies: this.graph.arcList
    */
    //## operation addBondBetween(int,Collection,int) 
    public void addBondBetween(int p_position1, Collection p_bondList, int p_position2) throws EmptyAtomException {
        //#[ operation addBondBetween(int,Collection,int) 
        Iterator iter = p_bondList.iterator();
        while (iter.hasNext()) {
        	if (!(iter.next() instanceof Bond)) {
        		throw new InvalidAtomListException();
        	}
        }
        
        graph.addArcBetween(p_position1,p_bondList,p_position2);
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: add a new bond between p_position1 and p_position2 to store p_bond
    Modifies: this.graph.arcList
    */
    //## operation addBondBetween(int,Bond,int) 
    public void addBondBetween(int p_position1, Bond p_bond, int p_position2) {
        //#[ operation addBondBetween(int,Bond,int) 
        graph.addArcBetween(p_position1,p_bond,p_position2);
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: get the stored node element at p_position
    Modifies:
    */
    //## operation getAtomAt(int) 
    public Object getAtomAt(int p_position) {
        //#[ operation getAtomAt(int) 
        Node node = (Node)(graph.getNodeAt(p_position));
        return (Collection)(node.getElement());
        //#]
    }
    
    /**
    Requires:
    Effects: return the stroed arc element between p_position1 and p_position2 
    Modifies:
    */
    //## operation getBondBetween(int,int) 
    public Object getBondBetween(int p_position1, int p_position2) {
        //#[ operation getBondBetween(int,int) 
        Arc arc = graph.getArcBetween(p_position1,p_position2);
        return (Collection)(arc.getElement());
        //#]
    }
    
    /**
    Requires:
    Effects: return a HashMap storing this functional group centers.
    Modifies:
    */
    //## operation getCentralNode() 
    protected HashMap getCentralNode() {
        //#[ operation getCentralNode() 
        return getGraph().getCentralNode();
        //#]
    }
    
    /**
    Requires:
    Effects: return the number of central nodes in this functional group.
    Modifies:
    */
    //## operation getCentralNodeNumber() 
    public int getCentralNodeNumber() {
        //#[ operation getCentralNodeNumber() 
        return getGraph().getCentralNodeNumber();
        //#]
    }
    
    //## operation getName() 
    public String getName() {
        //#[ operation getName() 
        return name;
        //#]
    }
    
    /**
    Requires:
    Effects: return this functional group's radical number
    Modifies:
    */
    //## operation getRadicalNumber() 
    public int getRadicalNumber() {
        //#[ operation getRadicalNumber() 
        int radicalNumber = 0;
        Iterator iter = getGraph().getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	radicalNumber += node.getFeElement().getOrder();
        }
        return radicalNumber;
        //#]
    }
    
    //## operation isParticipateInReaction(GraphComponent) 
    public boolean isParticipateInReaction(GraphComponent p_graphComponent) {
        //#[ operation isParticipateInReaction(GraphComponent) 
        if (!getGraph().contains(p_graphComponent)) throw new NotInGraphException();
        
        if (p_graphComponent instanceof Node) {
        	Node n = (Node)p_graphComponent;
        	return (n.getCentralID().intValue()>0);
        }
        else if (p_graphComponent instanceof Arc) {
        	Arc a = (Arc)p_graphComponent;
            for (Iterator i = a.getNeighbor(); i.hasNext(); ) {
            	Node neighbor = (Node)i.next();
            	if (neighbor.getCentralID().intValue()<0) return false;
            }
            return true;
        }
        else throw new InvalidGraphComponentException();
        
        
        
        
        
        //#]
    }
    
    //## operation isRadical() 
    public boolean isRadical() {
        //#[ operation isRadical() 
        return (getRadicalNumber()>0);
        //#]
    }
    
    //## operation isSubAtCentralNodes(FunctionalGroup) 
    public boolean isSubAtCentralNodes(FunctionalGroup p_functionalGroup) {
        //#[ operation isSubAtCentralNodes(FunctionalGroup) 
        return getGraph().isSubAtCentralNodes(p_functionalGroup.getGraph());
        //#]
    }
    
    /**
    Requires:
    Effects: if p_functional is a functional group, return if this functional group is matched with p_functional at the central nodes; if p_functional is a functional group collection, return if this functional group is matched with one of the p_functionals in the collection at the central nodes; all other cases, return false;
    Modifies:
    */
    //## operation isSubAtCentralNodes(Matchable) 
    public boolean isSubAtCentralNodes(Matchable p_functional) {
        //#[ operation isSubAtCentralNodes(Matchable) 
        if (this == p_functional) return false;
        if (p_functional instanceof FunctionalGroup) {
        	FunctionalGroup fg = (FunctionalGroup)p_functional;
        	return isSubAtCentralNodes(fg);
        }
        else if (p_functional instanceof FunctionalGroupCollection) {
        	if (((FunctionalGroupCollection)p_functional).contains(this)) return true;
        	Iterator iter = ((FunctionalGroupCollection)p_functional).getFunctionalGroups();
        	while (iter.hasNext()) {
        		FunctionalGroup fg = (FunctionalGroup)iter.next();
        		if (isSubAtCentralNodes(fg)) return true;
        	}
        	return false;
        }
        else {
        	return false;
        }
        //#]
    }
    
    //## operation make(String,Graph) 
    public static FunctionalGroup make(String p_name, Graph p_graph) {
        //#[ operation make(String,Graph) 
        /*FunctionalGroup fg = functionalGroupDictionary.getFunctionalGroup(p_uniqueID);
        if (fg == null) {
        	fg = new FunctionalGroup(p_name, p_uniqueID, p_graph);
        	functionalGroupDictionary.putFunctionalGroup(fg);
        }
        return fg;*/
        return new FunctionalGroup(p_name, p_graph);
        //#]
    }
    
    public static FunctionalGroup makeForbiddenStructureFG(String p_name, Graph p_graph) {
		// Read in the graph and store each line in its own string
    	String p_string = p_graph.toString();
		String[] indivRows = p_string.split("[\n]");
		
		// Initialize the radical, oxygen, and carbon counter
		int radCount = 0;
		int oCount = 0;
		int cCount = 0;
		
		for (int i=0; i<indivRows.length; i++) {
			StringTokenizer forbiddenST = new StringTokenizer(indivRows[i]);
			String counter = forbiddenST.nextToken();	// Skip over the counter, e.g. (1)
			
			// Determine if atom is carbon or oxygen
			String atom = forbiddenST.nextToken();
			if (atom.toUpperCase().startsWith("C")) ++cCount;
			else if (atom.toUpperCase().startsWith("O")) ++oCount;
			
			// Determine radical number
			String radical = forbiddenST.nextToken();
			if (radical.startsWith("1")) ++radCount;
			else if (radical.startsWith("2")) radCount += 2;
			else if (radical.startsWith("3")) radCount += 3;
			else if (radical.startsWith("4")) radCount += 4;
		}
    	return new FunctionalGroup(p_name, p_graph, radCount, oCount, cCount);
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return getGraph().repOk();
        //#]
    }
    
    protected void setName(String p_name) {
        name = p_name;
    }
    
    public Graph getGraph() {
        return graph;
    }
    
    public void setGraph(Graph p_Graph) {
        graph = p_Graph;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FunctionalGroup.java
*********************************************************************/


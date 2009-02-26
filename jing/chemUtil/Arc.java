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



package jing.chemUtil;
/**
 * Provides classes that are fundamental to the building and mantaining molecules and heirarchy tree structures.
 */

import java.util.*;
import jing.mathTool.*;
import jing.chemParser.*;

//## package jing::chemUtil 

//----------------------------------------------------------------------------
// jing\chemUtil\Arc.java                                                                  
//----------------------------------------------------------------------------

/**
 * Arc is the graph component connecting two nodes in a graph.  it can store special information defined by user. Arc's neighbor should be node, and each arc has a neighbor collection with the size equal to or less than 2 since the arc can at most connect two nodes.
 */
//## class Arc 
public class Arc extends GraphComponent {
    
    //private boolean backArc = false;		//## attribute backArc 
    
    
    // Constructors
    
    /**
     * constructor to store the pass-in element into Arc and put all the contents in the pass-in collection if the collection size is less than or equal to 2.  if the size is greater than 2, throw InvalidNeighborSizeException(). 
     */
    //## operation Arc(Object,Collection) 
    public  Arc(Object p_element, Collection p_neighbor) throws InvalidNeighborException {
        //#[ operation Arc(Object,Collection) 
        super(p_element);
        
        if (p_neighbor.size() > 2) throw new InvalidNeighborSizeException();
        
        Iterator iter = p_neighbor.iterator();
        while (iter.hasNext()) {
          try {
          	Node node = (Node)iter.next();
          	addNeighbor(node);
          }
          catch (ClassCastException e) {
            throw new InvalidNeighborException("Arc-Arc");
          }
        }
        
        p_neighbor = null;
        //#]
    }
    /**
    Store the pass-in element into this Arc
    
    */
    //## operation Arc(Object) 
    public  Arc(Object p_element) {
        //#[ operation Arc(Object) 
        super(p_element);
        //#]
    }
    public  Arc() {
    }
    
    /**
    If the pass-in GraphComponent is an instance of Node, add the pass-in GraphComponent into this.neighbor collection.   Otherwise, throw InvalidNeighborException.<br><br>
    <b>Modifies:</b><br> this.neighbor
    */
    //## operation addNeighbor(GraphComponent) 
    public void addNeighbor(GraphComponent p_graphComponent) throws InvalidNeighborException {
        //#[ operation addNeighbor(GraphComponent) 
        if (!(p_graphComponent instanceof Node)) {
        	throw new InvalidNeighborException();
        }
        super.addNeighbor(p_graphComponent);
        
        
        
        //#]
    }
    
    /**
    add the pass-in Node into this.neighbor collection.<br><br>
    <b>Modifies:</b><br> this.neighbor
    */
    //## operation addNeighbor(Node) 
    public void addNeighbor(Node p_node) {
        //#[ operation addNeighbor(Node) 
        super.addNeighbor(p_node);
        
        
        
        //#]
    }
    
    /**
     * This method is not recursive and compares only the arcs and not the neighbors. 
     * <li>Returns false if the passed in GraphCompnent is not an instance of Arc.</li>
     * <li>Returns true if the passed in Arc is the same as the this.</li>
     * <li>Returns true if the passed in Arc has a set of possible bonds in its elements (eq: {S,D,T}) and is a superset of this.element.</li>
     */
    //## operation contentSub(GraphComponent) 
    public boolean contentSub(GraphComponent p_graphComponent) {
        //#[ operation contentSub(GraphComponent) 
        if (this == p_graphComponent) return false;
        if (!(p_graphComponent instanceof Arc)) return false;
        
        Arc arc = (Arc)p_graphComponent;
        
        Object o1 = getElement();
        Object o2 = arc.getElement();
        
        return MathTool.isSub(o1, o2);
        //#]
    }
    
    /**
    If this arc connects p_node and another node, return that node; otherwise, return null;
    
    */
    //## operation getOtherNode(Node) 
    public Node getOtherNode(Node p_node) {
        //#[ operation getOtherNode(Node) 
        if (!neighborOk()) {
        	return null;
        }
        else {
        	Iterator iter = getNeighbor();
        	Node node1 = (Node)iter.next();
        	Node node2 = (Node)iter.next();	
        
        	if (node1 == p_node) return node2;
        	else if (node2 == p_node) return node1;
        	else return null;
        }
        //#]
    }
    
    
    
    /**
    If pass-in GraphComponent is an instance of Node, return super.isConnected(); otherwise, return false; since arc can't connect to arc.
    
    */
    //## operation isConnected(GraphComponent) 
    public boolean isConnected(GraphComponent p_graphComponent) {
        //#[ operation isConnected(GraphComponent) 
        if (!(p_graphComponent instanceof Node)) return false;
        else return super.isConnected(p_graphComponent);
        //#]
    }
    
    /**
    Return false since only node can be leaf.
    
    */
    //## operation isLeaf() 
    public boolean isLeaf() {
        //#[ operation isLeaf() 
        // an arc can't be a leaf of a graph
        return false;
        //#]
    }
    
    /**
     * Adds the nodes p_node1 and p_node2 as neighbors and then checks if the number of neighbors are greater than 2.
     * @param p_node1
     * @param p_node2
     */
    //## operation link(Node,Node) 
    public void link(Node p_node1, Node p_node2) {
        //#[ operation link(Node,Node) 
        addNeighbor(p_node1);
        addNeighbor(p_node2);
        repOk();
        return;
        //#]
    }
    
    /**
    Check if neighbor of this arc is okay.<br>
    <li>the neighbor collection size should be equal to 2.  if it is not 2, return false</li>
    <li>both the two neighbors should be nodes. if any of them is not a node, return false</li>
    
    */
    //## operation neighborOk() 
    private boolean neighborOk() {
        //#[ operation neighborOk() 
        // check if this arc has two neighbors
        if (neighbor.size() != 2) return false;
        
        // check if all the neighbors are Node
        Iterator iter = getNeighbor();
        while (iter.hasNext()) {
        	if (!(iter.next() instanceof Node)) return false;
        }
        
        return true;
        //#]
    }
    
   
    
    /**
    Check if the rep is ok in the following aspects:
    (1) neighborOk()?
    Modifies:
    */
    //## operation repOk() 
    private boolean repOk() {
        //#[ operation repOk() 
        return (neighborOk());
        //#]
    }
    
    /**
     * Writes the name of the bond eg: S, D, T etc.
     */
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        return ChemParser.writeBond(getElement());
 
    }
    
   
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\Arc.java
*********************************************************************/


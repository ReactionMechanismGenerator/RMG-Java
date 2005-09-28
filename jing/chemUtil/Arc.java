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



package jing.chemUtil;


import java.util.*;
import jing.mathTool.*;
import jing.chemParser.*;

//## package jing::chemUtil 

//----------------------------------------------------------------------------
// jing\chemUtil\Arc.java                                                                  
//----------------------------------------------------------------------------

/**
Arc is the graph component connecting two nodes in a graph.  it can store special information defined by user. Arc's neighbor should be node, and each arc has a neighbor collection with the size equal to or less than 2 since the arc can at most connect two nodes.
*/
//## class Arc 
public class Arc extends GraphComponent {
    
    protected boolean backArc = false;		//## attribute backArc 
    
    
    // Constructors
    
    /**
    Requires:
    Effects: constructor to store the pass-in element into Arc and put all the contents in the pass-in collection if the collection size is less than or equal to 2.  if the size is greater than 2, throw InvalidNeighborSizeException().
    Modifies: 
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
    Requires:
    Effects:constrctor.  store the pass-in element into this Arc
    Modifies:
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
    Requires:
    Effects:  if the pass-in GraphComponent is an instance of Node, add the pass-in GraphComponent into this.neighbor collection.   Otherwise, throw InvalidNeighborException.
    Modifies: this.neighbor
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
    Requires:
    Effects: use super's addNeighbor().
    Modifies:
    */
    //## operation addNeighbor(Node) 
    public void addNeighbor(Node p_node) {
        //#[ operation addNeighbor(Node) 
        super.addNeighbor(p_node);
        
        
        
        //#]
    }
    
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
    Requires:
    Effects: if this arc connects p_node and another node, return that node; otherwise, return null;
    Modifies:
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
    Requires:
    Effects: return super.isConnected()
    Modifies:
    */
    //## operation isConnected(Node) 
    public boolean isConnected(Node p_node) {
        //#[ operation isConnected(Node) 
        return super.isConnected(p_node);
        //#]
    }
    
    /**
    Requires:
    Effects: if pass-in GraphComponent is an instance of Node, return super.isConnected(); otherwise, return false; since arc can't connect to arc.
    Modifies:
    */
    //## operation isConnected(GraphComponent) 
    public boolean isConnected(GraphComponent p_graphComponent) {
        //#[ operation isConnected(GraphComponent) 
        if (!(p_graphComponent instanceof Node)) return false;
        else return super.isConnected(p_graphComponent);
        //#]
    }
    
    /**
    Requires:
    Effects: return false since only node can be leaf.
    Modifies:
    */
    //## operation isLeaf() 
    public boolean isLeaf() {
        //#[ operation isLeaf() 
        // an arc can't be a leaf of a graph
        return false;
        //#]
    }
    
    //## operation link(Node,Node) 
    public void link(Node p_node1, Node p_node2) {
        //#[ operation link(Node,Node) 
        addNeighbor(p_node1);
        addNeighbor(p_node2);
        return;
        //#]
    }
    
    /**
    Requires:
    Effects: check if neighbor of this arc is okay.
    (1) the neighbor collection size should be equal to 2.  if it is not 2, return false.
    (2) both the two neighbors should be nodes. if any of them is not a node, return false;
    Modifies:
    */
    //## operation neighborOk() 
    public boolean neighborOk() {
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
    
    //## operation rearrangeBenzene(String) 
    public void rearrangeBenzene(String p_startType) {
        //#[ operation rearrangeBenzene(String) 
        /*Bond bond = (Bond)getElement();
        if (!bond.isBenzene()) return;
        
        if (p_startType.compareToIgnoreCase("D") != 0 && p_startType.compareToIgnoreCase("S") != 0) {
        	throw new BenzeneRearrangeException();
        }
         
        Bond newBond = Bond.make(p_startType);
        setElement(newBond);
        setVisited(true);
        
        Iterator iter = getNeighbor();
        while (iter.hasNext()) {
        	Node n = (Node)iter.next();
        	if (!n.isVisited()) {
        		rearrangeBenzene(this);
        	}
        }
        
        return;
        */
        //#]
    }
    
    /**
    Requires:
    Effects: check if the rep is ok in the following aspects:
    (1) neighborOk()?
    Modifies:
    */
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return (neighborOk());
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        return ChemParser.writeBond(getElement());
        
        
        //#]
    }
    
    public boolean getBackArc() {
        return backArc;
    }
    
    public void setBackArc(boolean p_backArc) {
        backArc = p_backArc;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\Arc.java
*********************************************************************/


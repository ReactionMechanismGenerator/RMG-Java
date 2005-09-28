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


import jing.chem.*;
import java.util.*;
import jing.mathTool.*;

//## package jing::chemUtil 

//----------------------------------------------------------------------------
// jing\chemUtil\GraphComponent.java                                                                  
//----------------------------------------------------------------------------

/**
GraphComponenet is the super abstract class which represents a component in a graph.  The component could be an Arc or a Node, and user can put whatever they like as the element stroed in the component.  GraphComponent also save the information about one-step connectivity, which is the nearest neighbor of this GraphComponent.
*/
//## class GraphComponent 
abstract public class GraphComponent {
    
    /**
    The content this graph component stores.
    */
    protected Object element = null;		//## attribute element 
    
    /**
    A flag indicating if this graph component is in a cycle.
    */
    protected boolean inCycle = false;		//## attribute inCycle 
    
    /**
    Storing the correspondingly matched graph component when we comparing two graphs.
    */
    protected GraphComponent matchedGC = null;		//## attribute matchedGC 
    
    /**
    A flag indicating if this graph component has been visited when we are doing traversal/search in a graph.
    */
    protected boolean visited = false;		//## attribute visited 
    
    protected HashSet neighbor;
    
    // Constructors
    
    /**
    Requires:
    Effects: construct a new GraphComponent and store the pass-in object as the element of this GraphComponent.
    Modifies: 
    */
    //## operation GraphComponent(Object) 
    public  GraphComponent(Object p_element) {
        {
            neighbor=new HashSet();
        }
        //#[ operation GraphComponent(Object) 
        element = p_element;
        
        
        
        //#]
    }
    public  GraphComponent() {
        {
            neighbor=new HashSet();
        }
    }
    
    /**
    Requires:
    Effects: if this and p_GraphComponent is not the same instance, add the pass-in p_GraphComponent as the neighbor this, add this as the neighbor of p_GraphComponent, too.
    if this and the pass-in GraphComponents are the same object, throw InvalidNeighborException since we forbid one to be its own neighbor.
    Notice: this function has a good side-effect: it not only add p_GraphComponent as the neighbor of this, but also the reversed direction.  so if addNeighbor(p_GraphComponent) is called, it actually add two-way relation between this and p_GraphComponent.
    Modifies: this.neighbor, p_GraphComponent.neighbor.
    */
    //## operation addNeighbor(GraphComponent) 
    public void addNeighbor(GraphComponent p_GraphComponent) throws InvalidNeighborException {
        //#[ operation addNeighbor(GraphComponent) 
        // one component can't add itself as neighbor
        if (this == p_GraphComponent) {
        	throw new InvalidNeighborException("Add itself as neighbor!");
        }
        
        if (!this.neighbor.contains(p_GraphComponent)) {
        	neighbor.add(p_GraphComponent);
        }
        
        if (!p_GraphComponent.neighbor.contains(this)) {
        	p_GraphComponent.neighbor.add(this);
        }
        
        return;
        //#]
    }
    
    /**
    Requires:
    Effects: clear neighbor set.
    Modifies: neighbor set will be set empty.
    */
    //## operation clearNeighbor() 
    public void clearNeighbor() {
        //#[ operation clearNeighbor() 
        // clear
        neighbor.clear();
        //#]
    }
    
    //## operation contentSub(GraphComponent) 
    public abstract boolean contentSub(GraphComponent p_graphComponent);
    
    /**
    Requires:
    Effects: return the number of the neighbor of this graph component.
    Modifies:
    */
    //## operation getNeighborNumber() 
    public int getNeighborNumber() {
        //#[ operation getNeighborNumber() 
        return neighbor.size();
        //#]
    }
    
    /**
    Requires:
    Effects: identifiy and return all possible matched patterns starting from this GraphComponent and p_graphComponent.
    Modifies: visited status and matchedGC of nodes and arcs.
    */
    //## operation identifyAllMatchedSites(GraphComponent) 
    public HashSet identifyAllMatchedSites(GraphComponent p_graphComponent) {
        //#[ operation identifyAllMatchedSites(GraphComponent) 
        if (this == p_graphComponent) return null;
        
        // compare the this gc with the p_gc for one step
        if (!contentSub(p_graphComponent)) return null;
        
        // compare the neighbor
        Collection c1 = neighbor;
        Collection c2 = p_graphComponent.neighbor;
        int positionNum = c2.size();
        int candidateNum = c1.size();
        if (positionNum > candidateNum)	return null;
        
        boolean found = false;
        
        this.setMatchedGC(p_graphComponent);
        p_graphComponent.setMatchedGC(this);
        
        LinkedList matchedList = new LinkedList();
        Iterator iter2 = c2.iterator();
        while (iter2.hasNext()) {
        	found = false;
        	GraphComponent co2 = (GraphComponent)iter2.next();
        	// if a neighbor in c2 has not been visited, compare it with the neighbors not visited in c1
        	HashSet matchedAtThisPosition = new HashSet();
        
        	GraphComponent matched = co2.getMatchedGC();
        	if (matched == null) {
        		Iterator iter1 = c1.iterator();
        		GraphComponent co1 = null;
        		while (iter1.hasNext()) {
        			co1 = (GraphComponent)iter1.next();
        			if (co1.getMatchedGC() == null && co2.getMatchedGC() == null) {
        				HashSet ms = co1.identifyAllMatchedSites(co2);
        				if (ms != null) {
        					found = true;
        					matchedAtThisPosition.add(ms);
        				}
        			}
        		}
        		if (!found) {
        			this.setMatchedGC(null);
        			p_graphComponent.setMatchedGC(null);
        
        			return null;
        		}
        		else {
        			matchedList.add(matchedAtThisPosition);
        		}
        	}
        	else {
        		if (!c1.contains(matched)) {
        			this.setMatchedGC(null);
        			p_graphComponent.setMatchedGC(null);
        
         			return null;
        		}
        	}
        }
        
        this.setMatchedGC(null);
        p_graphComponent.setMatchedGC(null);
        
        HashSet result = new HashSet();
        
        if (matchedList.isEmpty()) {
        	if (this instanceof Node && p_graphComponent instanceof Node) {
        		MatchedSite ms = new MatchedSite();
        		Integer cid = ((Node)p_graphComponent).getCentralID();
        		if (!((Node)p_graphComponent).isCentralNode()) {
        			if (!ms.addPeriphery((Node)this)) return null;
        		}
        		else {
        			if (!ms.putCenter(cid, (Node)this)) return null;
        		}
        		result.add(ms);
        	}
        	else {
        		return null;
        	}
        }
        else {
        	int matchedSize = matchedList.size();
        	Collection comb = MathTool.expandDisjointly(matchedList.iterator());
        
        	if (comb.isEmpty()) return null;
        	for (Iterator iter = comb.iterator(); iter.hasNext(); ) {
        		Collection match = (Collection)iter.next();
        		if (match.size() != matchedSize) return null;
        		Collection oneMatchSet = MathTool.expand(match.iterator());
        		for (Iterator oneMatchIter = oneMatchSet.iterator(); oneMatchIter.hasNext();) {
        			Collection oneMatch = (Collection)oneMatchIter.next();
        			MatchedSite ms = MatchedSite.union(oneMatch);
        
        			if (ms != null) {
        				boolean add = true;
        				if (this instanceof Node && p_graphComponent instanceof Node) {
        					Integer cid = ((Node)p_graphComponent).getCentralID();
        					if (!((Node)p_graphComponent).isCentralNode()) {
        						if (!ms.addPeriphery((Node)this)) add = false;
        					}
        					else {
        						if (!ms.putCenter(cid, (Node)this)) add = false;
        					}
        				}
        				if (add) result.add(ms);
        			}
        		}
        	}
        }
        
        return result;
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: check if this and p_GraphComponent are neighbors.  
    (0) if (this == p_GraphComponent), return false;
    (1) if neighbor of this contains p_GraphComponent and neighbor of p_GraphComponent contains this, return true;
    (2) if neighbor of this contains p_GraphComponent and neighbor of p_GraphComponent does not contain this,  return false;
    (3) if neighbor of this does not contain p_GraphComponent and neighbor of p_GraphComponent contains this, return false;
    (4) if neighbor of this does not contain p_GraphComponent and neighbor of p_GraphComponent does not contain this, return false;
    Modifies: 
    */
    //## operation isConnected(GraphComponent) 
    public boolean isConnected(GraphComponent p_GraphComponent) {
        //#[ operation isConnected(GraphComponent) 
        // a component can't connect to itself
        if (this == p_GraphComponent) return false;
        
        return (neighbor.contains(p_GraphComponent) && p_GraphComponent.neighbor.contains(this));
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: check if this GraphComponent and the pass-in GraphComponent are logically equivalent, which means those two GraphComponents should store the same element and also have equivalent neighbors.
    Modifies: visited and matchedGC
    */
    //## operation isEquivalent(GraphComponent,Stack,Stack) 
    public boolean isEquivalent(GraphComponent p_graphComponent, Stack p_stack1, Stack p_stack2) {
        //#[ operation isEquivalent(GraphComponent,Stack,Stack) 
        if (this == p_graphComponent) return true;
        
        Object o1 = element;
        Object o2 = p_graphComponent.element;
        Collection c1 = neighbor;
        Collection c2 = p_graphComponent.neighbor;
        boolean found = false;
        
        // compare the element
        Collection co1,co2;
        // if o1 and o2 are both collection, we judge if o2 contains all elements in o1
        if ((o1 instanceof Collection) && (o2 instanceof Collection)) {
        	co1 = (Collection)o1;
        	co2 = (Collection)o2;
        	if ((co1.size() == co2.size()) && (co2.containsAll(co1))) {
        		found = true;
        	}
        }
        // if only o2 is collection, we judge if o2 contains o1 as the only element
        else if (!(o1 instanceof Collection) && (o2 instanceof Collection)) {
        	co2 = (Collection)o2;
        	if ((co2.size() == 1) && (co2.contains(o1))) {
        		found = true;
        	}
        }
        // if only o1 is collection, we judge if o1 contains o2 as the only element
        else if ((o1 instanceof Collection) && (!(o2 instanceof Collection))) {
        	co1 = (Collection)o1;
        	if ((co1.size() == 1) && (co1.contains(o2))) {
        		found = true;
        	}
        }
        // if o1 and o2 are both not collection, we judge if o1 equals o2
        else if (!(o1 instanceof Collection) && !(o2 instanceof Collection)) {
        	if (o2.equals(o1)) {
        		found = true;
        	}
        }
        if (found) {
        	setMatchedGC(p_graphComponent);
        	p_graphComponent.setMatchedGC(this);
        	p_stack1.push(this);
        	p_stack2.push(p_graphComponent);
        }
        else {
        	setMatchedGC(null);
        	p_graphComponent.setMatchedGC(null);
        	return false;
        }
        
        // compare the neighbor
        Iterator iter1 = c1.iterator();
        while (iter1.hasNext()) {
        	found = false;
        	GraphComponent gc1 = (GraphComponent)iter1.next();
        	GraphComponent matched = gc1.getMatchedGC();
        	// if a neighbor in c1 has not been visited, compare it with the neighbors not visited in c2
        	if (matched == null) {
        		Iterator iter2 = c2.iterator();
        		GraphComponent gc2 = null;
        		while (iter2.hasNext()) {
        			gc2 = (GraphComponent)iter2.next();
        			// if a non-visited neighbor in c2 has been found to be the father of this neighbor in c1
        			// set found to true and set visited for both two GC as true, break
        			if (gc1.getMatchedGC() == null && gc2.getMatchedGC() == null) {
        				if (gc1.isEquivalent(gc2,p_stack1,p_stack2)) {
        					found = true;
        					break;
        				}
        			}
        		}
        		if (!found) {
        			resetStack(p_stack1,this);
        			resetStack(p_stack2,p_graphComponent);
        			return false;
        		}
        	}
        	else {
        		if (!c2.contains(matched) || matched.getMatchedGC()!=gc1) {
        			resetStack(p_stack1,this);
        			resetStack(p_stack2,p_graphComponent);
         			return false;
        		}
        	}
        
        }
        
        return true;
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: return true iff this GC is equivalent to p_graphComponent, (checking centerID)
    Modifies: visited and matchedGC status
    */
    //## operation isEquivalentCenterMatched(GraphComponent,Stack,Stack) 
    public boolean isEquivalentCenterMatched(GraphComponent p_graphComponent, Stack p_stack1, Stack p_stack2) {
        //#[ operation isEquivalentCenterMatched(GraphComponent,Stack,Stack) 
        if (this == p_graphComponent) return true;
        
        if ((this instanceof Node) && (p_graphComponent instanceof Node)) {
        	if (((Node)this).getCentralID() != ((Node)p_graphComponent).getCentralID()) return false;
        }
        
        return isEquivalent(p_graphComponent,p_stack1,p_stack2);
        //#]
    }
    
    /**
    Requires: this graphcomponent and the pass-in graph component are not in the same graph.
    Effects: check if this GraphComponents is a SubGraph of the pass-in GraphComponent.  Here, SubGraph means two things:
    (1) the graph structure itself
    (2) the contents stored in the graph
    Notice: Node's central ID doesn't matter
    Modifies: this.visited.
    */
    //## operation isSub(GraphComponent,Stack,Stack) 
    public boolean isSub(GraphComponent p_graphComponent, Stack p_stack1, Stack p_stack2) {
        //#[ operation isSub(GraphComponent,Stack,Stack) 
        if (this == p_graphComponent) return false;
        
        // compare the this gc with the p_gc for one step
        boolean found = contentSub(p_graphComponent);
        if (!found) return false;
        else {
        	this.setMatchedGC(p_graphComponent);
        	p_graphComponent.setMatchedGC(this);
        	p_stack1.push(this);
        	p_stack2.push(p_graphComponent);
        }
        
        // compare the neighbor
        Collection c1 = neighbor;
        Collection c2 = p_graphComponent.neighbor;
        
        int foundNum = 0;
        Iterator iter2 = c2.iterator();
        while (iter2.hasNext()) {
        	foundNum = 0;
        	GraphComponent co2 = (GraphComponent)iter2.next();
        	// if a neighbor in c2 has not been visited, compare it with the neighbors not visited in c1
        	GraphComponent matched = co2.getMatchedGC();
        	if (matched == null) {
        		Iterator iter1 = c1.iterator();
        		GraphComponent co1 = null;
        		while (iter1.hasNext()) {
        			co1 = (GraphComponent)iter1.next();
        			// if a non-visited neighbor in c1 has been found to be the sub of this neighbor in c2
        			// set found to true and set visited for both two GC as true, break
        			if (co1.getMatchedGC() == null) {
        				if (co1.isSub(co2,p_stack1,p_stack2)) {
        					foundNum++;
        				}
        			}
        		}
        		if (foundNum == 0) {
        			resetStack(p_stack1,this);
        			resetStack(p_stack2,p_graphComponent);
        			return false;
        		}
        		else if (foundNum == 1) {
        		}
        
        	}
        	// if a neighbor of has been visited, check if p_graphComponent has a corresponding match neighbor
        	// if there is no match in p_chemGraph's neighbor, this and p_graphComponent are not matched, return false
        	else {
        		if (!c1.contains(matched) || matched.getMatchedGC()!=co2) {
        			resetStack(p_stack1,this);
        			resetStack(p_stack2,p_graphComponent);
         			return false;
        		}
        	}
        }
        
        return true;
        
        
        
        //#]
    }
    
    /**
    Requires: this graphcomponent and the pass-in graph component are not in the same graph.
    Effects: check if this GraphComponents is a SubGraph of the pass-in GraphComponent.  Here, SubGraph means two things:
    (1) the graph structure itself
    (2) the contents stored in the graph
    Notice: the central ID of the matched Nodes connected to this graph component will be set according to the matching grpah component connected to the pass-in graph component.
    Modifies: this.visited.
    */
    //## operation isSubAndSetCentralID(GraphComponent,Stack,Stack) 
    public boolean isSubAndSetCentralID(GraphComponent p_graphComponent, Stack p_stack1, Stack p_stack2) {
        //#[ operation isSubAndSetCentralID(GraphComponent,Stack,Stack) 
        if (this == p_graphComponent) return false;
        
        // compare the this gc with the p_gc for one step
        boolean found = contentSub(p_graphComponent);
        if (!found) return false;
        else {
        	this.setMatchedGC(p_graphComponent);
        	p_graphComponent.setMatchedGC(this);
        	p_stack1.push(this);
        	p_stack2.push(p_graphComponent);
        }
        
        // compare the neighbor
        Collection c1 = neighbor;
        Collection c2 = p_graphComponent.neighbor;
        
        Iterator iter2 = c2.iterator();
        while (iter2.hasNext()) {
        	found = false;
        	GraphComponent co2 = (GraphComponent)iter2.next();
        	// if a neighbor in c2 has not been visited, compare it with the neighbors not visited in c1
        	GraphComponent matched = co2.getMatchedGC();
        	if (matched == null) {
        		Iterator iter1 = c1.iterator();
        		GraphComponent co1 = null;
        		while (iter1.hasNext()) {
        			co1 = (GraphComponent)iter1.next();
        			// if a non-visited neighbor in c1 has been found to be the sub of this neighbor in c2
        			// set found to true and set visited for both two GC as true, break
        			if (co2.getMatchedGC() == null && co1.getMatchedGC() == null) {
        				if (co1.isSubAndSetCentralID(co2,p_stack1,p_stack2)) {
        					found = true;
        					break;
        				}
        			}
        		}
        		if (!found) {
        			resetStack(p_stack1,this);
        			resetStack(p_stack2,p_graphComponent);
        			return false;
        		}
        	}
        	// if a neighbor of has been visited, check if p_graphComponent has a corresponding match neighbor
        	// if there is no match in p_chemGraph's neighbor, this and p_graphComponent are not matched, return false
        	else {
        		if (!c1.contains(matched)) {
        			resetStack(p_stack1,this);
        			resetStack(p_stack2,p_graphComponent);
         			return false;
        		}
        	}
        }
        
        if (this instanceof Node && p_graphComponent instanceof Node) {
        	((Node)this).setCentralID(((Node)p_graphComponent).getCentralID());
        }
        
        return true;
        
        
        
        //#]
    }
    
    /**
    Requires: this graphcomponent and the pass-in graph component are not in the same graph.
    Effects: check if this GraphComponents is a SubGraph of the pass-in GraphComponent.  Here, SubGraph means two things:
    (1) the graph structure itself
    (2) the contents stored in the graph
    Notice: the central ID of the Nodes matching with each other should be matched.
    Modifies: this.visited.
    
    
    
    */
    //## operation isSubCentralMatched(GraphComponent,Stack,Stack) 
    public boolean isSubCentralMatched(GraphComponent p_graphComponent, Stack p_stack1, Stack p_stack2) {
        //#[ operation isSubCentralMatched(GraphComponent,Stack,Stack) 
        if (this == p_graphComponent) return false;
        
        if ((this instanceof Node) && (p_graphComponent instanceof Node)) {
        	if (((Node)this).getCentralID().intValue() != ((Node)p_graphComponent).getCentralID().intValue()) return false;
        }
        
        // compare the this gc with the p_gc for one step
        boolean found = contentSub(p_graphComponent);
        if (!found) return false;
        else {
        	this.setMatchedGC(p_graphComponent);
        	p_graphComponent.setMatchedGC(this);
        	p_stack1.push(this);
        	p_stack2.push(p_graphComponent);
        }
        
        // compare the neighbor
        Collection c1 = neighbor;
        Collection c2 = p_graphComponent.neighbor;
        
        int foundNum = 0;
        Iterator iter2 = c2.iterator();
        while (iter2.hasNext()) {
        	foundNum = 0;
        	GraphComponent co2 = (GraphComponent)iter2.next();
        	// if a neighbor in c2 has not been visited, compare it with the neighbors not visited in c1
        	GraphComponent matched = co2.getMatchedGC();
        	if (matched == null) {
        		Iterator iter1 = c1.iterator();
        		GraphComponent co1 = null;
        		while (iter1.hasNext()) {
        			co1 = (GraphComponent)iter1.next();
        			// if a non-visited neighbor in c1 has been found to be the sub of this neighbor in c2
        			// set found to true and set visited for both two GC as true, break
        			if (co2.getMatchedGC() == null && co1.getMatchedGC() == null) {
                    	if (co1.isSubCentralMatched(co2,p_stack1,p_stack2)) {
           					foundNum++;
        				}
                    }
        		}
        		if (foundNum == 0) {
        			resetStack(p_stack1,this);
        			resetStack(p_stack2,p_graphComponent);
        			return false;
        		}
        		else if (foundNum == 1) {
        		}
        
        	}
        	// if a neighbor of has been visited, check if p_graphComponent has a corresponding match neighbor
        	// if there is no match in p_chemGraph's neighbor, this and p_graphComponent are not matched, return false
        	else {
        		if (!c1.contains(matched) || matched.getMatchedGC()!=co2) {
        			resetStack(p_stack1,this);
        			resetStack(p_stack2,p_graphComponent);
         			return false;
        		}
        	}
        }
        
        return true;
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: return true iff p_gc1 and p_gc are two equivalent neighbors of this GC.
    Modifies:
    */
    //## operation isSymmetric(GraphComponent,GraphComponent) 
    public boolean isSymmetric(GraphComponent p_gc1, GraphComponent p_gc2) {
        //#[ operation isSymmetric(GraphComponent,GraphComponent) 
        if (p_gc1 == p_gc2) return true;
        
        if (!neighbor.contains(p_gc1) || !neighbor.contains(p_gc2)) return false;
        removeNeighbor(p_gc1);
        removeNeighbor(p_gc2);
        
        Stack s1 = new Stack();
        Stack s2 = new Stack();
        if (p_gc1.isEquivalent(p_gc2, s1, s2)) {
        	resetStack(s1);
        	resetStack(s2);
        	addNeighbor(p_gc1);
        	addNeighbor(p_gc2);
        	return true;
        }
        else {
        	resetStack(s1);
        	resetStack(s2);
        	addNeighbor(p_gc1);
        	addNeighbor(p_gc2);
        	return false;
        }
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: return the Visited status of this GraphComponent.
    Modifies:
    */
    //## operation isVisited() 
    public boolean isVisited() {
        //#[ operation isVisited() 
        return visited;
        //#]
    }
    
    /**
    Requires:
    Effects: abstract operation will be implemented by subclass of Node and Ard to check if rep of neighbor are ok.
    Modifies: 
    */
    //## operation neighborOk() 
    public abstract boolean neighborOk();
    
    //## operation printLable() 
    public String printLable() {
        //#[ operation printLable() 
        if (this instanceof Node) return ((Node)this).printLable();
        else if (this instanceof Arc) return ((Arc)this).printLable();
        else return null;
        //#]
    }
    
    /**
    Requires:
    Effects: if this and the pass-in component are connected, delete this connection, i.e., delete each other from their neighbor collection.
    Modifies: p_GraphComponent.neighbor, this.neighbor
    */
    //## operation removeNeighbor(GraphComponent) 
    public void removeNeighbor(GraphComponent p_GraphComponent) {
        //#[ operation removeNeighbor(GraphComponent) 
        if (isConnected(p_GraphComponent)) {
        	neighbor.remove(p_GraphComponent);
        	p_GraphComponent.neighbor.remove(this);
        }
        return;
        //#]
    }
    
    /**
    Requires:
    Effects:abstract operation should be implemented by subclass of Node and Arc to check if right rep is okay.
    Modifies:
    */
    //## operation repOk() 
    public abstract boolean repOk();
    
    //## operation resetStack(Stack,GraphComponent) 
    private void resetStack(Stack p_stack, GraphComponent p_graphComponent) {
        //#[ operation resetStack(Stack,GraphComponent) 
        if (p_stack.empty()) return;
        GraphComponent gc = null;
        do {
        	gc = (GraphComponent)p_stack.pop();
        	gc.setMatchedGC(null);
        } while (!gc.equals(p_graphComponent));
        return;
        //#]
    }
    
    //## operation resetStack(Stack) 
    public void resetStack(Stack p_stack) {
        //#[ operation resetStack(Stack) 
        while (!p_stack.empty()) {
        	GraphComponent gc = (GraphComponent)p_stack.pop();
        	gc.setMatchedGC(null);
        }
        return;
        //#]
    }
    
    //## operation setNeighborVisited(GraphComponent,boolean) 
    public void setNeighborVisited(GraphComponent p_firstNeighbor, boolean p_visit) {
        //#[ operation setNeighborVisited(GraphComponent,boolean) 
        Iterator iter = getNeighbor();
        while (iter.hasNext()) {
        	GraphComponent gc = (GraphComponent)iter.next();
        	if (!gc.equals(p_firstNeighbor)) {
        		gc.setVisited(p_visit);
        		gc.setNeighborVisited(this,p_visit);
        	}
        }
        return;
        //#]
    }
    
    /**
    Requires:
    Effects: retrun the toString() method of the stored element.
    Modifies: 
    */
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        return element.toString();
        //#]
    }
    
    public Object getElement() {
        return element;
    }
    
    public void setElement(Object p_element) {
        element = p_element;
    }
    
    public boolean getInCycle() {
        return inCycle;
    }
    
    public void setInCycle(boolean p_inCycle) {
        inCycle = p_inCycle;
    }
    
    public GraphComponent getMatchedGC() {
        return matchedGC;
    }
    
    public void setMatchedGC(GraphComponent p_matchedGC) {
        matchedGC = p_matchedGC;
    }
    
    public boolean getVisited() {
        return visited;
    }
    
    public void setVisited(boolean p_visited) {
        visited = p_visited;
    }
    
    public Iterator getNeighbor() {
        Iterator iter=neighbor.iterator();
        return iter;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\GraphComponent.java
*********************************************************************/


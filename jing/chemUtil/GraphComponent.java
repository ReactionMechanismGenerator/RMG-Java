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


import jing.chem.*;
import java.util.*;

import jing.mathTool.*;

//## package jing::chemUtil 

//----------------------------------------------------------------------------
// jing\chemUtil\GraphComponent.java                                                                  
//----------------------------------------------------------------------------

/**
* GraphComponenet is the super abstract class which represents a component in a graph.  The component could be an Arc or a Node, and user can put whatever they like as the element stroed in the component.  GraphComponent also save the information about one-step connectivity, which is the nearest neighbor of this GraphComponent.
*/
//## class GraphComponent 
abstract public class GraphComponent {
    
    
    protected Object element = null;		//## attribute element 
    
    
    protected boolean inCycle = false;		//## attribute inCycle 
    
    
    protected GraphComponent matchedGC = null;		//## attribute matchedGC 
    
    
    protected boolean visited = false;		//## attribute visited 
    
    
    protected LinkedHashSet neighbor = null;
    
    protected Integer centralID = new Integer(-1);
    
    // Constructors
    
    /**
     Construct a new GraphComponent and store the pass-in object as the element of this GraphComponent. 
    */
    //## operation GraphComponent(Object) 
    public  GraphComponent(Object p_element) {
        {
			if (this instanceof Arc)
				neighbor = new LinkedHashSet(3,1);
			else if (this instanceof Node) {
				if (p_element instanceof Atom) {
					int valency = (int)(((Atom)p_element).getValency()+1);
					neighbor = new LinkedHashSet(valency+1,1);
				}
				else
					neighbor = new LinkedHashSet(5,1);
			}
			else
				neighbor = new LinkedHashSet();
			
        }
        //#[ operation GraphComponent(Object) 
        element = p_element;
        
        
        
        //#]
    }
    public  GraphComponent() {
        {
			neighbor = new LinkedHashSet();
        }
    }
    
    /**
    If this and p_GraphComponent is not the same instance, add the pass-in p_GraphComponent as the neighbor this, add this as the neighbor of p_GraphComponent, too.
    if this and the pass-in GraphComponents are the same object, throw InvalidNeighborException since we forbid one to be its own neighbor.
    This function has a good side-effect: it not only add p_GraphComponent as the neighbor of this, but also the reversed direction.  so if addNeighbor(p_GraphComponent) is called, it actually add two-way relation between this and p_GraphComponent.
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
    Clear neighbor set. Neighbor set will be set empty.
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
    Return the number of the neighbor of this graph component.
    
    */
    //## operation getNeighborNumber() 
    public int getNeighborNumber() {
        //#[ operation getNeighborNumber() 
        return neighbor.size();
        //#]
    }
    
    /**
    Identifiy and return all possible matched patterns starting from this GraphComponent and p_graphComponent.<br>
    <b>Modifies</b><br>
    visited status and matchedGC of nodes and arcs.
    */
    //## operation identifyAllMatchedSites(GraphComponent) 
    LinkedList identifyAllMatchedSites(GraphComponent p_graphComponent) {
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
        boolean nonVisited = true;
        while (iter2.hasNext()) {
        	
        	GraphComponent co2 = (GraphComponent)iter2.next();
        
        	GraphComponent matched = co2.getMatchedGC();
        	
        	if (matched == null) {
        		nonVisited = false;
        		LinkedHashSet matchedAtThisSite = new LinkedHashSet();
        		Iterator iter1 = c1.iterator();
        		GraphComponent co1 = null;
        		while (iter1.hasNext()) {
        			co1 = (GraphComponent)iter1.next();
        			if (co1.getMatchedGC() == null && co2.getMatchedGC() == null) {
        				LinkedList listOfMatchedsites = co1.identifyAllMatchedSites(co2);
        				if (listOfMatchedsites != null) {
        					matchedAtThisSite.addAll(listOfMatchedsites);
        				}
        			}
        		}
        		if (matchedAtThisSite != null)
        			matchedList.add(matchedAtThisSite);
        		else {
        			this.setMatchedGC(null);
        			p_graphComponent.setMatchedGC(null);
        			return null;
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
        
        if (nonVisited) {
        	Integer cid = (p_graphComponent).getCentralID();
			MatchedSite ms = new MatchedSite();
			if (!(p_graphComponent).isCentralNode()) {
    			if (!ms.putPeriphery(p_graphComponent, this)) return null;
    		}
    		else {
    			if (!ms.putCenter(cid, this)) return null;
    		}
			matchedList.add(ms);
        	return matchedList;
        }
        else{
        	if (matchedList.isEmpty()) return null;
        	else {
        		int index  = 0;
            	MatchedSite matchedSite = new MatchedSite();
            	Integer cid = (p_graphComponent).getCentralID();
        		if (!(p_graphComponent).isCentralNode()) {
        			if (!matchedSite.putPeriphery(p_graphComponent, this)) return null;
        		}
        		else {
        			if (!matchedSite.putCenter(cid, this)) return null;
        		}
            	LinkedList result = uniteNextValidMs(matchedSite, matchedList, index);
            	return result;
        	}
        	
        }
        //#]
    }
    
    private LinkedList uniteNextValidMs(MatchedSite p_ms, LinkedList matchedList, int index) {
    	LinkedList result = new LinkedList();
		LinkedHashSet thisListOfMS = (LinkedHashSet)matchedList.get(index);
		Iterator thisListIterator = thisListOfMS.iterator();
		if (index == matchedList.size()-1){
			while (thisListIterator.hasNext()){
				MatchedSite ms = (MatchedSite)thisListIterator.next();
				MatchedSite mergedMS = MatchedSite.merge(p_ms, ms);
				if (mergedMS != null) result.add(mergedMS);
			}
		}
		else {
			while (thisListIterator.hasNext()){
				MatchedSite ms = (MatchedSite)thisListIterator.next();
				MatchedSite mergedMS = MatchedSite.merge(p_ms, ms);
				LinkedList unitedSoFar = new LinkedList();
				if (mergedMS != null)
					 unitedSoFar = uniteNextValidMs(mergedMS, matchedList, index +1);
				if (unitedSoFar != null) result.addAll(unitedSoFar);
			}
		}
		if (result.isEmpty())
			return null;
		else
			return result;
	}
    
    /* LinkedHashSet identifyAllMatchedSites(GraphComponent p_graphComponent) {
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
        	LinkedHashSet matchedAtThisPosition = new LinkedHashSet();
        
        	GraphComponent matched = co2.getMatchedGC();
        	if (matched == null) {
        		Iterator iter1 = c1.iterator();
        		GraphComponent co1 = null;
        		while (iter1.hasNext()) {
        			co1 = (GraphComponent)iter1.next();
        			if (co1.getMatchedGC() == null && co2.getMatchedGC() == null) {
        				LinkedHashSet ms = co1.identifyAllMatchedSites(co2);
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
        
        LinkedHashSet result = new LinkedHashSet();
        
        if (matchedList.isEmpty()) {
        	
        		MatchedSite ms = new MatchedSite();
        		Integer cid = (p_graphComponent).getCentralID();
        		if (!(p_graphComponent).isCentralNode()) {
        			if (!ms.addPeriphery(this)) return null;
        		}
        		else {
        			if (!ms.putCenter(cid, this)) return null;
        		}
        		result.add(ms);
        	
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
        				
        					Integer cid = (p_graphComponent).getCentralID();
        					if (!(p_graphComponent).isCentralNode()) {
        						if (!ms.addPeriphery(this)) add = false;
        					}
        					else {
        						if (!ms.putCenter(cid, this)) add = false;
        					}
        				if (add) result.add(ms);
        			}
        		}
        	}
        }
        
        return result;
        
        
        //#]
    }*/
    
    
   
	/**
    Check if this and p_GraphComponent are neighbors.<br>  
    <li> if (this == p_GraphComponent), return false </li>
    <li>if neighbor of this contains p_GraphComponent and neighbor of p_GraphComponent contains this, return true</li>
    <li>if neighbor of this contains p_GraphComponent and neighbor of p_GraphComponent does not contain this,  return false</li>
    <li>if neighbor of this does not contain p_GraphComponent and neighbor of p_GraphComponent contains this, return false</li>
    <li>if neighbor of this does not contain p_GraphComponent and neighbor of p_GraphComponent does not contain this, return false</li>
    
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
    Check if this GraphComponent and the pass-in GraphComponent are logically equivalent, which means those two GraphComponents should store the same element and also have equivalent neighbors.<br>
    <b>Modifies</b><br>
    It modifies visited and matchedGC
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
    Return true iff this GC is equivalent to p_graphComponent, also the central ID of the graph components should be the same.<br>
    <b>Modifies</b><br>
    It modifies visited and matchedGC status
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
    This graphcomponent and the pass-in graph component should not be in the same graph.<br>
    Check if this GraphComponents is a SubGraph of the pass-in GraphComponent.  Here, SubGraph means two things:<br>
    <li>the graph structure itself<\li>
    <li> the contents stored in the graph<\li>
    Node's central ID doesn't matter.<br>
    <b>Modifies</b><br>
    It modifies this.visited.
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
    This graphcomponent and the pass-in graph component should not be in the same graph.
    Check if this GraphComponents is a SubGraph of the pass-in GraphComponent.  Here, SubGraph means two things:
    <li>the graph structure itself<\li>
    <li>the contents stored in the graph<\li>
    The central ID of the matched Nodes connected to this graph component will be set according to the matching grpah component connected to the pass-in graph component.<br>
    <b>Modifies</b><br>
    It modifies this.visited.
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
    This graphcomponent and the pass-in graph component are not in the same graph.
    Check if this GraphComponents is a SubGraph of the pass-in GraphComponent.  Here, SubGraph means two things:
    <li>The graph structure itself</li> 
    <li>The contents stored in the graph</li> 
    The central ID of the Nodes matching with each other should be matched.<br>
    <b>Modifies</b><br>
     It modifies this.visited.
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
           					
           					boolean matchedRest = false;
           					matchedRest = isSubCentralMatched(p_graphComponent, p_stack1, p_stack2);
           					if (matchedRest) {
           						foundNum++;
           						co1.setMatchedGC(null); co2.setMatchedGC(null);
           						return true;
           					}
           					else {
           						p_graphComponent.setMatchedGC(this); this.setMatchedGC(p_graphComponent);
           						co1.setMatchedGC(null); co2.setMatchedGC(null);
           						continue;
           					}
        				}
                		
                    }
        		}
        		
        		p_graphComponent.setMatchedGC(null); this.setMatchedGC(null);
        		return false;
        		
        		//else if (foundNum == 1) {
        		//}
        
        	}
        	// if a neighbor of has been visited, check if p_graphComponent has a corresponding match neighbor
        	// if there is no match in p_chemGraph's neighbor, this and p_graphComponent are not matched, return false
        	else {
        		if (!c1.contains(matched) || matched.getMatchedGC()!=co2) {
        			p_graphComponent.setMatchedGC(null); this.setMatchedGC(null);
         			return false;
        		}
        	}
        }
        
        return true;
    	
    
    }
    
    /**
    Return true iff p_gc1 and p_gc are two equivalent neighbors of this GC.
    
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
    Return the Visited status of this GraphComponent.
    
    */
    //## operation isVisited() 
    public boolean isVisited() {
        //#[ operation isVisited() 
        return visited;
        //#]
    }
    
   
    
    //## operation printLable() 
    public String printLable() {
        //#[ operation printLable() 
        if (this instanceof Node) return ((Node)this).printLable();
        else if (this instanceof Arc) return ((Arc)this).printLable();
        else return null;
        //#]
    }
    
    /**
     *If this and the pass-in component are connected, delete this connection, i.e., delete each other from their neighbor collection.<br>
     <b>Modifies</b><br>
     It modifies p_GraphComponent.neighbor, this.neighbor
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
    private void resetStack(Stack p_stack) {
        //#[ operation resetStack(Stack) 
        while (!p_stack.empty()) {
        	GraphComponent gc = (GraphComponent)p_stack.pop();
        	gc.setMatchedGC(null);
        }
        return;
        //#]
    }
    
    //## operation isCentralNode()
    public boolean isCentralNode() {
        //#[ operation isCentralNode()
        return centralID.intValue()>0;
        //#]
    }
    
//  ## operation setCentralID(int)
    public void setCentralID(int p_centralID) {
        //#[ operation setCentralID(int)
        centralID = new Integer(p_centralID);


        //#]
    }
    
    public Integer getCentralID() {
        return centralID;
    }

    public void setCentralID(Integer p_centralID) {
        centralID = p_centralID;
    }

    
    //## operation setNeighborVisited(GraphComponent,boolean) 
    /*public void setNeighborVisited(GraphComponent p_firstNeighbor, boolean p_visit) {
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
    }*/
    
    /**
    Retrun the toString() method of the stored element.
    
    */
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        return element.toString();
        //#]
    }
    
    /**
    Returns the content this graph component stores.
    */
    public Object getElement() {
        return element;
    }
    
    public void setElement(Object p_element) {
        element = p_element;
    }
    
    /**
    Returns true if this graph component is in a cycle.
    */
    public boolean getInCycle() {
        return inCycle;
    }
    
    /**
    Sets the incycle of graph component as true.
    */
    public void setInCycle(boolean p_inCycle) {
    	
        inCycle = p_inCycle;
    }
    
    /**
     * returns the matched graph component. Usually useful when comparing two graphs during 
     * the test of equivalence or sub graph matching.
     * @return graphComponent
     */
    public GraphComponent getMatchedGC() {
        return matchedGC;
    }
    
    /**
     * set the p_matchedGC as the matched graph component of this. Used in equivalence and subgraph matching 
     * functions.
     * @param p_matchedGC
     */
    public void setMatchedGC(GraphComponent p_matchedGC) {
        matchedGC = p_matchedGC;
    }
    
    /**
     * has the graph component been visited
     * 
     */
    public boolean getVisited() {
        return visited;
    }
    
    /** 
     * changes the visited status of the graph component.
     * @param p_visited
     */
    public void setVisited(boolean p_visited) {
        visited = p_visited;
    }
    
    /**
     * returns the iterator of all the neighbors of the graph component.
     */
    public Iterator getNeighbor() {
        Iterator iter=neighbor.iterator();
        return iter;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\GraphComponent.java
*********************************************************************/


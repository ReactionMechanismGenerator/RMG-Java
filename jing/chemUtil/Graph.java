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

//## package jing::chemUtil

//----------------------------------------------------------------------------
// jing\chemUtil\Graph.java
//----------------------------------------------------------------------------

/**
A graph is a set of GraphComponents (nodes and edges).
*/
//## class Graph
public class Graph {

    /**
    number of nodes in this graph
    */
    protected static int MAXNODENUMBER = 1000;		//## attribute MAXNODENUMBER

    protected boolean acyclic;		//## attribute acyclic

    protected HashMap centralNode;		//## attribute centralNode

    protected LinkedList cycle;		//## attribute cycle

    protected int highestCentralID = 0;		//## attribute highestCentralID

    protected int highestNodeID = 0;		//## attribute highestNodeID

    protected ArrayList arcList;
    protected HashMap nodeList;

    // Constructors

    //## operation Graph()
    public  Graph() {
        {
            arcList=new ArrayList();
        }
        {
            nodeList=new HashMap();
        }
        //#[ operation Graph()
        centralNode = new HashMap();
        //#]
    }

    /**
    Requires:
    Effects: add a new arc with p_arcElement to connect those two nodes on those p_positions.  if there are no nodes on those p_positions, throw NotInGraphException.
    Modifies: this.nodeList, this.arcList
    */
    //## operation addArcBetween(int,Object,int)
    public Arc addArcBetween(int p_position1, Object p_arcElement, int p_position2) throws NotInGraphException {
        //#[ operation addArcBetween(int,Object,int)
        Node node1 = getNodeAt(p_position1);
        Node node2 = getNodeAt(p_position2);

        try {
        	return addArcBetween(node1,p_arcElement,node2);
        }
        catch (NotInGraphException e) {
        	throw new NotInGraphException();
        }



        //#]
    }

    /**
    Requires:
    Effects: if any of the two pass-in two nodes is null or not in the graph, throw NotInGraphException; otherwise, if there already is an arc between, throw PositionOccupiedException; otherwise, add an arc storing the pass-in p_arcElement and connecting the two pass-in nodes.
    Modifies: this.arcList, this.nodeList,p_node1, p_node2
    */
    //## operation addArcBetween(Node,Object,Node)
    public Arc addArcBetween(Node p_node1, Object p_arcElement, Node p_node2) {
        //#[ operation addArcBetween(Node,Object,Node)
        if (p_node1 == null || p_node2 == null) {
        	throw new NotInGraphException("null node");
        }

        if (!contains(p_node1) || !contains(p_node2)) {
        	throw new NotInGraphException("node");
        }

        Arc arc = getArcBetween(p_node1, p_node2);

        if (arc == null) {
        	arc = new Arc(p_arcElement);
         	arcList.add(arc);
        	connect(p_node1, arc, p_node2);
        	return arc;
        }
        else {
        	throw new PositionOccupiedException("arc");
        }



        //#]
    }

    /**
    Requires:
    Effects: saturate all the node atom's undefined valence by adding Hydrogen.
    Modifies: this.nodeList
    */
    //## operation addMissingHydrogen()
    public void addMissingHydrogen() {
        //#[ operation addMissingHydrogen()
        Atom H = Atom.make(ChemElement.make("H"), FreeElectron.make("0"));
        Bond S = Bond.make("S");
        HashMap addedH = new HashMap();

        Iterator iter = getNodeList();
        while (iter.hasNext()) {
                Node node = (Node)iter.next();
                Atom atom = (Atom)node.getElement();
                int val = (int)atom.getValency();

                int bondOrder = 0;
                Iterator neighbor_iter = node.getNeighbor();
                while (neighbor_iter.hasNext()) {
                        Arc arc = (Arc)neighbor_iter.next();
                        Bond bond = (Bond)arc.getElement();
                        bondOrder += bond.getOrder();
                }
                if (bondOrder > val) throw new InvalidConnectivityException();
                else if (bondOrder < val) {
                        addedH.put(node, new Integer(val-bondOrder));
                }
        }
        //Graph g = getGraph();
        iter = addedH.keySet().iterator();
        while (iter.hasNext()) {
                Node node = (Node)iter.next();
                int Hnum = ((Integer)addedH.get(node)).intValue();
                for (int i=0;i<Hnum; i++) {
                        Node newNode = addNode(H);
                        addArcBetween(node, S, newNode);
                }
                node.updateFgElement();
        }

        return;
        //#]
    }


    /**
    Requires:
    Effects: if p_element != null, add a new node with element = p_element; else do nothing and return null.
    Modifies: this.nodeList, this.highestNodeID
    */
    //## operation addNode(Object)
    public Node addNode(Object p_element) {
        //#[ operation addNode(Object)
        if (p_element != null) {
        	highestNodeID++;
        	return addNodeAt(highestNodeID,p_element);
        }
        else return null;
        //#]
    }

    /**
    Requires:
    Effects: if the node at the p_position in this graph has been ocuppied, throw PositionOcuppiedException; otherwise, add a new node storing p_nodeElement at the p_position of this graph
    Modifies: this.nodeList
    */
    //## operation addNodeAt(int,Object)
    public Node addNodeAt(int p_position, Object p_nodeElement) throws PositionOccupiedException {
        //#[ operation addNodeAt(int,Object)
        Node node = new Node(p_position,p_nodeElement);
        if (nodeList.put(node.getID(),node) != null) {
        	throw new PositionOccupiedException("node");
        }

        updateHighestNodeID(p_position);

        return node;



        //#]
    }

    /**
    Requires:
    Effects: add a new node with p_nodeElement and p_centralID at p_position of this graph.  If p_position is occupied or p_centralPosition is occupuied, throw PositionOccupiedException.
    Modifies: this.nodeList, this.centralNodeList, this.highestNodeID, this.highestCentralNodeID
    */
    //## operation addNodeAt(int,Object,int)
    public Node addNodeAt(int p_position, Object p_nodeElement, int p_centralPosition) {
        //#[ operation addNodeAt(int,Object,int)
        Node old = getNodeAt(p_position);
        if (old != null) {
        	throw new PositionOccupiedException("node: " + "node ID = " + old.getID().toString());
        }

        Node node = new Node(p_position,p_nodeElement);
        nodeList.put(node.getID(),node);
        updateHighestNodeID(p_position);

        Integer cenID = new Integer(p_centralPosition);
        node.setCentralID(cenID);
        if (p_centralPosition >= 0) {
        	old = getCentralNodeAt(p_centralPosition);
        	if (old != null) {
        		throw new PositionOccupiedException("central node: " + "node ID = " + old.getID().toString() + "Central ID = " + old.getCentralID().toString());
        	}
        	centralNode.put(cenID,node);
        	updateHighestCentralID(p_centralPosition);
        }

        return node;



        //#]
    }

    /**
    Requires:
    Effects: if all the graph components in this graph are visited, return true; otherwise, return false.
    Modifies:
    */
    //## operation allVisited()
    public boolean allVisited() {
        //#[ operation allVisited()
        Iterator iter1 = getArcList();
        while (iter1.hasNext()) {
        	if (!((Arc)iter1.next()).visited) return false;
        }

        Iterator iter2 = getNodeList();
        while (iter2.hasNext()) {
        	if (!((Node)iter2.next()).visited) return false;
        }

        return true;
        //#]
    }

    /**
    Requires:
    Effects: return true iff there is no center node
    Modifies:
    */
    //## operation centerIsEmpty()
    private boolean centerIsEmpty() {
        //#[ operation centerIsEmpty()
        if (centralNode == null) return true;
        return centralNode.isEmpty();
        //#]
    }

    /**
    Requires:
    Effects: clear node list, arc list, cycle list, and center node list.
    Modifies: this
    */
    //## operation clear()
    public void clear() {
        //#[ operation clear()
        clearNodeList();
        clearArcList();
        clearCycle();
        clearCentralNode();


        //#]
    }

    /**
    Requires:
    Effects: reset the centralID of node in the centralNode list  to -1, and clear centralNode list. reset highestCentralID to 0.
    Modifies: this.centralNode, and the centralID of the nodes in this graph.
    */
    //## operation clearCentralNode()
    public void clearCentralNode() {
        //#[ operation clearCentralNode()
        Iterator iter=getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	node.setCentralID(-1);
        }
        centralNode.clear();
        highestCentralID = 0;
        return;
        //#]
    }

    /**
    Requires:
    Effects: clear cycle list
    Modifies: this.cycle
    */
    //## operation clearCycle()
    public void clearCycle() {
        //#[ operation clearCycle()
        cycle.clear();
        //#]
    }

    /**
    Requires:
    Effects: clear the node list of this graph. reset highestNodeID to 0.
    Modifies: this.nodeList
    */
    //## operation clearNodeList()
    public void clearNodeList() {
        //#[ operation clearNodeList()
        nodeList.clear();
        highestNodeID = 0;


        //#]
    }

    /**
    Requires:
    Effects: put two graphs g1 and g2 into one overall graph, and return the combined graph.  The indexes of the nodes in g2 will be following the ones in g1.
    Modifies:
    */
    //## operation combine(Graph,Graph)
    public static Graph combine(Graph p_g1, Graph p_g2) {
        //#[ operation combine(Graph,Graph)
        Graph result;
        Graph temp;
        result = Graph.copy(p_g1);
        temp = Graph.copy(p_g2);

        int nodeID = result.getHighestNodeID();

        Iterator iter1 = temp.getNodeList();
        while (iter1.hasNext()) {
        	Node n = (Node)iter1.next();

        	int newID = nodeID + n.getID().intValue();
        	int newCentralID = -1;
        	int oldCentralID = n.getCentralID().intValue();

        	if (oldCentralID>0) {
        		Node present = result.getCentralNodeAt(oldCentralID);
        		if (present == null) {
        			newCentralID = oldCentralID;
        		}
        		else {
        			// this is for radical recombination case, bad hard-coded style
        			if (oldCentralID == 1) {
        				newCentralID = 2;
        			}
        			else {
        				throw new ReplaceCentralNodeException();
        			}
        		}
        	}

        	n.setID(newID);
        	n.setCentralID(newCentralID);

        	// add the ID-updated nodes into the combined graph
        	result.nodeList.put(n.getID(),n);
        	result.centralNode.put(n.getCentralID(),n);
        }

        Iterator iter2 = temp.getArcList();
        while (iter2.hasNext()) {
        	Arc a = (Arc)iter2.next();
        	result.arcList.add(a);
        }

        return result;

        //#]
    }

    /**
    Requires:
    Effects: if all pass-in arc and nodes exit in graph, link two nodes by the arc; otherwise, throw NotInGraphException.
    Modifies:
    */
    //## operation connect(Node,Arc,Node)
    public void connect(Node p_node1, Arc p_arc, Node p_node2) {
        //#[ operation connect(Node,Arc,Node)
        if (contains(p_node1) && contains(p_arc) && contains(p_node2)) {
        	p_arc.link(p_node1,p_node2);
        	return;
        }
        else throw new NotInGraphException();

        //#]
    }

    /**
    Requires:
    Effects: return true iff the nodeList contains p_node.
    Modifies:
    */
    //## operation contains(Node)
    public boolean contains(Node p_node) {
        //#[ operation contains(Node)
        return nodeList.containsValue(p_node);
        //#]
    }

    /**
    Requires:
    Effects: return true iff the arcList contains p_arc
    Modifies:
    */
    //## operation contains(Arc)
    public boolean contains(Arc p_arc) {
        //#[ operation contains(Arc)
        return arcList.contains(p_arc);
        //#]
    }

    /**
    Requires:
    Effects: return true iff this graph contains the p_graphComponent.
    Modifies:
    */
    //## operation contains(GraphComponent)
    public boolean contains(GraphComponent p_graphComponent) {
        //#[ operation contains(GraphComponent)
        if (p_graphComponent instanceof Node) return contains((Node)p_graphComponent);
        else if (p_graphComponent instanceof Arc) return contains((Arc)p_graphComponent);
        else throw new InvalidGraphComponentException();
        //#]
    }

    /**
    Requires:
    Effects: return a cloned object of the pass-in graph
    Modifies:
    */
    //## operation copy(Graph)
    public static Graph copy(Graph p_graph) throws InvalidNeighborException {
        //#[ operation copy(Graph)
        Graph result = new Graph();
        Iterator iter = p_graph.getNodeList();
        while (iter.hasNext()) {
        	Node n = (Node)iter.next();
        	Node newNode = result.addNodeAt(n.getID().intValue(),n.getElement(),n.getCentralID().intValue());
        }

        iter = p_graph.getArcList();
        while (iter.hasNext()) {
        	Arc a = (Arc)iter.next();
        	Iterator iter1 = a.getNeighbor();
        	if (!iter1.hasNext()) throw new InvalidNeighborException();
        	Node n1 = (Node)iter1.next();
        	if (!iter1.hasNext()) throw new InvalidNeighborException();
        	Node n2 = (Node)iter1.next();
        	if (iter1.hasNext()) throw new InvalidNeighborException();

        	result.addArcBetween(n1.getID().intValue(),a.getElement(),n2.getID().intValue());
        }
        return result;



        //#]
    }

    /**
    Requires:
    Effects: recursive function to identify possible cycle starting from p_node. Add identified cycle to this.cycle.
    Modifies: this.cycle.  visited status of nodes and arcs.
    */
    //## operation cycleIdentificationDFS(Node,LinkedList)
    public void cycleIdentificationDFS(Node p_node, LinkedList p_list) {
        //#[ operation cycleIdentificationDFS(Node,LinkedList)
        p_node.setVisited(true);
        p_list.add(p_node);

        Iterator iter = p_node.getNeighbor();
        while (iter.hasNext()) {
        	Arc arc = (Arc)iter.next();
        	if (!arc.getVisited()) {
        		arc.setVisited(true);
        		Node otherNode = arc.getOtherNode(p_node);
        		if (!otherNode.getVisited()) {
        			cycleIdentificationDFS(otherNode,p_list);
        		}
        		else {
        			int begin = p_list.indexOf(otherNode);
        			int end = p_list.indexOf(p_node);
        			if (begin == -1 || end == -1) throw new InvalidCycleDetectionException();
            		LinkedList found_cycle = new LinkedList();
            		for (int i = begin; i<=end; i++) {
            			Node node = (Node)p_list.get(i);
            			found_cycle.add(node);
                                    node.setInCycle(true);//svp
                                    if (i != end) {//svp added to include arcs in cycle list
                                      Arc new_arc = getArcBetween(node, (Node) p_list.get(i + 1));
                                      found_cycle.add(new_arc);
                                      new_arc.setInCycle(true);//svp
                                    }

            		}
                            Arc new_arc = getArcBetween( (Node) p_list.get(end),
                                                   (Node) p_list.get(begin));//svp
                           found_cycle.add(new_arc);//svp
                           new_arc.setInCycle(true);//svp

            		cycle.add(found_cycle);
        		}
        	}
        }
        p_list.remove(p_node);

        return;
        //#]
    }

    /**
    Requires: all the graph component's visited reset to false before the deepFirstSearch.
    Effects: mark all the graph components on the deep first search path as visited
    Modifies: visited status of all the graph components in this graph.
    */
    //## operation deepFirstSearch(GraphComponent)
    private static void deepFirstSearch(GraphComponent p_graphComponent) {
        //#[ operation deepFirstSearch(GraphComponent)
        p_graphComponent.setVisited(true);
        Iterator iter = p_graphComponent.getNeighbor();
        while (iter.hasNext()) {
        	GraphComponent gc = (GraphComponent)iter.next();
        	if (!gc.isVisited()) {
        		deepFirstSearch(gc);
        	}
        }
        //#]
    }

    /**
    Requires:
    Effects: remove all the nodes that is Leaf and the associated arc of the removed nodes.
    Modifies: this.nodeList, this.arcList
    */
    //## operation deleteLeaf(Graph)
    protected void deleteLeaf(Graph p_graph) {
        //#[ operation deleteLeaf(Graph)
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	if (node.isLeaf()) removeNode(node);
        }




        //#]
    }

    /**
    Requires:
    Effects: return true iff this graph and p_graph is equivalent at center nodes.  if no center nodes defined, return true iff this graph is equivalent with p_graph.
    Modifies: visited status of nodes and arcs.
    */
    //## operation equals(Object)
    public boolean equals(Object p_graph) {
        //#[ operation equals(Object)
        if (this == p_graph) return true;

        return isEquivalentAtCentralNodes(p_graph);



        //#]
    }

    /**
    Requires: there is no duplicated bond between two atoms.
    Effects: loop over the arc list to find out the arc whose neighbor set includes the position1 and position2.  if found, return this arc; otherwise, return null.
    Modifies:


    */
    //## operation getArcBetween(int,int)
    public Arc getArcBetween(int p_position1, int p_position2) throws NotInGraphException {
        //#[ operation getArcBetween(int,int)
        Node node1 = getNodeAt(p_position1);
        Node node2 = getNodeAt(p_position2);

        return getArcBetween(node1,node2);




        //#]
    }

    /**
    Requires:
    Effects: if there is an arc connecting the two nodes in this graph, return that arc, otherwise, return null.
    Modifies:
    */
    //## operation getArcBetween(Node,Node)
    public Arc getArcBetween(Node p_node1, Node p_node2) {
        //#[ operation getArcBetween(Node,Node)
        if (p_node1 == null || p_node2 == null) {
        	return null;
        }
        if (!contains(p_node1) || !contains(p_node2)) {
        	return null;
        }

        Iterator iter = p_node1.getNeighbor();
        while (iter.hasNext()) {
        	Arc arc = (Arc)(iter.next());
        	if (arc.isConnected(p_node2)) {
        		return arc;
        	}
        }

        return null;

        //#]
    }

    /**
    Requires:
    Effects: return the iterator over the arc collection
    Modifies:
    */
    //## operation getArcList()
    public Iterator getArcList() {
        //#[ operation getArcList()
        Iterator iter = arcList.iterator();
        return iter;
        //#]
    }

    /**
    Requires;
    Effects: return the number of arcs in this Graph
    Modifies:
    */
    //## operation getArcNumber()
    public int getArcNumber() {
        //#[ operation getArcNumber()
        return arcList.size();
        //#]
    }

    /**
    Requires:
    Effects: return the node with p_centralID as centralID in centralNode list
    Modifies:
    */
    //## operation getCentralNodeAt(int)
    public Node getCentralNodeAt(int p_centralID) {
        //#[ operation getCentralNodeAt(int)
        return getCentralNodeAt(new Integer(p_centralID));
        //#]
    }

    /**
    Requires:
    Effects: return the node with p_ID as centralID in centralNode list
    Modifies:
    */
    //## operation getCentralNodeAt(Integer)
    public Node getCentralNodeAt(Integer p_ID) {
        //#[ operation getCentralNodeAt(Integer)
        return (Node)centralNode.get(p_ID);
        //#]
    }

    /**
    Requires:
    Effects: return iterator of center node ID
    Modifies:
    */
    //## operation getCentralNodeKeys()
    public Iterator getCentralNodeKeys() {
        //#[ operation getCentralNodeKeys()
        return centralNode.keySet().iterator();
        //#]
    }

    /**
    Requires:
    Effects: return iterator over the cetnralNode list of this graph.
    Modifies:
    */
    //## operation getCentralNodeList()
    public Iterator getCentralNodeList() {
        //#[ operation getCentralNodeList()
        Iterator iter = centralNode.values().iterator();
        return iter;
        //#]
    }

    /**
    Requires:
    Effects: return number of central nodes in this graph.
    Modifies:
    */
    //## operation getCentralNodeNumber()
    public int getCentralNodeNumber() {
        //#[ operation getCentralNodeNumber()
        return centralNode.size();
        //#]
    }

    /**
    Requires:
    Effects: return the iterator over the cycle collection.
    Modifies:
    */
    //## operation getCycle()
    public Iterator getCycle() {
        //#[ operation getCycle()
        Iterator iter = cycle.iterator();
        return iter;
        //#]
    }

    /**
    Requies:
    Effects: return the highestCenteralID of this graph.
    Modifies:
    */
    //## operation getHighestCentralID()
    public int getHighestCentralID() {
        //#[ operation getHighestCentralID()
        return highestCentralID;
        //#]
    }

    /**
    Requires:
    Effects: return the node at the p_position in this graph, if there is no node at that position, or if at this position, there is a null node, throw NotInGraphException.
    Modifies:
    */
    //## operation getNodeAt(int)
    public Node getNodeAt(int p_position) throws NotInGraphException {
        //#[ operation getNodeAt(int)
        Integer id = new Integer(p_position);
        return getNodeAt(id);



        //#]
    }

    /**
    Requires:
    Effects: return the node at p_ID position in this graph, if there is no node at that position, or if at this position, there is a null node, throw NotInGraphException.
    Modifies:
    */
    //## operation getNodeAt(Integer)
    public Node getNodeAt(Integer p_ID) {
        //#[ operation getNodeAt(Integer)
        return (Node)nodeList.get(p_ID);


        //#]
    }

    /**
    Requires:
    Effects: return the iterator over the node collection
    Modifies:
    */
    //## operation getNodeList()
    public Iterator getNodeList() {
        //#[ operation getNodeList()
        Iterator iter = nodeList.values().iterator();
        return iter;
        //#]
    }

    /**
    Requires:
    Effects: return the number of nodes in this Graph
    Modifies:
    */
    //## operation getNodeNumber()
    public int getNodeNumber() {
        //#[ operation getNodeNumber()
        return nodeList.values().size();
        //#]
    }

    /**
    Requires:
    Effects: compute and return the hashcode of this graph.  hashcode = 10^6*(#node) + #arc
    Modifies:
    */
    //## operation hashCode()
    public int hashCode() {
        //#[ operation hashCode()
        return getNodeNumber() * 100000 + getArcNumber();
        //#]
    }

    /**
    Requires:
    Effects: identify and return all the possible matched pattens between this graph and p_graph.
    Note: the order of center node doesn't matter.
    Modifies: visited status of nodes and arcs in this graph.
    */
    //## operation identifyAllOrderedMatchedSites(Graph)
    public HashSet identifyAllOrderedMatchedSites(Graph p_graph) {
        //#[ operation identifyAllOrderedMatchedSites(Graph)
        HashSet allMatchedSites = new HashSet();

        clearCentralNode();
        resetMatchedGC();
        p_graph.resetMatchedGC();

        Iterator iter2 = p_graph.getCentralNodeList();
        if (iter2.hasNext()) {
        	Node node2 = (Node)iter2.next();
         	Iterator iter1 = getNodeList();
        	while (iter1.hasNext()) {
        		Node node1 = (Node)iter1.next();
          		HashSet matched = node1.identifyAllMatchedSites(node2);
        		if (matched!=null) allMatchedSites.addAll(matched);
        	}
        }

        return allMatchedSites;
        //#]
    }

    /**
    Requires:
    Effects: identify and return all the possible matched pattens between this graph and p_graph.
    Note: the order of center node does matter.
    Modifies: visited status of nodes and arcs in this graph.
    */
    //## operation identifyAllUnorderedMatchedSite(Graph)
    public HashSet identifyAllUnorderedMatchedSite(Graph p_graph) {
        //#[ operation identifyAllUnorderedMatchedSite(Graph)
        HashSet matchedSite = new HashSet();
        HashSet total = new HashSet();
        int size = p_graph.getCentralNodeNumber();
        Stack stack1 = new Stack();
        Stack stack2 = new Stack();

        Iterator iter2 = p_graph.getCentralNodeList();
        boolean found = false;
        if (iter2.hasNext()) {
        	Node node2 = (Node)iter2.next();
         	Iterator iter1 = getNodeList();
        	while (iter1.hasNext()) {
               	clearCentralNode();
        		Node node1 = (Node)iter1.next();
                resetMatchedGC();
                p_graph.resetMatchedGC();

        		found = node1.isSubAndSetCentralID(node2,stack1,stack2);
                if (found) {
                	refreshCentralNode();
                	// check if the combination of the id of central nodes has been found already
                	HashSet IdSet = new HashSet();
                	Iterator iter = getCentralNodeList();
                	while (iter.hasNext()) {
                		Integer id = ((Node)iter.next()).getID();
                		IdSet.add(id);
                	}
                	if (!total.contains(IdSet)) {
                		// if this is a new combination, add to matched site list
                		total.add(IdSet);
                		matchedSite.add((HashMap)(centralNode.clone()));
                	}
                	else {
                		// if this is an old combination, do not add to matched site list
                		IdSet.clear();
                		IdSet = null;
                	}
                }
        	}
        }
        total.clear();
        total = null;

        // there is no central node set in this graph, can't make comparison, return false
        if (matchedSite.size()==0) matchedSite = null;
        return matchedSite;
        //#]
    }

    /**
    Requires:
    Effects:
    (1) identify all the nodes/arcs involved in at least one cycle in this graph.
    (2) add those graph component involved in a ring into cycle collection.
    Modifies: the inCycle attributes of all the graph componenets, and this.cycle
    */
    //## operation identifyCycle()
    public void identifyCycle() throws InconnectedGraphException {
        //#[ operation identifyCycle()
        if (!isConnected()) throw new InconnectedGraphException();

        if ((getNodeNumber()-getArcNumber()) == 1) {
        	acyclic = true;
        	return;
        }

        cycle = new LinkedList();
        setVisited(false);

        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node n = (Node)iter.next();
        	if (!n.getVisited()) {
        		LinkedList list = new LinkedList();
        		cycleIdentificationDFS(n,list);
            }
        }

        if (cycle.size() == 0) acyclic = true;
        else acyclic = false;

        return;
        //#]
    }

    /**
    Requires:
    Effects: identify the type of each node in this graph as Cs, Cd, Ca, Ck, etc.
    Modifies:
    */
    //## operation identifyFgElement()
    public void identifyFgElement() {
        //#[ operation identifyFgElement()
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	node.identifyFgElement();
        }

        return;


        //#]
    }

    /**
    Requires:
    Effects: return true iff this graph is acyclic.
    Modifies:
    */
    //## operation isAcyclic()
    public boolean isAcyclic() {
        //#[ operation isAcyclic()
        if (cycle == null) identifyCycle();
        return acyclic;
        //#]
    }

    /**
    Requires:
    Effects: if all the graph components in this graph are connected to each other, i.e., there is no seperated parts in this graph, return true, otherwise, return false.
    Modifies: visited status of every graph component
    */
    //## operation isConnected()
    public boolean isConnected() {
        //#[ operation isConnected()
        if (isEmpty()) return true;
        setVisited(false);
        Iterator iter = getNodeList();
        Node n = (Node)iter.next();
        deepFirstSearch(n);
        if (!allVisited()) return false;
        else return true;



        //#]
    }

    /**
    Requires:
    Effects: return true iff there is no node and arc in this graph.
    Modifies:
    */
    //## operation isEmpty()
    public boolean isEmpty() {
        //#[ operation isEmpty()
        if (nodeList.isEmpty() && arcList.isEmpty()) return true;
        else return false;
        //#]
    }

    /**
    Requires:
    Effects: if this graph and the pass-in p_graph are equivalent to each other.  No check for centerID.
    Modifies:visited status of all the graph component
    */
    //## operation isEquivalent(Object)
    public boolean isEquivalent(Object p_graph) {
        //#[ operation isEquivalent(Object)
        if (this == p_graph) return true;

        if (!(p_graph instanceof Graph)) return false;
        Graph graph = (Graph)p_graph;

        if (isEmpty() && graph.isEmpty()) return true;

        resetMatchedGC();
        graph.resetMatchedGC();

        Iterator iter1 = getNodeList();
        if (iter1.hasNext()) {
        	GraphComponent gc1 = (GraphComponent)iter1.next();

        	Iterator iter2 = graph.getNodeList();
        	while (iter2.hasNext()) {
        		GraphComponent gc2 = (GraphComponent)iter2.next();
        		Stack stack1 = new Stack();
        		Stack stack2 = new Stack();
        		if (gc1.isEquivalent(gc2,stack1,stack2)) return true;
        	}
        }
        return false;



        //#]
    }

    /**
    Requires:
    Effects: if two graphs both have centers, compare equivalence at centers; if two graphs both have no centers, call isEquivalent to check the no-center equivalence; otherwise, return false.
    Modifies: visited status of nodes in both graphs
    */
    //## operation isEquivalentAtCentralNodes(Object)
    public boolean isEquivalentAtCentralNodes(Object p_graph) {
        //#[ operation isEquivalentAtCentralNodes(Object)
        if (this == p_graph) return true;

        if (!(p_graph instanceof Graph)) return false;
        Graph graph = (Graph)p_graph;

        if (isEmpty() && graph.isEmpty()) return true;

        boolean b1 = centerIsEmpty();
        boolean b2 = graph.centerIsEmpty();

        // if both graphs don't have center, use isEquivalent()
        if (b1 && b2) return isEquivalent(graph);
        // if any graph doesn't have center, but the other has, return false
        else if ((b1 && !b2) || (!b1 && b2)) return false;
        // if both graphs have center, compare center
        else {

        	resetMatchedGC();
        	graph.resetMatchedGC();

        	Iterator iter1 = getCentralNodeList();
        	if (!iter1.hasNext()) return false;
        	GraphComponent gc1 = (GraphComponent)(iter1.next());

        	Iterator iter2 = graph.getCentralNodeList();
        	while (iter2.hasNext()) {
        		GraphComponent gc2 = (GraphComponent)iter2.next();
        		Stack stack1 = new Stack();
        		Stack stack2 = new Stack();
        		if (gc1.isEquivalentCenterMatched(gc2,stack1,stack2)) return true;
        	}

        	return false;
        }


        //#]
    }

    /**
    Requires:
    Effects: return true iff this graph is a subgraph of p_graph, i.e., iff there is a node in this graph is a subNode of one node in p_graph.
    Modifies: Visited status of all the nodes in this and p_graph.
    */
    //## operation isSub(Graph)
    public boolean isSub(Graph p_graph) {
        //#[ operation isSub(Graph)
        for (Iterator iter1 = getNodeList(); iter1.hasNext();) {
        	Node node1 = (Node)iter1.next();
        	for (Iterator iter2 = p_graph.getNodeList(); iter2.hasNext(); ) {
        		Node node2 = (Node)iter2.next();
        		Stack stack1 = new Stack();
        		Stack stack2 = new Stack();
        		resetMatchedGC();
        		p_graph.resetMatchedGC();
        	    if (node1.isSub(node2,stack1,stack2)) return true;
        	}
        }
        return false;



        //#]
    }

    /**
    Requires:
    Effects: return true iff this graph is a subgraph of p_graph at central nodes, i.e., iff a central node in this graph is a subNode of one central node in p_graph.
    Modifies: visited status of all the nodes in this and p_graph.
    */
    //## operation isSubAtCentralNodes(Graph)
    public boolean isSubAtCentralNodes(Graph p_graph) {
        //#[ operation isSubAtCentralNodes(Graph)
        Iterator iter1 = p_graph.getCentralNodeList();

        Stack s1 = new Stack();
        Stack s2 = new Stack();

        if (iter1.hasNext()) {
        	Node node1 = (Node)iter1.next();
            // check if there is a matched central node in p_graph
        	Integer cid = node1.getCentralID();
        	Node node2 = getCentralNodeAt(cid);

        	if (node2==null) return false;

            resetMatchedGC();
            p_graph.resetMatchedGC();

        	return node2.isSubCentralMatched(node1,s1,s2);
        }
        // there is no central node set in this graph, can't make comparison, return false
        else {
        	return false;
        }



        //#]
    }

    /**
    Requires:
    Effects: if this graph is connected, do nothing; if the graph is not connected, partioning the graph into a list of graphs.
    Modifies: this
    */
    //## operation partition()
    public LinkedList partition() throws NotPartitionedException {
        //#[ operation partition()
        LinkedList result = new LinkedList();

        if (isEmpty()) {
        	return result;
        }

        setVisited(false);
        Iterator iter = getNodeList();
        deepFirstSearch((Node)iter.next());

        while (true) {
        	Graph newGraph = new Graph();
        	int nodeID = 0;
        	int centralID, highestCentralID=-1;
        	Iterator iter1 = nodeList.keySet().iterator();
        	while (iter1.hasNext()) {
        		Integer ID = (Integer)iter1.next();
        		Node n = (Node)nodeList.get(ID);
        		if (n.isVisited()) {
        			nodeID++;
        			centralID = n.getCentralID().intValue();
        			if (centralID >= 0) {
        				if (highestCentralID<centralID) highestCentralID = centralID;
        	 		}
        			iter1.remove();
        			if (centralID >= 0) centralNode.remove(n.getCentralID());
        			n.setID(nodeID);
        			n.setCentralID(centralID);
        			newGraph.nodeList.put(n.getID(),n);
        			if (centralID >= 0) newGraph.centralNode.put(n.getCentralID(),n);
        		}
        	}
        	newGraph.highestNodeID = nodeID;
        	newGraph.highestCentralID = highestCentralID;

        	Iterator iter2 = getArcList();
        	while (iter2.hasNext()) {
        		Arc a = (Arc)iter2.next();
        		if (a.isVisited()) {
        			newGraph.arcList.add(a);
        			iter2.remove();
         		}
        	}
        	result.add(newGraph);

        	if (nodeList.isEmpty()) break;
        	iter = getNodeList();
        	deepFirstSearch((Node)iter.next());
        }

        if (result.size() == 0) {
        	throw new NotPartitionedException();
        }
        else return result;
        //#]
    }

    /**
    Requires:
    Effects: clear central node list, and reform it according to centralID information of node in the graph
    Modifies: this.centralNodeList

    */
    //## operation refreshCentralNode()
    public void refreshCentralNode() {
        //#[ operation refreshCentralNode()
        centralNode.clear();
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Integer centralID = node.getCentralID();
        	if (centralID.intValue() > 0) centralNode.put(centralID, node);
        }




        //#]
    }

    /**
    Requires:
    Effects: refresh centralNode list according to the nodeID setting of each individual node in nodeList
    Modifies: this.centralNode
    */
    //## operation refreshHighestNodeID()
    public void refreshHighestNodeID() {
        //#[ operation refreshHighestNodeID()
        int hID = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	int ID = node.getID().intValue();
        	if (ID > hID) hID = ID;
        }
        highestNodeID = hID;


        //#]
    }

    /**
    Requires:
    Effects: if the pass-in arc is in this graph, remove it from this graph, also get rid of all the connectivity between this arc and all it neighbor nodes.  if it is not in the graph, throw NotInGraphException.
    Modifies: this
    */
    //## operation removeArc(Arc)
    public void removeArc(Arc p_arc) {
        //#[ operation removeArc(Arc)
        if (!contains(p_arc)) throw new NotInGraphException();

        // remove p_arc from all the node's neighbor collection if the node connect to p_arc
        Iterator iter = p_arc.getNeighbor();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
            node.neighbor.remove(p_arc);
        }

        arcList.remove(p_arc);
        p_arc=null;
        //#]
    }

    /**
    Requires:
    Effects: remove arc between the two nodes at p_position1 and p_position2, iff that arc exists.
    Modifies: this.arcList
    */
    //## operation removeArcAt(int,int)
    public void removeArcAt(int p_position1, int p_position2) {
        //#[ operation removeArcAt(int,int)
        try {
        	Arc arc = getArcBetween(p_position1,p_position2);
        	removeArc(arc);
        }
        catch (NotInGraphException e) {
        	return;
        }



        //#]
    }

    /**
    Requires:
    Effects: if p_node is in this graph, remove it, and remove all the arcs linked to it.  if p_node is not in this graph, throw NotInGraphException.
    Modifies: this.nodeList
    */
    //## operation removeNode(Node)
    public void removeNode(Node p_node) {
        //#[ operation removeNode(Node)
        if (!contains(p_node)) throw new NotInGraphException();

        Stack stack = new Stack();

        // find all arc linked to this removed node
        Iterator iter = p_node.getNeighbor();
        while (iter.hasNext()) {
        	Arc arc = (Arc)iter.next();
        	Node otherNode = p_node.getOtherNode(arc);
        	if (!contains(otherNode)) throw new InvalidConnectivityException();
        	otherNode.neighbor.remove(arc);
        	arcList.remove(arc);
        }
        // remove this node
        nodeList.remove(p_node.getID());

        if (p_node.getID().intValue() >= highestNodeID) {
        	refreshHighestNodeID();
        }

        p_node = null;






        //#]
    }

    /**
    Requires:
    Effects: remove the node at the p_position, iff it exists
    Modifies: this.nodeList
    */
    //## operation removeNodeAt(int)
    public void removeNodeAt(int p_position) {
        //#[ operation removeNodeAt(int)
        try {
        	Node node = getNodeAt(p_position);
        	removeNode(node);
        }
        catch (NotInGraphException e) {
        	return;
        }



        //#]
    }

    /**
    Requires:
    Effects: return true.  there is no restrict for graph invalidity
    Modifies:
    */
    //## operation repOk()
    public boolean repOk() {
        //#[ operation repOk()
        return true;
        //#]
    }

    /**
    Requires:
    Effects: reset all the matched GC for eac node and arc in the graph to null, to prepare for a new matching/equivalent testing
    Modifies: matchedGC of every node and arc in this graph
    */
    //## operation resetMatchedGC()
    public void resetMatchedGC() {
        //#[ operation resetMatchedGC()
        Iterator iter1 = getArcList();
        while (iter1.hasNext()) {
        	((GraphComponent)iter1.next()).setMatchedGC(null);
        }
        Iterator iter2 = getNodeList();
        while (iter2.hasNext()) {
        	((GraphComponent)iter2.next()).setMatchedGC(null);
        }



        //#]
    }

    /**
    Requires:
    Effects: if the pass-in node in this graph, set it as the ith central node of this graph, where i = p_position; if the node is not in graph, throw NotInGraphException.
    Modifies: centralID of the node in this graph, this.centralNode list
    */
    //## operation setCentralNode(int,Node)
    public void setCentralNode(int p_position, Node p_node) throws NotInGraphException {
        //#[ operation setCentralNode(int,Node)
        if (!contains(p_node)) {
        	throw new NotInGraphException();
        }

        p_node.setCentralID(p_position);
        Node old = (Node)centralNode.put(new Integer(p_position),p_node);
        updateHighestCentralID(p_position);
        if (old != null) old.setCentralID(-1);
        return;
        //#]
    }

    /**
    Requires:
    Effects: set all the current central nodes as normal ones, and reset the centralID and centralNode list according to the pass-in p_centralNode HashMap.
    Modifies: this.centralNode, centralID of all the old  and new central nodes
    */
    //## operation setCentralNodes(HashMap)
    public void setCentralNodes(HashMap p_centralNode) {
        //#[ operation setCentralNodes(HashMap)
        clearCentralNode();

        Iterator iter = p_centralNode.keySet().iterator();
        while (iter.hasNext()) {
        	Integer centralID = (Integer)iter.next();
        	Node node = (Node)(p_centralNode.get(centralID));
        	if (!contains(node)) throw new NotInGraphException();
        	node.setCentralID(centralID);
        	Node old = (Node)centralNode.put(centralID,node);
        	updateHighestCentralID(centralID.intValue());
        	if (old != null) old.setCentralID(-1);
        }
        //#]
    }

    /**
    Requires:
    Effects: set the visited status of all the graph components in this graph as the pass-in p_visited.
    Modifies: this.nodeList, this.arcList
    */
    //## operation setVisited(boolean)
    private void setVisited(boolean p_visited) {
        //#[ operation setVisited(boolean)
        Iterator iter1 = getArcList();
        while (iter1.hasNext()) {
        	((GraphComponent)iter1.next()).setVisited(p_visited);
        }
        Iterator iter2 = getNodeList();
        while (iter2.hasNext()) {
        	((GraphComponent)iter2.next()).setVisited(p_visited);
        }


        //#]
    }

    /**
    Requies:
    Effects: output a string of adjacency list of this graph
    Modifies:
    */
    //## operation toString()
    public String toString() {
        //#[ operation toString()
        String s = "";
        for (int i = 1; i<=highestNodeID; i++) {
        	Node n = getNodeAt(i);
        	if (n != null) s = s + n.toString() + '\n';
        }

        return s;
        //#]
    }

    /**
    Requies:
    Effects: output a string of adjacency list of this graph without outputting the centralID information of node
    Modifies:
    */
    //## operation toStringWithoutCentralID()
    public String toStringWithoutCentralID() {
        //#[ operation toStringWithoutCentralID()
        String s = "";
        for (int i = 1; i<=highestNodeID; i++) {
        	Node n = getNodeAt(i);
        	if (n != null) s = s + n.toStringWithoutCentralID() + '\n';
        }

        return s;
        //#]
    }

    /**
    Requies:
    Effects: output a string of adjacency list of this graph without outputting the centralID information of node and without outputting H node
    Modifies:
    */
    //## operation toStringWithoutCentralIDAndH()
    public String toStringWithoutCentralIDAndH() {
        //#[ operation toStringWithoutCentralIDAndH()
        String s = "";
        for (int i = 1; i<=highestNodeID; i++) {
        	Node n = getNodeAt(i);
        	if (n != null) {
        		Atom a = (Atom)n.getElement();
        		if (!a.isHydrogen()) {
        			s = s + n.toStringWithoutCentralIDAndH() + '\n';
        		}
        	}
        }

        return s;
        //#]
    }

    /**
    Requies:
    Effects: set highestCentralID to p_centralID iff present highestCentralID is less than p_centralID
    Modifies: this.highestCentralID
    */
    //## operation updateHighestCentralID(int)
    public void updateHighestCentralID(int p_centralID) {
        //#[ operation updateHighestCentralID(int)
        if (p_centralID > highestCentralID) highestCentralID = p_centralID;
        return;
        //#]
    }

    /**
    Requires:
    Effects: set highestNodeID = p_position iff p_position > highestNodeID
    Modifies: this.highestNodeID

    */
    //## operation updateHighestNodeID(int)
    public void updateHighestNodeID(int p_position) {
        //#[ operation updateHighestNodeID(int)
        if (p_position > highestNodeID) highestNodeID = p_position;
        return;
        //#]
    }

    protected static int getMAXNODENUMBER() {
        return MAXNODENUMBER;
    }

    public HashMap getCentralNode() {
        return centralNode;
    }

    public int getHighestNodeID() {
        return highestNodeID;
    }

    public void clearArcList() {
        arcList.clear();
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\Graph.java
*********************************************************************/


// //////////////////////////////////////////////////////////////////////////////
//
// RMG - Reaction Mechanism Generator
//
// Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
// RMG Team (rmg_dev@mit.edu)
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// //////////////////////////////////////////////////////////////////////////////
package jing.chemUtil;

import jing.chem.*;
import jing.mathTool.MathTool;
import java.util.*;

import jing.rxnSys.Logger;

// ## package jing::chemUtil
// ----------------------------------------------------------------------------
// jing\chemUtil\Graph.java
// ----------------------------------------------------------------------------
/**
 * A graph is a set of GraphComponents (nodes and edges).
 */
// ## class Graph
public class Graph {
    /**
     * number of nodes in this graph
     */
    private static int MAXNODENUMBER = 1000; // ## attribute MAXNODENUMBER
    public boolean acyclic; // ## attribute acyclic
    private LinkedHashMap centralNode; // ## attribute centralNode
    public LinkedList SSSRings; // ## attribute cycle
    private int highestCentralID = 0; // ## attribute highestCentralID
    private int highestNodeID = 0; // ## attribute highestNodeID
    private int lowestCentralID = 10000;
    private ArrayList arcList;
    private LinkedHashMap nodeList;
    private boolean[] isAromatic;

    // Constructors
    // ## operation Graph()
    public Graph() {
        {
            arcList = new ArrayList();
        }
        {
            nodeList = new LinkedHashMap();
        }
        centralNode = new LinkedHashMap();
    }

    /**
     * Add a new arc with p_arcElement to connect those two nodes on those p_positions. if there are no nodes on those
     * p_positions, throw NotInGraphException.<br>
     * <b>Modifies</b><br>
     * this.nodeList, this.arcList
     */
    // ## operation addArcBetween(int,Object,int)
    public Arc addArcBetween(int p_position1, Object p_arcElement,
            int p_position2) throws NotInGraphException {
        // #[ operation addArcBetween(int,Object,int)
        Node node1 = getNodeAt(p_position1);
        Node node2 = getNodeAt(p_position2);
        try {
            return addArcBetween(node1, p_arcElement, node2);
        } catch (NotInGraphException e) {
            throw new NotInGraphException();
        }
        // #]
    }

    /**
     * If any of the two pass-in two nodes is null or not in the graph, throw NotInGraphException; otherwise, if there
     * already is an arc between, throw PositionOccupiedException; otherwise, add an arc storing the pass-in
     * p_arcElement and connecting the two pass-in nodes.<br>
     * <b>Modifies</b><br>
     * this.arcList, this.nodeList,p_node1, p_node2
     */
    // ## operation addArcBetween(Node,Object,Node)
    public Arc addArcBetween(Node p_node1, Object p_arcElement, Node p_node2) {
        // #[ operation addArcBetween(Node,Object,Node)
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
        } else {
            Logger.error("Error in Graph.addArcBetween() method: "
                    + "RMG is about to add a bond between atoms "
                    + p_node1.getID() + " and " + p_node2.getID()
                    + " in the graph presented below.\n"
                    + "Error most likely occurs because a bond already exists "
                    + "between these two atoms.\nConsider adding graph to"
                    + " forbiddenGroups.txt file.\n" + this.toString());
            throw new PositionOccupiedException("arc");
        }
        // #]
    }

    /**
     * Saturate all the node atom's undefined valence by adding Hydrogen.<br>
     * <b>Modifies</b><br>
     * this.nodeList
     */
    // ## operation addMissingHydrogen()
    public void addMissingHydrogen() {
        // #[ operation addMissingHydrogen()
        Atom H = Atom.make(ChemElement.make("H"), FreeElectron.make("0"));
        Bond S = Bond.make("S");
        LinkedHashMap addedH = new LinkedHashMap();
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
            Node node = (Node) iter.next();
            Atom atom = (Atom) node.getElement();
            int val = (int) atom.getValency();
            double bondOrder = 0;
            Iterator neighbor_iter = node.getNeighbor();
            while (neighbor_iter.hasNext()) {
                Arc arc = (Arc) neighbor_iter.next();
                Bond bond = (Bond) arc.getElement();
                bondOrder += bond.getOrder();
            }
            if (bondOrder > val)
                throw new InvalidConnectivityException();
            else if (bondOrder < val) {
                addedH.put(node, new Integer((int) (val - bondOrder)));
            }
        }
        // Graph g = getGraph();
        iter = addedH.keySet().iterator();
        while (iter.hasNext()) {
            Node node = (Node) iter.next();
            int Hnum = ((Integer) addedH.get(node)).intValue();
            for (int i = 0; i < Hnum; i++) {
                Node newNode = addNode(H);
                addArcBetween(node, S, newNode);
            }
            node.updateFgElement();
        }
        return;
        // #]
    }

    /**
     * If p_element != null, add a new node with element = p_element; else do nothing and return null.<br>
     * <b>Modifies</b><br>
     * this.nodeList, this.highestNodeID
     */
    // ## operation addNode(Object)
    public Node addNode(Object p_element) {
        // #[ operation addNode(Object)
        if (p_element != null) {
            highestNodeID++;
            return addNodeAt(highestNodeID, p_element);
        } else
            return null;
        // #]
    }

    /**
     * If the node at the p_position in this graph has been ocuppied, throw PositionOcuppiedException; otherwise, add a
     * new node storing p_nodeElement at the p_position of this graph<br>
     * <b>Modifies</b><br>
     * this.nodeList
     */
    // ## operation addNodeAt(int,Object)
    public Node addNodeAt(int p_position, Object p_nodeElement)
            throws PositionOccupiedException {
        // #[ operation addNodeAt(int,Object)
        Node node = new Node(p_position, p_nodeElement);
        if (nodeList.put(node.getID(), node) != null) {
            throw new PositionOccupiedException("node");
        }
        updateHighestNodeID(p_position);
        return node;
    }

    /**
     * Add a new node with p_nodeElement and p_centralID at p_position of this graph. If p_position is occupied or
     * p_centralPosition is occupuied, throw PositionOccupiedException.<br>
     * <b>Modifies</b><br>
     * this.nodeList, this.centralNodeList, this.highestNodeID, this.highestCentralNodeID
     */
    // ## operation addNodeAt(int,Object,int)
    public Node addNodeAt(int p_position, Object p_nodeElement,
            int p_centralPosition) {
        // #[ operation addNodeAt(int,Object,int)
        Node old = getNodeAt(p_position);
        if (old != null) {
            throw new PositionOccupiedException("node: " + "node ID = "
                    + old.getID().toString());
        }
        Node node = new Node(p_position, p_nodeElement);
        nodeList.put(node.getID(), node);
        updateHighestNodeID(p_position);
        Integer cenID = new Integer(p_centralPosition);
        node.setCentralID(cenID);
        if (p_centralPosition >= 0) {
            old = getCentralNodeAt(p_centralPosition);
            if (old != null) {
                throw new PositionOccupiedException("central node: "
                        + "node ID = " + old.getID().toString()
                        + "Central ID = " + old.getCentralID().toString());
            }
            centralNode.put(cenID, node);
            updateHighestCentralID(p_centralPosition);
            updateLowestCentralID(p_centralPosition);
        }
        return node;
        // #]
    }

    /**
     * If all the graph components in this graph are visited, return true; otherwise, return false.
     */
    // ## operation allVisited()
    public boolean allVisited() {
        // #[ operation allVisited()
        Iterator iter1 = getArcList();
        while (iter1.hasNext()) {
            if (!((Arc) iter1.next()).visited)
                return false;
        }
        Iterator iter2 = getNodeList();
        while (iter2.hasNext()) {
            if (!((Node) iter2.next()).visited)
                return false;
        }
        return true;
        // #]
    }

    /**
     * Return true iff there is no center node
     */
    // ## operation centerIsEmpty()
    private boolean centerIsEmpty() {
        // #[ operation centerIsEmpty()
        if (centralNode == null)
            return true;
        return centralNode.isEmpty();
        // #]
    }

    /**
     * Clear node list, arc list, cycle list, and center node list.<br>
     * <b>Modifies</b><br>
     * this
     */
    // ## operation clear()
    public void clear() {
        // #[ operation clear()
        clearNodeList();
        clearArcList();
        clearCycle();
        clearCentralNode();
        // #]
    }

    /**
     * Reset the centralID of node in the centralNode list to -1, and clear centralNode list. reset highestCentralID to
     * 0.<br>
     * <b>Modifies</b><br>
     * this.centralNode, and the centralID of the nodes in this graph.
     */
    // ## operation clearCentralNode()
    public void clearCentralNode() {
        // #[ operation clearCentralNode()
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
            Node node = (Node) iter.next();
            node.setCentralID(-1);
        }
        centralNode.clear();
        highestCentralID = 0;
        lowestCentralID = 10000;
        return;
        // #]
    }

    /**
     * Clear cycle list<br>
     * <b>Modifies</b><br>
     * this.cycle
     */
    // ## operation clearCycle()
    public void clearCycle() {
        // #[ operation clearCycle()
        SSSRings.clear();
        // #]
    }

    /**
     * Clear the node list of this graph. reset highestNodeID to 0.<br>
     * <b>Modifies</b><br>
     * this.nodeList
     */
    // ## operation clearNodeList()
    public void clearNodeList() {
        // #[ operation clearNodeList()
        nodeList.clear();
        highestNodeID = 0;
        // #]
    }

    /**
     * Put two graphs g1 and g2 into one overall graph, and return the combined graph. The indexes of the nodes in g2
     * will be following the ones in g1.
     */
    // ## operation combine(Graph,Graph)
    public static Graph combine(Graph p_g1, Graph p_g2) {
        // #[ operation combine(Graph,Graph)
        Graph result;
        Graph temp;
        result = Graph.copy(p_g1);
        temp = Graph.copy(p_g2);
        int nodeID = result.getHighestNodeID();
        Iterator iter1 = temp.getNodeList();
        while (iter1.hasNext()) {
            Node n = (Node) iter1.next();
            int newID = nodeID + n.getID().intValue();
            int newCentralID = -1;
            int oldCentralID = n.getCentralID().intValue();
            if (oldCentralID > 0) {
                Node present = result.getCentralNodeAt(oldCentralID);
                if (present == null) {
                    newCentralID = oldCentralID;
                } else {
                    // this is for radical recombination case, bad hard-coded style
                    if (oldCentralID == 1) {
                        newCentralID = 2;
                    } else {
                        throw new ReplaceCentralNodeException();
                    }
                }
            }
            n.setID(newID);
            n.setCentralID(newCentralID);
            // add the ID-updated nodes into the combined graph
            result.nodeList.put(n.getID(), n);
            result.centralNode.put(n.getCentralID(), n);
        }
        Iterator iter2 = temp.getArcList();
        while (iter2.hasNext()) {
            Arc a = (Arc) iter2.next();
            result.arcList.add(a);
        }
        return result;
        // #]
    }

    /**
     * If all pass-in arc and nodes exit in graph, link two nodes by the arc; otherwise, throw NotInGraphException.
     */
    // ## operation connect(Node,Arc,Node)
    public void connect(Node p_node1, Arc p_arc, Node p_node2) {
        // #[ operation connect(Node,Arc,Node)
        if (contains(p_node1) && contains(p_arc) && contains(p_node2)) {
            p_arc.link(p_node1, p_node2);
            return;
        } else
            throw new NotInGraphException();
        // #]
    }

    /**
     * Return true iff the nodeList contains p_node.
     */
    // ## operation contains(Node)
    public boolean contains(Node p_node) {
        // #[ operation contains(Node)
        return nodeList.containsValue(p_node);
        // #]
    }

    /**
     * Return true iff the arcList contains p_arc
     */
    // ## operation contains(Arc)
    public boolean contains(Arc p_arc) {
        // #[ operation contains(Arc)
        return arcList.contains(p_arc);
        // #]
    }

    /**
     * Return true iff this graph contains the p_graphComponent.
     */
    // ## operation contains(GraphComponent)
    public boolean contains(GraphComponent p_graphComponent) {
        // #[ operation contains(GraphComponent)
        if (p_graphComponent instanceof Node)
            return contains((Node) p_graphComponent);
        else if (p_graphComponent instanceof Arc)
            return contains((Arc) p_graphComponent);
        else
            throw new InvalidGraphComponentException();
        // #]
    }

    /**
     * Return a cloned object of the pass-in graph
     */
    // ## operation copy(Graph)
    public static Graph copy(Graph p_graph) throws InvalidNeighborException {
        // #[ operation copy(Graph)
        Graph result = new Graph();
        Iterator iter = p_graph.getNodeList();
        while (iter.hasNext()) {
            Node n = (Node) iter.next();
            Node newNode = result.addNodeAt(n.getID().intValue(),
                    n.getElement(), n.getCentralID().intValue());
        }
        iter = p_graph.getArcList();
        while (iter.hasNext()) {
            Arc a = (Arc) iter.next();
            Iterator iter1 = a.getNeighbor();
            if (!iter1.hasNext())
                throw new InvalidNeighborException();
            Node n1 = (Node) iter1.next();
            if (!iter1.hasNext())
                throw new InvalidNeighborException();
            Node n2 = (Node) iter1.next();
            if (iter1.hasNext())
                throw new InvalidNeighborException();
            result.addArcBetween(n1.getID().intValue(), a.getElement(), n2
                    .getID().intValue());
        }
        return result;
        // #]
    }

    public static Graph copywithSSSR(Graph p_graph) throws InvalidNeighborException {
        // #[ operation copy(Graph)
        Graph result = new Graph();
        Iterator iter = p_graph.getNodeList();
        while (iter.hasNext()) {
            Node n = (Node) iter.next();
            Node newNode = result.addNodeAt(n.getID().intValue(),
                    n.getElement(), n.getCentralID().intValue());
        }
        iter = p_graph.getArcList();
        while (iter.hasNext()) {
            Arc a = (Arc) iter.next();
            Iterator iter1 = a.getNeighbor();
            if (!iter1.hasNext())
                throw new InvalidNeighborException();
            Node n1 = (Node) iter1.next();
            if (!iter1.hasNext())
                throw new InvalidNeighborException();
            Node n2 = (Node) iter1.next();
            if (iter1.hasNext())
                throw new InvalidNeighborException();
            result.addArcBetween(n1.getID().intValue(), a.getElement(), n2
                    .getID().intValue());
        }
        result.formSSSR();
        return result;
        // #]
    }

    /*
     * /** Recursive function to identify possible cycle starting from p_node. Add identified cycle to this.cycle.<br>
     * <b>Modifies</b><br> this.cycle. visited status of nodes and arcs.
     */
    // ## operation cycleIdentificationDFS(Node,LinkedList)
    /*
     * public void cycleIdentificationDFS(Node p_node, LinkedList p_list) { //#[ operation
     * cycleIdentificationDFS(Node,LinkedList) p_node.setVisited(true); p_list.add(p_node); Iterator iter =
     * p_node.getNeighbor(); while (iter.hasNext()) { Arc arc = (Arc)iter.next(); if (!arc.getVisited()) {
     * arc.setVisited(true); Node otherNode = arc.getOtherNode(p_node); if (!otherNode.getVisited()) {
     * cycleIdentificationDFS(otherNode,p_list); } else { int begin = p_list.indexOf(otherNode); int end =
     * p_list.indexOf(p_node); if (begin == -1 || end == -1) throw new InvalidCycleDetectionException(); LinkedList
     * found_cycle = new LinkedList(); for (int i = begin; i<=end; i++) { Node node = (Node)p_list.get(i);
     * found_cycle.add(node); node.setInCycle(true);//svp if (i != end) {//svp added to include arcs in cycle list Arc
     * new_arc = getArcBetween(node, (Node) p_list.get(i + 1)); found_cycle.add(new_arc); new_arc.setInCycle(true);//svp
     * } } Arc new_arc = getArcBetween( (Node) p_list.get(end), (Node) p_list.get(begin));//svp
     * found_cycle.add(new_arc);//svp new_arc.setInCycle(true);//svp cycle.add(found_cycle); } } }
     * p_list.remove(p_node); return; //#] }
     */
    /**
     * <b>Requires</b><br>
     * all the graph component's visited reset to false before the deepFirstSearch.<br>
     * <b>Effects</b><br>
     * mark all the graph components on the deep first search path as visited<br>
     * <b>Modifies</b><br>
     * visited status of all the graph components in this graph.
     */
    // ## operation deepFirstSearch(GraphComponent)
    private static void deepFirstSearch(GraphComponent p_graphComponent) {
        // #[ operation deepFirstSearch(GraphComponent)
        p_graphComponent.setVisited(true);
        Iterator iter = p_graphComponent.getNeighbor();
        while (iter.hasNext()) {
            GraphComponent gc = (GraphComponent) iter.next();
            if (!gc.isVisited()) {
                deepFirstSearch(gc);
            }
        }
        // #]
    }

    /**
     * Remove all the nodes that is Leaf and the associated arc of the removed nodes.<br>
     * <b>Modifies</b><br>
     * this.nodeList, this.arcList
     */
    // ## operation deleteLeaf(Graph)
    protected void deleteLeaf(Graph p_graph) {
        // #[ operation deleteLeaf(Graph)
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
            Node node = (Node) iter.next();
            if (node.isLeaf())
                removeNode(node);
        }
        // #]
    }

    public boolean strainedRing(Node node1,Node node2) { 
    boolean strained = false;
    for (int i = 0; i < SSSRings.size(); i++) {
        LinkedList cycle1 = (LinkedList) SSSRings.get(i);
        if (cycle1.contains(node1)) {
	    for (int j = 0; j < SSSRings.size(); j++) {
	    int atomCount = 0;
            LinkedList cycle2 = (LinkedList) SSSRings.get(j);
	    if (i!=j && cycle2.contains(node2)) {
		for (int k = 0; k < cycle2.size(); k++) {
		    GraphComponent gc = (GraphComponent) cycle2.get(k);
		    if (gc instanceof Node) {	
		       Node n = (Node) gc;
                       if(cycle1.contains(n))
			   atomCount += 1;
		       }	
		    }
		}
	    if(atomCount >= 2)
		strained = true;
	    }
	}
    }
    return strained;
    }

    public boolean inBiRing(Node node1) {
    int count = 0;
    for (int i = 0; i < SSSRings.size(); i++) {
        LinkedList cycle = (LinkedList) SSSRings.get(i);
        if(cycle.contains(node1)) {
	    count++;
            }
        }
    if(count >= 2) 
	return true;
    else
	return false;
    }

    public boolean sameRing(Node node1, Node node2) { 

    for (int i = 0; i < SSSRings.size(); i++) {
        LinkedList cycle = (LinkedList) SSSRings.get(i);
	if(cycle.contains(node1) && cycle.contains(node2)) {
	    return true;
	    }
    }

    return false;

    }

    public boolean hasExocyclicPi(LinkedList cycle) {
        boolean hasExoPi = false;
        for (int i = 0; i < cycle.size(); i = i + 2) {
            Node n = (Node) cycle.get(i);
            Iterator nodeNeighbor = n.getNeighbor();
            while (nodeNeighbor.hasNext()) {
                Arc arc = (Arc) nodeNeighbor.next();
                Bond b = (Bond) arc.getElement();
                if (!cycle.contains(arc) && b.isDouble())
                    return true;
            }
        }
        return hasExoPi;
    }

    public void getAromatic(int[] alreadyClassified) {
        isAromatic = new boolean[SSSRings.size()];
 
        LinkedList ringWithExoCyclicPi = new LinkedList();
        for (int i = 0; i < SSSRings.size(); i++) {
            LinkedList cycle = (LinkedList) SSSRings.get(i);
            boolean ringDoubleBonds = false;
            boolean quarternaryAtom = false;
            boolean saturatedCarbon = false;
            int numPiBonds = 0;
            // preliminary screening
            for (int j = 0; j < cycle.size(); j++) {
                GraphComponent gc = (GraphComponent) cycle.get(j);
                if (gc instanceof Node) {
                    Node n = (Node) gc;
                    Atom a = (Atom) n.getElement();
                    // UPDATE: for our purposes, we want the presence of one double bond (either in ring or outside ring
                    // (in the case of some naphthalene resonance isomers)) to be a necessary condition for aromaticity (radicals are
                    // tricky; phenyl should be aromatic, but not C3H2)
                    // we will do this check in the calling function from ChemGraph, as it will be faster as it won't
                    // get called as often
				    // UPDATE November 19th 2013 by AG Vandeputte
				    // Tried to resolve the naphtalene issue
				    // The code will loop trough all rings and only using exo-aromatic double bonds for the Huckel test
				    // hence each time a new aromatic ring is found we have to loop over all rings again => growing aromatic complex
				    // Before, double bonds were always added leading to fake hits which were then overwritten in Graph.java where a 
				    // a final test (check if double bond is inside) was performed. That test has been removed.
				    // So far code seems to work for benzene, napthalene and does not give a fake hit for InChI=1/C9H8/c1-2-5-9-7-3-6-8(9)4-1/h1-2,4-7H,3H2 

                    
                    // has more than 2 saturated carbon atoms, not aromatic
                    if (a.isCarbon() && n.getNeighborNumber() == 4) {

                     	    if (saturatedCarbon) {
                            isAromatic[i] = false;
                            alreadyClassified[i] = 1;
                            break;

                          } else
                          saturatedCarbon = true;
                    }
                    
                    // Added by AG Vandeputte for problems with phenyl radical, as all resonance isomers obey the Huckel theory they are all considered aromatic
                    // has a quarternary atom, not aromatic
                    if (n.getNeighborNumber() == 4) {
                        Iterator neighborNodes = n.getNeighboringNodes()
                                .iterator();
                        while (neighborNodes.hasNext()) {
                            Node neighborNode = (Node) neighborNodes.next();
                            if (((Atom) neighborNode.getElement()).isHydrogen()) {
                                quarternaryAtom = true;
                                alreadyClassified[i] = 1;
                                break;
                            }
                        }
                        if (quarternaryAtom) {
                            isAromatic[i] = false;
                            alreadyClassified[i] = 1;
                            break;
                        }
                    }
                } else {
                    Arc a = (Arc) gc;
                    if (((Bond) a.getElement()).isDouble())
                        ringDoubleBonds = true;
                }
            }
            if (!ringDoubleBonds) {
                isAromatic[i] = false;
                alreadyClassified[i] = 1;
            }
        } 


        int check=0;
        while(check < SSSRings.size()) {
            check=classifyAsAromatic(check, alreadyClassified);

        }
        for (int i = 0; i < SSSRings.size(); i++) {
            if(alreadyClassified[i] == 0) {
            isAromatic[i] = false;
            }

        } 

    }

    public int classifyAsAromatic(int j, int[] alreadyClassified) {
        if (alreadyClassified[j] == 1)
            return ++j;

        LinkedList cycle = (LinkedList) SSSRings.get(j);
        int aromaticExoPi = 0;
        int nonAromaticExoPi = 0;

        for (int i = 0; i < cycle.size(); i = i + 1) {

          GraphComponent gx = (GraphComponent) cycle.get(i);
          if (gx instanceof Node) {

            Node n = (Node) gx;

            Iterator neighbor = n.getNeighbor();
            while (neighbor.hasNext()) {
                Arc arc = (Arc) neighbor.next();
                Bond b = (Bond) arc.getElement();
                if ((!cycle.contains(arc) && b.isDouble()) || (!cycle.contains(arc) && b.isTriple())) {
                    boolean classifiedThisDouble = false;
                    // this is a exocycle, find out if it is part of any other cycle
                    for (int k = 0 ; k < SSSRings.size(); k++) {
                        if(k!=j) { 
                           LinkedList otherCycle = (LinkedList) SSSRings.get(k);
                           if (otherCycle.contains(arc)) {
                               // this cycle contains the arc, now find out if it is aromatic
                               // classifyAsAromatic(k, alreadyClassified);
                               if (isAromatic[k])
                                   aromaticExoPi++;
                               //else
                               //    nonAromaticExoPi++;
                               classifiedThisDouble = true;
                           }
                        }
                    }
                    if (!classifiedThisDouble)
                        nonAromaticExoPi++;
                }
            }
	  }
        }
        int numPiBonds = 0;
        numPiBonds = aromaticExoPi;
        // more thorough screening, for cycles with no exocyclic Pi bonds
        // if (aromaticExoPi == 0 && nonAromaticExoPi == 0 ){
        for (int k = 0; k < cycle.size(); k++) {
            GraphComponent gc = (GraphComponent) cycle.get(k);
            if (gc instanceof Node) {
                Node node = (Node) gc;
                Atom a = (Atom) node.getElement();
                // dont have to check for cationic carbon or boron
                // saturated heteroatoms contribute 2 pi bonds, but only O is the heteroatom
                if (a.isOxygen() && node.getNeighborNumber() == 2)
                    numPiBonds = numPiBonds + 2;
            } else {
                // if a double bond then 2 pi electrons and if a triple bond then 4 pi electrons
                Arc a = (Arc) gc;
                Bond b = (Bond) a.getElement();
                if (b.isDouble() || b.isTriple())
                    numPiBonds = numPiBonds + 2;
            }
        }
	// System.out.print("numPiBonds"+numPiBonds);
        if ( numPiBonds > 2 && (numPiBonds - 2) % 4 == 0) {
            alreadyClassified[j] = 1;
            isAromatic[j] = true;
            System.out.print("Ring is aromatic!");
            j=0;
        } 
	else {
            j++;
        }
        return j;
    }
     
	public void formSSSR() {
        // determine the number of SSSR
        if (SSSRings != null)
            return;
        int numSSSR = getArcNumber() - getNodeNumber() + 1;
        if (numSSSR == 0) {
            acyclic = true;
            return;
        } else {
            acyclic = false;
            SSSRings = new LinkedList();
        }
        // remove the nonRing Nodes
        Graph g = Graph.copy(this);
        boolean moreRemovals = true;
        LinkedHashSet removeNodes = new LinkedHashSet();
        while (moreRemovals) {
            moreRemovals = false;
            Iterator nodeIter = g.getNodeList();
            while (nodeIter.hasNext()) {
                Node n = (Node) nodeIter.next();
                if (n.getNeighborNumber() <= 1)
                    removeNodes.add(n);
            }
            Iterator removeNodesIter = removeNodes.iterator();
            while (removeNodesIter.hasNext()) {
                moreRemovals = true;
                Node n = (Node) removeNodesIter.next();
                g.removeNodeWithoutIDChange(n);
            }
            removeNodes.clear();
        }
        if (numSSSR == 1) {
            // get any node
            LinkedList ring = new LinkedList();
            Iterator nodeIter = g.getNodeList();
            Node n = (Node) nodeIter.next();
            ring.add(0, n);
            int i = 1;
            while (i < (g.getNodeNumber())) {
                Iterator neighborArc = n.getNeighbor();
                Arc arc = (Arc) neighborArc.next();
                Node nextNode = n.getOtherNode(arc);
                if (!ring.contains(nextNode))
                    ring.add(i, nextNode);
                else {
                    arc = (Arc) neighborArc.next();
                    nextNode = n.getOtherNode(arc);
                    ring.add(i, nextNode);
                }
                i++;
                n = nextNode;
            }
            SSSRings.add(ring);
        } else {
            // identify quarternary and tertiary atom centers
            LinkedHashSet quarternary = new LinkedHashSet();
            LinkedHashSet tertiary = new LinkedHashSet();
            Iterator nodeIter = g.getNodeList();
            while (nodeIter.hasNext()) {
                Node n = (Node) nodeIter.next();
                if (n.getNeighborNumber() == 3)
                    tertiary.add(n);
                if (n.getNeighborNumber() == 4)
                    quarternary.add(n);
            }
            Iterator quarternaryIter = quarternary.iterator();
            Iterator tertiaryIter = tertiary.iterator();
            Node n;
            int pathWidth = 0;
            int[] level = new int[4];
            boolean[] grow = new boolean[16];
            int[] numTerminalAtoms = new int[4];
            int nodeNum = g.getNodeNumber();
            while (!((!quarternaryIter.hasNext() && !tertiaryIter.hasNext()) || (SSSRings
                    .size() == numSSSR && allAtomsIncluded(SSSRings, g)))) {
                for (int i = 0; i < 16; i++)
                    grow[i] = true;
                Node[] paths = new Node[4 * nodeNum];
                LinkedList ringsFromNode = new LinkedList();
                if (quarternaryIter.hasNext()) {
                    n = (Node) quarternaryIter.next();
                    pathWidth = 4;
                    Iterator nextNodeIter = n.getNeighboringNodes().iterator();
                    for (int i = 0; i < pathWidth; i++) {
                        paths[i * nodeNum] = n;
                        paths[i * nodeNum + 1] = (Node) nextNodeIter.next();
                        level[i] = 1;
                        grow[i * 4] = true;
                        numTerminalAtoms[i] = 1;
                        checkPresenceOfRings(i, paths, level, pathWidth,
                                nodeNum, ringsFromNode, grow, numTerminalAtoms);
                    }
                } else {
                    n = (Node) tertiaryIter.next();
                    pathWidth = 3;
                    Iterator nextNodeIter = n.getNeighboringNodes().iterator();
                    for (int i = 0; i < pathWidth; i++) {
                        paths[i * nodeNum] = n;
                        paths[i * nodeNum + 1] = (Node) nextNodeIter.next();
                        level[i] = 1;
                        grow[i * 4] = true;
                        numTerminalAtoms[i] = 1;
                        // System.out.print( paths[i*nodeNum].getID()+ " "+ paths[i*nodeNum+1].getID()+" ");
                        checkPresenceOfRings(i, paths, level, pathWidth,
                                nodeNum, ringsFromNode, grow, numTerminalAtoms);
                    }
                }
                int maxlevel = 1;
                boolean allBranched = false;
                int prevlevel = 0;
                int ringsFound = SSSRings.size();
                while (maxlevel == prevlevel + 1) {
                    prevlevel = maxlevel;
                    for (int i = 0; i < pathWidth; i++) {
                        Node nlevel = paths[i * nodeNum + level[i]];
                        if (numTerminalAtoms[i] == 1 && grow[i * 4]) {
                            level[i]++;
                            if (maxlevel == prevlevel)
                                maxlevel++;
                            numTerminalAtoms[i] = 0;
                            Iterator nextNodeIter = nlevel
                                    .getNeighboringNodes().iterator();
                            while (nextNodeIter.hasNext()) {
                                Node nextNode = (Node) nextNodeIter.next();
                                if (nextNode != paths[i * nodeNum + level[i]
                                        - 2]) {
                                    paths[i * nodeNum + level[i]
                                            + numTerminalAtoms[i]] = nextNode;
                                    numTerminalAtoms[i]++;
                                }
                            }
                        }
                        checkPresenceOfRings(i, paths, level, pathWidth,
                                nodeNum, ringsFromNode, grow, numTerminalAtoms);
                    }
                }
                for (int i = 0; i < ringsFromNode.size(); i++) {
                    LinkedList cycle = (LinkedList) ringsFromNode.get(i);
                    if (!cycleAlreadyPresent(cycle, SSSRings))
                        SSSRings.add(cycle);
                }
                if (1 + pathWidth / 2 > ringsFromNode.size()) {
                    // find unfound rings
                    for (int i = 0; i < pathWidth; i++) {
                        for (int numTerm = 0; numTerm < numTerminalAtoms[i]; numTerm++) {
                            Node ni = paths[i * nodeNum + level[i] + numTerm];
                            Iterator niIterator = ni.getNeighboringNodes()
                                    .iterator();
                            while (niIterator.hasNext()) {
                                Node niNeighbor = (Node) niIterator.next();
                                if (niNeighbor != paths[i * nodeNum + level[i]
                                        - 1]) {
                                    for (int j = i + 1; j < pathWidth; j++) {
                                        for (int numTerm2 = 0; numTerm2 < numTerminalAtoms[j]; numTerm2++) {
                                            Node nj = paths[j * nodeNum
                                                    + level[j] + numTerm2];
                                            Iterator njIterator = nj
                                                    .getNeighboringNodes()
                                                    .iterator();
                                            while (njIterator.hasNext()) {
                                                Node njNeighbor = (Node) njIterator
                                                        .next();
                                                if (njNeighbor != paths[j
                                                        * nodeNum + level[j]
                                                        - 1]
                                                        && njNeighbor == niNeighbor
                                                        && (grow[i * 4
                                                                + numTerm] || grow[j
                                                                * 4 + numTerm2])) {
                                                    LinkedList cycle = new LinkedList();
                                                    for (int k = 0; k <= level[i] - 1; k++) {
                                                        cycle.add(paths[i
                                                                * nodeNum + k]);
                                                    }
                                                    cycle.add(ni);
                                                    cycle.add(niNeighbor);
                                                    cycle.add(nj);
                                                    for (int k = level[j] - 1; k > 0; k--)
                                                        cycle.add(paths[j
                                                                * nodeNum + k]);
                                                    if (!cycleAlreadyPresent(
                                                            cycle, SSSRings))
                                                        SSSRings.add(cycle);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (!acyclic) {
            LinkedList SSSRingsCopy = new LinkedList();
            // System.out.println(SSSRings.size()+" Rings found");
            for (int i = 0; i < SSSRings.size(); i++) {
                LinkedList cycle = (LinkedList) SSSRings.get(i);
                // System.out.println("Ring "+i+1+ " size= "+cycle.size());
                LinkedList actualCycle = new LinkedList();
                for (int j = 0; j < cycle.size(); j++) {
                    Node node1 = getNodeAt(((Node) cycle.get(j)).getID());
                    Node node2 = getNodeAt(((Node) cycle.get(0)).getID());
                    if (j != cycle.size() - 1) {
                        node2 = getNodeAt(((Node) cycle.get(j + 1)).getID());
                    }
                    actualCycle.add(node1);
                    actualCycle.add(getArcBetween(node1, node2));
                    // System.out.print(node1.getID()+ " ");
                }
                // System.out.println();
                cycle = null;
                SSSRingsCopy.add(i, actualCycle);
            }
            SSSRings = SSSRingsCopy;
        }
        g = null;
        Iterator iter = SSSRings.iterator();
        while (iter.hasNext()) {
            LinkedList ring = (LinkedList) iter.next();
            Iterator ringIter = ring.iterator();
            while (ringIter.hasNext()) {
                GraphComponent gc = (GraphComponent) ringIter.next();
                gc.setInCycle(true);
            }
        }
        return;
    }

    public boolean allAtomsIncluded(LinkedList SSSRings, Graph g) {
        Iterator nodeIter = g.getNodeList();
        while (nodeIter.hasNext()) {
            Node n = (Node) nodeIter.next();
            boolean nodeFound = false;
            for (int i = 0; i < SSSRings.size(); i++) {
                LinkedList cycle = (LinkedList) SSSRings.get(i);
                if (cycle.contains(n)) {
                    nodeFound = true;
                    break;
                }
            }
            if (!nodeFound)
                return false;
        }
        return true;
    }

    public void checkPresenceOfRings(int i, Node[] paths, int[] level,
            int pathWidth, int nodeNum, LinkedList ringsFromNode,
            boolean[] grow, int[] numTerminalAtoms) {
        // for ( i=0;i<pathWidth;i++){
        for (int numTerm = 0; numTerm < numTerminalAtoms[i]; numTerm++) {
            Node n = (Node) paths[i * nodeNum + level[i] + numTerm];
            if (!grow[i * 4 + numTerm])
                continue;
            Iterator neighborsIter = n.getNeighboringNodes().iterator();
            while (neighborsIter.hasNext()) {
                Node neighbor = (Node) neighborsIter.next();
                if (neighbor == paths[i * nodeNum + level[i] - 1])
                    continue;
                for (int j = 0; j < pathWidth; j++) {
                    if (i == j)
                        continue;
                    for (int numTerm2 = 0; numTerm2 < numTerminalAtoms[j]; numTerm2++) {
                        Node checkNode = (Node) paths[j * nodeNum + level[j]
                                + numTerm2];
                        if ((!grow[j * 4 + numTerm2] && numTerminalAtoms[j] >= 2))
                            continue;
                        if (neighbor == checkNode) {
                            // cycle found
                            grow[i * 4 + numTerm] = false;
                            grow[j * 4 + numTerm2] = false;
                            LinkedList cycle = new LinkedList();
                            for (int k = 0; k <= level[i] - 1; k++) {
                                cycle.add(paths[i * nodeNum + k]);
                            }
                            cycle.add(n);
                            cycle.add(neighbor);
                            for (int k = level[j] - 1; k >= 1; k--) {
                                cycle.add(paths[j * nodeNum + k]);
                            }
                            if (!cycleAlreadyPresent(cycle, ringsFromNode))
                                ringsFromNode.add(cycle);
                        }
                    }
                }
            }
        }
        // }
        return;
    }

    public boolean cycleAlreadyPresent(LinkedList cycle, LinkedList SSSRings) {
        if (SSSRings.size() == 0)
            return false;
        for (int i = 0; i < SSSRings.size(); i++) {
            LinkedList nextCycle = (LinkedList) SSSRings.get(i);
            if (MathTool.isListEqual(nextCycle, cycle))
                return true;
        }
        return false;
    }

    /**
     * Output a string of adjacency list of this graph Modifies:
     */
    // ## operation toString()
    public String toNewString() {
        // #[ operation toString()
        String s = "";
        for (int i = 1; i <= highestNodeID; i++) {
            Node n = getNodeAt(i);
            if (n != null)
                s = s + n.toString() + '\n';
        }
        return s;
        // #]
    }

    /**
     * Return true iff this graph and p_graph is equivalent at center nodes. if no center nodes defined, return true iff
     * this graph is equivalent with p_graph.<br>
     * <b>Modifies</b><br>
     * visited status of nodes and arcs.
     */
    // ## operation equals(Object)
    public boolean equals(Object p_graph) {
        // #[ operation equals(Object)
        if (this == p_graph)
            return true;
        return isEquivalentAtCentralNodes(p_graph);
        // #]
    }

    /**
     * <b>Requires</b><br>
     * there is no duplicated bond between two atoms.<br>
     * <b>Effects</b><br>
     * loop over the arc list to find out the arc whose neighbor set includes the position1 and position2. if found,
     * return this arc; otherwise, return null.
     */
    // ## operation getArcBetween(int,int)
    public Arc getArcBetween(int p_position1, int p_position2)
            throws NotInGraphException {
        // #[ operation getArcBetween(int,int)
        Node node1 = getNodeAt(p_position1);
        Node node2 = getNodeAt(p_position2);
        return getArcBetween(node1, node2);
        // #]
    }

    public boolean[] getIsAromatic() {
        return isAromatic;
    }

    /**
     * If there is an arc connecting the two nodes in this graph, return that arc, otherwise, return null.
     */
    // ## operation getArcBetween(Node,Node)
    public Arc getArcBetween(Node p_node1, Node p_node2) {
        // #[ operation getArcBetween(Node,Node)
        if (p_node1 == null || p_node2 == null) {
            return null;
        }
        if (!contains(p_node1) || !contains(p_node2)) {
            return null;
        }
        Iterator iter = p_node1.getNeighbor();
        while (iter.hasNext()) {
            Arc arc = (Arc) (iter.next());
            if (arc.isConnected(p_node2)) {
                return arc;
            }
        }
        return null;
        // #]
    }

    /**
     * Return the iterator over the arc collection
     */
    // ## operation getArcList()
    public Iterator getArcList() {
        // #[ operation getArcList()
        Iterator iter = arcList.iterator();
        return iter;
        // #]
    }

    /**
     * Return the number of arcs in this Graph
     */
    // ## operation getArcNumber()
    public int getArcNumber() {
        // #[ operation getArcNumber()
        return arcList.size();
        // #]
    }

    /**
     * Return the node with p_centralID as centralID in centralNode list
     */
    // ## operation getCentralNodeAt(int)
    public Node getCentralNodeAt(int p_centralID) {
        // #[ operation getCentralNodeAt(int)
        return getCentralNodeAt(new Integer(p_centralID));
        // #]
    }

    /**
     * Return the node with p_ID as centralID in centralNode list
     */
    // ## operation getCentralNodeAt(Integer)
    public Node getCentralNodeAt(Integer p_ID) {
        // #[ operation getCentralNodeAt(Integer)
        return (Node) centralNode.get(p_ID);
        // #]
    }

    /**
     * Return iterator of center node ID
     */
    // ## operation getCentralNodeKeys()
    public Iterator getCentralNodeKeys() {
        // #[ operation getCentralNodeKeys()
        return centralNode.keySet().iterator();
        // #]
    }

    /**
     * Return iterator over the cetnralNode list of this graph.
     */
    // ## operation getCentralNodeList()
    public Iterator getCentralNodeList() {
        // #[ operation getCentralNodeList()
        Iterator iter = centralNode.values().iterator();
        return iter;
        // #]
    }

    /**
     * Return number of central nodes in this graph.
     */
    // ## operation getCentralNodeNumber()
    public int getCentralNodeNumber() {
        // #[ operation getCentralNodeNumber()
        return centralNode.size();
        // #]
    }

    /**
     * Return the iterator over the cycle collection.
     */
    // ## operation getCycle()
    public LinkedList getCycle() {
        // #[ operation getCycle()
        // Iterator iter = SSSRings.iterator();
        return SSSRings;
        // #]
    }

    public void setCycle(LinkedList p_cycle) {
        SSSRings = p_cycle;
    }

    /**
     * REturn the number of cycles in the graph.
     */
    public int getCycleNumber() {
        if (SSSRings == null)
            formSSSR();
        if (acyclic)
            return 0;
        else
            return SSSRings.size();
    }

    /**
     * Return the highestCenteralID of this graph.
     */
    // ## operation getHighestCentralID()
    public int getHighestCentralID() {
        // #[ operation getHighestCentralID()
        return highestCentralID;
        // #]
    }

    public int getLowestCentralID() {
        return lowestCentralID;
    }

    /**
     * Return the node at the p_position in this graph, if there is no node at that position, or if at this position,
     * there is a null node, throw NotInGraphException.
     */
    // ## operation getNodeAt(int)
    public Node getNodeAt(int p_position) throws NotInGraphException {
        // #[ operation getNodeAt(int)
        Integer id = new Integer(p_position);
        return getNodeAt(id);
        // #]
    }

    /**
     * Return the node at p_ID position in this graph, if there is no node at that position, or if at this position,
     * there is a null node, throw NotInGraphException.
     */
    // ## operation getNodeAt(Integer)
    public Node getNodeAt(Integer p_ID) {
        // #[ operation getNodeAt(Integer)
        return (Node) nodeList.get(p_ID);
        // #]
    }

    // return a LinkedHashSet of the node IDs
    public LinkedHashSet getNodeIDs() {
        Iterator iter = nodeList.values().iterator();
        LinkedHashSet ids = new LinkedHashSet();
        while (iter.hasNext()) {
            Node n = (Node) iter.next();
            ids.add(n.getID());
        }
        return ids;
    }

    /**
     * Return the iterator over the node collection
     */
    // ## operation getNodeList()
    public Iterator getNodeList() {
        // #[ operation getNodeList()
        Iterator iter = nodeList.values().iterator();
        return iter;
        // #]
    }

    /**
     * Return the number of nodes in this Graph
     */
    // ## operation getNodeNumber()
    public int getNodeNumber() {
        // #[ operation getNodeNumber()
        return nodeList.values().size();
        // #]
    }

    /**
     * Compute and return the hashcode of this graph. hashcode = 10^6*(#node) + #arc
     */
    // ## operation hashCode()
    public int hashCode() {
// #[ operation hashCode()
        return getNodeNumber() * 100000 + getArcNumber();
    }

    /**
     * Identify and return all the possible matched pattens between this graph and p_graph. The order of center node
     * doesn't matter.<br>
     * <b>Modifies</b><br>
     * visited status of nodes and arcs in this graph.
     */
    // ## operation identifyAllOrderedMatchedSites(Graph)
    public LinkedHashSet identifyAllOrderedMatchedSites(Graph p_graph) {
        // #[ operation identifyAllOrderedMatchedSites(Graph)
        LinkedHashSet allMatchedSites = new LinkedHashSet();
        clearCentralNode();
        resetMatchedGC();
        p_graph.resetMatchedGC();
        Iterator iter2 = p_graph.getCentralNodeList();
        if (iter2.hasNext()) {
            Node node2 = (Node) iter2.next();
            Iterator iter1 = getNodeList();
            while (iter1.hasNext()) {
                Node node1 = (Node) iter1.next();
                LinkedList matched = new LinkedList();
                matched = node1.identifyAllMatchedSites(node2);
                if (matched != null)
                    allMatchedSites.addAll(matched);
            }
        }
        return allMatchedSites;
        // #]
    }

    /**
     * Identify and return all the possible matched pattens between this graph and p_graph. The order of center node
     * does matter.<br>
     * <b>Modifies</b><br>
     * visited status of nodes and arcs in this graph.
     */
    // ## operation identifyAllUnorderedMatchedSite(Graph)
    public LinkedHashSet identifyAllUnorderedMatchedSite(Graph p_graph) {
        // #[ operation identifyAllUnorderedMatchedSite(Graph)
        LinkedHashSet matchedSite = new LinkedHashSet();
        LinkedHashSet total = new LinkedHashSet();
        int size = p_graph.getCentralNodeNumber();
        Stack stack1 = new Stack();
        Stack stack2 = new Stack();
        Iterator iter2 = p_graph.getCentralNodeList();
        boolean found = false;
        if (iter2.hasNext()) {
            Node node2 = (Node) iter2.next();
            Iterator iter1 = getNodeList();
            while (iter1.hasNext()) {
                clearCentralNode();
                Node node1 = (Node) iter1.next();
                resetMatchedGC();
                p_graph.resetMatchedGC();
                found = node1.isSubAndSetCentralID(node2, stack1, stack2);
                if (found) {
                    refreshCentralNode();
                    // check if the combination of the id of central nodes has been found already
                    LinkedHashSet IdSet = new LinkedHashSet();
                    Iterator iter = getCentralNodeList();
                    while (iter.hasNext()) {
                        Integer id = ((Node) iter.next()).getID();
                        IdSet.add(id);
                    }
                    if (!total.contains(IdSet)) {
                        // if this is a new combination, add to matched site list
                        total.add(IdSet);
                        matchedSite.add((LinkedHashMap) (centralNode.clone()));
                    } else {
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
        if (matchedSite.size() == 0)
            matchedSite = null;
        return matchedSite;
        // #]
    }

    /**
     * <li>identify all the nodes/arcs involved in at least one cycle in this graph</li> <li>add those graph component
     * involved in a ring into cycle collection.</li> <b>Modifies</b><br>
     * the inCycle attributes of all the graph componenets, and this.cycle
     */
    // ## operation identifyCycle()
    /*
     * public void identifyCycle() throws InconnectedGraphException { //#[ operation identifyCycle() if (!isConnected())
     * throw new InconnectedGraphException(); if ((getNodeNumber()-getArcNumber()) == 1) { acyclic = true; return; }
     * cycle = new LinkedList(); setVisited(false); Iterator iter = getNodeList(); while (iter.hasNext()) { Node n =
     * (Node)iter.next(); if (!n.getVisited()) { LinkedList list = new LinkedList(); cycleIdentificationDFS(n,list); } }
     * if (SSSRings.size() == 0) acyclic = true; else acyclic = false; return; //#] }
     */
    /**
     * Identify the type of each node in this graph as Cs, Cd, Ca, Ck, etc.
     */
    // ## operation identifyFgElement()
    public void identifyFgElement() {
        // #[ operation identifyFgElement()
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
            Node node = (Node) iter.next();
            node.identifyFgElement();
        }
        return;
        // #]
    }

    /**
     * Return true iff this graph is acyclic.
     */
    // ## operation isAcyclic()
    public boolean isAcyclic() {
        // #[ operation isAcyclic()
        if (SSSRings == null)
            formSSSR();
        return acyclic;
        // #]
    }

    /**
     * If all the graph components in this graph are connected to each other, i.e., there is no seperated parts in this
     * graph, return true, otherwise, return false.<br>
     * <b>Modifies</b><br>
     * visited status of every graph component
     */
    // ## operation isConnected()
    public boolean isConnected() {
        // #[ operation isConnected()
        if (isEmpty())
            return true;
        setVisited(false);
        Iterator iter = getNodeList();
        Node n = (Node) iter.next();
        deepFirstSearch(n);
        if (!allVisited())
            return false;
        else
            return true;
        // #]
    }

    /**
     * Return true iff there is no node and arc in this graph.
     */
    // ## operation isEmpty()
    public boolean isEmpty() {
        // #[ operation isEmpty()
        if (nodeList.isEmpty() && arcList.isEmpty())
            return true;
        else
            return false;
        // #]
    }

    /**
     * If this graph and the pass-in p_graph are equivalent to each other. No check for centerID.<br>
     * <b>Modifies</b><br>
     * visited status of all the graph component
     */
    // ## operation isEquivalent(Object)
    public boolean isEquivalent(Object p_graph) {
        // #[ operation isEquivalent(Object)
        if (this == p_graph)
            return true;
        if (!(p_graph instanceof Graph))
            return false;
        Graph graph = (Graph) p_graph;
        if (isEmpty() && graph.isEmpty())
            return true;
        if (getNodeNumber() != graph.getNodeNumber()) {
            return false;
        }
        resetMatchedGC();
        graph.resetMatchedGC();
        Iterator iter1 = getNodeList();
        if (iter1.hasNext()) {
            GraphComponent gc1 = (GraphComponent) iter1.next();
            Iterator iter2 = graph.getNodeList();
            while (iter2.hasNext()) {
                GraphComponent gc2 = (GraphComponent) iter2.next();
                Stack stack1 = new Stack();
                Stack stack2 = new Stack();
                if (gc1.isEquivalent(gc2, stack1, stack2))
                    return true;
            }
        }
        return false;
        // #]
    }

    /**
     * If two graphs both have centers, compare equivalence at centers; if two graphs both have no centers, call
     * isEquivalent to check the no-center equivalence; otherwise, return false.<br>
     * <b>Modifies</b><br>
     * visited status of nodes in both graphs
     */
    // ## operation isEquivalentAtCentralNodes(Object)
    public boolean isEquivalentAtCentralNodes(Object p_graph) {
        // #[ operation isEquivalentAtCentralNodes(Object)
        if (this == p_graph)
            return true;
        if (!(p_graph instanceof Graph))
            return false;
        Graph graph = (Graph) p_graph;
        if (isEmpty() && graph.isEmpty())
            return true;
        boolean b1 = centerIsEmpty();
        boolean b2 = graph.centerIsEmpty();
        // if both graphs don't have center, use isEquivalent()
        if (b1 && b2)
            return isEquivalent(graph);
        // if any graph doesn't have center, but the other has, return false
        else if ((b1 && !b2) || (!b1 && b2))
            return false;
        // if both graphs have center, compare center
        else {
            resetMatchedGC();
            graph.resetMatchedGC();
            Iterator iter1 = getCentralNodeList();
            if (!iter1.hasNext())
                return false;
            GraphComponent gc1 = (GraphComponent) (iter1.next());
            Iterator iter2 = graph.getCentralNodeList();
            while (iter2.hasNext()) {
                GraphComponent gc2 = (GraphComponent) iter2.next();
                Stack stack1 = new Stack();
                Stack stack2 = new Stack();
                if (gc1.isEquivalentCenterMatched(gc2, stack1, stack2))
                    return true;
            }
            return false;
        }
        // #]
    }

    /**
     * Return true iff this graph is a subgraph of p_graph, i.e., iff there is a node in this graph is a subNode of one
     * node in p_graph.<br>
     * <b>Modifies</b><br>
     * Visited status of all the nodes in this and p_graph.
     */
    // ## operation isSub(Graph)
    public boolean isSub(Graph p_graph) {
        // #[ operation isSub(Graph)
        for (Iterator iter1 = getNodeList(); iter1.hasNext();) {
            Node node1 = (Node) iter1.next();
            for (Iterator iter2 = p_graph.getNodeList(); iter2.hasNext();) {
                Node node2 = (Node) iter2.next();
                Stack stack1 = new Stack();
                Stack stack2 = new Stack();
                resetMatchedGC();
                p_graph.resetMatchedGC();
                if (node1.isSub(node2, stack1, stack2))
                    return true;
            }
        }
        return false;
        // #]
    }

    /**
     * Return true iff this graph is a subgraph of p_graph at central nodes, i.e., iff a central node in this graph is a
     * subNode of one central node in p_graph.<br>
     * <b>Modifies</b><br>
     * visited status of all the nodes in this and p_graph.
     */
    // ## operation isSubAtCentralNodes(Graph)
    public boolean isSubAtCentralNodes(Graph p_graph) {
        // #[ operation isSubAtCentralNodes(Graph)
        Iterator iter1 = p_graph.getCentralNodeList();
        Stack s1 = new Stack();
        Stack s2 = new Stack();
        if (iter1.hasNext()) {
            Node node1 = (Node) iter1.next();
            // check if there is a matched central node in p_graph
            Integer cid = node1.getCentralID();
            Node node2 = getCentralNodeAt(cid);
            if (node2 == null)
                return false;
            resetMatchedGC();
            p_graph.resetMatchedGC();
            return node2.isSubCentralMatched(node1, s1, s2);
        }
        // there is no central node set in this graph, can't make comparison, return false
        else {
            return false;
        }
        // #]
    }

    /**
     * If this graph is connected, do nothing; if the graph is not connected, partioning the graph into a list of
     * graphs.<br>
     * <b>Modifies</b><br>
     * this
     */
    // ## operation partition()
    public LinkedList partition() throws NotPartitionedException {
        // #[ operation partition()
        LinkedList result = new LinkedList();
        if (isEmpty()) {
            return result;
        }
        setVisited(false);
        Iterator iter = getNodeList();
        deepFirstSearch((Node) iter.next());
        while (true) {
            Graph newGraph = new Graph();
            int nodeID = 0;
            int centralID, highestCentralID = -1, lowestCentralID = 10000;
            Iterator iter1 = nodeList.keySet().iterator();
            while (iter1.hasNext()) {
                Integer ID = (Integer) iter1.next();
                Node n = (Node) nodeList.get(ID);
                if (n.isVisited()) {
                    nodeID++;
                    centralID = n.getCentralID().intValue();
                    if (centralID >= 0) {
                        if (highestCentralID < centralID)
                            highestCentralID = centralID;
                        if (lowestCentralID > centralID)
                            lowestCentralID = centralID;
                    }
                    iter1.remove();
                    if (centralID >= 0)
                        centralNode.remove(n.getCentralID());
                    n.setID(nodeID);
                    n.setCentralID(centralID);
                    newGraph.nodeList.put(n.getID(), n);
                    if (centralID >= 0)
                        newGraph.centralNode.put(n.getCentralID(), n);
                }
            }
            newGraph.highestNodeID = nodeID;
            newGraph.highestCentralID = highestCentralID;
            newGraph.lowestCentralID = lowestCentralID;
            Iterator iter2 = getArcList();
            while (iter2.hasNext()) {
                Arc a = (Arc) iter2.next();
                if (a.isVisited()) {
                    newGraph.arcList.add(a);
                    iter2.remove();
                }
            }
            result.add(newGraph);
            if (nodeList.isEmpty())
                break;
            iter = getNodeList();
            deepFirstSearch((Node) iter.next());
        }
        if (result.size() == 0) {
            throw new NotPartitionedException();
        } else
            return result;
        // #]
    }

    // same as partition, above, except the node IDs are preserved
    // highestNodeID is NOT guaranteed to be stored correctly
    /**
     * If this graph is connected, do nothing; if the graph is not connected, partioning the graph into a list of
     * graphs.<br>
     * <b>Modifies</b><br>
     * this
     */
    // ## operation partition()
    public LinkedList partitionWithPreservedIDs()
            throws NotPartitionedException {
        // #[ operation partition()
        LinkedList result = new LinkedList();
        if (isEmpty()) {
            return result;
        }
        setVisited(false);
        Iterator iter = getNodeList();
        deepFirstSearch((Node) iter.next());
        while (true) {
            Graph newGraph = new Graph();
            int centralID, highestCentralID = -1, lowestCentralID = 10000;
            Iterator iter1 = nodeList.keySet().iterator();
            while (iter1.hasNext()) {
                Integer ID = (Integer) iter1.next();
                Node n = (Node) nodeList.get(ID);
                if (n.isVisited()) {
                    centralID = n.getCentralID().intValue();
                    if (centralID >= 0) {
                        if (highestCentralID < centralID)
                            highestCentralID = centralID;
                        if (lowestCentralID > centralID)
                            lowestCentralID = centralID;
                    }
                    iter1.remove();
                    if (centralID >= 0)
                        centralNode.remove(n.getCentralID());
                    n.setID(ID);
                    n.setCentralID(centralID);
                    newGraph.nodeList.put(n.getID(), n);
                    if (centralID >= 0)
                        newGraph.centralNode.put(n.getCentralID(), n);
                }
            }
            newGraph.highestNodeID = this.highestNodeID;// this will be correct for one of the resulting "pieces" but
// not for the others; in any case
            newGraph.highestCentralID = highestCentralID;
            newGraph.lowestCentralID = lowestCentralID;
            Iterator iter2 = getArcList();
            while (iter2.hasNext()) {
                Arc a = (Arc) iter2.next();
                if (a.isVisited()) {
                    newGraph.arcList.add(a);
                    iter2.remove();
                }
            }
            result.add(newGraph);
            if (nodeList.isEmpty())
                break;
            iter = getNodeList();
            deepFirstSearch((Node) iter.next());
        }
        if (result.size() == 0) {
            throw new NotPartitionedException();
        } else
            return result;
        // #]
    }

    /**
     * Clear central node list, and reform it according to centralID information of node in the graph<br>
     * <b>Modifies</b><br>
     * this.centralNodeList
     */
    // ## operation refreshCentralNode()
    public void refreshCentralNode() {
        // #[ operation refreshCentralNode()
        centralNode.clear();
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
            Node node = (Node) iter.next();
            Integer centralID = node.getCentralID();
            if (centralID.intValue() > 0)
                centralNode.put(centralID, node);
        }
        // #]
    }

    /**
     * Refresh centralNode list according to the nodeID setting of each individual node in nodeList<br>
     * <b>Modifies</b><br>
     * this.centralNode
     */
    // ## operation refreshHighestNodeID()
    public void refreshHighestNodeID() {
        // #[ operation refreshHighestNodeID()
        int hID = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
            Node node = (Node) iter.next();
            int ID = node.getID().intValue();
            if (ID > hID)
                hID = ID;
        }
        highestNodeID = hID;
        // #]
    }

    /**
     * If the pass-in arc is in this graph, remove it from this graph, also get rid of all the connectivity between this
     * arc and all it neighbor nodes. if it is not in the graph, throw NotInGraphException.<br>
     * <b>Modifies</b><br>
     * this
     */
    // ## operation removeArc(Arc)
    public void removeArc(Arc p_arc) {
        // #[ operation removeArc(Arc)
        if (!contains(p_arc))
            throw new NotInGraphException();
        // remove p_arc from all the node's neighbor collection if the node connect to p_arc
        Iterator iter = p_arc.getNeighbor();
        while (iter.hasNext()) {
            Node node = (Node) iter.next();
            node.neighbor.remove(p_arc);
        }
        arcList.remove(p_arc);
        p_arc = null;
        // #]
    }

    /**
     * Remove arc between the two nodes at p_position1 and p_position2, iff that arc exists.<br>
     * <b>Modifies</b><br>
     * this.arcList
     */
    // ## operation removeArcAt(int,int)
    public void removeArcAt(int p_position1, int p_position2) {
        // #[ operation removeArcAt(int,int)
        try {
            Arc arc = getArcBetween(p_position1, p_position2);
            removeArc(arc);
        } catch (NotInGraphException e) {
            return;
        }
        // #]
    }

    /**
     * If p_node is in this graph, remove it, and remove all the arcs linked to it. if p_node is not in this graph,
     * throw NotInGraphException.<br>
     * <b>Modifies</b><br>
     * this.nodeList
     */
    // ## operation removeNode(Node)
    public void removeNode(Node p_node) {
        // #[ operation removeNode(Node)
        if (!contains(p_node))
            throw new NotInGraphException();
        Stack stack = new Stack();
        // find all arc linked to this removed node
        Iterator iter = p_node.getNeighbor();
        while (iter.hasNext()) {
            Arc arc = (Arc) iter.next();
            Node otherNode = p_node.getOtherNode(arc);
            if (!contains(otherNode))
                throw new InvalidConnectivityException();
            otherNode.neighbor.remove(arc);
            arcList.remove(arc);
        }
        // remove this node
        nodeList.remove(p_node.getID());
        if (p_node.getID().intValue() >= highestNodeID) {
            refreshHighestNodeID();
        }
        p_node = null;
        // #]
    }

    /**
     * If p_node is in this graph, remove it, and remove all the arcs linked to it. if p_node is not in this graph,
     * throw NotInGraphException.<br>
     * <b>Modifies</b><br>
     * this.nodeList
     */
    // ## operation removeNode(Node)
    public void removeNodeWithoutIDChange(Node p_node) {
        // #[ operation removeNode(Node)
        if (!contains(p_node))
            throw new NotInGraphException();
        Stack stack = new Stack();
        // find all arc linked to this removed node
        Iterator iter = p_node.getNeighbor();
        while (iter.hasNext()) {
            Arc arc = (Arc) iter.next();
            Node otherNode = p_node.getOtherNode(arc);
            if (!contains(otherNode))
                throw new InvalidConnectivityException();
            otherNode.neighbor.remove(arc);
            arcList.remove(arc);
        }
        // remove this node
        nodeList.remove(p_node.getID());
        p_node = null;
        // #]
    }

    /**
     * Remove the node at the p_position, iff it exists<br>
     * <b>Modifies</b><br>
     * this.nodeList
     */
    // ## operation removeNodeAt(int)
    public void removeNodeAt(int p_position) {
        // #[ operation removeNodeAt(int)
        try {
            Node node = getNodeAt(p_position);
            removeNode(node);
        } catch (NotInGraphException e) {
            return;
        }
        // #]
    }

    /**
     * Return true. there is no restrict for graph invalidity
     */
    // ## operation repOk()
    public boolean repOk() {
        // #[ operation repOk()
        return true;
        // #]
    }

    /**
     * Reset all the matched GC for eac node and arc in the graph to null, to prepare for a new matching/equivalent
     * testing<br>
     * <b>Modifies</b><br>
     * matchedGC of every node and arc in this graph
     */
    // ## operation resetMatchedGC()
    public void resetMatchedGC() {
        // #[ operation resetMatchedGC()
        Iterator iter1 = getArcList();
        while (iter1.hasNext()) {
            ((GraphComponent) iter1.next()).setMatchedGC(null);
        }
        Iterator iter2 = getNodeList();
        while (iter2.hasNext()) {
            ((GraphComponent) iter2.next()).setMatchedGC(null);
        }
        // #]
    }

    /**
     * If the pass-in node in this graph, set it as the ith central node of this graph, where i = p_position; if the
     * node is not in graph, throw NotInGraphException.<br>
     * <b>Modifies</b><br>
     * centralID of the node in this graph, this.centralNode list
     */
    // ## operation setCentralNode(int,Node)
    public void setCentralNode(int p_position, Node p_node)
            throws NotInGraphException {
        // #[ operation setCentralNode(int,Node)
        if (!contains(p_node)) {
            throw new NotInGraphException();
        }
        p_node.setCentralID(p_position);
        Node old = (Node) centralNode.put(new Integer(p_position), p_node);
        updateHighestCentralID(p_position);
        updateLowestCentralID(p_position);
        if (old != null)
            old.setCentralID(-1);
        return;
        // #]
    }

    /**
     * Set all the current central nodes as normal ones, and reset the centralID and centralNode list according to the
     * pass-in p_centralNode LinkedHashMap.<br>
     * <b>Modifies</b><br>
     * this.centralNode, centralID of all the old and new central nodes
     */
    // ## operation setCentralNodes(LinkedHashMap)
    public void setCentralNodes(LinkedHashMap p_centralNode) {
        // #[ operation setCentralNodes(LinkedHashMap)
        clearCentralNode();
        Iterator iter = p_centralNode.keySet().iterator();
        while (iter.hasNext()) {
            Integer centralID = (Integer) iter.next();
            Node node = (Node) (p_centralNode.get(centralID));
            if (!contains(node))
                throw new NotInGraphException();
            node.setCentralID(centralID);
            Node old = (Node) centralNode.put(centralID, node);
            updateHighestCentralID(centralID.intValue());
            updateLowestCentralID(centralID.intValue());
            if (old != null)
                old.setCentralID(-1);
        }
        // #]
    }

   /**
   * Finds the shortest number of bonds between 2 nodes. Must enter a starting distance of 0.
   */
    public int minimumDistance(Node node1, Node node2) {
            LinkedHashSet pathlist = new LinkedHashSet();
            return node1.minimumNumBonds(node2, 0, pathlist);
            // not sure if we need to reset the visited status of all the nodes and arcs?
    }

    public LinkedHashSet minimumPath(Node node1, Node node2) {
	    LinkedHashSet pathlist = new LinkedHashSet();
	    return node1.minimumPath(node2, pathlist);
    }

    public int countCyclicsAlongMinPathInSameRing (Node node1, Node node2) {
	    LinkedHashSet pathlist = new LinkedHashSet();
	    LinkedHashSet minpath = node1.minimumPath(node2, pathlist);
	    Iterator it = minpath.iterator();
            int ncyclics = 0;
 	    int debug = 0;
            while(it.hasNext()) {
		debug += 1;
		if (sameRing(node1, (Node) it.next())) {
		   ncyclics +=1;
		   }
		}	
	    return ncyclics;
    }

    /**
     * Set the visited status of all the graph components in this graph as the pass-in p_visited.<br>
     * <b>Modifies</b><br>
     * this.nodeList, this.arcList
     */
    // ## operation setVisited(boolean)
    private void setVisited(boolean p_visited) {
        // #[ operation setVisited(boolean)
        Iterator iter1 = getArcList();
        while (iter1.hasNext()) {
            ((GraphComponent) iter1.next()).setVisited(p_visited);
        }
        Iterator iter2 = getNodeList();
        while (iter2.hasNext()) {
            ((GraphComponent) iter2.next()).setVisited(p_visited);
        }
        // #]
    }

    
    
    /**
     * Output a string of adjacency list of this graph Modifies:
     */
    // ## operation toString()
    public String toString() {
        // #[ operation toString()
        String s = "";
        for (int i = 1; i <= highestNodeID; i++) {
            Node n = getNodeAt(i);
            if (n != null)
                s = s + n.toString() + '\n';
        }
        return s;
        // #]
    }

    /**
     * Output a string of adjacency list of this graph without outputting the centralID information of node Modifies:
     */
    // ## operation toStringWithoutCentralID()
    public String toStringWithoutCentralID() {
        // #[ operation toStringWithoutCentralID()
        String s = "";
        for (int i = 1; i <= highestNodeID; i++) {
            Node n = getNodeAt(i);
            if (n != null)
                s = s + n.toStringWithoutCentralID() + '\n';
        }
        return s;
        // #]
    }

    /**
     * Output a string of adjacency list of this graph without outputting the centralID information of node and without
     * outputting H node Modifies:
     */
    // ## operation toStringWithoutCentralIDAndH()
    public String toStringWithoutCentralIDAndH() {
        // #[ operation toStringWithoutCentralIDAndH()
        String s = "";
        for (int i = 1; i <= highestNodeID; i++) {
            Node n = getNodeAt(i);
            if (n != null) {
                Atom a = (Atom) n.getElement();
                if (!a.isHydrogen()) {
                    s = s + n.toStringWithoutCentralIDAndH() + '\n';
                }
            }
        }
        return s;
        // #]
    }

    /**
     * Set highestCentralID to p_centralID iff present highestCentralID is less than p_centralID<br>
     * <b>Modifies</b><br>
     * this.highestCentralID
     */
    // ## operation updateHighestCentralID(int)
    public void updateHighestCentralID(int p_centralID) {
        // #[ operation updateHighestCentralID(int)
        if (p_centralID > highestCentralID)
            highestCentralID = p_centralID;
        return;
        // #]
    }

    public void updateLowestCentralID(int p_centralID) {
        if (p_centralID < lowestCentralID)
            lowestCentralID = p_centralID;
    }

    /**
     * Set highestNodeID = p_position iff p_position > highestNodeID<br>
     * <b>Modifies</b><br>
     * this.highestNodeID
     */
    // ## operation updateHighestNodeID(int)
    public void updateHighestNodeID(int p_position) {
        // #[ operation updateHighestNodeID(int)
        if (p_position > highestNodeID)
            highestNodeID = p_position;
        return;
        // #]
    }

    protected static int getMAXNODENUMBER() {
        return MAXNODENUMBER;
    }

    public LinkedHashMap getCentralNode() {
        return centralNode;
    }

    public int getHighestNodeID() {
        return highestNodeID;
    }

    public void clearArcList() {
        arcList.clear();
    }

    /**
     * returns a set of Node's of the ChemGraph that belong to two or more fused rings. A fused ring system is a set of
     * rings that share at least one common bond. An example is decalin, with two fused cyclohexane rings.
     * 
     * @return
     */
    public Set<Node> getFusedRingAtoms() {
        if (SSSRings == null) {
            return null;
        } else if (acyclic) {
            return null;
        } else if (getCycleNumber() == 1) {// monocyclic case
            return null;
        } else {
            Set<Node> fusedRingAtoms = new LinkedHashSet<Node>();
            // First retrieve lists with all the Nodes that belong to rings:
            List<Set<Node>> ringNodesList = new LinkedList<Set<Node>>();
            Iterator iterCycle = SSSRings.iterator(); // loop over all SSSRs
            while (iterCycle.hasNext()) {
                Set<Node> nodesSet = new LinkedHashSet<Node>();
                LinkedList cycle = (LinkedList) iterCycle.next();
                Iterator iter = cycle.iterator();
                while (iter.hasNext()) {
                    GraphComponent node = (GraphComponent) iter.next();
                    if (node instanceof Node) {
                        nodesSet.add((Node) node);
                    }
                }
                ringNodesList.add(nodesSet);
            }
            // Now look for Nodes that belong to more than one Set of Nodes = fused ring atoms
            for (int index = 0; index < ringNodesList.size() - 1; index++) {// end at the one before last
                Set<Node> set = ringNodesList.get(index);
                for (Node node : set) {
                    if (ringNodesList.get(index + 1).contains(node)) {
                        fusedRingAtoms.add(node);
                    }
                }
            }
            if (fusedRingAtoms.isEmpty()) {
                return null;
            } else {
                return fusedRingAtoms;
            }
        }
    }

    /**
     * Returns the atoms belonging to one or more rings. This is retrieved by iterating over all nodes belonging to a
     * cycle and incrementing a Set of nodes corresponding to the atoms in the molecule. Since the java.util.Set type
     * only allows unique elements, adding an already existing reference to a node will not be permitted. <BR>
     * <BR>
     * Keep in mind that just adding the number of atoms in each of the SSSRings fails (overestimates) for fused ring
     * systems, where a single atom can belong to multiple rings at the time.
     * 
     * @return
     */
    public Set<Node> getRingAtoms() {
        Set<Node> ringAtoms = new LinkedHashSet<Node>();
        List<Set<Node>> cycleNodes = getCycleNodes();
        for (Set set : cycleNodes)
            for (Iterator iter = set.iterator(); iter.hasNext();) {
                Node node = (Node) iter.next();
                ringAtoms.add(node);
            }
        return ringAtoms;
    }

    /**
     * Returns a list with sets of Nodes that belong to SSS Rings.
     * 
     * @return
     */
    public List<Set<Node>> getCycleNodes() {
        if (SSSRings == null) {
            return null;
        } else if (acyclic) {
            return null;
        } else {
            // First retrieve lists with all the Nodes that belong to rings:
            List<Set<Node>> ringNodesList = new LinkedList<Set<Node>>();
            Iterator iterCycle = SSSRings.iterator(); // loop over all SSSRs
            while (iterCycle.hasNext()) {
                Set<Node> nodesSet = new LinkedHashSet<Node>();
                LinkedList cycle = (LinkedList) iterCycle.next();
                Iterator iter = cycle.iterator();
                while (iter.hasNext()) {
                    GraphComponent node = (GraphComponent) iter.next();
                    if (node instanceof Node) {
                        nodesSet.add((Node) node);
                    }
                }
                ringNodesList.add(nodesSet);
            }
            return ringNodesList;
        }
    }

    /**
     * Iterates over the SSSR's and searches for rings that solely consist of Cb atoms.<BR>
     * <BR>
     * The algorithms iterates over all nodes of each ring and then iterates over all neighbouring bonds.<BR>
     * <BR>
     * If at least 2 "B" bonds have been found, the atom must be an aromatic atom.<BR>
     * <BR>
     * If all atoms in the ring are aromatic atoms, then an aromatic ring is found, and true is returned.<BR>
     * <BR>
     * Efficiency is increased by first checking the number of nodes in each ring. If this number is different from 6,
     * then skip the ring. This assumes that aromatic rings always consist of 6 atoms. In cases of naphthalene, where a
     * 6 + 4 aromatic system exists, there will be at least one 6 membered aromatic ring so this algorithm will not fail
     * for fused aromatic rings.
     * 
     * @return
     */
    public boolean containsAromaticRing() {
        if (SSSRings == null || acyclic) {
            return false;
        } else {
            Iterator iterCycle = SSSRings.iterator(); // loop over all SSSRs
            while (iterCycle.hasNext()) {
                Set<Node> nodesSet = new LinkedHashSet<Node>();
                List cycle = (LinkedList) iterCycle.next();
                Iterator iter = cycle.iterator();
                int numberOfNodes = 0;
                while (iter.hasNext()) {// count number of node elements
                    if ((GraphComponent) iter.next() instanceof Node)
                        numberOfNodes++;
                }
                if (numberOfNodes == 6) {// Assume aromatic rings always count 6 atoms:
                    iter = cycle.iterator();
                    boolean isAromaticRing = true;
                    while (iter.hasNext() && isAromaticRing) {
                        GraphComponent node = (GraphComponent) iter.next();
                        if (node instanceof Node) {
                            node = (Node) node;
                            Iterator neighbours = node.getNeighbor();
                            int counter = 0;
                            while (neighbours.hasNext() && counter < 2) {
                                GraphComponent neighbour = (GraphComponent) neighbours
                                        .next();
                                if (neighbour instanceof Arc) {
                                    neighbour = (Arc) neighbour;
                                    if (((Bond) neighbour.getElement())
                                            .getName().equals("B")) {
                                        counter++;
                                    }
                                }
                            }
                            if (counter < 2) {// this is not an aromatic atom
                                isAromaticRing = false;
                            }
                        }
                    }
                    if (isAromaticRing) {// as soon as one aromatic ring is found, return true
                        return true;
                    }
                }
            }
        }
        return false;
    }
}
/*********************************************************************
 * File Path : RMG\RMG\jing\chemUtil\Graph.java
 *********************************************************************/

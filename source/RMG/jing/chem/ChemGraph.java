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


import java.io.*;
import java.util.*;
import jing.chemUtil.*;
import jing.chemParser.*;
import jing.chemUtil.Arc;
import jing.chemUtil.Node;
import jing.chemUtil.Graph;
import jing.param.Global;
import jing.param.Temperature;
import jing.rxnSys.ReactionModelGenerator;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\ChemGraph.java
//----------------------------------------------------------------------------

//## class ChemGraph
public class ChemGraph implements Matchable {

    protected static int MAX_OXYGEN_NUM = 10; //20 Modified by AJ		//## attribute MAX_OXYGEN_NUM
	protected static  int MAX_CARBON_NUM = 30;//100 Modified by AJ       //SS
	protected static int MAX_CYCLE_NUM = 10;		//SS (no fused rings); gmagoon: to turn fused rings off, set this to 1 (will exclude multiple non-fused rings, as well, I think)

	/**
    Maximal radical number allowed in a ChemGraph.
    */
    protected static int MAX_RADICAL_NUM = 10;		//## attribute MAX_RADICAL_NUM
    
    protected static int MAX_SILICON_NUM = 10;
    protected static int MAX_SULFUR_NUM = 10;
    
    protected static int MAX_HEAVYATOM_NUM = 100;
    protected static String repOkString = null;
    
    public static boolean useQM = false;//gmagoon 6/15/09: flag for thermo estimation using quantum results; there may be a better place for this (Global?) but for now, this should work
    public static boolean useQMonCyclicsOnly=false;
    /**
    Chemical Formula of a ChemGraph.
    */
    protected String chemicalFormula = null;		//## attribute chemicalFormula

    /**
    The overall forbidden structure.  When any new ChemGraph instance is generated, RMG check if it has any of the forbidden structure.  If it has, it wont be generated.
    */
    protected static HashSet forbiddenStructure = new HashSet();		//## attribute forbiddenStructure

    protected int internalRotor = -1;		//## attribute internalRotor

    /**
    A collection of all the possible symmetry Axis in a ChemGraph.
    For example: the C=C=O skeleton in (CH3)2C=C=O
    */
    //protected HashSet symmetryAxis = null;		//## attribute symmetryAxis

    /**
    Symmetry number of a ChemGraph.  Used in calculating entropy.
    */
    protected int symmetryNumber = -1;		//## attribute symmetryNumber

    /**
    It is a unique string representation for this chem structure.  The method to generating it has not implemented yet.
    */
    protected String uniqueString;		//## attribute uniqueString

    protected Graph graph;
    protected Species species;
    protected ThermoData thermoData;
    protected AbramData abramData;
    protected UnifacData unifacData;
    protected GeneralGAPP thermoGAPP;
    protected GeneralSolvationGAPP SolvationGAPP;
    protected GeneralAbramGAPP abramGAPP;
    protected GeneralUnifacGAPP unifacGAPP;
    protected ThermoData solvthermoData;
    protected TransportData transportData;
    protected GATransportP transportGAPP;    

    protected boolean fromprimarythermolibrary = false;
    protected boolean isAromatic = false;
    protected String InChI;
    protected String InChIKey;
    // Constructors

    private  ChemGraph() {
    }
	
    private  ChemGraph(Graph p_graph) throws ForbiddenStructureException {
        graph = p_graph;
       // isAromatic = isAromatic();

        if (isForbiddenStructure(p_graph,getRadicalNumber(),getOxygenNumber(),getCarbonNumber()) || getRadicalNumber() > MAX_RADICAL_NUM || getOxygenNumber() > MAX_OXYGEN_NUM || getCycleNumber() > MAX_CYCLE_NUM) {
		//if (getRadicalNumber() > MAX_RADICAL_NUM || getOxygenNumber() > MAX_OXYGEN_NUM || getCycleNumber() > MAX_CYCLE_NUM) {		        
			graph = null;
			String message = p_graph.toString() + " is forbidden by "+whichForbiddenStructures(p_graph) +"and not allowed.";
        	throw new ForbiddenStructureException(message);
        }
    }

    /*private boolean isAromatic() {
		//first see if it is cyclic
    	if (graph.getCycleNumber() == 0)
    		return false;
    	
    	//if cyclic then iterate over all the cycles
    	Iterator cyclesIter = graph.getCycle();
    	while (cyclesIter.hasNext()){
    		LinkedList cycle = (LinkedList)cyclesIter.next();
    		boolean hasdoublebond = false;
    		boolean quarternaryAtom = false;
    		boolean monoRadical = false;
    		int unsaturatedCarbon =0 ;
    		for (int i=0; i<=cycle.size(); i++){
    			GraphComponent gc = (GraphComponent)cycle.get(i);
    			
    			
    			if (gc instanceof Arc){
    				if ( ((Bond)((Arc)gc).getElement()).isDouble())
    					hasdoublebond = true;
    			}
    			else {
    				
    				Atom a = (Atom)((Node)gc).getElement();
    				
//    				is it unsaturate Carbon
    				if (a.isCarbon() && !a.isRadical())
    					unsaturatedCarbon++;
    				
    				//it is a monoradical
    				if (a.freeElectron.order == 1){
    					monoRadical = true;
    					return false;
    				}
    				
    				//you have to check for SO2RR radical presence 
    				///what is a quarternary atom...i will check
    			}
    		}
    		if (!hasdoublebond || unsaturatedCarbon > 1)
    			return false;
    		
    		//first check if the ring has exocyclic pi bonds
    		for (int i=0; i<=cycle.size();i++){
    			GraphComponent gc = (GraphComponent)cycle.get(i);
    			if (gc instanceof Node){
    				Iterator neighbors = gc.getNeighbor();
    				while (neighbors.hasNext()){
    					Arc neighborArc = (Arc)neighbors.next();
    					
    				}
    				
    			}
    		}
    	}
		
	}*/
	
	
	/**
    Requires:
    Effects: saturate all the node atom's undefined valence by adding Hydrogen.
    Modifies: this.graph.nodeList
    */
    //## operation addMissingHydrogen()
    public void addMissingHydrogen() {
        //#[ operation addMissingHydrogen()
        Atom H = Atom.make(ChemElement.make("H"), FreeElectron.make("0"));
        Bond S = Bond.make("S");
        LinkedHashMap addedH = new LinkedHashMap();

        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();
        	int val = (int)atom.getValency();

        	double bondOrder = 0;
        	Iterator neighbor_iter = node.getNeighbor();
        	while (neighbor_iter.hasNext()) {
        		Arc arc = (Arc)neighbor_iter.next();
        		Bond bond = (Bond)arc.getElement();
        		bondOrder += bond.getOrder();
        	}
//        	if (bondOrder > val) throw new InvalidConnectivityException();
//        	else if (bondOrder < val) {
//        		addedH.put(node, new Integer(val-bondOrder));
//        	}
        	if (bondOrder < val) {
        		addedH.put(node, new Integer(val-(int)(bondOrder+1.0e-8)));
        	}
        }
        Graph g = getGraph();
        iter = addedH.keySet().iterator();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	int Hnum = ((Integer)addedH.get(node)).intValue();
        	for (int i=0;i<Hnum; i++) {
        		Node newNode = g.addNode(H);
        		g.addArcBetween(node, S, newNode);
        	}
        	node.updateFgElement();
        }

        return;
        //#]
    }

    public static ChemGraph saturate(ChemGraph p_chemGraph) {
       //#[ operation saturate(ChemGraph)
       int max_radNum_molecule = ChemGraph.getMAX_RADICAL_NUM();
       int max_radNum_atom = Math.min(8,max_radNum_molecule);

       ChemGraph result = null;
       try {
               result = ChemGraph.copy(p_chemGraph);
       }
       catch (Exception e) {
               System.out.println(e.getMessage());
               System.exit(0);
       }

       FreeElectron satuated = FreeElectron.make("0");
       Atom H = Atom.make(ChemElement.make("H"),satuated);
       Bond S = Bond.make("S");
       Graph g = result.getGraph();
       int nn = g.getHighestNodeID();
       for (int i = 0 ; i < nn; i++) {
               Node node = g.getNodeAt(i);
               if (node != null) {
                       Atom atom = (Atom)node.getElement();
                       int HNum = atom.getRadicalNumber();
                       if (atom.isRadical()) {
                               Atom newAtom = new Atom(atom.getChemElement(),satuated);
                               node.setElement(newAtom);
                               node.updateFeElement();
                               for (int j = 0; j < HNum; j++) {
                                       Node n = g.addNode(H);
                                       g.addArcBetween(node,S,n);
                               }
                               node.updateFgElement();
                       }
               }
       }

       return result;
       /*
       int max_radNum_molecule = ChemGraph.getMAX_RADICAL_NUM();
       int max_radNum_atom = Math.min(8,max_radNum_molecule);
       int [] idArray = new int[max_radNum_molecule];
       Atom []  atomArray = new Atom[max_radNum_molecule];
       Node [][] newnode = new Node[max_radNum_molecule][max_radNum_atom];

       int radicalSite = 0;
       Iterator iter = p_chemGraph.getNodeList();
       FreeElectron satuated = FreeElectron.make("0");
       while (iter.hasNext()) {
               Node node = (Node)iter.next();
               Atom atom = (Atom)node.getElement();
               if (atom.isRadical()) {
                       radicalSite ++;
                       // save the old radical atom
                       idArray[radicalSite-1] = node.getID().intValue();
                       atomArray[radicalSite-1] = atom;
                       // new a satuated atom and replace the old one
                       Atom newAtom = new Atom(atom.getChemElement(),satuated);
                       node.setElement(newAtom);
                       node.updateFeElement();
               }
       }

       // add H to satuate chem graph
       Atom H = Atom.make(ChemElement.make("H"),satuated);
       Bond S = Bond.make("S");
       for (int i=0;i<radicalSite;i++) {
               Node node = p_chemGraph.getNodeAt(idArray[i]);
               Atom atom = atomArray[i];
               int HNum = atom.getRadicalNumber();
               for (int j=0;j<HNum;j++) {
                       newnode[i][j] = g.addNode(H);
                       g.addArcBetween(node,S,newnode[i][j]);
               }
               node.updateFgElement();
       }
       */
       //#]
   }

    /**
    Requires: acyclic ChemGraph
    Effects: calculate and return the symmetry number centered at p_node atom
    Modifies:
    */
    //## operation calculateAtomSymmetryNumber(Node)
    public int calculateAtomSymmetryNumber(Node p_node) {
        //#[ operation calculateAtomSymmetryNumber(Node)
        // note: acyclic structure!!!!!!!!!!!!!
        int sn = 1;

        // if no neighbor or only one neighbor, sigma = 1, return 1;
        int neighborNumber = p_node.getNeighborNumber();
        if (neighborNumber < 2) return sn;

        Atom atom = (Atom)p_node.getElement();
        Iterator neighbor_iter = p_node.getNeighbor();
        FGElement fge = (FGElement)p_node.getFgElement();

        if (!atom.isRadical()) {
        // satuated atom symmetric number calculation
        	if (fge.equals(FGElement.make("Cs")) || fge.equals(FGElement.make("Sis"))) {
        	// Cs:
               	Arc a1 = (Arc)neighbor_iter.next();
               	Arc a2 = (Arc)neighbor_iter.next();
               	Arc a3 = (Arc)neighbor_iter.next();
               	Arc a4 = (Arc)neighbor_iter.next();
               	if (p_node.isSymmetric(a1,a2)) {
               		if (p_node.isSymmetric(a1,a3)) {
          				if (p_node.isSymmetric(a1,a4)) {
        				// AAAA
        					sn *= 12;
               			}
               			else {
               			// AAAB
         					sn *= 3;
               			}
               		}
               		else {
          				if (p_node.isSymmetric(a1,a4)) {
        				// AAAB
               				sn *= 3;
               			}
               			else if (p_node.isSymmetric(a3,a4)) {
               			// AABB
               				sn *= 2;
               			}
               			else {
        				// AABC
               			}
               		}
               	}
               	else {
               		if (p_node.isSymmetric(a1,a3)) {
               			if (p_node.isSymmetric(a1,a4)) {
               			// AAAB
               				sn *= 3;
               			}
               			else if (p_node.isSymmetric(a2,a4)) {
               			// AABB
               				sn *= 2;
               			}
               			else {
               			// AABC
               			}
               		}
               		else if (p_node.isSymmetric(a2,a3)) {
               			if (p_node.isSymmetric(a1,a4)) {
               			// AABB
               				sn *= 2;
               			}
               			else if (p_node.isSymmetric(a2,a4)) {
               			// AAAB
               				sn *= 3;
               			}
               			else {
               			// AABC
               			}
               		}
               		else {
               			// AABC or ABCD
               		}
               	}
        	}
        	else if (fge.equals(FGElement.make("Os")) || fge.equals(FGElement.make("Ss"))) {
        	// Os:
               	Arc a1 = (Arc)neighbor_iter.next();
               	Arc a2 = (Arc)neighbor_iter.next();
               	if (p_node.isSymmetric(a1,a2)) {
               		sn *= 2;
               	}
        	}
        	else if (fge.equals(FGElement.make("Cdd")) || fge.equals(FGElement.make("Sidd"))) {
        	// Cdd:
               	Arc a1 = (Arc)neighbor_iter.next();
               	Arc a2 = (Arc)neighbor_iter.next();
               	if (p_node.isSymmetric(a1,a2)) {
               		sn *= 2;
                }
          	}
        }
        else {
        // radical symmetric number calculation
        	if (fge.equals(FGElement.make("Cs")) || fge.equals(FGElement.make("Sis"))) {
        	// only consider Cs. and Cs..
        		FreeElectron fe = atom.getFreeElectron();
            	if (fe.getOrder() == 1) {
        		// mono-radical Cs.
                	Arc a1 = (Arc)neighbor_iter.next();
                	Arc a2 = (Arc)neighbor_iter.next();
                	Arc a3 = (Arc)neighbor_iter.next();
                	if (p_node.isSymmetric(a1,a2)) {
                		if (p_node.isSymmetric(a1,a3))
                			sn *= 6;
                		else
                			sn *= 2;
                	}
                	else {
                		if (p_node.isSymmetric(a1,a3) || p_node.isSymmetric(a2,a3))
                			sn *= 2;
                	}
        		}
        		else if (fe.getOrder() == 2) {
        		// bi-radical Cs..
                	Arc a1 = (Arc)neighbor_iter.next();
                	Arc a2 = (Arc)neighbor_iter.next();
                	if (p_node.isSymmetric(a1,a2))
                		sn *= 2;
        		}
        	}
        }

        return sn;
        //#]
    }

    //takes the fragments and corresponding rotor nodes for each side of the rotor
    public int calculateRotorSymmetryNumber(Node p_node1, Node p_node2) {
	//first calculate the symmetry number for each fragment
	int frag1sn=calculateRotorFragmentSymmetryNumber(p_node1);//will return 1, 2, or 3
	int frag2sn=calculateRotorFragmentSymmetryNumber(p_node2);//will return 1, 2, or 3
	if(frag1sn==3 || frag2sn==3){
	    if(frag1sn==3){
		if(frag2sn==3) return 3; //example: ethane
		else if (frag2sn==2) return 6; //example: toluene
		else return 3; //frag2sn==1; example: methanol
	    }
	    else{//frag2sn==3 (and frag1sn!=3)
		if (frag2sn==2) return 6;//see above
		else return 3; //frag1sn==1; see above
	    }
	}
	else if(frag1sn==2 || frag2sn==2){
	    return 2; //see "full" code below
	    //if(frag1sn==2){
		
		//if(frag2sn==2) return 2; //example: biphenyl
		//else return 2; //frag2sn==1; example: phenol
	    //}
	    //else{//frag2sn==2 (and frag1sn=1)
		//return 2;//see above
	    //}
	}
	//otherwise, both frag1sn and frag2sn equal 1

	return 1;
    }

    //returns 1, 2, or 3 based on the rotor fragment and rotor node passed in
    //code is based off of calculateAtomSymmetryNumber as the idea is similar, but we must be able to handle incomplete fragments, triple bonds, and aromatic structures
    //code does not handle radicals at present
    public int calculateRotorFragmentSymmetryNumber(Node p_node) {

	Atom atom = (Atom)p_node.getElement();
        Iterator neighbor_iter = p_node.getNeighbor();
        FGElement fge = (FGElement)p_node.getFgElement();
	if(atom.isRadical()){
	    System.out.println("Error: calculateRotorFragmentSymmetryNumber() does not support radical sites");
	    System.exit(0);
	}
        // if no neighbor or only one neighbor, sigma = 1, return 1;
        int neighborNumber = p_node.getNeighborNumber();
        if (neighborNumber < 2){
	    if(fge.equals(FGElement.make("Ct"))){//triple bond case
		//find the end of the cumulated triple bond system
		Node n1 = p_node;
                LinkedHashSet axis = new LinkedHashSet();
		Arc arc = (Arc)neighbor_iter.next();
                axis.add(arc);
                p_node = getToEndOfCumulatedTripleBondSystem(arc,p_node,axis);
		//use the new, non-Ct node for subsequent calculations
		fge = (FGElement)p_node.getFgElement();
		neighbor_iter = p_node.getNeighbor();
		fge = (FGElement)p_node.getFgElement();
	    }
	    else{//e.g. peroxides, alcohols, etc.
		return 1;
	    }
	}

	if (fge.equals(FGElement.make("Cs")) || fge.equals(FGElement.make("Sis"))) {
		Arc a1 = (Arc)neighbor_iter.next();
		Arc a2 = (Arc)neighbor_iter.next();
		Arc a3 = (Arc)neighbor_iter.next();
		if (p_node.isSymmetric(a1,a2)) {
			if (p_node.isSymmetric(a1,a3))
				return 3;//AAA*//e.g. methyl rotor
		}
	}
	else if(fge.equals(FGElement.make("Cb")) || fge.equals(FGElement.make("Cbf"))){
	    	Arc a1 = (Arc)neighbor_iter.next();
		Arc a2 = (Arc)neighbor_iter.next();
		if (p_node.isSymmetric(a1,a2)) {
		    return 2;//AA*//e.g. phenyl group
		}
	}

        return 1;
    }

    /**
    Requires: acyclic ChemGraph
    Effects: calculate and return the symmetry number by all the possible symmetry axis in this ChemGraph.
    Modifies:
    */
    //## operation calculateAxisSymmetryNumber()
    public int calculateAxisSymmetryNumber() {
        //#[ operation calculateAxisSymmetryNumber()
        int sn = 1;
        // note: acyclic structure!!!!!!!!!!!!!
        LinkedHashSet symmetryAxis = new LinkedHashSet();

        Iterator iter = getArcList();
        while (iter.hasNext()) {
        	Arc arc = (Arc)iter.next();
        	Bond bond = (Bond)arc.getElement();
        	if (bond.isDouble()&&!arc.getInCycle()) {//2/22/10 gmagoon: added check to make sure the arc is not in a cycle; hopefully this should prevent infinite recursion errors in getToEndOfAxis caused by cyclic species like cyclic, cumulenic C3 while still preserving the ability to estimate axis symmetry number for "cyclic" cases where the axis is not part of the cycle
        		Iterator neighbor_iter = arc.getNeighbor();
        		Node n1 = (Node)neighbor_iter.next();
        		Node n2 = (Node)neighbor_iter.next();

//        		FGElement Cd = FGElement.make("Cd");
        		FGElement Cdd = FGElement.make("Cdd");
//        		FGElement Sid = FGElement.make("Sid");
        		FGElement Sidd = FGElement.make("Sidd");

        		FGElement fge1 = (FGElement)n1.getFgElement();
        		FGElement fge2 = (FGElement)n2.getFgElement();

                LinkedHashSet axis = new LinkedHashSet();
                axis.add(arc);
                if (fge1.equals(Cdd) || fge1.equals(Sidd)) n1 = getToEndOfAxis(arc,n1,axis);
                if (fge2.equals(Cdd) || fge2.equals(Sidd)) n2 = getToEndOfAxis(arc,n2,axis);

               	Atom atom1 = (Atom)n1.getElement();
               	Atom atom2 = (Atom)n2.getElement();

               	if (atom1.isRadical() || atom2.isRadical()) return sn;

               	Bond D = Bond.make("D");

               	if (!symmetryAxis.contains(axis)) {
               		symmetryAxis.add(axis);
               		boolean l1 = n1.isLeaf();
               		boolean l2 = n2.isLeaf();
               		if (!l1 && !l2) {
        				Iterator i1 = n1.getNeighbor();
        				Iterator i2 = n2.getNeighbor();
        				Arc a1 = (Arc)i1.next();
        				if (((Bond)a1.getElement()).equals(D)) a1 = (Arc)i1.next();
        				Arc a2 = (Arc)i1.next();
        				if (((Bond)a2.getElement()).equals(D)) a2 = (Arc)i1.next();

        				Arc a3 = (Arc)i2.next();
        				if (((Bond)a3.getElement()).equals(D)) a3 = (Arc)i2.next();
        				Arc a4 = (Arc)i2.next();
        				if (((Bond)a4.getElement()).equals(D)) a4 = (Arc)i2.next();

        				if (n1.isSymmetric(a1,a2) && n2.isSymmetric(a3,a4)) {
        					sn *= 2;
        				}
        			}
        			else if (!l1 && l2) {
        				Iterator i = n1.getNeighbor();
        				Arc a1 = (Arc)i.next();
        				if (((Bond)a1.getElement()).equals(D)) a1 = (Arc)i.next();
        				Arc a2 = (Arc)i.next();
        				if (((Bond)a2.getElement()).equals(D)) a2 = (Arc)i.next();

        				if (n1.isSymmetric(a1,a2)) {
        					sn *= 2;
        				}
        			}
        			else if (l1 && !l2) {
        				Iterator i = n2.getNeighbor();
        				Arc a1 = (Arc)i.next();
        				if (((Bond)a1.getElement()).equals(D)) a1 = (Arc)i.next();
        				Arc a2 = (Arc)i.next();
        				if (((Bond)a2.getElement()).equals(D)) a2 = (Arc)i.next();

        				if (n2.isSymmetric(a1,a2)) {
        					sn *= 2;
        				}
        			}
        		}
        	}
        }

        return sn;
        //#]
    }

    /**
    Requires: acyclic ChemGraph
    Effects: calculate and return the symmetry number centered at p_arc bond
    Modifies:
    */
    //## operation calculateBondSymmetryNumber(Arc)
    public int calculateBondSymmetryNumber(Arc p_arc) {
        //#[ operation calculateBondSymmetryNumber(Arc)
        // note: acyclic structure!!!!!!!!!!!!!

        int sn = 1;

        // if no neighbor or only one neighbor, sigma = 1, return 1;
        int neighborNumber = p_arc.getNeighborNumber();
        if (neighborNumber != 2) throw new InvalidNeighborException("bond has " + neighborNumber + " neighbor!");

        Bond bond = (Bond)p_arc.getElement();
        Iterator neighbor_iter = p_arc.getNeighbor();
        Node n1 = (Node)neighbor_iter.next();
        Node n2 = (Node)neighbor_iter.next();

        if (bond.isSingle() || bond.isDouble() || bond.isTriple()) {
        	if (p_arc.isSymmetric(n1,n2)) {
        		boolean opt = checkOpticalIsomer(p_arc);
        		if (!opt) sn *= 2;
        	}
        }

        return sn;
        //#]
    }

    /**
    Requires:
    Effects: return Cp(T)
    Modifies:
    */
    //## operation calculateCp(Temperature)
    public double calculateCp(Temperature p_temperature) {
        //#[ operation calculateCp(Temperature)
        return getThermoData().calculateCp(p_temperature);
        //#]
    }

    /**
Requires:
Effects: calculate and return the symmetry number of the cyclic portion of this ChemGraph
Modifies:
svp
*/
//## operation calculateCyclicSymmetryNumber()
public int calculateCyclicSymmetryNumber(){
  //#[ operation calculateCyclicSymmetryNumber()
 int sn = 1;
 LinkedList ring_structures = new LinkedList();//list of ring structures
 LinkedList cycle_list = new LinkedList();//list of cycles
 Iterator cycle_iter = getGraph().getCycle().iterator();
 while (cycle_iter.hasNext()){
   LinkedList current_cycle = (LinkedList)cycle_iter.next();
   cycle_list.add(current_cycle);
 }
 //if 2 graph components share at least one cycle, they belong to the same ring structure
 for (int i = 0; i <= cycle_list.size()-1; i++){
   LinkedList current_ring = (LinkedList)cycle_list.get(i);
   if (ring_structures.isEmpty()){
     ring_structures.add(current_ring);
   }
   else{
     int same_gc = 0;
     Iterator ring_structure_iter = ring_structures.iterator();
     Iterator current_ring_iter = current_ring.iterator();
     while (ring_structure_iter.hasNext()){
       LinkedList ring_structure = (LinkedList) ring_structure_iter.next();
       while (current_ring_iter.hasNext()) {
         GraphComponent gc = (GraphComponent) current_ring_iter.next();
         if (ring_structure.contains(gc)) {
           Iterator current_iter = current_ring.iterator();
           while (current_iter.hasNext()){
             GraphComponent current_gc = (GraphComponent)current_iter.next();
             if (!ring_structure.contains(current_gc)){
               ring_structure.add(current_gc);
               ring_structures.set(ring_structures.indexOf(ring_structure), ring_structure);
             }
           }
           same_gc++;
         }
       }
     }
     if (same_gc == 0){
       ring_structures.add(current_ring);
     }
   }
   if (i != cycle_list.size()-1){
     for (int j = 1; j <= cycle_list.size()-1; j++){
       LinkedList next_cycle = (LinkedList)cycle_list.get(j);
       Iterator ring_structures_iter = ring_structures.iterator();
       Iterator next_cycle_iter = next_cycle.iterator();
       while (ring_structures_iter.hasNext()){
         LinkedList current_ring_structure = (LinkedList)ring_structures_iter.next();
         while (next_cycle_iter.hasNext()){
           GraphComponent gc = (GraphComponent)next_cycle_iter.next();
           if (current_ring_structure.contains(gc)){
             Iterator ring_iter = next_cycle.iterator();
             while(ring_iter.hasNext()){
               GraphComponent current_gc = (GraphComponent)ring_iter.next();
               if (!current_ring_structure.contains(current_gc)){
                 current_ring_structure.add(current_gc);
                 ring_structures.set(ring_structures.indexOf(current_ring_structure),current_ring_structure);
               }
             }
             break;
           }
         }
       }
     }
   }
 }
 Iterator iter = ring_structures.iterator();
 while (iter.hasNext()){
   LinkedList current_ring_structure = (LinkedList)iter.next();
   Iterator gc_iter = current_ring_structure.iterator();
   LinkedList node_list = new LinkedList(); //list of all cyclic nodes
   LinkedList arc_list = new LinkedList(); //list of all cyclic arcs
   while (gc_iter.hasNext()){
     GraphComponent gc = (GraphComponent)gc_iter.next();
     gc.setVisited(false);
     if (gc instanceof Node){
       node_list.add(gc);
     }
     else {
       arc_list.add(gc);
     }
   }
   //find all sets of equal nodes
   LinkedList equal_node_list = new LinkedList();
   for (int i = 0; i <= node_list.size() - 2; i++) {
     Node current = (Node) node_list.get(i);
     if (!current.isVisited()) {
       LinkedList equal_list = new LinkedList(); //list of equivalent nodes
       current.setVisited(true);
       equal_list.add(current);
       for (int j = i + 1; j <= node_list.size() - 1; j++) {
         Node next = (Node) node_list.get(j);
         Iterator list_iter = equal_list.iterator();
         while (list_iter.hasNext()) {
           current = (Node) list_iter.next();
           if (isSymmetric(current, next)) {
             equal_list.add(next);
             next.setVisited(true);
             break;
           }
         }
       }
       equal_node_list.add(equal_list); //add list of equivalent nodes to list of sets of equivalent nodes
     }
   }
   //find all sets of equal arcs
   LinkedList equal_arc_list = new LinkedList();
   for (int i = 0; i <= arc_list.size() - 2; i++) {
     Arc current = (Arc) arc_list.get(i);
     if (!current.isVisited()) {
       LinkedList equal_list = new LinkedList(); //list of equivalent arcs
       current.setVisited(true);
       equal_list.add(current);
       for (int j = i + 1; j <= arc_list.size() - 1; j++) {
         Arc next = (Arc) arc_list.get(j);
         Iterator list_iter = equal_list.iterator();
         while (list_iter.hasNext()) {
           current = (Arc) list_iter.next();
           if (isSymmetric(current, next)) {
             equal_list.add(next);
             next.setVisited(true);
             break;
           }
         }
       }
       equal_arc_list.add(equal_list); //add list of equivalent arcs to list of sets of equivalent arcs
     }
   }
   //find largest set of equal nodes
   int node_sn = 1;
   Iterator node_list_iter = equal_node_list.iterator();
   while (node_list_iter.hasNext()) {
     LinkedList current = (LinkedList) node_list_iter.next();
     if (current.size() > node_sn) {
       node_sn = current.size(); //node symmetry number = size of largest set of equivalent nodes
     }
   }
   //find largest set of equal arcs
   int arc_sn = 1;
   Iterator arc_list_iter = equal_arc_list.iterator();
   while (arc_list_iter.hasNext()) {
     LinkedList current = (LinkedList) arc_list_iter.next();
     if (current.size() > arc_sn) {
       arc_sn = current.size(); //arc symmetry number = size of largest set of equivalent arcs
     }
   }
   if (node_sn == node_list.size() && arc_sn == arc_list.size()) { //all nodes equal and all arcs equal
     sn *= node_sn;
     sn *= 2;
     Node first_node = (Node)node_list.getFirst();
     FGElement fge = (FGElement)first_node.getFgElement();
     if (fge.equals(FGElement.make("Cs"))){
       LinkedList acyclic_neighbor = new LinkedList();
       Iterator neighbor_iter = first_node.getNeighbor();
       while (neighbor_iter.hasNext()){
         Arc arc = (Arc)neighbor_iter.next();
         if (!arc.getInCycle()){
           acyclic_neighbor.add(arc);
         }
       }
       if (acyclic_neighbor.size() == 2){
         Arc a1 = (Arc) acyclic_neighbor.getFirst();
         Arc a2 = (Arc) acyclic_neighbor.getLast();
         if (!first_node.isSymmetric(a1, a2)) {
           sn /= 2;
         }
       }
     }
   }
   else {
     if (node_sn >= arc_sn) {
       sn *= node_sn;
     }
     else {
       sn *= arc_sn;
     }
   }
   //if (sn >= 2 && sn%2 == 0){//added by Sally for non-planar PAH's
     //sn = correctSymmetryNumber(sn);
   //}
 }
 graph.setCycle(null);
 graph.formSSSR();
return sn;
//#]
}


    /**
    Requires:
    Effects: return G(T)
    Modifies:
    */
    //## operation calculateG(Temperature)
    public double calculateG(Temperature p_temperature) {
        //#[ operation calculateG(Temperature)
        return getThermoData().calculateG(p_temperature);
        //#]
    }

    /**
    Requires:
    Effects:return H(T)
    Modifies:
    */
    //## operation calculateH(Temperature)
    public double calculateH(Temperature p_temperature) {
        //#[ operation calculateH(Temperature)
        return getThermoData().calculateH(p_temperature);
        //#]
    }

    //## operation calculateInternalRotor()
    public void calculateInternalRotor() {
        //#[ operation calculateInternalRotor()
        // add more check for resonance axis!!
        int rotor = 0;
        Graph g = getGraph();
        for (Iterator iter = g.getArcList(); iter.hasNext();) {
        	Arc a = (Arc)iter.next();
        	Bond bond = (Bond)a.getElement();
        	if (bond.isSingle() && !a.getInCycle()) {
        		Iterator atomIter = a.getNeighbor();
        		Node n1 = (Node)atomIter.next();
        		Node n2 = (Node)atomIter.next();
        		if (!n1.isLeaf() && !n2.isLeaf()) {
        			rotor++;
        		}
        	}
        }
        internalRotor = rotor;

        //#]
    }

    //uses same algorithm as calculate internal rotor, but stores the two atoms involved in the rotor and all the atoms on one side of the rotor (specifically, the side corresponding to the 2nd atom of the rotor)
    //the information is returned in a LinkedHashMap where the Key is an array of [atom0 atom1 atom2 atom3] (with atom1 and atom2 being the "rotor atoms" and the others being the dihedral atoms) and the value is a Collection of atoms IDs associated with atom2
    //note that at this stage, there is no check to make sure that the dihedral atoms are not collinear
    public LinkedHashMap getInternalRotorInformation(){
	LinkedHashMap rotorInfo = new LinkedHashMap();
	Graph g = getGraph();
        for (Iterator iter = g.getArcList(); iter.hasNext();) {
        	Arc a = (Arc)iter.next();
        	Bond bond = (Bond)a.getElement();
        	if (bond.isSingle() && !a.getInCycle()) {
        		Iterator atomIter = a.getNeighbor();
        		Node n1 = (Node)atomIter.next();
        		Node n2 = (Node)atomIter.next();
        		if (!n1.isLeaf() && !n2.isLeaf()) {
        			//rotor++;
				//above here is the rotor identification algorithm; below here is the code that stores the necessary information about the rotor
				Graph f=Graph.copy(g);//copy the graph so we don't modify the original
				f.removeArc(f.getArcBetween(n1.getID(), n2.getID()));//this should separate the graph into disconnected pieces (unless it is part of a cycle; if it is part of a cycle, however, this section of code shouldn't be reached)
				LinkedList pieces = f.partitionWithPreservedIDs();//partition into the two separate graphs
				Graph sideA = (Graph)pieces.getFirst();
				Graph sideB = (Graph)pieces.getLast();
				//look for the piece that has node2
				if(sideA.getNodeIDs().contains(n2.getID())){
				    Node atom1 = sideB.getNodeAt(n1.getID());
				    Node atom2 = sideA.getNodeAt(n2.getID());
				    Node dihedral1 = (Node)atom1.getNeighboringNodes().iterator().next();//get a neighboring node
				    Node dihedral2 = (Node)atom2.getNeighboringNodes().iterator().next();//get a neighboring node
				    int rotorSym = calculateRotorSymmetryNumber(atom1,atom2);
				    int[] rotorAtoms = {dihedral1.getID(), n1.getID(), n2.getID(), dihedral2.getID(), rotorSym};
				    rotorInfo.put(rotorAtoms, sideA.getNodeIDs());
				}
				else if (sideB.getNodeIDs().contains(n2.getID())){
				    Node atom1 = sideA.getNodeAt(n1.getID());
				    Node atom2 = sideB.getNodeAt(n2.getID());
				    Node dihedral1 = (Node)atom1.getNeighboringNodes().iterator().next();//get a neighboring node
				    Node dihedral2 = (Node)atom2.getNeighboringNodes().iterator().next();//get a neighboring node
				    int rotorSym = calculateRotorSymmetryNumber(atom1,atom2);
				    int[] rotorAtoms = {dihedral1.getID(), n1.getID(), n2.getID(), dihedral2.getID(), rotorSym};
				    rotorInfo.put(rotorAtoms, sideB.getNodeIDs());
				}
				else{
				    System.out.println("Error in getInternalRotorInformation(): Cannot find node "+ n2.getID()+" after splitting from "+ n1.getID() +" in the following graph:\n"+ g.toString());
				    System.exit(0);
				}
        		}
        	}
        }

	return rotorInfo;
    }

//	## operation isLinear() 
    public boolean isLinear() {
        //#[ operation isLinear() 
        // only check for linearity in molecules with at least two atoms
        if (getAtomNumber() == 1) return false;
        
        // cyclic molecules are not linear
        if (!isAcyclic()) return false;
        
        // biatomic molecules are always linear
        if (getAtomNumber() == 2) return true;
        
        // molecules with only double bonds are linear (e.g. CO2)
        boolean allDouble = true;
        Iterator iter = getArcList();
        while (iter.hasNext()) {
        	Arc arc = (Arc)iter.next();
        	Bond bond = (Bond)arc.getElement();
        	if (!bond.isDouble()) allDouble = false;
        }
        if (allDouble) return true;
        
        // molecule with alternating single and triple bonds are linear (e.g. acetylene)
        boolean alternatingSingleTriple = true;
        Iterator node_iter = getNodeList();
        while (node_iter.hasNext()) {
        	Node node = (Node)node_iter.next();
        	int neighborNumber = node.getNeighborNumber();
        	if (neighborNumber == 2) {
        		Iterator neighbor_iter = node.getNeighbor();
        		Arc a1 = (Arc)neighbor_iter.next();
        		Bond b1 = (Bond)a1.getElement();
        		Arc a2 = (Arc)neighbor_iter.next();
        		Bond b2 = (Bond)a2.getElement();
        		if (! ((b1.isTriple() && b2.isSingle()) || (b1.isSingle() && b2.isTriple())))
        			alternatingSingleTriple = false;
        	}
        	else if (neighborNumber > 2)
        		alternatingSingleTriple = false;
        }
            
            if (alternatingSingleTriple) return true;
        
        // if none of the above are true, it's nonlinear
        return false;
        
        //#]
    }
	
    /**
    Requires:
    Effects: return S(T)
    Modifies:
    */
    //## operation calculateS(Temperature)
    public double calculateS(Temperature p_temperature) {
        //#[ operation calculateS(Temperature)
        return getThermoData().calculateS(p_temperature);




        //#]
    }

    //## operation calculateSymmetryNumber()
   public int calculateSymmetryNumber() {
       //#[ operation calculateSymmetryNumber()
       try {
         getGraph().formSSSR();//svp
               int sn = 1;
               Iterator iter = getNodeList();
               while (iter.hasNext()) {
                       Node node = (Node)iter.next();
                       if (!node.getInCycle()){//svp
                        sn *= calculateAtomSymmetryNumber(node);
                        }

               }
               iter = getArcList();
               while (iter.hasNext()) {
                       Arc arc = (Arc)iter.next();
                       if (!arc.getInCycle()){//svp
                       sn *= calculateBondSymmetryNumber(arc);
                     }

               }
               sn *= calculateAxisSymmetryNumber();
               if (!isAcyclic()) {//svp
                sn *= calculateCyclicSymmetryNumber();
              }

               symmetryNumber = sn;
               return sn;
       }
       catch (ClassCastException e) {
               throw new InvalidChemGraphException();
       }

       //#]
   }

    /**
    Requies:
    Effects: check if p_arc bond is the single bond between two Oxygens, which is considered as optical isomers
    Modifies:
    */
    //## operation checkOpticalIsomer(Arc)
    public boolean checkOpticalIsomer(Arc p_arc) {
        //#[ operation checkOpticalIsomer(Arc)
        // check if the p_arc is -O-O-
        Bond b = (Bond)p_arc.getElement();
        if (b.isSingle()) {
        	Iterator neighbor_iter = p_arc.getNeighbor();
        	Node n1 = (Node)neighbor_iter.next();
        	Node n2 = (Node)neighbor_iter.next();
        	Atom a1 = (Atom)n1.getElement();
        	Atom a2 = (Atom)n2.getElement();
        	FGElement fge1 = (FGElement)n1.getFgElement();
        	FGElement fge2 = (FGElement)n2.getFgElement();
        	if (!a1.isRadical() && fge1.equals(FGElement.make("Os")) && !a2.isRadical() && fge2.equals(FGElement.make("Os"))) {
        		return true;
        	} 
//        	else if (!a1.isRadical() && fge1.equals(FGElement.make("Sis")) && !a2.isRadical() && fge2.equals(FGElement.make("Sis"))) {
//        		return true;
//        	}
        }

        return false;


        //#]
    }

    /**
    Requires:
    Effects: clear the central node list
    Modifies:
    */
    //## operation clearCentralNode()
    public void clearCentralNode() {
        //#[ operation clearCentralNode()
        getGraph().clearCentralNode();
        //#]
    }

    /**
    Requires:
    Effects: return a new instance identical to this ChemGraph
    Modifies:
    */
    //## operation copy(ChemGraph)
    public static ChemGraph copy(ChemGraph p_chemGraph) throws ForbiddenStructureException {
        Graph g = Graph.copy(p_chemGraph.getGraph());

        ChemGraph cg = new ChemGraph(g);
        cg.uniqueString = p_chemGraph.getUniqueString();
        cg.chemicalFormula = p_chemGraph.getChemicalFormula();
        cg.species = p_chemGraph.getSpecies();
        cg.symmetryNumber = p_chemGraph.symmetryNumber;
        cg.thermoData = p_chemGraph.thermoData;
        cg.thermoGAPP = p_chemGraph.thermoGAPP;
        cg.InChI = p_chemGraph.InChI;
        cg.internalRotor = p_chemGraph.internalRotor;
        cg.solvthermoData = p_chemGraph.solvthermoData;

        /*HashSet oldSymmetryAxis = p_chemGraph.getSymmetryAxis();
        if (oldSymmetryAxis != null) {
        	cg.symmetryAxis = new HashSet();
        	for (Iterator iAxis = oldSymmetryAxis.iterator(); iAxis.hasNext(); ) {
        		HashSet newAxis = new HashSet();
        		HashSet oldAxis = (HashSet)iAxis.next();
        		for (Iterator iArc = oldAxis.iterator(); iArc.hasNext(); ) {
        			Arc arc = (Arc)iArc.next();
        			Iterator iNode = arc.getNeighbor();
         			int n1 = ((Node)iNode.next()).getID().intValue();
           			int n2 = ((Node)iNode.next()).getID().intValue();
           			Arc newArc = cg.getArcBetween(n1,n2);
           			newAxis.add(newArc);
           		}
        		cg.symmetryAxis.add(newAxis);
           	}
        }*/

        return cg;
    }

    /**
    Requires:
    Effects: return true iff two chemgraph have equivalent graph structures
    Modifies:
    */
    //## operation equals(Object)
    public boolean equals(Object p_chemGraph) {
        //#[ operation equals(Object)
        if (this == p_chemGraph) return true;

        return isEquivalent((ChemGraph)p_chemGraph);
        //#]
    }

    /**
    Requires:
    Effects: generate the chemical formula of this chem graph and return it.  if the graph is not initialized, return null.
    Modifies:
    */
    //## operation generateChemicalFormula()
    public String generateChemicalFormula() {
        //#[ operation generateChemicalFormula()
        if (getGraph() == null) return null;

        /*int cap = ChemElementDictionary.size();
        Vector type = new Vector(cap);
        int[] number = new int[cap];
        */
        int C_number = 0;
        int H_number = 0;
        int O_number = 0;
        int radical = 0;
        // Added by MRH on 18-Jun-2009
        //	Hardcoding Si and S into RMG-java
        int Si_number = 0;
        int S_number = 0;

        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();
        	radical += atom.getRadicalNumber();

        	if (atom.isCarbon()) {
        		C_number++;
        	}
        	else if (atom.isHydrogen()) {
        		H_number++;
        	}
        	else if (atom.isOxygen()) {
        		O_number++;
        	}
            // Added by MRH on 18-Jun-2009
            //	Hardcoding Si and S into RMG-java
        	else if (atom.isSilicon()) {
        		Si_number++;
        	}
        	else if (atom.isSulfur()) {
        		S_number++;
        	}
        	else {
        		throw new InvalidChemNodeElementException();
        	}
        }

        String s = "";
        if (C_number>0) {
        	s = s + "C";
        	if (C_number >1) {
        		s = s + String.valueOf(C_number);
        	}
        }
        if (H_number>0) {
        	s = s + "H";
        	if (H_number >1) {
        		s = s + String.valueOf(H_number);
        	}
        }
        if (O_number>0) {
        	s = s + "O";
        	if (O_number >1) {
        		s = s + String.valueOf(O_number);
        	}
        }
        // Added by MRH on 18-Jun-2009
        //	Hardcoding Si and S into RMG-java
        if (Si_number>0) {
        	s += "Si";
        	if (Si_number>1) {
        		s += String.valueOf(Si_number);
        	}
        }
        if (S_number>0) {
        	s += "S";
        	if (S_number>1) {
        		s += String.valueOf(S_number);
        	}
        }

        chemicalFormula = s;
        if (radical == 1) {
        	chemicalFormula = chemicalFormula + "J";
        }
        else if (radical == 2) {
        	chemicalFormula = chemicalFormula + "JJ";
        }
        else if (radical == 3) {
        	chemicalFormula = chemicalFormula + "JJJ";
        }

        return chemicalFormula;


        //#]
    }
    
    //6/9/09 gmagoon: modified to also generate InChIKey; see comments in corresponding Species function
    public String [] generateInChI() {
    	String [] result = Species.generateInChI(this);
        InChI = result[0];
        InChIKey = result[1];
    	return result;
    }
    
    public String generateMolFileString() {
    	String mfs = Species.generateMolFileString(this);
    	return mfs;
    }
    
    /**
    Requires:
    Effects: if the thermoGAPP is not set, set default GAPP.  Use it to calculate the thermoData of this chem graph.  if there is any exception during this process, throw FailGenerateThermoDataException.
    Modifies: this.thermoData
    */
    //## operation generateThermoData()
    public ThermoData generateThermoData() throws FailGenerateThermoDataException {
        //#[ operation generateThermoData()
        // use GAPP to generate Thermo data
        try {
                if (useQM){
                    if(useQMonCyclicsOnly && this.isAcyclic()) thermoGAPP=GATP.getINSTANCE();//use GroupAdditivity for acyclic compounds if this option is set
                    else  thermoGAPP=QMTP.getINSTANCE();
                }
                else if (thermoGAPP == null) setDefaultThermoGAPP();
        	thermoData = thermoGAPP.generateThermoData(this);

            //thermoData = thermoGAPP.generateAbramData(this);
        	return thermoData;
        }
		catch (MultipleGroupFoundException e) {
			throw e;
		}
        catch (Exception e) {
		e.printStackTrace();
        	throw new FailGenerateThermoDataException();
        }
        //#]
    }
    
    public TransportData generateTransportData() {
    	if (transportGAPP == null) setDefaultTransportGAPP();
    	transportData = transportGAPP.generateTransportData(this);
    	return transportData;
    }
	
    // Amrit Jalan 05/09/2009
	public ThermoData generateSolvThermoData() throws FailGenerateThermoDataException {
        // use GAPP to generate Thermo data
        try {
        	if (SolvationGAPP == null) setDefaultSolvationGAPP();
            solvthermoData = SolvationGAPP.generateSolvThermoData(this);
        	return solvthermoData;
        }
        catch (Exception e) {
        	throw new FailGenerateThermoDataException();
        }
    }
	
    public AbramData generateAbramData() throws FailGenerateThermoDataException {
        // use GAPP to generate Thermo data
        try {
        	if (abramGAPP == null) setDefaultAbramGAPP();
        	abramData = abramGAPP.generateAbramData(this);
        	return abramData;
        }
        catch (Exception e) {
        	throw new FailGenerateThermoDataException();
        }
    }
	
	public UnifacData generateUnifacData() throws FailGenerateThermoDataException {
        try {
        	if (unifacGAPP == null) setDefaultUnifacGAPP();
        	unifacData = unifacGAPP.generateUnifacData(this);
        	return unifacData;
        }
        catch (Exception e) {
        	throw new FailGenerateThermoDataException();
        }
    }
	
    /**
    Requires:
    Effects: return the Arc between two positions in this ChemGraph
    Modifies:
    */
    public Arc getArcBetween(int p_position1, int p_position2) {
        return getGraph().getArcBetween(p_position1,p_position2);
    }

    /**
    Requires:
    Effects: return an arc iterator of this ChemGraph
    Modifies:
    */
    //## operation getArcList()
    public Iterator getArcList() {
        //#[ operation getArcList()
        return getGraph().getArcList();
        //#]
    }

    /**
    Requires:
    Effects: return the atom at the p_position in this chem graph; if p_position is empty, return null;
    Modifies:
    */
    //## operation getAtomAt(int)
    public Atom getAtomAt(int p_position) throws EmptyAtomException {
        //#[ operation getAtomAt(int)
        try {
        	return (Atom)(getNodeAt(p_position).getElement());
        }
        catch (NotInGraphException e) {
        	return null;
        }



        //#]
    }

    /**
    Requires:
    Effects: return the total atom number in this ChemGraph.
    Modifies:
    */
    //## operation getAtomNumber()
    public int getAtomNumber() {
        //#[ operation getAtomNumber()
        return getGraph().getNodeNumber();
        //#]
    }

    /**
    Requires:
    Effects: if there is a bond connecting p_position1 and p_position2, return that bond; otherwise, return null.
    Modifies:
    */
    //## operation getBondBetween(int,int)
    public Bond getBondBetween(int p_position1, int p_position2) throws EmptyAtomException {
        //#[ operation getBondBetween(int,int)
        try {
        	return (Bond)(getArcBetween(p_position1,p_position2).getElement());
        }
        catch (ClassCastException e) {
        	return null;
        }
        //#]
    }

    //## operation getCarbonNumber()
    public int getCarbonNumber() {
        //#[ operation getCarbonNumber()
        int cNum = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();

        	if (atom.isCarbon()) {
        		cNum++;
        	}
        }
        return cNum;
        //#]
    }

    /**
    Requires:
    Effects: return the hashMap of centralNode in this ChemGraph
    Modifies:
    */
    //## operation getCentralNode()
    public LinkedHashMap getCentralNode() {
        //#[ operation getCentralNode()
        return getGraph().getCentralNode();
        //#]
    }

    /**
    Requires:
    Effects: return the node whose centralID equals p_position
    Modifies:
    */
    //## operation getCentralNodeAt(int)
    public Node getCentralNodeAt(int p_position) {
        //#[ operation getCentralNodeAt(int)
        return getGraph().getCentralNodeAt(p_position);


        //#]
    }

    /**
    Requires:
    Effects: return the number of the central nodes in this chem graph
    Modifies:
    */
    //## operation getCentralNodeNumber()
    public int getCentralNodeNumber() {
        //#[ operation getCentralNodeNumber()
        return getGraph().getCentralNodeNumber();
        //#]
    }

    //## operation getChemicalFormula()
    public String getChemicalFormula() {
        //#[ operation getChemicalFormula()
        if (chemicalFormula == null || chemicalFormula.length() == 0) generateChemicalFormula();
        return chemicalFormula;
        //#]
    }

    /**
    Requires:
    Effects: return the iterator loop over the graph's cycle list.
    Modifies:
    */
    //## operation getCycle()
    public Iterator getCycle() {
        //#[ operation getCycle()
        return getGraph().getCycle().iterator();
        //#]
    }

	public int getCycleNumber(){
		return getGraph().getCycleNumber();
	}
    //## operation getHydrogenNumber()
    public int getHydrogenNumber() {
        //#[ operation getHydrogenNumber()
        int hNum = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();

        	if (atom.isHydrogen()) {
        		hNum++;
        	}
        }
        return hNum;
        //#]
    }
    
    public String getInChI() {
    	if (InChI == null || InChI.length() == 0) generateInChI();
    	return InChI;
    }
    
    public String getInChIKey() {
    	if (InChIKey == null || InChIKey.length() == 0) generateInChI();
    	return InChIKey;
    }
    
    //gmagoon 7/25/09: same as above function, except it will replace an existing InChI; this is useful when ChemGraph changes (e.g. during HBI process)
    public String getInChIAnew() {
    	generateInChI();
    	return InChI;
    }
    
    //gmagoon 7/25/09: same as above function, except it will replace an existing InChIKey; this is useful when ChemGraph changes (e.g. during HBI process)
    public String getInChIKeyAnew() {
    	generateInChI();
    	return InChIKey;
    }

    //gmagoon 9/1/09: these functions (getModifiedInChIAnew and getModifiedInChIKeyAnew) do the same as above, except return a modified version of the InChI and InChIKey with an indication of the multiplicity based on the number of radicals if the number of radicals is 2 or greater
    public String getModifiedInChIAnew(){
        generateInChI();
        String newInChI=null;
        int radicalNumber = this.getUnpairedRadicalNumber();
       // System.out.println("Radical number:"+radicalNumber);//for debugging purposes
        if (radicalNumber == 2){
            newInChI = InChI.concat("/mult3");
        }
        else if (radicalNumber == 3){
            newInChI = InChI.concat("/mult4");
        }
        else if (radicalNumber == 4){
            newInChI = InChI.concat("/mult5");
        }
        else if (radicalNumber == 5){
            newInChI = InChI.concat("/mult6");
        }
        else if (radicalNumber > 5){
            System.out.println("Unexpectedly high multiplicity ("+ (radicalNumber+1)+ ") for " + InChI);
            System.exit(0);
        }
        else{
            newInChI = InChI;
        }
        
        return newInChI;
    }
    
    public String getModifiedInChIKeyAnew(){
        generateInChI();
        String newInChIKey=null;
        int radicalNumber = this.getUnpairedRadicalNumber();
       // System.out.println("Radical number:"+radicalNumber);//for debugging purposes
        if (radicalNumber == 2){
            newInChIKey = InChIKey.concat("mult3");
        }
        else if (radicalNumber == 3){
            newInChIKey = InChIKey.concat("mult4");
        }
        else if (radicalNumber == 4){
            newInChIKey = InChIKey.concat("mult5");
        }
        else if (radicalNumber == 5){
            newInChIKey = InChIKey.concat("mult6");
        }
        else if (radicalNumber > 5){
            System.out.println("Unexpectedly high multiplicity ("+ (radicalNumber+1)+ ") for " + InChI);
            System.exit(0);
        }
        else{
            newInChIKey= InChIKey;
        }
        
        return newInChIKey;
    }
    
    /**
    Requires:
    Effects: add the weight of all the atoms in this graph to calculate the molecular weight of this chem graph
    Modifies:
    */
    //## operation getMolecularWeight()
    public double getMolecularWeight() {
        //#[ operation getMolecularWeight()
        double MW = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	MW += ((Atom)(node.getElement())).getWeight();
        }
        return MW;
        //#]
    }

    /**
    Requires:
    Effects: return the "legal" name of this ChemGraph.  The order of priority
    (1) uniqueString
    (2) species.name
    (3) chemicalFormula
    Modifies:
    */
    //## operation getName()
    public String getName() {
        //#[ operation getName()
        if (uniqueString != null && uniqueString.length() > 0) return uniqueString;

        String name = getSpecies().getName();
        if (name != null && name.length() > 0) return name;

        return getChemicalFormula();
        //#]
    }

    /**
    Requires:
    Effects: return the node whose ID equals p_position.
    Modifies:
    */
    //## operation getNodeAt(int)
    public Node getNodeAt(int p_position) {
        //#[ operation getNodeAt(int)
        return getGraph().getNodeAt(p_position);
        //#]
    }

    /**
    Requires:
    Effects: return the node whose ID equals p_ID.
    Modifies:
    */
    //## operation getNodeAt(Integer)
    public Node getNodeAt(Integer p_ID) {
        //#[ operation getNodeAt(Integer)
        return getGraph().getNodeAt(p_ID);
        //#]
    }

    /**
    Requires:
    Effects: return an iterator over the node collection of this graph
    Modifies:
    */
    //## operation getNodeList()
    public Iterator getNodeList() {
        //#[ operation getNodeList()
        return getGraph().getNodeList();
        //#]
    }

	public int getDoubleBonds(){
		int dBond = 0;
		Iterator iter = getArcList();
		while (iter.hasNext()){
			Arc arc = (Arc)iter.next();
			Bond bond = (Bond)arc.getElement();
			if (bond.order==2){
				dBond++;
			}
		}
		return dBond;
	}
	
	public int getSingleBonds(){
		int sBond = 0;
		Iterator iter = getArcList();
		while (iter.hasNext()){
			Arc arc = (Arc)iter.next();
			Bond bond = (Bond)arc.getElement();
			if (bond.order==1){
				sBond++;
			}
		}
		return sBond;
	}
	
	public int getTripleBonds(){
		int tBond = 0;
		Iterator iter = getArcList();
		while (iter.hasNext()){
			Arc arc = (Arc)iter.next();
			Bond bond = (Bond)arc.getElement();
			if (bond.order==3){
				tBond++;
			}
		}
		return tBond;
	}
	
    //## operation getOxygenNumber()
    public int getOxygenNumber() {
        //#[ operation getOxygenNumber()
        int oNum = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();

        	if (atom.isOxygen()) {
        		oNum++;
        	}
        }
        return oNum;
        //#]
    }


    /**
    Requires:
    Effects: return a collection of all the radical sites
    Modifies:
    */
    //## operation getRadicalNode()
    public LinkedHashSet getRadicalNode() {
        //#[ operation getRadicalNode()
        LinkedHashSet radicalNode = new LinkedHashSet();
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node n = (Node)iter.next();
        	Atom a = (Atom)n.getElement();
        	if (a.isRadical()) {
        		radicalNode.add(n);
        	}
        }
        return radicalNode;
        //#]
    }

    /**
    Requires:
    Effects: calculate the total radical number in this chem graph.
    Modifies:
    */
    //## operation getRadicalNumber()
    public int getRadicalNumber() {
        //#[ operation getRadicalNumber()
        int radicalNumber = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Object element = ((Node)(iter.next())).getElement();
        	radicalNumber += ((Atom)element).getRadicalNumber();
        }
        return radicalNumber;
        //#]
    }

    //gmagoon 4/30/10: modified version of getRadicalNumber; same as getRadicalNumber, except it will not count 2S as radical
    public int getUnpairedRadicalNumber() {
        int radicalNumber = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Object element = ((Node)(iter.next())).getElement();
        	radicalNumber += ((Atom)element).getUnpairedRadicalNumber();
        }
        return radicalNumber;
    }
    
    public int getSiliconNumber() {
        int siNum = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();
        	if (atom.isSilicon()) {
        		siNum++;
        	}
        }
        return siNum;
    }
    
    public int getSulfurNumber() {
        int sNum = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();
        	if (atom.isSulfur()) {
        		sNum++;
        	}
        }
        return sNum;
    }

    //## operation getSymmetryNumber()
    public int getSymmetryNumber() {
        //#[ operation getSymmetryNumber()
        if (symmetryNumber < 0) calculateSymmetryNumber();
        return symmetryNumber;
        //#]
    }

    //## operation getThermoData()
    public ThermoData getThermoData() {
        //#[ operation getThermoData()
        if (thermoData == null)
        { 
            generateThermoData();
            // thermoData is brand new gas phase estimate
            if (Species.useSolvation) {
                if (solvthermoData==null) generateSolvThermoData();
                // solvthermoData is estimate of correction
                thermoData.plus(solvthermoData);
                thermoData.comments+=" corrected for solvation";
                // thermoData is corrected
            }
        }
        return thermoData;
        //#]
    }
    
    public TransportData getTransportData() {
        if (transportData == null) generateTransportData();
        return transportData;
    }

        //## operation getThermoData()
    public ThermoData getSolvationData() {
        //#[ operation getThermoData()
        if (solvthermoData == null) generateSolvThermoData();
        return solvthermoData;
        //#]
    }
    
    public AbramData getAbramData() {
        //#[ operation getThermoData()
        if (abramData == null) generateAbramData();
        return abramData;
        //#]
    }

        public UnifacData getUnifacData() {
        //#[ operation getThermoData()
        if (unifacData == null) generateUnifacData();
        return unifacData;
        //#]
    }

    /**
	 Added by: Amrit Jalan
	 Effects: calculate the raduis of the chemGraph using UNIFAC Ri values. (UNITSof radius = m)
	 */
    public double getRadius() {
		
        double ri;
		
        if (getCarbonNumber() == 0 && getOxygenNumber() == 0){    // Which means we ar dealing with HJ or H2
			double ri3;
			ri3 = 21*8.867/88;                              // 8.867 Ang^3 is the volume of a single Hydrogen Atom
            if (getHydrogenNumber() == 1){                        // i.e. we are dealing with the Hydrogen radical
                ri = Math.pow(ri3,0.333) * Math.pow(10,-10);
                return ri;
            }
            if (getHydrogenNumber() == 2){                        // i.e. we are dealing with the Hydrogen molecule
                ri3 = 2*ri3;                                      // Assumption: volume of H2 molecule ~ 2 * Volume of H atom
                ri = Math.pow(ri3,0.333) * Math.pow(10,-10);
                return ri;                
            }
        }
		
        double Ri=getUnifacData().R;
        ri=3.18*Math.pow(Ri,0.333)*Math.pow(10,-10);   // From Koojiman Ind. Eng. Chem. Res 2002, 41 3326-3328
        return ri;
		
    }
	
	/**
	 Added by: Amrit Jalan
	 Effects: calculate the diffusivity of the chemGraph using radii, solvent viscosity and Stokes Einstein. (UNITS m2/sec)
	 */
    public double getDiffusivity() {
		double speRad=getRadius();
		double solventViscosity = 0.473*Math.pow(10,-3);
		double diffusivity;
		Temperature sysTemp = ReactionModelGenerator.getTemp4BestKinetics();
		double denom = 132*solventViscosity*speRad/7;
		diffusivity = 1.381*500* Math.pow(10,-23)/denom;  //sysTemp.getK()
		return diffusivity;
    }
    
    /**
    Requires:
    Effects: find out the end of C=C=C... pattern
    Modifies:
    */
    //## operation getToEndOfAxis(Arc,Node,HashSet)
    private static final Node getToEndOfAxis(Arc p_beginArc, Node p_beginNode, HashSet p_axis) {
        //#[ operation getToEndOfAxis(Arc,Node,HashSet)
        Arc nextArc = null;
        Iterator iter = p_beginNode.getNeighbor();
        while (iter.hasNext()) {
        	nextArc = (Arc)iter.next();
        	if (nextArc != p_beginArc) break;
        }

        p_axis.add(nextArc);

        Node nextNode = nextArc.getOtherNode(p_beginNode);
        FGElement fge = (FGElement)nextNode.getFgElement();
        FGElement Cdd = FGElement.make("Cdd");
        FGElement Sidd = FGElement.make("Sidd");

        if (!fge.equals(Cdd) & !fge.equals(Sidd)) return nextNode;
        else {
        	return getToEndOfAxis(nextArc,nextNode,p_axis);
        }
        //#]
    }

    //find the end of the Ct-Ct-Ct... pattern
    //based on similar function getToEndOfAxis
    private static final Node getToEndOfCumulatedTripleBondSystem(Arc p_beginArc, Node p_beginNode, HashSet p_axis) {
        //#[ operation getToEndOfAxis(Arc,Node,HashSet)
        Arc nextArc = null;
        Iterator iter = p_beginNode.getNeighbor();
        while (iter.hasNext()) {
        	nextArc = (Arc)iter.next();
        	if (nextArc != p_beginArc) break;
        }

        p_axis.add(nextArc);

        Node nextNode = nextArc.getOtherNode(p_beginNode);
        FGElement fge = (FGElement)nextNode.getFgElement();
        FGElement Ct = FGElement.make("Ct");

        if (!fge.equals(Ct)) return nextNode;
        else {
        	return getToEndOfCumulatedTripleBondSystem(nextArc,nextNode,p_axis);
        }
    }

    /**
    Requires:
    Effects: check if the element of every node is an atom, and if the element of every node is a bond.
    Modifies:
    */
    //## operation graphContentsOk(Graph)
    public static boolean graphContentsOk(Graph p_graph) {
        //#[ operation graphContentsOk(Graph)
        Iterator iter = p_graph.getNodeList();
        while (iter.hasNext()) {
        	Object atom = ((Node)iter.next()).getElement();
        	if (!(atom instanceof Atom)) return false;
        }
        iter = p_graph.getArcList();
        while (iter.hasNext()) {
        	Object bond = ((Arc)iter.next()).getElement();
        	if (!(bond instanceof Bond)) return false;
        }
        return true;
        //#]
    }

   /**
    Requires:
    Effects: return chemicalFormula's hashcode.  i.e., all the isomers have the same hashcode
    Modifies:
    */
    //## operation hashCode()
    public int hashCode() {
        //#[ operation hashCode()
        if (chemicalFormula == null) generateChemicalFormula();
        return chemicalFormula.hashCode()+getTripleBonds()*300 + getSingleBonds()*1 + getDoubleBonds()*20;
        //#]
    }

    /**
    Requires:
    Effects: check all the possible reacted sites in this chemgraph according to the pass-in functional group or functional group collection.  If there are any matches, return all the matches in a linked list; otheriwse, return an empty list.
    Modifies:
    */
    //## operation identifyReactionMatchedSite(Matchable)
    public LinkedHashSet identifyReactionMatchedSite(Matchable p_functionalGroup) {
        //#[ operation identifyReactionMatchedSite(Matchable)
        if (p_functionalGroup instanceof FunctionalGroup) {
        	FunctionalGroup fg = (FunctionalGroup)p_functionalGroup;
        	//boolean thisIsRadical = this.isRadical();
        	//boolean fgIsRadical = fg.isRadical();
        	//if (thisIsRadical == fgIsRadical) {
        		return getGraph().identifyAllOrderedMatchedSites(fg.getGraph());
        	//}
        	//else {
        		//return new LinkedHashSet();
        	//}
        }
        else if (p_functionalGroup instanceof FunctionalGroupCollection) {
        	LinkedHashSet result = new LinkedHashSet();
        	FunctionalGroupCollection fgc = (FunctionalGroupCollection)p_functionalGroup;
        	Iterator iter = fgc.getFunctionalGroups();
        	while (iter.hasNext()) {
        		FunctionalGroup fg = (FunctionalGroup)iter.next();
        		//boolean thisIsRadical = this.isRadical();
        		//boolean fgIsRadical = fg.isRadical();
        		//if (thisIsRadical == fgIsRadical) {
           			LinkedHashSet site = getGraph().identifyAllOrderedMatchedSites(fg.getGraph());
        			result.addAll(site);
        		//}
        	}
        	return result;
        }
        else {
        	throw new InvalidFunctionalGroupException();
        }

        //#]
    }

    /**
    Requires:
    Effects: get a collection of all the thermo matched site, this for corrections in GAPP.
    Modifies:
    */
    //## operation identifyThermoMatchedSite(FunctionalGroup)
    public LinkedHashSet identifyThermoMatchedSite(FunctionalGroup p_functionalGroup) {
        //#[ operation identifyThermoMatchedSite(FunctionalGroup)
        return getGraph().identifyAllUnorderedMatchedSite(p_functionalGroup.getGraph());
        //#]
    }

    /**
    Requires:
    Effects: return true if there is no cycle in the chem graph
    Modifies:
    */
    //## operation isAcyclic()
    public boolean isAcyclic() {
        //#[ operation isAcyclic()
        return getGraph().isAcyclic();
        //#]
    }

    /**
    Requires:
    Effects: return true iff this and p_chemGraph are equivalent chemgraphs.
    Modifies:
    */
    //## operation isEquivalent(ChemGraph)
    public boolean isEquivalent(ChemGraph p_chemGraph) {
        //#[ operation isEquivalent(ChemGraph)
        if (chemicalFormula == null) generateChemicalFormula();
        if (p_chemGraph.chemicalFormula == null) p_chemGraph.generateChemicalFormula();

        if (!getChemicalFormula().equals(p_chemGraph.getChemicalFormula())) return false;

        if (!getGraph().isEquivalent(p_chemGraph.getGraph())) return false;

        return true;
        //#]
    }

    /**
    Requires:
    Effects: return true iff this chemGraph contains forbidden structure.
    Modifies:
    */
	
    public static boolean isForbiddenStructure(Graph p_graph, int radNumber, int oNumber, int cNumber) {

		/*Iterator iter = p_graph.getNodeList();
		while (iter.hasNext()){
			Node n = (Node)iter.next();
			Atom atom = (Atom)n.getElement();
			if (atom.isOxygen()) {
				Iterator neighborarc = n.getNeighbor();
				if (neighborarc.hasNext()){
					Arc arc1 = (Arc)neighborarc.next();
					if (neighborarc.hasNext()){
						Node neighbornode1 = n.getOtherNode(arc1);
						Arc arc2 = (Arc)neighborarc.next();
						Node neighbornode2 = n.getOtherNode(arc2);
						if (((Atom)neighbornode1.getElement()).isOxygen() && ((Atom)neighbornode2.getElement()).isOxygen())
							return true;
					}
				}
			}
		}
		return false;*/
        
        for (Iterator iter = forbiddenStructure.iterator(); iter.hasNext(); ) {
        	FunctionalGroup fg = (FunctionalGroup)iter.next();
        	if (radNumber >= fg.rad_count && oNumber >= fg.O_count && cNumber >= fg.C_count) {
	        	Graph g = fg.getGraph();
	        	if (p_graph.isSub(g)) {
	        		return true;
	        	}
        	}
        }
        return false;
    }
	
	
	// Which forbidden structure forbade this chemgraph?
	// returns the names of forbidden structures.
    public static String whichForbiddenStructures(Graph p_graph) {
		String forbidden_by = "";
		for (Iterator iter = forbiddenStructure.iterator(); iter.hasNext(); ) {
        	FunctionalGroup fg = (FunctionalGroup)iter.next();
	        Graph g = fg.getGraph();
	        if (p_graph.isSub(g)) {
	        	forbidden_by += fg.getName() +", ";
	        }
        }
        if (forbidden_by=="") return "no forbidden structures, ";
		return forbidden_by;
    }

    /**
    Requires:
    Effects: return true iff this chemgraph contains radical site
    Modifies:
    */
    //## operation isRadical()
    public boolean isRadical() {
        //#[ operation isRadical()
        return (getRadicalNumber() > 0);
        //#]
    }

    /**
    Requires:
    Effects: if p_functionalGroup is a FunctionalGroup, return if this chemgraph is matched with it at the central nodes; if p_functionalGroup is a FunctionalGroupCollection, return if this chemgraph is matched with any of the functionalgroup in the collection at the central node.  for all other case, return false.
    Modifies:
    */
    //## operation isSubAtCentralNodes(Matchable)
    public boolean isSubAtCentralNodes(Matchable p_functional) {
        //#[ operation isSubAtCentralNodes(Matchable)
        if (this == p_functional) return false;
        if (p_functional instanceof FunctionalGroup) {
        	return isSubAtCentralNodes((FunctionalGroup)p_functional);
        }
        else if (p_functional instanceof FunctionalGroupCollection) {
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

    /**
   Requires: both graph components belong to the same graph
   Effects: return true if two graph components are equal
   Modifies:
   svp
   */
   //## operation isSymmetric(GraphComponent, GraphComponent)
   public boolean isSymmetric(GraphComponent p_gc1, GraphComponent p_gc2) {
   //#[ operation isSymmetric(GraphComponent, GraphComponent)
   Stack s1 = new Stack();
   Stack s2 = new Stack();
   if (p_gc1.isEquivalent(p_gc2, s1, s2)) {
     resetStack(s1);
       resetStack(s2);
       getGraph().resetMatchedGC();
       return true;
     }
     else {
       resetStack(s1);
       resetStack(s2);
       getGraph().resetMatchedGC();
       return false;
     }
   //#]
 }


    /**
    Requires:
    Effects: return if this chem graph is matched with p_functionalGroup at the central nodes.
    Modifies:
    */
    //## operation isSubAtCentralNodes(FunctionalGroup)
    public boolean isSubAtCentralNodes(FunctionalGroup p_functionalGroup) {
        //#[ operation isSubAtCentralNodes(FunctionalGroup)
        return getGraph().isSubAtCentralNodes(p_functionalGroup.getGraph());



        //#]
    }

    /**
    Requires:
    Effects: factory method for ChemGraph. make a new instance of ChemGraph.  Do such things for it:
    (1) add missing H
    (2) generate chemical formula
    (3) calculate symmetry number
    (4) calculate thermal properties
    Modifies:
    */
    //## operation make(Graph)
    /*public static ChemGraph make(Graph p_graph) throws InvalidChemGraphException, ForbiddenStructureException {
        //#[ operation make(Graph)
		double pT = System.currentTimeMillis();
		Global.makeChemG++;
		Species sp = SpeciesDictionary.getInstance().getSpeciesFromGraph(p_graph);
		ChemGraph cg = null;
		if (sp !=null){
			if (sp.hasResonanceIsomers()){
				Iterator rIter = sp.getResonanceIsomers();
				while(rIter.hasNext()){
					ChemGraph chemGraph = (ChemGraph)rIter.next();
					if (chemGraph.graph.equals(p_graph)){
						chemGraph.graph = p_graph;
						return chemGraph;
						
					}
						
				}
			}
			else {
				sp.chemGraph.graph = p_graph;
				return sp.chemGraph;
			}
		}
		if (cg == null){
			try {
	        	cg = new ChemGraph(p_graph);
	        	cg.addMissingHydrogen();

	        	if (cg.repOk()){
	        		cg.generateChemicalFormula();
	        		//cg.calculateSymmetryNumber();
					Global.symmetryCG+=(System.currentTimeMillis()-pT)/1000/60;
					pT = System.currentTimeMillis();
	        		//cg.generateThermoData();
					Global.thermoCG += (System.currentTimeMillis()-pT)/1000/60;
					pT = System.currentTimeMillis();
	        		//cg.calculateInternalRotor();
					
	        	}
	        	else {
	        		throw new InvalidChemGraphException();
	        	}
	        }
	        catch (ForbiddenStructureException e) {
	        	throw new ForbiddenStructureException(e.getMessage());
	        }
		}
		Global.IRCG+=(System.currentTimeMillis()-pT)/1000/60;
		return cg;

        //#]
    }*/

	 public static ChemGraph make(Graph p_graph, boolean hasHydrogen) throws InvalidChemGraphException, ForbiddenStructureException {
	        //#[ operation make(Graph)
			double pT = System.currentTimeMillis();
			
			ChemGraph cg  = null;
			if (cg == null){
				try {
		        	cg = new ChemGraph(p_graph);
					if (!hasHydrogen)
						cg.addMissingHydrogen();

		        	if (cg.repOk()){
		        		cg.generateChemicalFormula();
						
						
		        	}
		        	else {
		        		throw new InvalidChemGraphException();
		        	}
					//cgd.putSpecies(cg);
		        }
		        catch (ForbiddenStructureException e) {
		        	throw new ForbiddenStructureException(e.getMessage());
		        }
			}
			
			return cg;

	        //#]
	    }
	
	 public static ChemGraph make(Graph p_graph) throws InvalidChemGraphException, ForbiddenStructureException {
	        //#[ operation make(Graph)
			double pT = System.currentTimeMillis();
			
			//ChemGraphDictionary cgd = ChemGraphDictionary.getInstance();
			ChemGraph cg  = null;//= cgd.getChemGraphFromGraph(p_graph);
			if (cg == null){
				try {
					
		        	cg = new ChemGraph(p_graph);
					cg.addMissingHydrogen();

		        	if (cg.repOk()){
		        		cg.generateChemicalFormula();
		        	}
		        	else {
		        		System.err.println(getRepOkString());
		        		throw new InvalidChemGraphException();
		        	}
					//cgd.putSpecies(cg);
		        }
		        catch (ForbiddenStructureException e) {
		        	throw new ForbiddenStructureException(e.getMessage());
		        }
			}
			
			return cg;

	        //#]
	    }
	
	
	public static void addForbiddenStructure(FunctionalGroup fg) {
		forbiddenStructure.add(fg);
	}
	
    /**
    Requires:
    Effects: read in forbidden structure for ChemGraph
    Modifies: this.forbiddenStructure
    */
    //## operation readForbiddenStructure()
    public static void readForbiddenStructure() throws IOException {
        //#[ operation readForbiddenStructure()
        try {
        	String forbiddenStructureFile = System.getProperty("jing.chem.ChemGraph.forbiddenStructureFile");
        	if (forbiddenStructureFile == null) {
        		System.out.println("undefined system property: jing.chem.ChemGraph.forbiddenStructureFile!");
        		System.out.println("No forbidden structure file defined!");
				throw new IOException("Undefined system property: jing.chem.ChemGraph.forbiddenStructureFile");
        		//return;
        	}

        	FileReader in = new FileReader(forbiddenStructureFile);
        	BufferedReader data = new BufferedReader(in);

        	// step 1: read in structure
        	String line = ChemParser.readMeaningfulLine(data);
        	read: while (line != null) {
        		StringTokenizer token = new StringTokenizer(line);
        		String fgname = token.nextToken();
        		Graph fgGraph = null;
        		try {
        			fgGraph = ChemParser.readFGGraph(data);
        		}
        		catch (InvalidGraphFormatException e) {
        			throw new InvalidFunctionalGroupException(fgname + ": " + e.getMessage());
        		}
        		if (fgGraph == null) throw new InvalidFunctionalGroupException(fgname);
        		//FunctionalGroup fg = FunctionalGroup.make(fgname, fgGraph);
        		FunctionalGroup fg = FunctionalGroup.makeForbiddenStructureFG(fgname, fgGraph);

        		forbiddenStructure.add(fg);

        		line = ChemParser.readMeaningfulLine(data);
        	}

            in.close();
        	return;
        }
        catch (Exception e) {
        	throw new IOException(e.getMessage());
        }




        //#]
    }

    /**
    Requires:
    Effects: reset centralNode list in this ChemGraph according to node's CentralID information
    Modifies: this.graph.centralNode
    */
    //## operation refreshCentralNode()
    public void refreshCentralNode() {
        //#[ operation refreshCentralNode()
        getGraph().refreshCentralNode();
        //#]
    }

    /**
    Requires:
    Effects: check four aspects:
    (1) if graph.repOk() defined in Graph class
    (2) if graph is connected
    (3) if the graph contents are atoms and bonds
    (4) if the valency are okay for all atom
    (5) if the radical number is in the limit
    */
    //## operation repOk()
    public boolean repOk() {
        //#[ operation repOk()
        // check if the graph is connected
        if (!getGraph().repOk()){
        	setRepOkString("Something is wrong with the following chemgraph: " + this.toString());
        	return false;
        }

        // a chemical species should be a connected graph
        if (!getGraph().isConnected()) {
        	setRepOkString("The following chemgraph is not a connected graph: " + this.toString());
        	return false;
        }

        // check if the elements stored in graph are atomd/bonds
        if (!graphContentsOk(getGraph())) {
        	setRepOkString("The following chemgraph contains contents other than atoms or bonds: " + this.toString());
        	return false;
        }

        // check if the valency of every atom is satuated
        //if (!valencyOk()) return false;

        // check if the radical number greater than MAX_RADICAL_NUM
        if (getRadicalNumber() > MAX_RADICAL_NUM) {
        	setRepOkString("The following chemgraph exceeds the maximum radicals allowed in a species: " + this.toString());
        	return false;
        }

        // check if the oxygen atom number is too large
        if (getOxygenNumber() > MAX_OXYGEN_NUM) {
        	setRepOkString("The following chemgraph exceeds the maximum oxygens allowed (" + MAX_OXYGEN_NUM + ") in a species:\n" + this.toString());
        	return false;
        }
		if (getCarbonNumber() > MAX_CARBON_NUM) {
        	setRepOkString("The following chemgraph exceeds the maximum carbons allowed (" + MAX_CARBON_NUM + ") in a species:\n" + this.toString());
        	return false;
        }
		if (getSulfurNumber() > MAX_SULFUR_NUM) {
        	setRepOkString("The following chemgraph exceeds the maximum sulfurs allowed (" + MAX_SULFUR_NUM + ") in a species:\n" + this.toString());
        	return false;
        }
		if (getSiliconNumber() > MAX_SILICON_NUM) {
        	setRepOkString("The following chemgraph exceeds the maximum silicons allowed (" + MAX_SILICON_NUM + ") in a species:\n" + this.toString());
        	return false;
        }
		
		if (getHeavyAtomNumber() > MAX_HEAVYATOM_NUM) {
        	setRepOkString("The following chemgraph exceeds the maximum heavy atoms allowed (" + MAX_HEAVYATOM_NUM + ") in a species:\n" + this.toString());
        	return false;
        }
		
        return true;



        //#]
    }

    /**
    Requires:
    Effects: reset reacting site as the pass-in p_site.
    Modifies: this.graph.centralNode
    */
    //## operation resetReactedSite(HashMap)
    public void resetReactedSite(HashMap p_site) throws SiteNotInSpeciesException {
        //#[ operation resetReactedSite(HashMap)
        setCentralNode(p_site);
        //#]
    }

    //## operation resetStack(Stack)
//svp
  public void resetStack(Stack p_stack) {
         //#[ operation resetStack(Stack)
         while (!p_stack.empty()) {
                 GraphComponent gc = (GraphComponent)p_stack.pop();
                 gc.setMatchedGC(null);
         }
         return;
         //#]
     }


    /**
    Requires:
    Effects: reset the only center to the p_node in this chem graph for thermo calculation
    Modifies: this.graph.centralNode, and  centralIDs in its associated nodes.
    */
    //## operation resetThermoSite(Node)
    public void resetThermoSite(Node p_node) {
        //#[ operation resetThermoSite(Node)
        getGraph().clearCentralNode();
        getGraph().setCentralNode(1,p_node);
        return;
        //#]
    }

    /**
    Requires:
    Effects: reset centreNode list as the pass-in p_site.
    Modifies: this.graph.centralNode
    */
    //## operation setCentralNode(HashMap)
    protected void setCentralNode(HashMap p_site) {
        //#[ operation setCentralNode(HashMap)
        try {
        	Graph g = getGraph();
        	g.clearCentralNode();
        	g.setCentralNodes(p_site);
        }
        catch (NotInGraphException e) {
        	throw new SiteNotInSpeciesException();
        }
        //#]
    }

    /**
    Requires:
    Effects: set GTPP as the thermoGAPP of this chem graph
    Modifies: thermoGAPP
    */
    //## operation setDefaultThermoGAPP()
    public void setDefaultThermoGAPP() {
        //#[ operation setDefaultThermoGAPP()
        thermoGAPP = GATP.getINSTANCE();
        return;
        //#]
    }
    
    public void setDefaultTransportGAPP() {
        transportGAPP = GATransportP.getINSTANCE();
        return;
    }

    public void setDefaultSolvationGAPP() {
        //#[ operation setDefaultThermoGAPP()
        SolvationGAPP = GATP_Solvation.getINSTANCE();
        return;
        //#]
    }

    
    public void setDefaultAbramGAPP() {
        //#[ operation setDefaultThermoGAPP()
        abramGAPP = GATP_Abraham.getINSTANCE();
        return;
        //#]
    }

        public void setDefaultUnifacGAPP() {
        //#[ operation setDefaultThermoGAPP()
        unifacGAPP = GATP_Unifac.getINSTANCE();
        return;
        //#]
    }



    /**
    Requires:
    Effects: return a string of this chemgraph.  the string includes two parts:
    (1) chemical formula
    (2) the string for the graph
    Modifies:
    */
    //## operation toString()
    public String toString() {
        //#[ operation toString()
        String s = "ChemFormula: " + getChemicalFormula() + '\n';
        s = s + getGraph().toStringWithoutCentralID();
        return s;
        //#]
    }

	 public String toString(int i) {
	        //#[ operation toString()
	        String s ="";// "ChemFormula: " + getChemicalFormula() + '\n';
	        s = s + getGraph().toStringWithoutCentralID();
	        return s;
	        //#]
	    }

    /**
    Requires:
    Effects: return a short string for this ChemGraph.  it includes two parts:
    (1) chemical formula
    (2) the graph string without H
    Modifies:
    */
    //## operation toStringWithoutH()
    public String toStringWithoutH() {
        //#[ operation toStringWithoutH()
        String s = "ChemFormula: " + getChemicalFormula() + '\n';
        if (getHeavyAtomNumber() == 0)
        	s += getGraph().toStringWithoutCentralID();
        else
        	s = s + getGraph().toStringWithoutCentralIDAndH();
        return s;
        //#]
    }

	/**
    Requires:
    Effects: return a short string for this ChemGraph.  it includes two parts:
    (1) the graph string without H
    Modifies:
    */
    //## operation toStringWithoutH()
    public String toStringWithoutH(int i) {
        //#[ operation toStringWithoutH()
        String s = "";//= "ChemFormula: " + getChemicalFormula() + '\n';
        if (getHeavyAtomNumber() == 0)
        	s += getGraph().toStringWithoutCentralID();
        else
        	s = s + getGraph().toStringWithoutCentralIDAndH();
        return s;
        //#]
    }
    /**
    Requires:
    Effects: check all the atom to see if it satisfies:
    Val of atom = radical number + ion number + sum of bond orders.
    if it is satisfied, return true, otherwise return false.
    Modifies:
    */
    //## operation valencyOk()
    public boolean valencyOk() throws InvalidNodeElementException {
        //#[ operation valencyOk()
        Iterator node_iter = graph.getNodeList();
        while (node_iter.hasNext()) {
        	Node node = (Node)node_iter.next();
        	Atom atom = (Atom)node.getElement();
        	double val = atom.getValency();
        	Iterator arc_iter = node.getNeighbor();
        	while (arc_iter.hasNext()) {
        		Bond bond = (Bond)((Arc)(arc_iter.next())).getElement();
        		val -= bond.getOrder();
        	}
        	if (Math.abs(val) >= 1e-10) return false;
        }
        return true;


        //#]
    }

    public int getMAX_OXYGEN_NUM() {
        return MAX_OXYGEN_NUM;
    }

    public static int getMAX_RADICAL_NUM() {
        return MAX_RADICAL_NUM;
    }

    public static HashSet getForbiddenStructure() {
        return forbiddenStructure;
    }

    public int getInternalRotor() {
		if (internalRotor <0 ){
			calculateInternalRotor();
		}
        return internalRotor;
    }

    /*public HashSet getSymmetryAxis() {
        return symmetryAxis;
    }*/

    protected String getUniqueString() {
        return uniqueString;
    }

    public Graph getGraph() {
        return graph;
    }
	
	public void setGraph(Graph p_graph){
		graph = p_graph;
	}

    public Species getSpecies() {
        return species;
    }

    public void setSpecies(Species p_Species) {
        species = p_Species;
    }

    public GeneralGAPP getThermoGAPP() {
        return thermoGAPP;
    }

    public GeneralAbramGAPP getAbramGAPP() {
        return abramGAPP;
    }

    public void setThermoGAPP(GeneralGAPP p_GeneralGAPP) {
        thermoGAPP = p_GeneralGAPP;
    }

    public void setAbramGAPP(GeneralAbramGAPP p_GeneralAbramGAPP) {
        abramGAPP = p_GeneralAbramGAPP;
    }
    
    public static void setMaxCarbonNumber(int maxCNumber) {
    	MAX_CARBON_NUM = maxCNumber;
    }
    
    public static void setMaxOxygenNumber(int maxONumber) {
    	MAX_OXYGEN_NUM = maxONumber;
    }
    
    public static void setMaxRadicalNumber(int maxRadNumber) {
    	MAX_RADICAL_NUM = maxRadNumber;
    }
    
    public static void setMaxSulfurNumber(int maxSNumber) {
    	MAX_SULFUR_NUM = maxSNumber;
    }
    
    public static void setMaxSiliconNumber(int maxSiNumber) {
    	MAX_SILICON_NUM = maxSiNumber;
    }
    
    public static void setMaxHeavyAtomNumber(int maxHANumber) {
    	MAX_HEAVYATOM_NUM = maxHANumber;
    }
    
    public int getHeavyAtomNumber() {
    	return getCarbonNumber() + getOxygenNumber() + getSulfurNumber() + getSiliconNumber();
    }
    
    public void setRepOkString(String s) {
    	repOkString = s;
    }
    
    public static String getRepOkString() {
    	return repOkString;
    }
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\ChemGraph.java
*********************************************************************/


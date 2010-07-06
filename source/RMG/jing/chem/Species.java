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

import jing.mathTool.*;
//import jing.mathTool.Queue;
import jing.chemUtil.*;
import jing.chemParser.*;
import jing.chemUtil.Node;
import jing.param.Global;
import jing.param.Temperature;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\Species.java
//----------------------------------------------------------------------------

//## class Species
public class Species {

	protected boolean GATPFitExecuted = false;
    protected int ID;		//## attribute ID

    protected static int TOTAL_NUMBER = 0;		//## attribute TOTAL_NUMBER
    protected static boolean addID = true;
    protected ChemGraph chemGraph;		//## attribute chemGraph

    /**
    The initial value is the parameter for N2, now it is treated as a default value if no further detailed information provided.
    unit: cm-1
    */
    protected double deltaEDown = 461;		//## attribute deltaEDown

    protected String name = null;		//## attribute name

    protected LinkedHashSet resonanceIsomers = new LinkedHashSet();		//## attribute resonanceIsomers

    protected boolean therfitExecuted = false;		//## attribute therfitExecuted
    
    protected String InChI = null;	//## attribute InChI

    protected TransportData chemkinTransData;
	protected NASAThermoData nasaThermoData;
	protected String nasaThermoSource;
    protected ThreeFrequencyModel threeFrequencyModel;
	//protected WilhoitThermoData wilhoitThermoData;

	/**
	 * The spectroscopic data for the species (vibrational frequencies,
	 * rotational frequencies, symmetry number, and hindered frequency-barrier 
	 * pairs. Will eventually replace ThreeFrequencyModel.
	 */
	protected SpectroscopicData spectroscopicData;
	
    // Flag to tag certain species as library only... i.e. we won't try them against RMG templates.
    // They will only react as defined in the primary kinetic library.   GJB
    protected boolean IsReactive = true; 
    
    // flag to tag certain species as having a constant concentration 
    // they won't be consumed by reactions (will lead to mass creation!)
    protected boolean constantConcentration = false;
    
    // Flag which specifies whether to generate InChIs
	public static boolean useInChI = false;
    public static boolean useSolvation = false;
    
    protected HashSet paths;
    
    // Constructors

    //## operation Species()
    private  Species() {
        initRelations();
        //#[ operation Species()
        //#]
    }
    //## operation Species(String,ChemGraph)
    private  Species(int id, String p_name, ChemGraph p_chemGraph) {
        initRelations();
		ID = id;
        //#[ operation Species(String,ChemGraph)
        name = p_name;
        chemGraph = p_chemGraph;
        generateResonanceIsomers();
        if (!constantConcentration) {
            findStablestThermoData();
        } else {
            //findSolvationData();
        }
        calculateTransportParameters();
        selectDeltaEDown();
		generateNASAThermoData();
       // generateSpectroscopicData(); // only get it if you need it!!!
		/*
		 * MRH 9MAR2010:
		 * Commenting InChI generation when making a new species.
		 * 
		 * Presently, the only time an InChI is needed in any RMG "default package"
		 * 	class is in the post-processing (e.g. write the InChI with the
		 * 	RMG_Dictionary.txt file, write the InChI with the chem.inp file).
		 * 	Generating an InChI for every species made slows down RMG and is
		 * 	unnecessary.  RMG will now only generate InChIs for the core species.
		 */
//        if (useInChI) InChI = p_chemGraph.getInChI();
        //#]
    }

    //## operation addResonanceIsomer(ChemGraph)
    public boolean addResonanceIsomer(ChemGraph p_resonanceIsomer) {
        //#[ operation addResonanceIsomer(ChemGraph)
        if (resonanceIsomers == null) resonanceIsomers = new LinkedHashSet();

        p_resonanceIsomer.setSpecies(this);
        return resonanceIsomers.add(p_resonanceIsomer);



        //#]
    }

    //## operation calculateCp(Temperature)
    public double calculateCp(Temperature p_temperature) {
        //#[ operation calculateCp(Temperature)
        return getThermoData().calculateCp(p_temperature);
        //#]
    }

    //## operation calculateG(Temperature)
    public double calculateG(Temperature p_temperature) {
        //#[ operation calculateG(Temperature)
        //return getThermoData().calculateG(p_temperature);
		return nasaThermoData.calculateFreeEnergy(p_temperature);
        //#]
    }

    //## operation calculateGLowerBound(Temperature)
    //svp
      public double calculateGLowerBound(Temperature p_temperature) {
        //#[ operation calculateGLowerBound(Temperature)
        return getThermoData().calculateGLowerBound(p_temperature);
        //#]
      }

    //## operation calculateGUpperBound(Temperature)
    //svp
      public double calculateGUpperBound(Temperature p_temperature) {
        //#[ operation calculateGUpperBound(Temperature)
        return getThermoData().calculateGUpperBound(p_temperature);
        //#]
      }


    //## operation calculateH(Temperature)
    public double calculateH(Temperature p_temperature) {
        //#[ operation calculateH(Temperature)
        //return getThermoData().calculateH(p_temperature);
        return nasaThermoData.calculateEnthalpy(p_temperature);
		
		//#]
    }

    //## operation calculateLJParameters()
    public void calculateTransportParameters() {
    	if (hasResonanceIsomers()){
			Iterator cgIter = getResonanceIsomers();
			/*
			 * TransportData will be found for all resonance structures for a
			 * 	given species.  The purpose of this while loop is to determine
			 * 	if one of the resonance structures exists in a primary transport
			 * 	library.  So, the while loop will continue until:
			 * 		(a) all resonance structures are exhausted
			 * 		(b) the chemgraph is found in a primary transport library
			 * If (a), we want the transportdata of the most stable chemgraph
			 * 	(for consistency)
			 * If (b), we want the found transport data
			 */
			while (cgIter.hasNext()) {
				ChemGraph cg = (ChemGraph)cgIter.next();
				chemkinTransData = cg.getTransportData();
				if (chemkinTransData.fromPrimaryTransportLibrary) break;
			}
			if (!chemkinTransData.fromPrimaryTransportLibrary)
				chemkinTransData = getChemGraph().getTransportData();
		}
    	else chemkinTransData = getChemGraph().getTransportData();
    	
        //int cNum = getChemGraph().getCarbonNumber();
//		int cNum = getChemGraph().getHeavyAtomNumber();
//
//        selectLJParametersForSpecialMolecule();
//        if (cNum == 1) chemkinTransData = new TransportData(3.758, 148.6);
//        else if (cNum == 2) chemkinTransData = new TransportData(4.443, 110.7);
//        else if (cNum == 3) chemkinTransData = new TransportData(5.118, 237.1);
//        else if (cNum == 4) chemkinTransData = new TransportData(4.687, 531.4);
//        else if (cNum == 5) chemkinTransData = new TransportData(5.784, 341.1);
//        else chemkinTransData = new TransportData(5.949, 399.3);

        return;
    }

    //## operation calculateS(Temperature)
    public double calculateS(Temperature p_temperature) {
        //#[ operation calculateS(Temperature)
        //return getThermoData().calculateS(p_temperature);
        return nasaThermoData.calculateEntropy(p_temperature);
		//#]
    }


    //## operation doDelocalization(ChemGraph,Stack)
    private ChemGraph doDelocalization(ChemGraph p_chemGraph, Stack p_path) {
        //#[ operation doDelocalization(ChemGraph,Stack)
        if (p_path.isEmpty() || p_path.size() != 3) throw new InvalidDelocalizationPathException();

        Graph graph = Graph.copy(p_chemGraph.getGraph());

        // n1-a1-n2-a2-n3.
        Node node1 = graph.getNodeAt((Integer)p_path.pop());
        Node node2 = graph.getNodeAt((Integer)p_path.pop());
        Node node3 = graph.getNodeAt((Integer)p_path.pop());
        Arc arc1 = graph.getArcBetween(node1,node2);
        Arc arc2 = graph.getArcBetween(node2,node3);

        // deal with node1
        Atom atom1 = (Atom)node1.getElement();
        Atom newAtom1 = (Atom)atom1.changeRadical(1,null);
        node1.setElement(newAtom1);

        // deal with arc1
        Bond bond1 = (Bond)arc1.getElement();
        Bond newBond1 = bond1.changeBond(-1);
        arc1.setElement(newBond1);

        // deal with node2, actually do nothing

        // deal with arc2
        Bond bond2 = (Bond)arc2.getElement();
        Bond newBond2 = bond2.changeBond(1);
        arc2.setElement(newBond2);

        // deal with node3
        Atom atom3 = (Atom)node3.getElement();
        Atom newAtom3 = (Atom)atom3.changeRadical(-1,null);
        node3.setElement(newAtom3);

        node1.updateFgElement();
        node1.updateFeElement();
        node2.updateFgElement();
        node2.updateFeElement();
        node3.updateFgElement();
        node3.updateFeElement();

        p_path = null;

        try {
        	ChemGraph newIsomer = ChemGraph.make(graph);

        	if (addResonanceIsomer(newIsomer)) return newIsomer;
        	else return null;
        }
        catch (ForbiddenStructureException e) {
        	return null;
        }
        //#]
    }

    //## operation findAllDelocalizationPaths(Node)
    private HashSet findAllDelocalizationPaths(Node p_radical) {
        //#[ operation findAllDelocalizationPaths(Node)
        HashSet allPaths = new HashSet();
        Atom atom = (Atom)p_radical.getElement();
        if (!atom.isRadical()) return allPaths;

        Iterator iter = p_radical.getNeighbor();
        while (iter.hasNext()) {
        	Arc arc1 = (Arc)iter.next();
        	Bond bond1 = (Bond)arc1.getElement();
        	if (bond1.isSingle() || bond1.isDouble()) {
        		Node node1 = arc1.getOtherNode(p_radical);
        		Iterator iter2 = node1.getNeighbor();
        		while (iter2.hasNext()) {
        			Arc arc2 = (Arc)iter2.next();
        			if (arc2 != arc1) {
        				Bond bond2 = (Bond)arc2.getElement();
        				if (bond2.isDouble() || bond2.isTriple()) {
        					Node node2 = arc2.getOtherNode(node1);
        					Stack path = new Stack();
        					path.push(p_radical.getID());
        					path.push(node1.getID());
        					path.push(node2.getID());
        					allPaths.add(path);
        				}
        			}
        		}

        	}
        }

        return allPaths;
        //#]
    }

    //## operation findStablestThermoData()
    public void findStablestThermoData() {
        //#[ operation findStablestThermoData()
        double H = chemGraph.getThermoData().getH298();
        ChemGraph stablest = chemGraph;
        if (resonanceIsomers != null) {
        	Iterator iter = resonanceIsomers.iterator();
        	while (iter.hasNext()) {
        		ChemGraph g = (ChemGraph)iter.next();
        		
        		double newH = g.getThermoData().getH298();
        		if (g.fromprimarythermolibrary){
        			H = newH;
        			stablest = g;
                    chemGraph = stablest; // presumably this should be before the return? rwest
        			return;
        		}
        			
        		if (newH < H) {
        			H = newH;
        			stablest = g;
        		}
        	}
        }

        chemGraph = stablest;
        //#]
    }
    
	public void generateNASAThermoData() {
        //nasaThermoData = Therfit.generateNASAThermoData(this);
        nasaThermoData = GATPFit.generateNASAThermoData(this);
        nasaThermoSource = getThermoData().source;
		GATPFitExecuted = (nasaThermoData != null);
    }


    /**
    Requires:
    Effects: generate all the possible resonance isomers for the primary chemGraph
    Modifies: this.resonanceIsomers
    */
    //## operation generateResonanceIsomers()
    public void generateResonanceIsomers() {
        //#[ operation generateResonanceIsomers()
        if (chemGraph == null) throw new NullPointerException();

        // check if the representation of chemGraph is correct
        if (!chemGraph.repOk()) {
        	resonanceIsomers = null;
        	throw new InvalidChemGraphException();
        }

        // generate RI for radical
        generateResonanceIsomersFromRadicalCenter();

        // generaate RI for O2, removed, don't allow .o-o.
        //generateResonanceIsomersForOxygen();

        if (resonanceIsomers.size() == 1) {
        	ChemGraph cg = (ChemGraph)(resonanceIsomers.iterator().next());
        	if (cg == chemGraph) resonanceIsomers.clear();
        	else addResonanceIsomer(chemGraph);
        }
        
        if (chemGraph.getRadicalNumber() >= 2){
        	//find if there are radicals next to each other and in that case 
        	//increase the bond order by 1
        	Iterator radicalnodeIter = chemGraph.getRadicalNode().iterator();
        	while (radicalnodeIter.hasNext()){
        		Node n1 = (Node)radicalnodeIter.next();
        		Iterator arcs = n1.getNeighbor();
        		while (arcs.hasNext()){
        			Arc arc1 = (Arc)arcs.next();
        			Bond bond1 = (Bond)arc1.getElement();
        			Node n2 = arc1.getOtherNode(n1);
        			Atom a2 = (Atom)n2.getElement();
        			if (a2.isRadical() && !bond1.isTriple()){
        			   Graph newG = Graph.copy(chemGraph.getGraph());

        			   Node newn1 = newG.getNodeAt(n1.getID());
        			   Atom newa1 = (Atom)newn1.getElement();
        			   newa1.changeRadical(-1, null);
        			   newn1.setElement(newa1);

        			   Node newn2 = newG.getNodeAt(n2.getID());
        			   Atom newa2 = (Atom)newn2.getElement();
        			   newa2.changeRadical(-1,null);
        			   newn2.setElement(newa2);

        			   Arc newarc1 = newG.getArcBetween(n2.getID(), n1.getID());
        			   Bond newb1 = (Bond)newarc1.getElement();
        			   newb1.changeBond(1);
        			   newarc1.setElement(newb1);
        			   ChemGraph newchemG;
        			   try {
						  newchemG = ChemGraph.make(newG);
						   addResonanceIsomer(newchemG);

					   } catch (InvalidChemGraphException e) {
						  // TODO Auto-generated catch block
						  e.printStackTrace();
					   } catch (ForbiddenStructureException e) {
						  // TODO Auto-generated catch block
						  e.printStackTrace();
					   }
        			}
        			
        		}
        	}
        }
        /*Graph g = Graph.copy(chemGraph.getGraph());
        // generate node-electron stucture
        int nodeNumber = g.getNodeNumber();
        int [] electronOnNode = new int[nodeNumber+1];

        Iterator nodeIter = g.getNodeList();
        while (nodeIter.hasNext()) {
        	Node node = (Node)nodeIter.next();
        	Atom atom = (Atom)node.getElement();
        	if (atom.isRadical()) {
        		int ID = node.getID().intValue();
        		electronOnNode[ID] = atom.getRadicalNumber();
        		Iterator arcIter = node.getNeighbor();
        		arcLoop:while (arcIter.hasNext()) {
        			Arc arc = (Arc)arcIter.next();
        			Bond bond = (Bond)arc.getElement();
        			if (bond.isBenzene() || bond.isTriple()) {
        				electronOnNode[ID] = 1;
        				break arcLoop;
        			}
        			else if (bond.isDouble()) {
        				electronOnNode[node.getID().intValue()-1] += 1;
        			}
        		}
        	}
        }
        */

        // combine them accordingly
        //resonanceIsomer = combineElectronBetweenNodes(g);

        return;




        //#]
    }

    //## operation generateResonanceIsomersForOxygen()
    public void generateResonanceIsomersForOxygen() {
        //#[ operation generateResonanceIsomersForOxygen()
        // form a O2 graph
        Graph g1 = new Graph();
        Node n1 = g1.addNodeAt(1,Atom.make(ChemElement.make("O"), FreeElectron.make("0")));
        Node n2 = g1.addNodeAt(2,Atom.make(ChemElement.make("O"), FreeElectron.make("0")));
        g1.addArcBetween(n1,Bond.make("D"),n2);

        // form a .o-o. graph
        Graph g2 = new Graph();
        Node n3 = g2.addNodeAt(1,Atom.make(ChemElement.make("O"), FreeElectron.make("1")));
        Node n4 = g2.addNodeAt(2,Atom.make(ChemElement.make("O"), FreeElectron.make("1")));
        g2.addArcBetween(n3,Bond.make("S"),n4);

        try {
        	if (chemGraph.getGraph().isEquivalent(g1)) {
        		ChemGraph cg = ChemGraph.make(g2);
        		addResonanceIsomer(cg);
        		}
        	else if (chemGraph.getGraph().isEquivalent(g2)) {
        		ChemGraph cg = ChemGraph.make(g1,true);
        		addResonanceIsomer(cg);
        	}
        	return;
        }
        catch (ForbiddenStructureException e) {
        	return;
        }



        //#]
    }

    //## operation generateResonanceIsomersFromRadicalCenter()
    private void generateResonanceIsomersFromRadicalCenter() {
        //#[ operation generateResonanceIsomersFromRadicalCenter()
    	/*
    	 * Updated on 27-May-2009 by Michael Harper
    	 * 		(in consultation w/Richard West and Josh Allen)
    	 * Problem: Running RMG produced a "Exception in thread "main" ...
    	 * 		java.lang.ArrayIndexOutOfBoundsException at ...
    	 * 		jing.mathTool.Queue.enqueue(Queue.java)", etc.
    	 * The out-of-bounds index was the undoChemGraph index.  The size
    	 * 		of the undoChemGraph was set to 4 times the number
    	 * 		of atoms in the chemgraph.  For some species, this
    	 * 		value was exceeded, causing RMG to fail.
    	 * Solution: Change undoChemGraph from Queue() to LinkedList().
    	 * 		Note: RMG-defined Queue.java is obsolete
    	 * Examples: Try running a condition.txt file with the following
    	 * 		sets of species.  Before, an exception would be caught.
    	 * Test Case 1:
    	 * 		1 C 0 {2,T}
    	 * 		2 C 0 {1,T} {3,S}
    	 * 		3 C 0 {2,S} {4,T}
    	 * 		4 C 0 {3,T}
    	 * 
    	 * 		1 C 1 {2,D}
    	 * 		2 C 0 {1,D} {3,D}
    	 * 		3 C 1 {2,D}
    	 * Test Case 2:
    	 * 		1 C 2 {2,D}
    	 * 		2 O 0 {1,D}
    	 * 
    	 * 		1 C 1 {2,D}
    	 * 		2 C 0 {1,D} {3,D}
    	 * 		3 C 0 {2,D} {4,D}
    	 * 		4 O 0 {3,D}
    	 */
    	
    	
        // only radical is considered here
        if (chemGraph.getRadicalNumber() <= 0) return;

        addResonanceIsomer(chemGraph);
        
        LinkedList undoChemGraph = new LinkedList();
        undoChemGraph.add(chemGraph);
//        Queue undoChemGraph = new Queue(4*chemGraph.getAtomNumber());
//        undoChemGraph.enqueue(chemGraph);

        HashSet processedChemGraph = new HashSet();
        while (!undoChemGraph.isEmpty()) {
        	ChemGraph cg = (ChemGraph)undoChemGraph.remove();
//        	ChemGraph cg = (ChemGraph)undoChemGraph.dequeue();
        	HashSet radicalNode = cg.getRadicalNode();
        	Iterator radicalIter = radicalNode.iterator();
        	while (radicalIter.hasNext()) {
        		Node radical = (Node)radicalIter.next();
        		int radicalNumber = ((Atom)radical.getElement()).getRadicalNumber();
        		if (radicalNumber > 0) {
        			HashSet allPath = findAllDelocalizationPaths(radical);
        			Iterator pathIter = allPath.iterator();
        			while (pathIter.hasNext()) {
        				Stack path = (Stack)pathIter.next();
        				ChemGraph newCG = doDelocalization(cg, path);
        				if (newCG!=null && !processedChemGraph.contains(newCG)) {
        					undoChemGraph.add(newCG);
//        					undoChemGraph.enqueue(newCG);
        				}
        			}
        		}
        		
        	}
        	processedChemGraph.add(cg);
        }
        /*for (Iterator iter = getResonanceIsomers(); iter.hasNext(); ){
        	ChemGraph cg = (ChemGraph)iter.next();
        	makeSingletAndTriplet(cg);
        	
        }*/
        
        //#]
    }

    private void makeSingletAndTriplet(ChemGraph cg) {
    	HashSet radicalNode = cg.getRadicalNode();
    	Iterator radicalIter = radicalNode.iterator();
    	while (radicalIter.hasNext()) {
    		Node radical = (Node)radicalIter.next();
    		int radicalNumber = ((Atom)radical.getElement()).getRadicalNumber();
    		if (radicalNumber == 2 && radical.getFeElement().spin == null) {
    			// make singlet
    			Graph graph = Graph.copy(cg.getGraph());
    			int nodeID = radical.getID();
    			Node newNode = graph.getNodeAt(nodeID);
    			newNode.getFeElement().spin = "S";
    			
    			try {
    				ChemGraph newIsomer = ChemGraph.make(graph);
    				
    				addResonanceIsomer(newIsomer);
    			}
    			catch (ForbiddenStructureException e) {
    				
    			}
    			
//    			 make triplet
    			Graph graph2 = Graph.copy(cg.getGraph());
    			
    			Node newNode2 = graph.getNodeAt(nodeID);
    			newNode.getFeElement().spin = "T";
    			try {
    				ChemGraph newIsomer = ChemGraph.make(graph2);
    				
    				addResonanceIsomer(newIsomer);
    			}
    			catch (ForbiddenStructureException e) {
    				
    			}
    		}
    		else if (radicalNumber == 2 && radical.getFeElement().spin.equals("T")) {
    			// make singlet
    			Graph graph = Graph.copy(cg.getGraph());
    			int nodeID = radical.getID();
    			Node newNode = graph.getNodeAt(nodeID);
    			newNode.getFeElement().spin = "S";
    			
    			try {
    				ChemGraph newIsomer = ChemGraph.make(graph);
    				
    				addResonanceIsomer(newIsomer);
    			}
    			catch (ForbiddenStructureException e) {
    				
    			}
    			
//    			 make nonSinglet
    			Graph graph2 = Graph.copy(cg.getGraph());
    			
    			Node newNode2 = graph.getNodeAt(nodeID);
    			newNode.getFeElement().spin = null;
    			try {
    				ChemGraph newIsomer = ChemGraph.make(graph2);
    				
    				addResonanceIsomer(newIsomer);
    			}
    			catch (ForbiddenStructureException e) {
    				
    			}
    		}
    		else if (radicalNumber == 2 && radical.getFeElement().spin.equals("S")) {
    			// make triplet
    			Graph graph = Graph.copy(cg.getGraph());
    			int nodeID = radical.getID();
    			Node newNode = graph.getNodeAt(nodeID);
    			newNode.getFeElement().spin = "T";
    			try {
    				ChemGraph newIsomer = ChemGraph.make(graph);
    				
    				addResonanceIsomer(newIsomer);
    			}
    			catch (ForbiddenStructureException e) {
    				
    			}
    			
    			
//   			 make nonSinglet
    			Graph graph2 = Graph.copy(cg.getGraph());
    			
    			Node newNode2 = graph.getNodeAt(nodeID);
    			newNode.getFeElement().spin = null;
    			try {
    				ChemGraph newIsomer = ChemGraph.make(graph2);
    				
    				addResonanceIsomer(newIsomer);
    			}
    			catch (ForbiddenStructureException e) {
    				
    			}
    		}
    		
    	}
		
	}
	
	public void generateSpectroscopicData() {
		if (SpectroscopicData.mode == SpectroscopicData.Mode.THREEFREQUENCY) {
			generateThreeFrequencyModel();
			spectroscopicData = null;
		}
		else if (SpectroscopicData.mode == SpectroscopicData.Mode.FREQUENCYGROUPS) {
			spectroscopicData = FrequencyGroups.getINSTANCE().generateFreqData(this);
			threeFrequencyModel = null;
		}
		else {
			spectroscopicData = null;
			threeFrequencyModel = null;
		}
	}
	
	public void generateThreeFrequencyModel() {
        
		threeFrequencyModel = null;

		// Do nothing if molecule is triatomic or smaller
		if (isTriatomicOrSmaller()) 
			return;

        // Create three frequency model for this species
        threeFrequencyModel = Therfit.generateThreeFrequencyModel(this);

		therfitExecuted = (threeFrequencyModel != null);
    }

	  public String getChemkinName() {
	        //#[ operation getChemkinName()
		  if (addID){
			  String chemkinName = getName() + "(" + getID() + ")";
			  /* Updated by MRH on 1-Jun-2008
			  	If statement used to check if chemkinName length was greater than 16
			  	I've changed it to length 10.  Chemkin format dictates that any 
			  	"A+B(+m)=C+D(+m)" be of length 52 or less.  Assuming we have 2 reactants
			  	and 2 products for a pressure-dependent network, we already use:
			  		- 2 characters (for the + symbol)
			  		- 1 character (for the = symbol)
			  		- 8 characters (4 for each of the (+m) symbols)
			  	This leaves 52 - 11 - 41 characters for 4 species.  Thus, I've changed
			  	the criteria from length 16 to length 10.
			  */
			  
		        if (chemkinName.length() > 10) chemkinName = "SPC(" + getID() + ")";
		        return chemkinName;
		  }
		  else
			  return getName();

	        //#]
	    }

    //## operation getInternalRotor()
    public int getInternalRotor() {
        //#[ operation getInternalRotor()
        return getChemGraph().getInternalRotor();
        //#]
    }

    //## operation getMolecularWeight()
    public double getMolecularWeight() {
        //#[ operation getMolecularWeight()
        return getChemGraph().getMolecularWeight();
        //#]
    }

	   //## operation getNasaThermoData()
    public NASAThermoData getNasaThermoData() {
        //#[ operation getNasaThermoData()
        //if (nasaThermoData==null && !therfitExecuted) generateNASAThermoData();
		if (nasaThermoData==null) generateNASAThermoData();
        return nasaThermoData;
        //#]
    }

    //## operation getResonanceIsomers()
    public Iterator getResonanceIsomers() {
        //#[ operation getResonanceIsomers()
        return resonanceIsomers.iterator();
        //#]
    }

	public HashSet getResonanceIsomersHashSet(){
		return resonanceIsomers;
	}

    /**
    Requires:
    Effects: return the thermo data of the stablest resonance isomer.
    Modifies:
    */
    //## operation getThermoData()
    public ThermoData getThermoData() {
        //#[ operation getThermoData()
        return chemGraph.getThermoData();
        //#]
    }

    public ThermoData getSolvationDAta(){
        return chemGraph.getSolvationData();
        //#]

    }

    //## operation getThreeFrequencyMode()
    public ThreeFrequencyModel getThreeFrequencyMode() {
        //#[ operation getThreeFrequencyMode()
        if (threeFrequencyModel==null && !therfitExecuted) generateThreeFrequencyModel();

        return threeFrequencyModel;
        //#]
    }

    //## operation hasResonanceIsomers()
    public boolean hasResonanceIsomers() {
        //#[ operation hasResonanceIsomers()
        if (resonanceIsomers == null) return false;
        else return (resonanceIsomers.size() > 0);
        //#]
    }

    //## operation hasThreeFrequencyModel()
    public boolean hasThreeFrequencyModel() {
        //#[ operation hasThreeFrequencyModel()
        return (getThreeFrequencyModel() != null);
        //#]
    }
	
	public boolean hasSpectroscopicData() {
        return (spectroscopicData != null || threeFrequencyModel != null);
    }

    //## operation isRadical()
    public boolean isRadical() {
        //#[ operation isRadical()
        return chemGraph.isRadical();
        //#]
    }

    //## operation isTriatomicOrSmaller()
    public boolean isTriatomicOrSmaller() {
        //#[ operation isTriatomicOrSmaller()
        return (getChemGraph().getAtomNumber()<=3);
        //#]
    }

	public boolean isMonatomic() {
        return (getChemGraph().getAtomNumber() == 1);
    }

    //## operation make(String,ChemGraph)
    public static Species make(String p_name, ChemGraph p_chemGraph) {
        //#[ operation make(String,ChemGraph)
		double pT = System.currentTimeMillis();
        SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
        
        // first try to get it from the dictionary (which now uses a cache to speed it up)
        Species spe = (Species)(dictionary.getSpecies(p_chemGraph));
        
        // if it wasn't there then it's unique and we need to add it
        if (spe == null) {
        	
        	String name = p_name;
        	if (name == null || name.length()==0) {
        		name = p_chemGraph.getChemicalFormula();
        	}
			int id= ++TOTAL_NUMBER;
			
        	spe = new Species(id,name,p_chemGraph);
        	//spe.ID =
        	dictionary.putSpecies(spe, true);
        	
			// DEBUG: Tell console I made this species
			System.out.println("Created new species: " + spe.getName() + "(" + spe.getID() + ")");

        }
        else {
			if (spe.chemGraph.equals(p_chemGraph)){
				//spe.chemGraph.graph = p_chemGraph.graph;
				//p_chemGraph = spe.chemGraph;

				p_chemGraph.thermoData = spe.chemGraph.thermoData;
				p_chemGraph.symmetryNumber = spe.chemGraph.symmetryNumber;
				p_chemGraph.internalRotor = spe.chemGraph.internalRotor;
                p_chemGraph.solvthermoData = spe.chemGraph.solvthermoData;
			}
			else if (spe.hasResonanceIsomers()){
				Iterator cgIter = spe.getResonanceIsomers();
				while(cgIter.hasNext()){
					ChemGraph cg = (ChemGraph)cgIter.next();
					if (cg.equals(p_chemGraph)){
						p_chemGraph.thermoData = spe.chemGraph.thermoData;
						p_chemGraph.symmetryNumber = spe.chemGraph.symmetryNumber;
						p_chemGraph.internalRotor = spe.chemGraph.internalRotor;
                        p_chemGraph.solvthermoData = spe.chemGraph.solvthermoData;
						break;
					}
				}
			}
			else {
				System.out.println("Cannot make species which has a chemgraph: "+p_chemGraph.toString());
				System.exit(0);
			}
        }
		p_chemGraph.setSpecies(spe);
		Global.makeSpecies += (System.currentTimeMillis()-pT)/1000/60;
        
        /*
        // added by rwest 2009/05/07 to see how many species are considered and in what order
        // N.B. this file is not cleared at the start of a run; results are just appended
        String restartFileContent="";
		try{
			File consideredSpecies = new File ("Restart/consideredSpecies.txt");
			FileWriter fw = new FileWriter(consideredSpecies, true);
			restartFileContent = restartFileContent + spe.getChemkinName() + " \n ";  // name and number
		//	restartFileContent = restartFileContent + spe.toString(1) + "\n\n";  // full chemgraph
			fw.write(restartFileContent);
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart consideredSpecies file");
	    	System.exit(0);
		}
        */
        
        return spe;
    
    }

//	## operation make(String,ChemGraph)
    /*public static Species make(String p_name, Graph p_graph) throws InvalidChemGraphException, ForbiddenStructureException {
        //#[ operation make(String,ChemGraph)
		double pT = System.currentTimeMillis();
        SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
        Species spe = dictionary.getSpeciesFromGraph(p_graph);
		
        if (spe == null) {
        	ChemGraph cg = ChemGraph.make(p_graph);
			spe = make(null, cg);
			cg.setSpecies(spe);
			
        }
        
        return spe;
        //#]
    }*/
	
	/*public static Species make(String p_name, ChemGraph p_chemGraph, int id) {
        //#[ operation make(String,ChemGraph)
        SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
        Species spe = (Species)(dictionary.getSpecies(p_chemGraph));

        if (spe == null) {
        	String name = p_name;
        	if (name == null || name.length()==0) {
        		name = p_chemGraph.getChemicalFormula();
        	}
        	spe = new Species(id, name,p_chemGraph);
			if (id > TOTAL_NUMBER) TOTAL_NUMBER=id;
        	dictionary.putSpecies(spe, false);
 
        }
        p_chemGraph.setSpecies(spe);
        return spe;
        //#]
    }*/
    
    /**
     * Species.make
     * 	As of 4-Sept-2009, this function is only called when making species
     * 		read in from "Restart" files.  The Restart files are generated
     * 		after an integration step in RMG, meaning all stored species (core
     * 		and edge) should be unique.  Therefore, we do not need to check
     * 		each new species against the SpeciesDictionary.
     * MRH 4-Sept-2009
     */
	public static Species make(String p_name, ChemGraph p_chemGraph, int id) {
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		
		String name = p_name;
		Species spe = new Species(id, name, p_chemGraph);
		if (id > TOTAL_NUMBER) TOTAL_NUMBER=id;
		dictionary.putSpecies(spe, false);
        p_chemGraph.setSpecies(spe);
        
     // DEBUG: Tell console I made this species
		System.out.println("Created new species: " + spe.getName() + "(" + spe.getID() + ")");
        
        return spe;
    }

    //## operation repOk()
    public boolean repOk() {
        //#[ operation repOk()
        if (name == null || name.length() == 0) return false;

        if (ID < 0 || ID > TOTAL_NUMBER) return false;

        if (chemGraph != null && !chemGraph.repOk()) return false;

        Iterator iter = resonanceIsomers.iterator();
        while (iter.hasNext()) {
        	ChemGraph cg = (ChemGraph)iter.next();
        	if (cg == null || !cg.repOk()) return false;
        }

        return true;
        //#]
    }

    //## operation selectDeltaEDown()
    public void selectDeltaEDown() {
        //#[ operation selectDeltaEDown()
        String name = getName();

        if (name.equals("AR")) deltaEDown = 374.0;
        else if (name.equals("H2")) deltaEDown = 224.0;
        else if (name.equals("O2")) deltaEDown = 517.0;
        else if (name.equals("N2")) deltaEDown = 461.0;
        else if (name.equals("He")) deltaEDown = 291.0;
        else if (name.equals("CH4")) deltaEDown = 1285.0;
        else if (name.equals("HO2")) deltaEDown = 975.0;
        else if (name.equals("H2O2")) deltaEDown = 975.0;
        else if (name.equals("CO2")) deltaEDown = 417.0;
        else if (name.equals("CO")) deltaEDown = 283.0;

        return;



        //#]
    }

    //## operation selectLJParametersForSpecialMolecule()
    public void selectLJParametersForSpecialMolecule() {
        //#[ operation selectLJParametersForSpecialMolecule()
        String name = getName();
        if (name.equals("H2")) chemkinTransData = new TransportData(2.8327, 59.7);
        else if (name.equals("O2")) chemkinTransData = new TransportData(3.467, 106.7);
        else if (name.equals("H2O")) chemkinTransData = new TransportData(2.641,809.1);
        else if (name.equals("H2O2")) chemkinTransData = new TransportData(4.196,289.3);
        else if (name.equals("CO2")) chemkinTransData = new TransportData(3.941,195.2);
        else if (name.equals("CO")) chemkinTransData = new TransportData(3.690,91.7);
        return;
        //#]
    }
	public String toChemkinString() {
		//#[ operation toChemkinString()
		return getChemkinName();
		//#]
	}

    //## operation toString()
    public String toString() {
        //#[ operation toString()
        String s = "Species " + String.valueOf(ID) + '\t';
        s = s + "Name: " + getName() + '\n';
        ChemGraph primary = getChemGraph();
        s = s + primary.toString();

        int index=0;
        for (Iterator iter = resonanceIsomers.iterator(); iter.hasNext(); ) {
        	index++;
        	ChemGraph isomer = (ChemGraph)iter.next();
        	if (isomer!=primary) {
        		s = s + '\n';
        		s = s + "isomer" + String.valueOf(index) + ":\n";
        		s = s + isomer.toString();
        	}
        }
        return s;
    }

//		## operation toString()
	    public String toString(int i) {
	        //#[ operation toString()
			String s ="";
			ChemGraph primary = getChemGraph();
	        s = s + primary.toString(i);

	        int index=0;
	        /*for (Iterator iter = resonanceIsomers.iterator(); iter.hasNext(); ) {
	        	index++;
	        	ChemGraph isomer = (ChemGraph)iter.next();
	        	if (isomer!=primary) {
	        		s = s + '\n';
	        		s = s + "isomer" + String.valueOf(index) + ":\n";
	        		s = s + isomer.toStringWithoutH(i);
	        	}
	        }*/
	        return s;
    }

    //## operation toStringWithoutH()
    public String toStringWithoutH() {
        //#[ operation toStringWithoutH()
        String s = "Species " + String.valueOf(ID) + '\t';
        s = s + "Name: " + getName() + '\n';
        ChemGraph primary = getChemGraph();
        s = s + primary.toStringWithoutH();

        int index=0;
        for (Iterator iter = resonanceIsomers.iterator(); iter.hasNext(); ) {
        	index++;
        	ChemGraph isomer = (ChemGraph)iter.next();
        	if (isomer!=primary) {
        		s = s + '\n';
        		s = s + "isomer" + String.valueOf(index) + ":\n";
        		s = s + isomer.toStringWithoutH();
        	}
        }
        return s;




        //#]
    }

	public String toStringWithoutH(int i) {
        //#[ operation toStringWithoutH()
        String s ="";
		ChemGraph primary = getChemGraph();
        s = s + primary.toStringWithoutH(i);

        int index=0;
        /*for (Iterator iter = resonanceIsomers.iterator(); iter.hasNext(); ) {
        	index++;
        	ChemGraph isomer = (ChemGraph)iter.next();
        	if (isomer!=primary) {
        		s = s + '\n';
        		s = s + "isomer" + String.valueOf(index) + ":\n";
        		s = s + isomer.toStringWithoutH(i);
        	}
        }*/
        return s;




        //#]
    }
	
	// Convert ChemGraph (p_chemGraph) to String (InChIstring)
        // 6/9/09 gmagoon: updated to also read InChIKey
        //output is string array with first ([0]) element being the InChI string (with "InChI="); the second ([1]) element will be the InChIKey (WITHOUT "InChIKey=")
        // the function can be easily changed to accomodate other versions of the InChI program by changing the command line options (optionsArgument); for example, if InChIKey generation is turned off, the second element of the string array will be an empty string;
	public static String [] generateInChI(ChemGraph p_chemGraph) {
		File inchi = new File("InChI");
		inchi.mkdir();
		
		String [] result = new String[2];
		String InChIstring = "";
                String InChIKeystring = "";
                String line = "";
                String cTable = generateMolFileString(p_chemGraph);
		
		//
                File molFile = null;
		String workingDirectory = System.getProperty("RMG.workingDirectory");
		String inchiDirectory = "InChI";

		// Write the cTable to species.mol file
        try {
        	molFile = new File(inchiDirectory + "/species.mol");
        	FileWriter fw = new FileWriter(molFile);
        	fw.write(cTable);
        	fw.close();
        } catch (IOException e) {
        	String err = "Error writing species.mol file for InChI generation: ";
        	err += e.toString();
        	System.out.println(err);
        }
        
        // Call cINChI-1 executable file
//        String optionsArgument = "";
//        if (getOs().toLowerCase().contains("windows"))
//        	optionsArgument = "/DoNotAddH /FixedH /Key";//6/9/09 gmagoon: added fixed H so tautomers are considered separately; this is apparently not an option for version 1.02 (standard inchi); also added Key option to generate inchikey...this is only available in version 1.02beta or later; therefore, this keyword combination will only work for version 1.02beta (as of now)
//        else if (getOs().toLowerCase().contains("linux"))
//        	optionsArgument = "-DoNotAddH -FixedH -Key";
//        else if (getOs().toLowerCase().contains("mac"))
//        	optionsArgument = "-DoNotAddH -FixedH -Key";
        //6/9/09 gmagoon: I have reorganized the function in order to get around difficulty when using multiple options...if I just try to pass all options as a string (e.g. "-DoNotAddH -FixedH -Key") it is not correctly interpreted by the InChI program (note that leading "/" is not shown in the resulting "unrecognized option" warning);
        //I would have liked to have specified optionsArgument as an array of strings, but I was getting "illegal start of expression" errors in NetBeans IDE
        int exitValue = -1;
        //while (exitValue != 0) {
        
        /*
         * MRH 9MAR2010:
         * Addressing Issue#13: InChI not working consistently on Linux / Mac
         * 
         * Switching how the input/error/output of the "InChI" Process is handled.
         * MRH has copied the format JWA has used for the fame executable.
         * InChI now appears to run very consistently on Linux (in the limited types
         * of jobs I have ran - RMG, AdjList2InChI), in addition to running 
         * consistently on Windows.
         * 
         * Hopefully RW can test the new implementation on Mac
         * 
         * NOTE: For Windows, the InChI.waitFor() must follow the .close()
         * 	methods for the stdout and stderr BufferedReaders.  Otherwise,
         * 	the RMG simulation gets hung up at the .waitFor() method and
         * 	is perfecting happy waiting ... causing an infinite runtime.
         * 
         * NOTE: For Linux, 
         */
            try {
                if (getOs().toLowerCase().contains("windows")){
                    String[] command = {workingDirectory + "/bin/cInChI-1",
                        "species.mol",
                        "species.txt",
                        "/DoNotAddH", "/FixedH", "/Key"};//6/9/09 gmagoon: added fixed H so tautomers are considered separately; this is apparently not an option for version 1.02 (standard inchi); also added Key option to generate inchikey...this is only available in version 1.02beta or later; therefore, this keyword combination will only work for version 1.02beta (as of now)
                    File runningDir = new File("InChI");
                    Process InChI = Runtime.getRuntime().exec(command, null, runningDir);

                    BufferedReader stdout = new BufferedReader(new InputStreamReader(InChI.getInputStream()));
                    BufferedReader stderr = new BufferedReader(new InputStreamReader(InChI.getErrorStream()));

        			// Clean up i/o streams
        			stdout.close();
        			stderr.close();
        			
                    exitValue = InChI.waitFor();                        
                }
                else if (getOs().toLowerCase().contains("linux")){
                    String[] command = {workingDirectory + "/bin/cInChI-1",
                        "species.mol",
                        "species.txt",
                        "-DoNotAddH", "-FixedH", "-Key"};
                    File runningDir = new File("InChI");
                    Process InChI = Runtime.getRuntime().exec(command, null, runningDir);
                    
                    BufferedReader stdout = new BufferedReader(new InputStreamReader(InChI.getInputStream()));
                    BufferedReader stderr = new BufferedReader(new InputStreamReader(InChI.getErrorStream()));

                    exitValue = InChI.waitFor();

        			// Clean up i/o streams
        			stdout.close();
        			stderr.close();
                    
                }
                else if (getOs().toLowerCase().contains("mac")){
                    String[] command = {workingDirectory + "/bin/cInChI-1",
                        "species.mol",
                        "species.txt",
                        "-DoNotAddH", "-FixedH", "-Key"};
                    File runningDir = new File("InChI");
                    Process InChI = Runtime.getRuntime().exec(command, null, runningDir);
                    
                    BufferedReader stdout = new BufferedReader(new InputStreamReader(InChI.getInputStream()));
                    BufferedReader stderr = new BufferedReader(new InputStreamReader(InChI.getErrorStream()));
        			
                    exitValue = InChI.waitFor();

        			// Clean up i/o streams
        			stdout.close();
        			stderr.close();
                }
            }
            catch (Exception e) {
                String err = "Error running cINChI-1: ";
                err += e.toString();
                System.out.println(err);
            }
        //}
		
		// Read in the output of the cINChI-1 executable file (species.txt)
            /*
             * MRH 9MAR2010:
             * This portion of code (reading in the species.txt file and searching
             * 	for the InChI and InChIKey) is now obsolete (see code above).
             */
        FileReader in = null;
		try {
			in = new FileReader(inchiDirectory + "/species.txt");
		} catch (FileNotFoundException e) {
			String err = "Error reading species.txt file in generating InChI for species " + p_chemGraph.chemicalFormula + " : ";
			err += e.toString();
			System.out.println(err);
		}
        
		if (in != null) {
			BufferedReader reader = new BufferedReader(in);
	        line = ChemParser.readMeaningfulLine(reader);
	        read: while (line != null) {
	        	if (line.startsWith("InChI=")) {//changed from InChI to InChI= (to distinguish fro InChIKey
	        		InChIstring = line;
	        	}
	                else if (line.startsWith("InChIKey=")) {//changed from "InChI" to "InChI=" (to distinguish from "InChIKey="
	        		InChIKeystring = line.replace("InChIKey=", "");//read in the InChIKey without the preceding "InChIKey="
	        		break;
	        	}
	        	line = ChemParser.readMeaningfulLine(reader);
	            }
	        result[0]=InChIstring;
	        result[1]=InChIKeystring;
		}
		else {
			result[0] = "";
			result[1] = "";
		}
        return result;
	}
	
        //convert a chemgraph into a string that represents a 2D molefile with all atom positions initialized to zero
        //gmagoon 6/2/09: I separated this out from generateInChI so it could be easily be used elsewhere
        public static String generateMolFileString(ChemGraph p_chemGraph) {
    		// Convert chemGraph to string
		int randNum = 1;
		String p_string = p_chemGraph.toString(randNum);
		StringTokenizer st = new StringTokenizer(p_string);
		
		// Extract the necessary information from the adjacency list
		//		- Element symbol
		//		- Radical number
		//		- Bonds (strength and connectivity)
		
		// Define two running counters, atomCount & atomBondCount
		int atomCount = 1;
		int atomBondCount = 0;
		
		// Assume molecule has <= 100 atoms
		// Assume no atom can have > 6 unique bonds
		int maxAtoms = 100;
		int maxBonds = 6;
		String[] elementSymbol = new String[maxAtoms];
		String[] radical = new String[maxAtoms];
		String[][] atomConnect = new String[maxAtoms][maxBonds];
		String[][] bondOrder = new String[maxAtoms][maxBonds];
		
		// Atom-by-atom, extract the necessary information
		String line = st.nextToken();
		while (st.hasMoreTokens()) {
			// If the next token is the next integer in the series 
			if (line.equals(Integer.toString(atomCount))) {
				// Grab element symbol and radical number
				elementSymbol[atomCount-1] = st.nextToken();
				radical[atomCount-1] = st.nextToken();
				// If bonds exist ....
				if (st.hasMoreTokens()) {
					// Grab the first {bond}
					line = st.nextToken();
					int atomPartnerCount = 0;
					while (line.endsWith("}")) {
						++atomBondCount;
						String insideBraces = ChemParser.removeBrace(line);
						int commaPos = insideBraces.indexOf(",");
						// Store bond order and atom connectivity
						atomConnect[atomCount-1][atomPartnerCount] = insideBraces.substring(0,commaPos);
						bondOrder[atomCount-1][atomPartnerCount] = insideBraces.substring(commaPos+1);
						if (st.hasMoreTokens())
							line = st.nextToken();
						else
							line = "";
						++atomPartnerCount; 
					}
				} else
					line = "";
				++atomCount;
			} else
				line = st.nextToken();
		}
		
		int mRad = 1;
		int radCount = 0;
		// Assume the species does not contain > 4 radicals
		int maxRad = 4;
		int[] radLocation = new int[maxRad];
		int[] radType = new int[maxRad];
		for (int numRads=0; numRads<(atomCount-1); numRads++) {
			if (!radical[numRads].equals("0") & !radical[numRads].equals("null")) {
				// If a radical exists on the atom, store its location
				radLocation[radCount] = numRads+1;
				// Convert RMG's radical number definition to a .mol file's
				// 		radical number definition
				if (radical[numRads].equals("2S"))
					radType[radCount] = 1;
				else if (radical[numRads].equals("1") | radical[numRads].equals("3"))
					radType[radCount] = 2;
				else if (radical[numRads].equals("2") | radical[numRads].equals("2T"))
					radType[radCount] = 3;
				++radCount;
				// If at least one radical exist, .mol file will have 2 "M" lines
				mRad = 2;
			}
		}
		
		// Convert the information into a connection table (as defined by a .mol file)
		
		// Create the "Header Block"
		String cTable = "\n\n\n";
		// Create the "Counts Line"
		if (atomCount-1 < 10) {
			cTable += "  " + (atomCount-1);
		} else
			cTable += " " + (atomCount-1);
		if (atomBondCount/2 < 10) {
			cTable += "  " + (atomBondCount/2) + "  0  0  0  0  0  0  0  0  " + mRad + " V2000";
		} else 
			cTable += " " + (atomBondCount/2) + "  0  0  0  0  0  0  0  0  " + mRad + " V2000";
		// Create the "Atom Block"
		for (int numAtoms=0; numAtoms<atomCount-1; numAtoms++) {
			// Assume no 3-d information available
			//		Set all atoms x,y,z-coordinates at (0.0,0.0,0.0)
			cTable += "\n    0.0       0.0       0.0    " + elementSymbol[numAtoms] + "  " + " 0  0  0  0  0  0  0  0  0  0  0  0";
		}
		
		// Create the "Bond Block"
		// Convert Chemical Bond from S, D, T, etc. to 1, 2, 3, etc.
		int[] bondStrength = new int[atomBondCount/2];
		int bondCount = 0;
		for (int i=0; i<maxAtoms; i++) {	// Assuming <= 50 atoms in a species
			for (int j=0; j<maxBonds; j++) {	// Assuming no atom has > 6 bonds
				
				if (atomConnect[i][j] != null) {
					int secondAtom = Integer.parseInt(atomConnect[i][j]);
					
					if (secondAtom > i) {	// Do not want to double count bonds
						// Convert Chemical Bond from S, D, T, etc. to 1, 2, 3, etc.
						if (bondOrder[i][j].equals("S")) {
							bondStrength[bondCount] = 1;
						} else if (bondOrder[i][j].equals("D")) {
							bondStrength[bondCount] = 2;
						} else if (bondOrder[i][j].equals("T")) {
							bondStrength[bondCount] = 3;
						} else if (bondOrder[i][j].equals("B")) {
							bondStrength[bondCount] = 4;
						}
						
						cTable += "\n";
						
						if (i+1 < 10 && secondAtom < 10) {
							cTable += "  " + (i+1) + "  " + secondAtom + "  " + bondStrength[bondCount] + "  0  0  0  0";
						} else if (i+1 < 10 && secondAtom >= 10) {
							cTable += "  " + (i+1) + " " + secondAtom + "  " + bondStrength[bondCount] + "  0  0  0  0";
						} else if (i+1 >= 10 && secondAtom < 10) {
							cTable += " " + (i+1) + "  " + secondAtom + "  " + bondStrength[bondCount] + "  0  0  0  0";
						} else 
							cTable += " " + (i+1) + " " + secondAtom + "  " + bondStrength[bondCount] + "  0  0  0  0";
						
						++bondCount;
					}
				
				} else	// There are no more bonds for this atom
					break;
			}
		}
		
		// Create the "Properties Block"
		for (int i=0; i<mRad-1; i++) {
			cTable += "\nM  RAD  " + radCount;
			for (int j=0; j<radCount; j++) {
				if (radLocation[j] < 10) {
					cTable += "   " + radLocation[j] + "   " + radType[j]; 
				} else
					cTable += "  " + radLocation[j] + "   " + radType[j];
			}
		}
		cTable += "\nM  END";
                
                return cTable;
    }
    
    /**
     * Converts a single InChI to its RMG adjacency lists
     * @param p_inchi: String containing the InChI
     * @return
     */
	public static String inchi2AdjList (String p_inchi) {
		// Convert InChI to .mol file
		inchi2mol(p_inchi);
		// Convert .mol file to adjacency list
        String adjList = mol2AdjList(p_inchi,"InChI/temp.mol");
        return adjList;
	}
	
	/**
	 * Converts a single InChI to its .mol file format.  The .mol file
	 * 	is saved to the $RMG/InChI directory
	 * @param p_inchi: String containing the InChI
	 */
	public static void inchi2mol(String p_inchi) {
		String workingDirectory = System.getenv("RMG");
		File inchi = new File("InChI");
		inchi.mkdir();
		String inchiDirectory = "InChI";
		File inchiFile = null;
		
		// Save InChI string in inchi.txt file
        try {
        	inchiFile = new File(inchiDirectory + "/inchi.txt");
        	FileWriter fw = new FileWriter(inchiFile);
        	/*
        	 * The following change is suggested by G.Magoon at:
        	 * http://github.com/gmagoon/RMG-Java/commit/f652decc547d0a5107435d9bc4ca64d11a9bbe7c
        	 * This fix allows the conversion of InChI --> .mol for InChI v. 1.02beta
        	 * MRH will check if all other InChI executable operations are still functional,
        	 * 	and for which versions
        	 * 
        	 * Windows: InChI --> .mol and .mol --> InChI work for InChI v. 1.01 &
        	 * 	v. 1.02beta
        	 * Linux: .mol --> InChI works for both v. 1.01 and 1.02beta;
        	 * 	InChI --> .mol works only for 1.01 (1.02beta states:
        	 * 		"Option InChI2Struct is currently not supported (use classic interface)
        	 */
        	//fw.write(p_inchi);
        	fw.write(p_inchi+"\n");
        	fw.close();
        } catch (IOException e) {
        	String err = "Error writing inchi.txt file for InChI-to-molFile conversion: ";
        	err += e.toString();
        	System.out.println(err);
        }
		
        // Call cINChI-1 executable file
        /*
         * MRH 11MAR2010:
         * Making the InChI executable work on Linux
         * 
         * NOTE: The Windows and Linux/Mac codes are slightly different (the position of
         * 	the InChI.waitFor() method.  Getting InChI to work on Windows & Linux is very
         * 	sensitive to the position of this line.
         * 
         * Still need to test on Mac; MRH assumes it will follow Linux.
         */
        String[] optionsArgument = new String[2];
        int exitValue = -1;
        if (getOs().toLowerCase().equals("windows")) {
        	optionsArgument[0] = "/InChI2Struct";
        	optionsArgument[1] = "/OutputSDF";
            String[] command1 = {workingDirectory + "/bin/cInChI-1",
                    "inchi.txt",
                    "temp.txt",
                    optionsArgument[0]};
            String[] command2 = {workingDirectory + "/bin/cInChI-1",
            		"temp.txt",
            		"temp.mol",
            		optionsArgument[1]};
            try {
	            File runningDir = new File("InChI");
	            Process InChI = Runtime.getRuntime().exec(command1, null, runningDir);
	            
	            BufferedReader stdout = new BufferedReader(new InputStreamReader(InChI.getInputStream()));
	            BufferedReader stderr = new BufferedReader(new InputStreamReader(InChI.getErrorStream()));
	            // Clean up i/o streams
	            stdout.close();
	            stderr.close();  			
	            exitValue = InChI.waitFor();
	            
	            InChI = Runtime.getRuntime().exec(command2, null, runningDir);
	            
	            stdout = new BufferedReader(new InputStreamReader(InChI.getInputStream()));
	            stderr = new BufferedReader(new InputStreamReader(InChI.getErrorStream()));
	            // Clean up i/o streams
	            stdout.close();
	            stderr.close();  			
	            exitValue = InChI.waitFor();	            
            }
            catch (Exception e) {
                String err = "Error running cInChI-1 while converting InChI to .mol file: ";
                err += e.toString();
                System.out.println(err);
            }
        } else if (getOs().toLowerCase().equals("linux")) {
        	optionsArgument[0] = "-InChI2Struct";
        	optionsArgument[1] = "-OutputSDF";
            String[] command1 = {workingDirectory + "/bin/cInChI-1",
                    "inchi.txt",
                    "temp.txt",
                    optionsArgument[0]};
            String[] command2 = {workingDirectory + "/bin/cInChI-1",
            		"temp.txt",
            		"temp.mol",
            		optionsArgument[1]};
            try {
	            File runningDir = new File("InChI");
	            Process InChI = Runtime.getRuntime().exec(command1, null, runningDir);
	            
	            BufferedReader stdout = new BufferedReader(new InputStreamReader(InChI.getInputStream()));
	            BufferedReader stderr = new BufferedReader(new InputStreamReader(InChI.getErrorStream()));
	            exitValue = InChI.waitFor();
	            // Clean up i/o streams
	            stdout.close();
	            stderr.close();  			
	            
	            InChI = Runtime.getRuntime().exec(command2, null, runningDir);
	            
	            stdout = new BufferedReader(new InputStreamReader(InChI.getInputStream()));
	            stderr = new BufferedReader(new InputStreamReader(InChI.getErrorStream()));
	            exitValue = InChI.waitFor();	
	            // Clean up i/o streams
	            stdout.close();
	            stderr.close();  		            
            }
            catch (Exception e) {
                String err = "Error running cInChI-1 while converting InChI to .mol file: ";
                err += e.toString();
                System.out.println(err);
            }
        } else if (getOs().toLowerCase().equals("mac")) {
        	optionsArgument[0] = "-InChI2Struct";
        	optionsArgument[1] = "-OutputSDF";
            String[] command1 = {workingDirectory + "/bin/cInChI-1",
                    "inchi.txt",
                    "temp.txt",
                    optionsArgument[0]};
            String[] command2 = {workingDirectory + "/bin/cInChI-1",
            		"temp.txt",
            		"temp.mol",
            		optionsArgument[1]};
            try {
	            File runningDir = new File("InChI");
	            Process InChI = Runtime.getRuntime().exec(command1, null, runningDir);
	            
	            BufferedReader stdout = new BufferedReader(new InputStreamReader(InChI.getInputStream()));
	            BufferedReader stderr = new BufferedReader(new InputStreamReader(InChI.getErrorStream()));
	            exitValue = InChI.waitFor();
	            // Clean up i/o streams
	            stdout.close();
	            stderr.close();  			
	            
	            InChI = Runtime.getRuntime().exec(command2, null, runningDir);
	            
	            stdout = new BufferedReader(new InputStreamReader(InChI.getInputStream()));
	            stderr = new BufferedReader(new InputStreamReader(InChI.getErrorStream()));
	            exitValue = InChI.waitFor();	
	            // Clean up i/o streams
	            stdout.close();
	            stderr.close();  		            
            }
            catch (Exception e) {
                String err = "Error running cInChI-1 while converting InChI to .mol file: ";
                err += e.toString();
                System.out.println(err);
            }
        }
	}
	
	/**
	 * Convert a .mol file to a RMG adjacency list
	 * @param inchi: String containing the InChI
	 * @param filePath: Location of the .mol file
	 * @return
	 */
	public static String mol2AdjList (String inchi, String filePath) {
		// Read in the .mol file
        FileReader in = null;
		try {
			in = new FileReader(filePath);
		} catch (FileNotFoundException e) {
			String err = "Error reading .mol file: " + e.toString();
			System.out.println(err);
		}
        
		BufferedReader reader = new BufferedReader(in);
		String line = ChemParser.readMeaningfulLine(reader);
		
		while (!line.toUpperCase().endsWith("V2000")) {
			line = ChemParser.readMeaningfulLine(reader);
		}
		
		String molFile = "";
		while (line != null) {
			molFile += line + "\r";
			line = ChemParser.readMeaningfulLine(reader);
		}
		
		// Determine how many lines are in the .mol file
		String[] molFileLines = molFile.split("[\r]",0);
		int numOfLines = molFileLines.length;
		StringTokenizer st = new StringTokenizer(molFileLines[numOfLines-1]);
		if (st.nextToken().equals("$$$$")) --numOfLines;
		
		// Extract the information in the first line (Count Line) of the .mol file
		st = new StringTokenizer(molFileLines[0]);
		int numOfAtoms = Integer.parseInt(st.nextToken());
		int numOfBonds = Integer.parseInt(st.nextToken());
		// Next few are irrelevant for RMG (as of 10-Feb-2009)
		int numOfAtomLists = Integer.parseInt(st.nextToken());
		String obsoleteString1 = st.nextToken();
		String chiralFlag = st.nextToken();
		int stextEntries = Integer.parseInt(st.nextToken());
		String obsoleteString2 = st.nextToken();
		String obsoleteString3 = st.nextToken();
		String obsoleteString4 = st.nextToken();
		String obsoleteString5 = st.nextToken();
		// Extract the number of M lines		
		int numOfMLines = Integer.parseInt(st.nextToken());
		
		// Construct each individual line of the adjacency list
		String[] adjListElement = new String[numOfAtoms];
		String[] adjListRadical = new String[numOfAtoms];
		String[] adjListConnectivity = new String[numOfAtoms];
		
		for (int i=0; i<numOfAtoms; i++) {
			adjListConnectivity[i] = "";
			adjListRadical[i] = "0 ";
		}
		
		//	Extract the element symbol
		for (int i=1; i<numOfAtoms+1; i++) {
			st = new StringTokenizer(molFileLines[i]);
			// These 3-d geometries may be helpful in the future
			double x_coord = Double.parseDouble(st.nextToken());
			double y_coord = Double.parseDouble(st.nextToken());
			double z_coord = Double.parseDouble(st.nextToken());
			// Extract the element symbol
			adjListElement[i-1] = " " + st.nextToken() + " ";
			
			/*
			 * Added by MRH on 4-Sept-2009
			 * 	RMG did not convert single heavy atom radicals properly (e.g. CH3, OH, H).  The
			 * 		problem was the structure of the .mol file.  I created mol2AdjList assuming
			 * 		the radical information would be present in the "M RAD" lines located at the
			 * 		bottom of the .mol file.  For InChI-v.1.01, this is not the case for single
			 * 		heavy atom radicals.  For the time being, I am hardcoding in the fix.
			 */
			if (numOfAtoms == 1) {
				String totalLine = molFileLines[i];
				int length = totalLine.length();
				// This length-15,length-12 may depend on the version of InChI
				int valencyFlag = Integer.parseInt(totalLine.substring(length-15,length-12).trim());
				if (valencyFlag != 0) {
					// One Hydrogen
					if (valencyFlag == 1) {
						if (adjListElement[i-1].equals(" C ")) adjListRadical[i-1] = "3 ";
						else if (adjListElement[i-1].equals(" Si ")) adjListRadical[i-1] = "3 ";
						else if (adjListElement[i-1].equals(" N ")) adjListRadical[i-1] = "2 ";
						else if (adjListElement[i-1].equals(" O ")) adjListRadical[i-1] = "1 ";
						else if (adjListElement[i-1].equals(" S ")) adjListRadical[i-1] = "1 ";
					// Two Hydrogens
					} else if (valencyFlag == 2) {
						if (adjListElement[i-1].equals(" C ")) adjListRadical[i-1] = "2 ";
						else if (adjListElement[i-1].equals(" Si ")) adjListRadical[i-1] = "2 ";
						else if (adjListElement[i-1].equals(" N ")) adjListRadical[i-1] = "1 ";
					// Three Hydrogens
					} else if (valencyFlag == 3) {
						if (adjListElement[i-1].equals(" C ")) adjListRadical[i-1] = "1 ";
						else if (adjListElement[i-1].equals(" Si ")) adjListRadical[i-1] = "1 ";
					// Zero Hydrogens
					} else if (valencyFlag == 15) {
						if (adjListElement[i-1].equals(" C ")) adjListRadical[i-1] = "4 ";
						else if (adjListElement[i-1].equals(" Si ")) adjListRadical[i-1] = "4 ";
						else if (adjListElement[i-1].equals(" N ")) adjListRadical[i-1] = "3 ";
						else if (adjListElement[i-1].equals(" O ")) adjListRadical[i-1] = "2 ";
						else if (adjListElement[i-1].equals(" S ")) adjListRadical[i-1] = "2 ";
						else if (adjListElement[i-1].equals(" H ")) adjListRadical[i-1] = "1 ";
					}
				}
			}
		}
		
		//	Extract the connectivity
		int Counter = numOfAtoms+1;
		while (!molFileLines[Counter].startsWith("M")) {
			st = new StringTokenizer(molFileLines[Counter]);
			// Extract the two atoms associated with this connection
			int atom1 = Integer.parseInt(st.nextToken());
			int atom2 = Integer.parseInt(st.nextToken());
			// Place the connected atoms in the braces
			adjListConnectivity[atom1-1] += "{" + atom2;
			adjListConnectivity[atom2-1] += "{" + atom1;
			// Extract the type of connection
			int connection = Integer.parseInt(st.nextToken());
			// Convert the connection and place in braces
			if (connection == 1) {
				adjListConnectivity[atom1-1] += ",S} ";
				adjListConnectivity[atom2-1] += ",S} ";
			} else if (connection == 2) {
				adjListConnectivity[atom1-1] += ",D} ";
				adjListConnectivity[atom2-1] += ",D} ";
			} else if (connection == 3) {
				adjListConnectivity[atom1-1] += ",T} ";
				adjListConnectivity[atom2-1] += ",T} ";
			} else if (connection == 4) {
				adjListConnectivity[atom1-1] += ",B} ";
				adjListConnectivity[atom2-1] += ",B} ";
			}
			++Counter;
		}
		
		// Determine the position and type of the radicals
		for (int i=numOfLines-numOfMLines; i<numOfLines-1; i++) {
			st = new StringTokenizer(molFileLines[numOfLines-numOfMLines]);
			// The following variables hold no meaning
			String M = st.nextToken();
			String RAD = st.nextToken();
                        if (RAD.equals("RAD")){
			//if (RAD.equals("RAD")||RAD.equals("CHG")) {
				// Extract radical information
				int numOfRads = Integer.parseInt(st.nextToken());
				for (int j=0; j<numOfRads; j++) {
					int atom = Integer.parseInt(st.nextToken());
					int radType = Integer.parseInt(st.nextToken());
                                       // if(RAD.equals("RAD")){
                                            if (radType == 1)
                                                    adjListRadical[atom-1] = "2S ";
                                            else if (radType == 2)
                                                    adjListRadical[atom-1] = "1 ";
                                            else if (radType== 3)
                                                    adjListRadical[atom-1] = "2T ";
                                            else
                                                    adjListRadical[atom-1] = "3 ";
                                      //  }
                                       // else if(RAD.equals("CHG")){//gmagoon 9/14/09: adding support for CHG, which is produced in InChI generated mol files for biradicals; there still could be issues for biradicals with radical sites on adjacent atoms, so this should be manually checked
                                       //      if (radType == 1 || radType == -1){
                                      //          System.out.println("Assuming CHG indicates biradical for " + inchi +  " If this is correct, two of these messages for the same InChI should be generated.");
                                      //           adjListRadical[atom-1] = "1 ";
                                      //       }
                                      //       else {
                                      //          System.out.println("Ignoring unknown M flag " + RAD + " for " + inchi);
                                      //       }  
                                      //  }
                                        
				}
                                
			}
			/*
			 * If the M line flag is not equal to "RAD", RMG cannot store the information
			 * 	(e.g. CHG for InChI=1/C20H40O2/c1-3-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19(4-2)20(21)22/h19H,3-18H2,1-2H3,(H,21,22)/p-1/t19-/m1/s1/fC20H39O2/q-1)
			 * 	Inform the user the field is being ignored, but continue constructing
			 * 	the adjacency list.
			 */
                        //update gmagoon 9/14/09: CHG is now allowed: see above; update 2: I returned to original case due to difficulties with handling adjacent biradicals
			else {
				System.out.println("Ignoring unknown M flag " + RAD + " for " + inchi);
			}      
		}
		
		// Construct the entire adjacency list from its individual lines
		String adjList = "";
		for (int i=0; i<numOfAtoms-1; i++) {
			adjList += (i+1) + adjListElement[i] + adjListRadical[i] + adjListConnectivity[i] + "\r";
		}
		adjList += numOfAtoms + adjListElement[numOfAtoms-1] + adjListRadical[numOfAtoms-1] + adjListConnectivity[numOfAtoms-1];
		
		return adjList;
	}

    public int getID() {
        return ID;
    }

    public static int getTOTAL_NUMBER() {
        return TOTAL_NUMBER;
    }

    public ChemGraph getChemGraph() {
        return chemGraph;
    }

    public void setChemGraph(ChemGraph p_chemGraph) {
        chemGraph = p_chemGraph;
    }

    public double getDeltaEDown() {
        return deltaEDown;
    }

    public String getName() {
        return name;
    }

    public void setName(String p_name) {
        name = p_name;
    }

    public boolean getTherfitExecuted() {
        return therfitExecuted;
    }

    public TransportData getChemkinTransportData() {
        return chemkinTransData;
    }

    public TransportData newChemkinTransportData() {
        chemkinTransData = new TransportData();
        return chemkinTransData;
    }
    
    public ThreeFrequencyModel getThreeFrequencyModel() {
        return threeFrequencyModel;
    }

    public SpectroscopicData getSpectroscopicData() {
		// Generate data if needed
		if (spectroscopicData == null && threeFrequencyModel == null)
			generateSpectroscopicData();
		// Return data in appropriate form
		if (spectroscopicData != null)
			return spectroscopicData;
		else if (threeFrequencyModel != null)
			return new SpectroscopicData(threeFrequencyModel);
		else if (isTriatomicOrSmaller())
			return new SpectroscopicData();
		else return null;
    }

	protected void initRelations() {
        chemkinTransData = newChemkinTransportData();
    }

    public boolean isReactive() {
    	return IsReactive;
    }	
    
    public void setReactivity(boolean reactive) {
    	IsReactive = reactive;
    }
    
    public boolean isConstantConcentration() {
        return constantConcentration;
    }
    
    public void setConstantConcentration(boolean hasconstantconcentration){
        constantConcentration = hasconstantconcentration;
    }
    
	public void addPdepPaths(HashSet pdepReactionSet) {
		
		if (paths == null)
			paths = pdepReactionSet;
		else
			paths.addAll(pdepReactionSet);
	}
	
	public HashSet getPdepPaths(){
		return paths;
	}
	public static void setAddID(boolean p_addID){
		addID = p_addID;
	}
	
	/**
	 * Checks to see if the species is an isomer of another species by comparing
	 * the numbers of each atom type in each species.
	 * @param species The species to check the current one against
	 * @return true if they are isomers, false otherwise
	 */
	public boolean isIsomerOf(Species species) {
//		ChemGraph cg1 = getChemGraph();
//		ChemGraph cg2 = species.getChemGraph();
		boolean areIsomers = true;
//		if (cg1.getCarbonNumber() != cg2.getCarbonNumber())
//			areIsomers = false;
//		if (cg1.getOxygenNumber() != cg2.getOxygenNumber())
//			areIsomers = false;
//		if (cg1.getHydrogenNumber() != cg2.getHydrogenNumber())
//			areIsomers = false;
		String inchi1 = getInChI();
		String inchi2 = species.getInChI();
		String[] inchi1Layers = inchi1.split("/");
		String[] inchi2Layers = inchi2.split("/");
		if (inchi1Layers[1] != inchi2Layers[1]) areIsomers = false;
		// Should add other atom types in the future!
		return areIsomers;
	}
	
    public String getInChI() {
    	if (InChI == null) {
    		String[] inchiANDinchikey = generateInChI(getChemGraph());
    		InChI = inchiANDinchikey[0];
    	}
        return InChI;
    }
    
    public void setInChI(String inchi) {
    	InChI = inchi;
    }
    
    public static String getOs() {
    	  String os = "";
    	  if (System.getProperty("os.name").toLowerCase().contains("windows")) {
    	    os = "windows";
    	  } else if (System.getProperty("os.name").toLowerCase().contains("linux")) {
    	    os = "linux";
    	  } else if (System.getProperty("os.name").toLowerCase().contains("mac")) {
    	    os = "mac";
    	  }
    	  return os;
    }
    
    public String getNasaThermoSource() {
    	if (nasaThermoSource == null)
    		nasaThermoSource = "Estimated by RMG using Group Additivity";
    	return nasaThermoSource;
    }
    
    public boolean equals(Species species) {
    	return this.getChemGraph().equals(species.getChemGraph());
    }
	
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\Species.java
*********************************************************************/


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
import jing.chemUtil.HierarchyTree;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\ThermoGAGroupLibrary.java
//----------------------------------------------------------------------------

/**
There are six libraries:
(1) group
(2) radical
(3) ring correction
(4) other correction
(5) gauche correction
(6) 1,5 correction
In each library, the key should be functional group (name + adjList), and the value should be ThermoGAValue.
for (2), (3), (4), we scan the library to find match between chemgraph and functional group each time. search time O(n), where n is the library size.
for (1), we first match chemgraph with a tree structure to find out the proper functional group, and then access the library by the key functional group, so the search tiem is O(1) + O(logN), where N is the tree size.
*/
//## class ThermoGAGroupLibrary
public class ThermoGAGroupLibrary {

    protected static ThermoGAGroupLibrary INSTANCE = new ThermoGAGroupLibrary();		//## attribute INSTANCE

    protected HashMap groupDictionary;		//## attribute groupDictionary

    protected HashMap groupLibrary;		//## attribute groupLibrary

    /**
    Note: this kind of tree is different with the kinetics tree.  In kinetics tree, tree nodes are FunctionalGroup or FunctionalGroupCollection.  In thermo tree, tree nodes are Nodes with connectivity,
    */
    protected HierarchyTree groupTree;		//## attribute groupTree
    protected HashMap otherDictionary;		//## attribute otherDictionary
    protected HashMap otherLibrary;		//## attribute otherLibrary
    protected HierarchyTree otherTree;		//## attribute otherTree
    protected HashMap radicalDictionary;		//## attribute radicalDictionary
    protected HashMap radicalLibrary;		//## attribute radicalLibrary
    protected HierarchyTree radicalTree;		//## attribute radicalTree
    protected HierarchyTree ringTree;
    protected HashMap ringDictionary;
    protected HashMap ringLibrary;		//## attribute ringLibrary    
    protected HashMap gaucheDictionary;		
    protected HashMap gaucheLibrary;		
    protected HierarchyTree gaucheTree;	
    protected HashMap oneFiveDictionary;		
    protected HashMap oneFiveLibrary;		
    protected HierarchyTree oneFiveTree;	
    protected HashMap abramDictionary;
    protected HashMap abramLibrary;
    protected HierarchyTree abramTree;
    protected HashMap unifacDictionary;
    protected HashMap unifacLibrary;
    protected HierarchyTree unifacTree;

    //protected HashMap solventDictionary;
    //protected HashMap solventLibrary;
    // Constructors

    //## operation ThermoGAGroupLibrary()
    private  ThermoGAGroupLibrary() {
        groupTree = new HierarchyTree();
        groupDictionary = new HashMap();
        groupLibrary = new HashMap();

        radicalTree = new HierarchyTree();
        radicalDictionary = new HashMap();
        radicalLibrary = new HashMap();

        ringLibrary = new HashMap();
        ringDictionary = new HashMap();
        ringTree = new HierarchyTree();

        otherLibrary = new HashMap();
        otherDictionary = new HashMap();
        otherTree = new HierarchyTree();
        
        gaucheLibrary = new HashMap();
        gaucheDictionary = new HashMap();
        gaucheTree = new HierarchyTree();
        
        oneFiveLibrary = new HashMap();
        oneFiveDictionary = new HashMap();
        oneFiveTree = new HierarchyTree();

        abramLibrary= new HashMap();
        abramDictionary=new HashMap();
        abramTree=new HierarchyTree();

        unifacLibrary= new HashMap();
        unifacDictionary=new HashMap();
        unifacTree=new HierarchyTree();

      //  solventDictionary=new HashMap();
       // solventLibrary=new HashMap();

        String directory = System.getProperty("jing.chem.ThermoGAGroupLibrary.pathName");
        if (directory == null) {
        	System.out.println("undefined system property: jing.chem.ThermoGAGroupLibrary.pathName, exit!");
        	System.exit(0);
        }

        String separator = System.getProperty("file.separator");
        if (!directory.endsWith(separator)) directory = directory + separator;
			
			System.out.println("\nReading thermo database from "+directory);

        String gDictionary = directory + "Group_Dictionary.txt";
        String gTree = directory + "Group_Tree.txt";
        String gLibrary = directory + "Group_Library.txt";
        String rDictionary = directory + "Radical_Dictionary.txt";
        String rTree = directory + "Radical_Tree.txt";
        String rLibrary = directory + "Radical_Library.txt";
        String ringDictionary = directory + "Ring_Dictionary.txt";
        String ringTree = directory + "Ring_Tree.txt";
        String ringLibrary = directory + "Ring_Library.txt";
        String otherLibrary = directory + "Other_Library.txt";
        String otherDictionary = directory + "Other_Dictionary.txt";
        String otherTree = directory + "Other_Tree.txt";
        
        String gauDictionary = directory + "Gauche_Dictionary.txt";
        String gauTree = directory + "Gauche_Tree.txt";
        String gauLibrary = directory + "Gauche_Library.txt";
        
        String one5Dictionary = directory + "15_Dictionary.txt";
        String one5Tree = directory + "15_Tree.txt";
        String one5Library = directory + "15_Library.txt";
	
        String AbDictionary=directory+"Abraham_Dictionary.txt";
        String AbTree=directory+"Abraham_Tree.txt";
        String AbLibrary=directory+"Abraham_Library.txt";

        String UnDictionary=directory+"Unifac_Dictionary.txt";
        String UnTree=directory+"Unifac_Tree.txt";
        String UnLibrary=directory+"Unifac_Library.txt";

        //String solventdict=directory+"Solvent_Dictionary.txt";
        //String solventlib=directory+"Solvent_Library.txt";

        read(gDictionary,gTree,gLibrary,rDictionary,rTree,rLibrary,ringDictionary,ringTree,ringLibrary,otherDictionary,otherLibrary,otherTree,gauDictionary,gauTree,gauLibrary,one5Dictionary,one5Tree,one5Library,AbDictionary,AbTree,AbLibrary,UnDictionary,UnTree,UnLibrary);

    }

    //## operation findCorrectionInLibrary(ChemGraph,HashMap)
    private ThermoData findCorrectionInLibrary(ChemGraph p_chemGraph, HashMap p_library) {
        //#[ operation findCorrectionInLibrary(ChemGraph,HashMap)
        p_chemGraph.clearCentralNode();
        ThermoData result=new ThermoData();
        int redundance;

        Iterator iter = p_library.keySet().iterator();
        while (iter.hasNext()) {
        	redundance = 0;
        	FunctionalGroup f = (FunctionalGroup)iter.next();
        	HashSet gv = p_chemGraph.identifyThermoMatchedSite(f);
        	if (gv != null) {
        		redundance = gv.size();
        		if (redundance > 0) {
        			ThermoGAValue ga = (ThermoGAValue)p_library.get(f);
        			if (ga != null) {
        				ThermoData temp = new ThermoData(ga);
        				temp.multiply(redundance);
        				result.plus(temp);
        				temp = null;
        			}
        		}
        	}
        }
        p_chemGraph.getGraph().resetMatchedGC();
        return result;
        //#]
    }

    /**
    Requires: the central node of p_chemGraph has been set to the thermo center atom.
    Effects: find a matched thermo functional group in the group tree for the pass-in p_chemGraph, return this functional group's thermo value.  If no leaf is found, throw  GroupNotFoundException
    Modifies:
    */
    //## operation findGAGroup(ChemGraph)
    public ThermoGAValue findGAGroup(ChemGraph p_chemGraph) throws GroupNotFoundException, MultipleGroupFoundException, InvalidCenterTypeException {
        //#[ operation findGAGroup(ChemGraph)
        if (p_chemGraph == null) return null;

        Stack stack = groupTree.findMatchedPath(p_chemGraph);
        p_chemGraph.getGraph().resetMatchedGC();//2/13/09 gmagoon: resetting the matched GC value...for some reason, thermoLibrary.findGAGroup (within getGAGroup in GATP.java) ended up modifiying the central node so that it was matched; this ended up wreaking havoc with subsequent symmetry number calculations; ideally, I would probably want to fix the code so that it didn't end up modifying the matchedGC from the null value after it is done with it, but I do not immediately see how to due so, and debugging proved extremely difficult; I have also tried to put this elsewhere in this class where it might be appropriate
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
        	Matchable fg = (Matchable)node.getElement();
        	ThermoGAValue ga = (ThermoGAValue)groupLibrary.get(fg);
        	if (ga != null) //{
                //System.out.println("Group found: " + fg.getName());
                return ga;
            //}
        }
        
        return null;



        //#]
    }

    //## operation findOtherCorrection(ChemGraph)
    public ThermoGAValue findOtherCorrection(ChemGraph p_chemGraph) {
        //#[ operation findOtherCorrection(ChemGraph)
        if (p_chemGraph == null) return null;

        Stack stack = otherTree.findMatchedPath(p_chemGraph);
        p_chemGraph.getGraph().resetMatchedGC();
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
        	FunctionalGroup fg = (FunctionalGroup)node.getElement();
        	ThermoGAValue ga = (ThermoGAValue)otherLibrary.get(fg);
        	if (ga != null) return ga;
        }
     
        return null;




        //#]
    }

    //## operation findRadicalGroup(ChemGraph)
    public ThermoGAValue findRadicalGroup(ChemGraph p_chemGraph) throws InvalidThermoCenterException {
        //#[ operation findRadicalGroup(ChemGraph)
        if (p_chemGraph == null) return null;

        Stack stack = radicalTree.findMatchedPath(p_chemGraph);
        p_chemGraph.getGraph().resetMatchedGC();
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
        	Matchable fg = (Matchable)node.getElement();
        	ThermoGAValue ga = (ThermoGAValue)radicalLibrary.get(fg);
        	if (ga != null) return ga;
        }

        return null;




        //#]
    }

    //## operation findRingCorrection(ChemGraph)
	public ThermoGAValue findRingCorrection(ChemGraph p_chemGraph) {
	    // end pey
	        //#[ operation findRingCorrection(ChemGraph)

	        // begin pey -- read tree instead of library
	        // return findCorrectionInLibrary(p_chemGraph,ringLibrary);

	        if (p_chemGraph == null) return null;

	        // initialize
	        int deepest = -1;
	        Stack deepestStack = new Stack();
	        deepestStack = null;

	        // iterate through nodes in chemgraph that are in a cycle
	        Iterator iter = p_chemGraph.getNodeList();
	        while (iter.hasNext()) {
	          Node node = (Node) iter.next();
	          Atom atom = (Atom)node.getElement();
	          if (!(atom.getType().equals("H"))) {
	          // waiting on Jing to fix the inCycle issue for radicals that get saturated
	          // if (node.getInCycle()) {
	            // make the current node the central atom
	            p_chemGraph.resetThermoSite(node);
	            // find the match in the thermo tree
	            Stack stack = ringTree.findMatchedPath(p_chemGraph);
	            // check if it's the deepest match
	            if (!stack.empty()) {
	              HierarchyTreeNode htn = (HierarchyTreeNode) stack.peek();
	              if (htn.getDepth() > deepest) {
	                deepestStack = stack;
	                deepest = htn.getDepth();
	              }
	            }

	          }
	        }

	        if (deepestStack == null) return null;

	        while (!deepestStack.empty()) {
	                HierarchyTreeNode node = (HierarchyTreeNode)deepestStack.pop();
	                FunctionalGroup fg = (FunctionalGroup)node.getElement();
	                ThermoGAValue ga = (ThermoGAValue)ringLibrary.get(fg);
	                if (ga != null) return ga;
	        }
                p_chemGraph.getGraph().resetMatchedGC();
	        return null;
	        // end pey

	        //#]
	    }

    //2/5/09 gmagoon: new functions for gauche and 1,5-interactions
        /**
    Requires: the central node of p_chemGraph has been set to the thermo center atom.
    Effects: find a matched thermo functional group in the group tree for the pass-in p_chemGraph, return this functional group's thermo value.  If no leaf is found, throw  GroupNotFoundException
    Modifies:
    */
    public ThermoGAValue findGaucheGroup(ChemGraph p_chemGraph) throws MultipleGroupFoundException, InvalidCenterTypeException {
        //#[ operation findGAGroup(ChemGraph)
        if (p_chemGraph == null) return null;

        Stack stack = gaucheTree.findMatchedPath(p_chemGraph);
        p_chemGraph.getGraph().resetMatchedGC();
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
        	Matchable fg = (Matchable)node.getElement();
        	ThermoGAValue ga = (ThermoGAValue)gaucheLibrary.get(fg);
        	if (ga != null) return ga;
        }

        return null;
    }
    
            /**
    Requires: the central node of p_chemGraph has been set to the thermo center atom.
    Effects: find a matched thermo functional group in the group tree for the pass-in p_chemGraph, return this functional group's thermo value.  If no leaf is found, throw  GroupNotFoundException
    Modifies:
    */
    public ThermoGAValue find15Group(ChemGraph p_chemGraph) throws MultipleGroupFoundException, InvalidCenterTypeException {
        //#[ operation findGAGroup(ChemGraph)
        if (p_chemGraph == null) return null;

        Stack stack = oneFiveTree.findMatchedPath(p_chemGraph);
        p_chemGraph.getGraph().resetMatchedGC();
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
        	Matchable fg = (Matchable)node.getElement();
        	ThermoGAValue ga = (ThermoGAValue)oneFiveLibrary.get(fg);
        	if (ga != null) return ga;
        }

        return null;
    }

    public AbrahamGAValue findAbrahamGroup(ChemGraph p_chemGraph) throws MultipleGroupFoundException, InvalidCenterTypeException {
        //#[ operation findGAGroup(ChemGraph)
        if (p_chemGraph == null) return null;

        Stack stack = abramTree.findMatchedPath(p_chemGraph);
        p_chemGraph.getGraph().resetMatchedGC();
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
        	Matchable fg = (Matchable)node.getElement();
        	AbrahamGAValue ga = (AbrahamGAValue)abramLibrary.get(fg);
        	if (ga != null) {
            //System.out.println("Platts Group found: " + fg.getName());
                return ga;
            }
             
        }

        return null;
    }

        public UnifacGAValue findUnifacGroup(ChemGraph p_chemGraph) throws MultipleGroupFoundException, InvalidCenterTypeException {
        //#[ operation findGAGroup(ChemGraph)
        if (p_chemGraph == null) return null;

        Stack stack = unifacTree.findMatchedPath(p_chemGraph);
        p_chemGraph.getGraph().resetMatchedGC();
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
        	Matchable fg = (Matchable)node.getElement();
        	UnifacGAValue ga = (UnifacGAValue)unifacLibrary.get(fg);
        	if (ga != null) {
                //System.out.println("Unifac Group found: " + fg.getName());
                return ga;
            }

        }

        return null;
    }


    //## operation read(String,String,String,String,String,String,String,String,String)
	public void read(String p_groupDictionary, String p_groupTree, String p_groupLibrary, String p_radicalDictionary, String p_radicalTree, String p_radicalLibrary, String p_ringDictionary, String p_ringTree, String p_ringLibrary, String p_otherDictionary, String p_otherLibrary, String p_otherTree, String p_gaucheDictionary, String p_gaucheTree, String p_gaucheLibrary, String p_15Dictionary, String p_15Tree, String p_15Library,String p_abramDictionary,String p_abramTree,String p_abramLibrary,String p_unifacDictionary,String p_unifacTree,String p_unifacLibrary) { //,String p_solventDictionary,String p_solventLibrary) {

	        	// step 1: read in GA Groups
					System.out.println("Reading thermochemistry groups");
                    // read thermo functional Group dictionary
                    readGroupDictionary(p_groupDictionary);
                    // read thermo functional Group tree structure
                    readGroupTree(p_groupTree);
                    // read group values
                    readGroupLibrary(p_groupLibrary);

	        	// step 2: read in Radical Corrections
					System.out.println("Reading radical correction groups");
                    // read radical dictionary
                    readRadicalDictionary(p_radicalDictionary);
                    // read radical tree
                    readRadicalTree(p_radicalTree);
                    // read radical value
                    readRadicalLibrary(p_radicalLibrary);

	        	// step 3: read in Ring Correction
					System.out.println("Reading ring correction groups");
	                readRingDictionary(p_ringDictionary);
	                readRingTree(p_ringTree);
                    readRingLibrary(p_ringLibrary);

	        	// step 4: read in Other Correction
					System.out.println("Reading other correction groups");
                    readOtherDictionary(p_otherDictionary);
                    readOtherLibrary(p_otherLibrary);
                    readOtherTree(p_otherTree);

                // step 5: read in Gauche and 15 Correction libraries
					System.out.println("Reading gauche and 1/5 correction groups");
                    readGaucheDictionary(p_gaucheDictionary);
                    readGaucheTree(p_gaucheTree);
                    readGaucheLibrary(p_gaucheLibrary);
                    read15Dictionary(p_15Dictionary);
                    read15Tree(p_15Tree);
                    read15Library(p_15Library);

		if (Species.useSolvation) {
			// Definitions of Platts dictionary, library and tree for Abraham Model Implementation
			System.out.println("Reading Abraham solvation groups");
			readAbrahamDictionary(p_abramDictionary);
			readAbrahamTree(p_abramTree);
			readAbrahamLibrary(p_abramLibrary);
			System.out.println("Reading UNIFAC solvation groups");
			readUnifacDictionary(p_unifacDictionary);
			readUnifacTree(p_unifacTree);
			readUnifacLibrary(p_unifacLibrary);
		}


	    }


    //## operation readGroupDictionary(String)
    public void readGroupDictionary(String p_fileName) {
        //#[ operation readGroupDictionary(String)
        try {
        	groupDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read group dictionary!");
        	System.exit(0);
        }
        //#]
    }

    //## operation readGroupLibrary(String)
    public void readGroupLibrary(String p_fileName) {
        try {
        	groupLibrary = readStandardLibrary(p_fileName, groupDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read Group Library!");
        	System.exit(0);
        }
    }

    //## operation readGroupTree(String)
    public void readGroupTree(String p_fileName) {
        try {
        	groupTree = readStandardTree(p_fileName,groupDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read thermo group tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }
    }

    //## operation readOtherDictionary(String)
    public void readOtherDictionary(String p_fileName) {
        try {
        	otherDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read other dictionary!");
        	System.exit(0);
        }
    }

    //## operation readOtherLibrary(String)
    public void readOtherLibrary(String p_fileName) {
        try {
        	otherLibrary = readStandardLibrary(p_fileName, otherDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read Other Library!");
        	System.exit(0);
        }
    }

    //## operation readOtherTree(String)
    public void readOtherTree(String p_fileName) {
        try {
        	otherTree = readStandardTree(p_fileName,otherDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read thermo Other tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }
    }
    
    //2/5/09 gmagoon: new functions for gauche and 1,5 correction reading (based on analogs for regular values, e.g. readGroupDictionary)
    public void readGaucheDictionary(String p_fileName) {
        try {
        	gaucheDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read gauche dictionary!");
        	System.exit(0);
        }
    }

    public void readGaucheLibrary(String p_fileName) {
        try {
        	gaucheLibrary = readStandardLibrary(p_fileName, gaucheDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read gauche library!");
        	System.exit(0);
        }
    }

    public void readGaucheTree(String p_fileName) {
        try {
        	gaucheTree = readStandardTree(p_fileName,gaucheDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read gauche tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }
    }

	public void readAbrahamDictionary(String p_fileName) {
        try {
        	abramDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read Abraham dictionary!");
        	System.exit(0);
        }
    }

	public void readUnifacDictionary(String p_fileName) {
        try {
        	unifacDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read Unifac dictionary!");
        	System.exit(0);
        }
    }

	public void readAbrahamLibrary(String p_fileName) {
        try {
        	abramLibrary = readAbramLibrary(p_fileName, abramDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read Abraham library!");
        	System.exit(0);
        }
    }

	public void readUnifacLibrary(String p_fileName) {
        try {
        	unifacLibrary = readUNIFACLibrary(p_fileName, unifacDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read Unifac library!");
        	System.exit(0);
        }
    }

	public void readAbrahamTree(String p_fileName) {
        try {
        	abramTree = readStandardTree(p_fileName,abramDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read Abraham tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }
    }

	public void readUnifacTree(String p_fileName) {
        try {
        	unifacTree = readStandardTree(p_fileName,unifacDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read Unifac tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }
    }

        //public void readSolventDictionary(String p_fileName) {
        //try {
        //	solventDictionary = readStandard/Dictionary(p_fileName);
        //	return;
        //}
        //catch (Exception e) {
        //	System.err.println("Error in read Solvent dictionary!");
        //	System.exit(0);
        //}
        //#]
    //}

      //      public void readSolventLibrary(String p_fileName) {
      //  try {
      //  	solventLibrary = readStandardLibrary(p_fileName, solventDictionary);
      //  	return;
      //  }
      //  catch (Exception e) {
      //  	System.err.println("Can't read solvent library!");
      //  	System.exit(0);
      //  }

        //#]
    //}

	public void read15Dictionary(String p_fileName) {
        try {
        	oneFiveDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read 1,5 dictionary!");
        	System.exit(0);
        }
    }

    public void read15Library(String p_fileName) {
        try {
        	oneFiveLibrary = readStandardLibrary(p_fileName, oneFiveDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read 1,5 library!");
        	System.exit(0);
        }
    }

    public void read15Tree(String p_fileName) {
        try {
        	oneFiveTree = readStandardTree(p_fileName,oneFiveDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read 1,5 tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }
    }

    //## operation readRadicalDictionary(String)
    public void readRadicalDictionary(String p_fileName) {
        try {
        	radicalDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read radical dictionary!\n" + e.getMessage());
        	System.exit(0);
        }
    }

    //## operation readRadicalLibrary(String)
    public void readRadicalLibrary(String p_fileName) {
        try {
        	radicalLibrary = readStandardLibrary(p_fileName, radicalDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read radical Library!");
        	System.exit(0);
        }
    }

    //## operation readRadicalTree(String)
    public void readRadicalTree(String p_fileName) {
         try {
        	radicalTree = readStandardTree(p_fileName,radicalDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read thermo group tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }
    }

    //## operation readRingLibrary(String)
	public void readRingDictionary(String p_fileName) {
        try {
                ringDictionary = readStandardDictionary(p_fileName);
                return;
        }
        catch (Exception e) {
                System.err.println("Error in read ring dictionary!\n" + e.getMessage());
                System.exit(0);
        }
    }

    //## operation readRingTree(String)
    public void readRingTree(String p_fileName) {
         try {
                ringTree = readStandardTree(p_fileName,ringDictionary,0);
        }
        catch (Exception e) {
                System.err.println("Can't read ring tree file!");
                System.err.println("Error: " + e.getMessage());
                System.exit(0);
        }
    }

// end pey

    //## operation readRingLibrary(String)
    public void readRingLibrary(String p_fileName) {
        try {
			ringLibrary = readStandardLibrary(p_fileName, ringDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read Ring Correction Library!");
			System.err.println("Error: " + e);
        	System.exit(0);
        }
    }


    //## operation readStandardCorrectionLibrary(String,HashMap)
    protected void readStandardCorrectionLibrary(String p_fileName, HashMap p_library) throws IOException {
        try {
                FileReader in = new FileReader(p_fileName);
                BufferedReader data = new BufferedReader(in);

                String line = ChemParser.readMeaningfulLine(data);
                 while (line != null) {
                        // step 1: read in index and name
                        StringTokenizer token = new StringTokenizer(line);
                        int index = Integer.parseInt(token.nextToken());
                        String name = token.nextToken();
                        if (p_library == ringLibrary) {
                                String fomula = token.nextToken();
                                String sigma = token.nextToken();
                        }

                        // setp 2: read in thermoGAValue
                        String thermo="";
                for (int i=0;i<12;i++) {
                                thermo = thermo.concat(token.nextToken());
                                thermo = thermo.concat(" ");
                        }
                        ThermoGAValue gaValue = ChemParser.parseThermoGAValue(thermo);
                        String comments = "";
                        while (token.hasMoreTokens()) {
                               comments = comments + " " + token.nextToken();
                       }
                        ThermoGAValue newGaValue = new ThermoGAValue(name,gaValue,comments);

                        // step3: read in graph of the functional group
                        Graph g = ChemParser.readFGGraph(data);
                        if (g == null) throw new NullGraphException();
                        FunctionalGroup fg = FunctionalGroup.make(name, g);

                // step4: put in library
                        Object previous = p_library.put(fg, newGaValue);
                        if (previous != null) {
                                throw new ReplaceThermoGAValueException();
                        }
                        line = ChemParser.readMeaningfulLine(data);
                }

                in.close();
                return;
        }
        catch (IOException e) {
                throw new IOException();
        }
    }


    //## operation readStandardDictionary(String)
    public HashMap readStandardDictionary(String p_fileName) throws FileNotFoundException, IOException {
        //#[ operation readStandardDictionary(String)
        try {
                FileReader in = new FileReader(p_fileName);
                BufferedReader data = new BufferedReader(in);
            HashMap dictionary = new HashMap();
            HashMap unRead = new HashMap();

                 String line = ChemParser.readMeaningfulLine(data);

                read: while (line != null) {
                        StringTokenizer st = new StringTokenizer(line);
                        String fgname = st.nextToken();
                    data.mark(10000);
                        line = ChemParser.readMeaningfulLine(data);
                        if (line == null) break read;
                        line = line.trim();
                        String prefix = line.substring(0,5);
                        if (prefix.compareToIgnoreCase("union") == 0) {
                                HashSet union = ChemParser.readUnion(line);
                                 unRead.put(fgname,union);
                        }
                        else {
                                data.reset();
                                Graph fgGraph = null;
                                try {
                                        fgGraph = ChemParser.readFGGraph(data);
                                }
                                catch (Exception e) {
                                        throw new InvalidFunctionalGroupException(fgname + ": " + e.getMessage());
                                }
                                if (fgGraph == null) throw new InvalidFunctionalGroupException(fgname);

                                FunctionalGroup fg = FunctionalGroup.make(fgname, fgGraph);
                                Object old = dictionary.get(fgname);
                                if (old == null) {
                                        dictionary.put(fgname,fg);
                                }
                                else {
                                        FunctionalGroup oldFG = (FunctionalGroup)old;
                                        if (!oldFG.equals(fg)) throw new ReplaceFunctionalGroupException(fgname);
                                }
						}
//System.out.println(line);
                        line = ChemParser.readMeaningfulLine(data);
                }

                while (!unRead.isEmpty()) {
                        String fgname = (String)(unRead.keySet().iterator().next());
                        ChemParser.findUnion(fgname,unRead,dictionary);
                }

                in.close();
                return dictionary;
        }
        catch (FileNotFoundException e) {
                throw new FileNotFoundException(p_fileName);
        }
        catch (IOException e) {
                throw new IOException(p_fileName + ": " + e.getMessage());
        }
    }


    //## operation readStandardLibrary(String,HashMap)
    protected HashMap readStandardLibrary(String p_fileName, HashMap p_dictionary) throws IOException {
        try {
                FileReader in = new FileReader(p_fileName);
                BufferedReader data = new BufferedReader(in);
            HashMap library = new HashMap();

            String line = ChemParser.readMeaningfulLine(data);
                while (line != null) {
//System.out.println(line);//
                        // step 1: read in index and name
                        StringTokenizer token = new StringTokenizer(line);
                        String index = token.nextToken(); //1/6/09 gmagoon changed index from integer to string, so that if/when ChemGreen/RMGVE adds a decimal after the entry number (after editing thermo library), RMG will still be able to read it
                        String name = token.nextToken();

                        // step 2: find this functional group in dictionary by name
                        Matchable fg = (Matchable)p_dictionary.get(name);
                        if (fg == null) {
                                throw new FunctionalGroupNotFoundException();
                                //System.out.println(name);
                        }

                        // step 3: read in thermoGAValue
                        String thermo = token.nextToken();
                        // if there is a set of real thermo numbers, read them in and put the thermo data into library
                        try {
                                double H = Double.parseDouble(thermo);
                                thermo = thermo.concat(" ");
                                for (int i=0;i<11;i++) {
                                        thermo = thermo.concat(token.nextToken());
                                        thermo = thermo.concat(" ");
                                }
                                ThermoGAValue gaValue = ChemParser.parseThermoGAValue(thermo);
                                String comments = "";
                                while (token.hasMoreTokens()) {
                                   comments = comments + " " + token.nextToken();
                           }
                                ThermoGAValue newGaValue=new ThermoGAValue(name,gaValue,comments);

                        // step4: put in library
                                Object previous = library.put(fg, newGaValue);
                                if (previous != null) {
                                        throw new ReplaceThermoGAValueException();
                                }

                        }
                        // if there is a referenced name, put the name into library
                        catch (NumberFormatException e) {
                                Object o = p_dictionary.get(thermo);
                                if (o == null) {
                                        //throw new FunctionalGroupNotFoundException(thermo);
                                        System.out.print(index);
                                        System.out.println(": " + thermo);
                                }
                                Object previous = library.put(fg, thermo);
                                if (previous != null) {
                                        throw new ReplaceThermoGAValueException();
                                }
                        }

                        line = ChemParser.readMeaningfulLine(data);
                }

                // scan the library to give the ones having referenced name the real thermo data
                Iterator iter = library.keySet().iterator();
                while (iter.hasNext()) {
                        Matchable fg = (Matchable)iter.next();
                        Object gaValue = library.get(fg);
                        String path = "";
                        if (gaValue instanceof String) {
                                do {
                                        String name = (String)gaValue;
                                        path = path + "->" + name;
                                    gaValue = library.get((Matchable)p_dictionary.get(name));
                                } while (gaValue instanceof String);

                                if (gaValue == null || !(gaValue instanceof ThermoGAValue)) {
                                        throw new InvalidReferenceThermoGAValueException();
                                }

                                ThermoGAValue newGaValue = new ThermoGAValue(fg.getName(),(ThermoGAValue)gaValue, "Use the value of " + path);
                                library.put(fg,newGaValue);
                         }
                }

            in.close();
                return library;
        }
        catch (IOException e) {
                throw new IOException();
        }
    }

    //## operation readAbramLibrary(String,HashMap)
    protected HashMap readAbramLibrary(String p_fileName, HashMap p_dictionary) throws IOException {
        //#[ operation readStandardLibrary(String,HashMap)
        try {
                FileReader in = new FileReader(p_fileName);
                BufferedReader data = new BufferedReader(in);
            HashMap library = new HashMap();

            String line = ChemParser.readMeaningfulLine(data);
                while (line != null) {

                        // step 1: read in index and name
                        StringTokenizer token = new StringTokenizer(line);
                        String index = token.nextToken(); //1/6/09 gmagoon changed index from integer to string, so that if/when ChemGreen/RMGVE adds a decimal after the entry number (after editing thermo library), RMG will still be able to read it
                        String name = token.nextToken();

                        // step 2: find this functional group in dictionary by name
                        Matchable fg = (Matchable)p_dictionary.get(name);
                        if (fg == null) {
                                throw new FunctionalGroupNotFoundException();
                                //System.out.println(name);
                        }

                        // step 3: read in AbrahamGAValue
                        String thermo = token.nextToken();
                        // if there is a set of real thermo numbers, read them in and put the thermo data into library
                        try {
                                double H = Double.parseDouble(thermo);
                                thermo = thermo.concat(" ");
                                for (int i=0;i<4;i++) {
                                        thermo = thermo.concat(token.nextToken());
                                        thermo = thermo.concat(" ");
                                }
                                AbrahamGAValue gaValue = ChemParser.parseAbrahamGAValue(thermo);
                                String comments = "";
                                while (token.hasMoreTokens()) {
                                   comments = comments + " " + token.nextToken();
                           }
                                AbrahamGAValue newGaValue=new AbrahamGAValue(gaValue);

                        // step4: put in library
                                Object previous = library.put(fg, newGaValue);
                                if (previous != null) {
                                        throw new ReplaceThermoGAValueException();
                                }

                        }
                        // if there is a referenced name, put the name into library
                        catch (NumberFormatException e) {
                                Object o = p_dictionary.get(thermo);
                                if (o == null) {
                                        //throw new FunctionalGroupNotFoundException(thermo);
                                        System.out.print(index);
                                        System.out.println(": " + thermo);
                                }
                                Object previous = library.put(fg, thermo);
                                if (previous != null) {
                                      throw new ReplaceThermoGAValueException();
                                }
                        }

                        line = ChemParser.readMeaningfulLine(data);
                }

                // scan the library to give the ones having referenced name the real thermo data
                Iterator iter = library.keySet().iterator();
                while (iter.hasNext()) {
                        Matchable fg = (Matchable)iter.next();
                        Object gaValue = library.get(fg);
                        String path = "";
                        if (gaValue instanceof String) {
                                do {
                                        String name = (String)gaValue;
                                        path = path + "->" + name;
                                    gaValue = library.get((Matchable)p_dictionary.get(name));
                                } while (gaValue instanceof String);

                                if (gaValue == null || !(gaValue instanceof AbrahamGAValue)) {
                                        throw new InvalidReferenceThermoGAValueException();
                                }

                                AbrahamGAValue newGaValue = new AbrahamGAValue((AbrahamGAValue)gaValue);
                                library.put(fg,newGaValue);
                         }
                }

            in.close();
                return library;
        }
        catch (IOException e) {
                throw new IOException();
        }
    }

    //## operation readAbramLibrary(String,HashMap)
    protected HashMap readUNIFACLibrary(String p_fileName, HashMap p_dictionary) throws IOException {
        try {
                FileReader in = new FileReader(p_fileName);
                BufferedReader data = new BufferedReader(in);
            HashMap library = new HashMap();

            String line = ChemParser.readMeaningfulLine(data);
                while (line != null) {

                        // step 1: read in index and name
                        StringTokenizer token = new StringTokenizer(line);
                        String index = token.nextToken(); //1/6/09 gmagoon changed index from integer to string, so that if/when ChemGreen/RMGVE adds a decimal after the entry number (after editing thermo library), RMG will still be able to read it
                        String name = token.nextToken();

                        // step 2: find this functional group in dictionary by name
                        Matchable fg = (Matchable)p_dictionary.get(name);
                        if (fg == null) {
                                throw new FunctionalGroupNotFoundException();
                                //System.out.println(name);
                        }

                        // step 3: read in AbrahamGAValue
                        String thermo = token.nextToken();
                        // if there is a set of real thermo numbers, read them in and put the thermo data into library
                        try {
                                double H = Double.parseDouble(thermo);
                                thermo = thermo.concat(" ");
                                for (int i=0;i<1;i++) {
                                        thermo = thermo.concat(token.nextToken());
                                        thermo = thermo.concat(" ");
                                }
                                UnifacGAValue gaValue = ChemParser.parseUnifacGAValue(thermo);
                                String comments = "";
                                while (token.hasMoreTokens()) {
                                   comments = comments + " " + token.nextToken();
                           }
                                UnifacGAValue newGaValue=new UnifacGAValue(gaValue);

                        // step4: put in library
                                Object previous = library.put(fg, newGaValue);
                                if (previous != null) {
                                        throw new ReplaceThermoGAValueException();
                                }

                        }
                        // if there is a referenced name, put the name into library
                        catch (NumberFormatException e) {
                                Object o = p_dictionary.get(thermo);
                                if (o == null) {
                                        //throw new FunctionalGroupNotFoundException(thermo);
                                        System.out.print(index);
                                        System.out.println(": " + thermo);
                                }
                                Object previous = library.put(fg, thermo);
                                if (previous != null) {
                                      throw new ReplaceThermoGAValueException();
                                }
                        }

                        line = ChemParser.readMeaningfulLine(data);
                }

                // scan the library to give the ones having referenced name the real thermo data
                Iterator iter = library.keySet().iterator();
                while (iter.hasNext()) {
                        Matchable fg = (Matchable)iter.next();
                        Object gaValue = library.get(fg);
                        String path = "";
                        if (gaValue instanceof String) {
                                do {
                                        String name = (String)gaValue;
                                        path = path + "->" + name;
                                    gaValue = library.get((Matchable)p_dictionary.get(name));
                                } while (gaValue instanceof String);

                                if (gaValue == null || !(gaValue instanceof UnifacGAValue)) {
                                        throw new InvalidReferenceThermoGAValueException();
                                }

                                UnifacGAValue newGaValue = new UnifacGAValue((UnifacGAValue)gaValue);
                                library.put(fg,newGaValue);
                         }
                }

            in.close();
                return library;
        }
        catch (IOException e) {
                throw new IOException();
        }
    }








    //## operation readStandardTree(String,HashMap,int)
    public HierarchyTree readStandardTree(String p_fileName, HashMap p_dictionary, int p_level) throws IOException {
        //#[ operation readStandardTree(String,HashMap,int)
        try {
        	FileReader in = new FileReader(p_fileName);
        	BufferedReader data = new BufferedReader(in);

        	HierarchyTree tree = ChemParser.readHierarchyTree(data,p_dictionary,p_level);

        	in.close();

        	return tree;
        }
        catch (IOException e) {
        	throw new IOException(p_fileName);
        }
    }

    protected static ThermoGAGroupLibrary getINSTANCE() {
        return INSTANCE;
    }

    public HashMap getGroupDictionary() {
        return groupDictionary;
    }

    public HashMap getGroupLibrary() {
        return groupLibrary;
    }

    public void setGroupLibrary(HashMap p_groupLibrary) {
        groupLibrary = p_groupLibrary;
    }

    protected HierarchyTree getGroupTree() {
        return groupTree;
    }

    public HashMap getOtherDictionary() {
        return otherDictionary;
    }

    public void setOtherDictionary(HashMap p_otherDictionary) {
        otherDictionary = p_otherDictionary;
    }

    public HashMap getOtherLibrary() {
        return otherLibrary;
    }

    public HierarchyTree getOtherTree() {
        return otherTree;
    }

    public void setOtherTree(HierarchyTree p_otherTree) {
        otherTree = p_otherTree;
    }

    public HashMap getRadicalDictionary() {
        return radicalDictionary;
    }

    public void setRadicalDictionary(HashMap p_radicalDictionary) {
        radicalDictionary = p_radicalDictionary;
    }

    protected HashMap getRadicalLibrary() {
        return radicalLibrary;
    }

    public HierarchyTree getRadicalTree() {
        return radicalTree;
    }

    public void setRadicalTree(HierarchyTree p_radicalTree) {
        radicalTree = p_radicalTree;
    }

    protected HashMap getRingLibrary() {
        return ringLibrary;
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\ThermoGAGroupLibrary.java
*********************************************************************/


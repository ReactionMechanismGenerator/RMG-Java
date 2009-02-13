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

	   // begin pey
    protected HierarchyTree ringTree;
    protected HashMap ringDictionary;
 // end pey

    protected HashMap ringLibrary;		//## attribute ringLibrary
    
    protected HashMap gaucheDictionary;		
    protected HashMap gaucheLibrary;		
    protected HierarchyTree gaucheTree;	
    
    protected HashMap oneFiveDictionary;		
    protected HashMap oneFiveLibrary;		
    protected HierarchyTree oneFiveTree;	


    // Constructors

    //## operation ThermoGAGroupLibrary()
    private  ThermoGAGroupLibrary() {
        //#[ operation ThermoGAGroupLibrary()
        groupTree = new HierarchyTree();
        groupDictionary = new HashMap();
        groupLibrary = new HashMap();
        radicalTree = new HierarchyTree();
        radicalDictionary = new HashMap();
        radicalLibrary = new HashMap();
        ringLibrary = new HashMap();
//		 begin pey
        ringDictionary = new HashMap();
        ringTree = new HierarchyTree();
        // end pey
        otherLibrary = new HashMap();
        otherDictionary = new HashMap();
        otherTree = new HierarchyTree();
        
        gaucheLibrary = new HashMap();
        gaucheDictionary = new HashMap();
        gaucheTree = new HierarchyTree();
        
        oneFiveLibrary = new HashMap();
        oneFiveDictionary = new HashMap();
        oneFiveTree = new HierarchyTree();

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
//		 begin pey
        //String ring = directory + "Ring_Corrections.txt";
        String ringDictionary = directory + "Ring_Dictionary.txt";
        String ringTree = directory + "Ring_Tree.txt";
        String ringLibrary = directory + "Ring_Library.txt";
        // end pey
        String otherLibrary = directory + "Other_Library_Dictionary.txt";
        String otherTree = directory + "Other_Tree.txt";
        
        String gauDictionary = directory + "Gauche_Dictionary.txt";
        String gauTree = directory + "Gauche_Tree.txt";
        String gauLibrary = directory + "Gauche_Library.txt";
        
        String one5Dictionary = directory + "15_Dictionary.txt";
        String one5Tree = directory + "15_Tree.txt";
        String one5Library = directory + "15_Library.txt";
	
        read(gDictionary,gTree,gLibrary,rDictionary,rTree,rLibrary,ringDictionary,ringTree,ringLibrary,otherLibrary,otherTree,gauDictionary,gauTree,gauLibrary,one5Dictionary,one5Tree,one5Library);



        //#]
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
        	if (ga != null) return ga;
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

    //## operation read(String,String,String,String,String,String,String,String,String)
	public void read(String p_groupDictionary, String p_groupTree, String p_groupLibrary, String p_radicalDictionary, String p_radicalTree, String p_radicalLibrary, String p_ringDictionary, String p_ringTree, String p_ringLibrary, String p_otherLibrary, String p_otherTree, String p_gaucheDictionary, String p_gaucheTree, String p_gaucheLibrary, String p_15Dictionary, String p_15Tree, String p_15Library) {
	    // end pey
	        //#[ operation read(String,String,String,String,String,String,String,String,String)
	        // try {
	        	// step 1: read in GA Groups
	        	// read thermo functional Group dictionary
	         	readGroupDictionary(p_groupDictionary);
	        	// read thermo functional Group tree structure
	        	readGroupTree(p_groupTree);
	        	// read group values
	        	readGroupLibrary(p_groupLibrary);

	        	// step 2: read in Radical Corrections
	        	// read radical dictionary
	        	readRadicalDictionary(p_radicalDictionary);
	        	// read radical tree
	        	readRadicalTree(p_radicalTree);
	        	// read radical value
	        	readRadicalLibrary(p_radicalLibrary);

	        	// step 3: read in Ring Correction
	                // begin pey
	                readRingDictionary(p_ringDictionary);
	                readRingTree(p_ringTree);
	        	readRingLibrary(p_ringLibrary);
	                // System.out.println("tree height = " + ringTree.height());
	                // end pey

	        	// step 4: read in Other Correction
	        	readOtherLibrary(p_otherLibrary);
	        	readOtherTree(p_otherTree);

                        // step 5: read in Gauche and 15 Correction libraries
                        readGaucheDictionary(p_gaucheDictionary);
                        readGaucheTree(p_gaucheTree);
	        	readGaucheLibrary(p_gaucheLibrary);
                        read15Dictionary(p_15Dictionary);
                        read15Tree(p_15Tree);
	        	read15Library(p_15Library);
	        /*}
	        catch (Exception e) {
	        	throw new ThermoIOException(e.getMessage());
	        }
	        */



	        //#]
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
        //#[ operation readGroupLibrary(String)
        try {
        	groupLibrary = readStandardLibrary(p_fileName, groupDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read Group Library!");
        	System.exit(0);
        }




        //#]
    }

    //## operation readGroupTree(String)
    public void readGroupTree(String p_fileName) {
        //#[ operation readGroupTree(String)
        try {
        	groupTree = readStandardTree(p_fileName,groupDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read thermo group tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }



        //#]
    }

    //## operation readOtherLibrary(String)
    public void readOtherLibrary(String p_fileName) {
        //#[ operation readOtherLibrary(String)
        try {
        	readStandardCorrectionLibrary(p_fileName, otherLibrary);
        	Iterator iter = otherLibrary.keySet().iterator();
        	while (iter.hasNext()) {
        		FunctionalGroup fg = (FunctionalGroup)iter.next();
        		otherDictionary.put(fg.name, fg);
        	}

        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read other correction Library!");
        	System.err.println("Error: " + e.getClass().getName() + "  :  " + e.getMessage());
        	System.exit(0);
        }




        //#]
    }

    //## operation readOtherTree(String)
    public void readOtherTree(String p_fileName) {
        //#[ operation readOtherTree(String)
        try {
        	otherTree = readStandardTree(p_fileName,otherDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read thermo Other tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }



        //#]
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
        //#]
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




        //#]
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



        //#]
    }
    
        public void read15Dictionary(String p_fileName) {
        try {
        	oneFiveDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read 1,5 dictionary!");
        	System.exit(0);
        }
        //#]
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




        //#]
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



        //#]
    }

    //## operation readRadicalDictionary(String)
    public void readRadicalDictionary(String p_fileName) {
        //#[ operation readRadicalDictionary(String)
        try {
        	radicalDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read radical dictionary!\n" + e.getMessage());
        	System.exit(0);
        }
        //#]
    }

    //## operation readRadicalLibrary(String)
    public void readRadicalLibrary(String p_fileName) {
        //#[ operation readRadicalLibrary(String)
        try {
        	radicalLibrary = readStandardLibrary(p_fileName, radicalDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read radical Library!");
        	System.exit(0);
        }




        //#]
    }

    //## operation readRadicalTree(String)
    public void readRadicalTree(String p_fileName) {
        //#[ operation readRadicalTree(String)
         try {
        	radicalTree = readStandardTree(p_fileName,radicalDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read thermo group tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }



        //#]
    }

    //## operation readRingLibrary(String)
	public void readRingDictionary(String p_fileName) {
        //#[ operation readRingDictionary(String)
        try {
                ringDictionary = readStandardDictionary(p_fileName);
                return;
        }
        catch (Exception e) {
                System.err.println("Error in read ring dictionary!\n" + e.getMessage());
                System.exit(0);
        }
        //#]
    }

    //## operation readRingTree(String)
    public void readRingTree(String p_fileName) {
        //#[ operation readRingTree(String)
         try {
                ringTree = readStandardTree(p_fileName,ringDictionary,0);
        }
        catch (Exception e) {
                System.err.println("Can't read ring tree file!");
                System.err.println("Error: " + e.getMessage());
                System.exit(0);
        }
        //#]
    }

// end pey

    //## operation readRingLibrary(String)
    public void readRingLibrary(String p_fileName) {
        //#[ operation readRingLibrary(String)
        try {
                // begin pey
        	// readStandardCorrectionLibrary(p_fileName, ringLibrary);
                ringLibrary = readStandardLibrary(p_fileName, ringDictionary);
                // end pey
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read Ring Correction Library!");
                // begin pey
                System.err.println("Error: " + e);
        	System.exit(0);
                // end pey
        }

        //#]
    }


    //## operation readStandardCorrectionLibrary(String,HashMap)
    protected void readStandardCorrectionLibrary(String p_fileName, HashMap p_library) throws IOException {
        //#[ operation readStandardCorrectionLibrary(String,HashMap)
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







        //#]
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
        //#]
    }


    //## operation readStandardLibrary(String,HashMap)
    protected HashMap readStandardLibrary(String p_fileName, HashMap p_dictionary) throws IOException {
        //#[ operation readStandardLibrary(String,HashMap)
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







        //#]
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



        //#]
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


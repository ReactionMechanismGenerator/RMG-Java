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

public class TransportGALibrary {

    protected static TransportGALibrary INSTANCE = new TransportGALibrary();
    protected HashMap groupDictionary;
    protected HashMap groupLibrary;
    protected HierarchyTree groupTree;
    // begin pey
    protected HierarchyTree ringTree;
    protected HashMap ringDictionary;
    // end pey
    protected HashMap ringLibrary;
    

    private TransportGALibrary() {
    	groupTree = new HierarchyTree();
    	groupDictionary = new HashMap();
    	groupLibrary = new HashMap();
    	ringLibrary = new HashMap();
    	//		 begin pey
    	ringDictionary = new HashMap();
    	ringTree = new HierarchyTree();
    	// end pey

        String directory = System.getProperty("jing.chem.LJDatabase.pathName");
        if (directory == null) {
        	System.out.println("undefined system property: jing.chem.LJDatabase.pathName, exit!");
        	System.exit(0);
        }

        String separator = System.getProperty("file.separator");
        if (!directory.endsWith(separator)) directory = directory + separator;

        System.out.println("\nReading Lennard-Jones database from "+directory);

        String gDictionary = directory + "nonring_Dictionary.txt";
        String gTree = directory + "nonring_Tree.txt";
        String gLibrary = directory + "nonring_Library.txt";
        //		 begin pey
        String ringDictionary = directory + "ring_Dictionary.txt";
        String ringTree = directory + "ring_Tree.txt";
        String ringLibrary = directory + "ring_Library.txt";
        // end pey

        read(gDictionary,gTree,gLibrary,ringDictionary,ringTree,ringLibrary);
    }

    /**
    Requires: the central node of p_chemGraph has been set to the transport center atom.
    Effects: find a matched transport functional group in the group tree for the pass-in 
    p_chemGraph, return this functional group's transport value.  If no leaf is found, 
    throw GroupNotFoundException
    Modifies:
    */
    public LJGroupData findGroup(ChemGraph p_chemGraph) throws GroupNotFoundException, MultipleGroupFoundException, InvalidCenterTypeException {
        if (p_chemGraph == null) return null;

        Stack stack = groupTree.findMatchedPath(p_chemGraph);
        p_chemGraph.getGraph().resetMatchedGC();//2/13/09 gmagoon: resetting the matched GC value...for some reason, thermoLibrary.findGAGroup (within getGAGroup in GATP.java) ended up modifiying the central node so that it was matched; this ended up wreaking havoc with subsequent symmetry number calculations; ideally, I would probably want to fix the code so that it didn't end up modifying the matchedGC from the null value after it is done with it, but I do not immediately see how to due so, and debugging proved extremely difficult; I have also tried to put this elsewhere in this class where it might be appropriate
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
        	Matchable fg = (Matchable)node.getElement();
        	LJGroupData ga = (LJGroupData)groupLibrary.get(fg);
        	if (ga != null) //{
                //System.out.println("Group found: " + fg.getName());
                return ga;
            //}
        }
        
        return null;
    }

    public LJGroupData findRingGroup(ChemGraph p_chemGraph) throws GroupNotFoundException, MultipleGroupFoundException, InvalidCenterTypeException {
        if (p_chemGraph == null) return null;

        Stack stack = ringTree.findMatchedPath(p_chemGraph);
        p_chemGraph.getGraph().resetMatchedGC();//2/13/09 gmagoon: resetting the matched GC value...for some reason, thermoLibrary.findGAGroup (within getGAGroup in GATP.java) ended up modifiying the central node so that it was matched; this ended up wreaking havoc with subsequent symmetry number calculations; ideally, I would probably want to fix the code so that it didn't end up modifying the matchedGC from the null value after it is done with it, but I do not immediately see how to due so, and debugging proved extremely difficult; I have also tried to put this elsewhere in this class where it might be appropriate
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
        	Matchable fg = (Matchable)node.getElement();
        	LJGroupData ga = (LJGroupData)ringLibrary.get(fg);
        	if (ga != null) //{
                //System.out.println("Group found: " + fg.getName());
                return ga;
            //}
        }
        
        return null;
    }

    public void read(String p_groupDictionary, String p_groupTree, String p_groupLibrary, String p_ringDictionary, String p_ringTree, String p_ringLibrary) { //,String p_solventDictionary,String p_solventLibrary) {
    	// step 1: read in GA Groups
    	// read thermo functional Group dictionary
    	readGroupDictionary(p_groupDictionary);
    	// read thermo functional Group tree structure
    	readGroupTree(p_groupTree);
    	// read group values
    	readGroupLibrary(p_groupLibrary);
    	
    	// step 2: read in Ring Correction
    	// begin pey
    	readRingDictionary(p_ringDictionary);
    	readRingTree(p_ringTree);
    	readRingLibrary(p_ringLibrary);
    	// end pey
	}

    public void readGroupDictionary(String p_fileName) {
        try {
        	groupDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read nonring_Dictionary ("+p_fileName+") for transport properties!");
        	System.exit(0);
        }
    }

    public void readGroupLibrary(String p_fileName) {
        try {
        	groupLibrary = readStandardLibrary(p_fileName, groupDictionary);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read nonring_Library ("+p_fileName+") for transport properties!");
        	System.exit(0);
        }
    }

    public void readGroupTree(String p_fileName) {
        try {
        	groupTree = readStandardTree(p_fileName,groupDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read nonring_Tree ("+p_fileName+") for transport properties!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }
    }


	public void readRingDictionary(String p_fileName) {
        try {
        	ringDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
                System.err.println("Can't read ring_Dictionary ("+p_fileName+") for transport properties!\n" + e.getMessage());
                System.exit(0);
        }
    }

    public void readRingTree(String p_fileName) {
    	try {
    		ringTree = readStandardTree(p_fileName,ringDictionary,0);
    	}
    	catch (Exception e) {
                System.err.println("Can't read ring_Tree ("+p_fileName+") for transport properties!");
                System.err.println("Error: " + e.getMessage());
                System.exit(0);
        }
    }

// end pey

    public void readRingLibrary(String p_fileName) {
        try {
        	// begin pey
        	ringLibrary = readStandardLibrary(p_fileName, ringDictionary);
        	// end pey
        	return;
        }
        catch (Exception e) {
        	System.err.println("Can't read ring_Library for transport properties!");
        	// begin pey
        	System.err.println("Error: " + e);
        	System.exit(0);
        	// end pey
        }
    }

    public HashMap readStandardDictionary(String p_fileName) throws FileNotFoundException, IOException {
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


    protected HashMap readStandardLibrary(String p_fileName, HashMap p_dictionary) throws IOException {
        try {
        	FileReader in = new FileReader(p_fileName);
        	BufferedReader data = new BufferedReader(in);
        	HashMap library = new HashMap();
        	
            String line = ChemParser.readMeaningfulLine(data);
            while (line != null) {
            	// step 1: read in index and name
            	StringTokenizer token = new StringTokenizer(line);
            	//String index = token.nextToken(); //1/6/09 gmagoon changed index from integer to string, so that if/when ChemGreen/RMGVE adds a decimal after the entry number (after editing thermo library), RMG will still be able to read it
            	String name = token.nextToken();

            	// step 2: find this functional group in dictionary by name
            	Matchable fg = (Matchable)p_dictionary.get(name);
            	if (fg == null) {
            		throw new FunctionalGroupNotFoundException();
            	}

            	// step 3: read in LJGroupData
            	String dTc = token.nextToken();

            	try{
            		double H = Double.parseDouble(dTc);//dummy line to check whether first entry is a number; if not, it is a reference, and will raise exception that will be caught below
            		String dPc = token.nextToken();
            		String dVc = token.nextToken();
            		String dTb = token.nextToken();
            		int shapeIndex = Integer.parseInt(token.nextToken());
            		String comments = "";
            		while (token.hasMoreTokens()) {
            			comments = comments + " " + token.nextToken();
            		}
            		LJGroupData gaValue = new LJGroupData(Double.parseDouble(dTc), Double.parseDouble(dPc), Double.parseDouble(dVc), Double.parseDouble(dTb), shapeIndex, comments);
            		LJGroupData newGaValue=new LJGroupData(name,gaValue,comments);

            		// step4: put in library
            		Object previous = library.put(fg, newGaValue);
            		if (previous != null) {
            			throw new ReplaceThermoGAValueException();
            		}
            	}
            	// if there is a referenced name, put the name into library
            	catch (NumberFormatException e) {
            		Object o = p_dictionary.get(dTc);
            		if (o == null) {
            			System.out.print(name);
            			System.out.println(": " + dTc);
            		}
            		Object previous = library.put(fg, dTc);
            		if (previous != null) {
            			throw new ReplaceThermoGAValueException();
            		}
            	}
            	line = ChemParser.readMeaningfulLine(data);
            }


            // scan the library to give the ones having referenced name the real transport data
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
            		if (gaValue == null || !(gaValue instanceof LJGroupData)) {
            			throw new InvalidReferenceThermoGAValueException();
            		}
            		
            		LJGroupData newGaValue = new LJGroupData(fg.getName(),(LJGroupData)gaValue, "Use the value of " + path);
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

    public HierarchyTree readStandardTree(String p_fileName, HashMap p_dictionary, int p_level) throws IOException {
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

    protected static TransportGALibrary getINSTANCE() {
        return INSTANCE;
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\LJDatabase.java
*********************************************************************/


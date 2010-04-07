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

public class FrequencyDatabase {

    protected static FrequencyDatabase INSTANCE = new FrequencyDatabase();		//## attribute INSTANCE

    protected HashMap freqDictionary;		//## attribute freqDictionary

    protected HierarchyTree freqTree;		//## attribute freqTree


    // Constructors

    private  FrequencyDatabase() {
        freqTree = new HierarchyTree();
        freqDictionary = new HashMap();
        
		String directory = System.getProperty("jing.chem.FrequencyDatabase.pathName");
        if (directory == null) {
        	System.out.println("undefined system property: jing.chem.FrequencyDatabase.pathName, exit!");
        	System.exit(0);
        }

        String separator = System.getProperty("file.separator");
        if (!directory.endsWith(separator)) directory = directory + separator;
        System.out.println("\nReading frequency database from "+directory);

        String gDictionary = directory + "Dictionary.txt";
        String gTree = directory + "Tree.txt";

        read(gDictionary,gTree);

    }


    //gmagoon 11/17/08: findFreqGroupName based on findGAGroup from ThermoGAGroupLibrary
    /**
    Requires: the central node of p_chemGraph has been set to the thermo center atom.
    Effects: find a matched frequency functional group in the group tree for the pass-in p_chemGraph, return this functional group's name.  If no leaf is found, return a null string
    Modifies:
    */
    public String findFreqGroupName(ChemGraph p_chemGraph) throws MultipleGroupFoundException, InvalidCenterTypeException {
        if (p_chemGraph == null) return null;

        Stack stack = freqTree.findMatchedPath(p_chemGraph);
        if (stack == null) return null;

        while (!stack.empty()) {
        	HierarchyTreeNode node = (HierarchyTreeNode)stack.pop(); //gets the "deepest" node of the stack
        	Matchable fg = (Matchable)node.getElement();
        	String fgname = fg.getName();//get the group name from the dictionary/tree
        	if (fgname != null) return fgname;
        }

        return null;



        //#]
    }



    //## operation read(String,String)
	public void read(String p_groupDictionary, String p_groupTree) {

            // read thermo functional Group dictionary
            readGroupDictionary(p_groupDictionary);
            // read thermo functional Group tree structure
            readGroupTree(p_groupTree);

        }


    //## operation readGroupDictionary(String)
    public void readGroupDictionary(String p_fileName) {
        //#[ operation readGroupDictionary(String)
        try {
        	freqDictionary = readStandardDictionary(p_fileName);
        	return;
        }
        catch (Exception e) {
        	System.err.println("Error in read freq dictionary!");
                System.err.println("Error: " + e.getMessage());
        	System.exit(0);
        }
        //#]
    }


    //## operation readGroupTree(String)
    public void readGroupTree(String p_fileName) {
        //#[ operation readGroupTree(String)
        try {
        	freqTree = readStandardTree(p_fileName,freqDictionary,0);
        }
        catch (Exception e) {
        	System.err.println("Can't read freq group tree file!");
        	System.err.println("Error: " + e.getMessage());
        	System.exit(0);
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

    protected static FrequencyDatabase getINSTANCE() {
        return INSTANCE;
    }

    public HashMap getGroupDictionary() {
        return freqDictionary;
    }


    protected HierarchyTree getGroupTree() {
        return freqTree;
    }
}

/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FrequencyDatabase.java
*********************************************************************/


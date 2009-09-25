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
//import jing.chem.*;
//import jing.LJGroupData;

//----------------------------------------------------------------------------
// based on jing\chem\ThermoGAGroupLibrary.java
//----------------------------------------------------------------------------


//## class LJDatabase
public class LJDatabase {

    protected static LJDatabase INSTANCE = new LJDatabase();		//## attribute INSTANCE

    protected HashMap groupDictionary;		//## attribute groupDictionary

    protected HashMap groupLibrary;		//## attribute groupLibrary

   protected HierarchyTree groupTree;		//## attribute groupTree
	   // begin pey
    protected HierarchyTree ringTree;
    protected HashMap ringDictionary;
 // end pey

    protected HashMap ringLibrary;		//## attribute ringLibrary
    

    //## operation LJDatabase()
    private  LJDatabase() {
        //#[ operation ThermoGAGroupLibrary()
        groupTree = new HierarchyTree();
        groupDictionary = new HashMap();
        groupLibrary = new HashMap();
        ringLibrary = new HashMap();
//		 begin pey
        ringDictionary = new HashMap();
        ringTree = new HierarchyTree();
        // end pey

        String directory = System.getProperty("LJDatabase.pathName");
        if (directory == null) {
        	System.out.println("undefined system property: LJDatabase.pathName, exit!");
        	System.exit(0);
        }

        String separator = System.getProperty("file.separator");
        if (!directory.endsWith(separator)) directory = directory + separator;
			
			System.out.println("\nReading LJ database from "+directory);

        String gDictionary = directory + "nonring"+separator+ "Dictionary.txt";
        String gTree = directory + "nonring"+separator+ "Tree.txt";
        String gLibrary = directory + "nonring"+separator+ "Library.txt";
//		 begin pey
        String ringDictionary = directory + "ring"+separator+ "Dictionary.txt";
        String ringTree = directory + "ring"+separator+ "Tree.txt";
        String ringLibrary = directory + "ring"+separator+ "Library.txt";
        // end pey


        read(gDictionary,gTree,gLibrary,ringDictionary,ringTree,ringLibrary);



        //#]
    }

    //## operation findCorrectionInLibrary(ChemGraph,HashMap)
//    private LJGroupData findCorrectionInLibrary(ChemGraph p_chemGraph, HashMap p_library) {
//        //#[ operation findCorrectionInLibrary(ChemGraph,HashMap)
//        p_chemGraph.clearCentralNode();
//        LJGroupData result=new LJGroupData();
//        int redundance;
//
//        Iterator iter = p_library.keySet().iterator();
//        while (iter.hasNext()) {
//        	redundance = 0;
//        	FunctionalGroup f = (FunctionalGroup)iter.next();
//        	HashSet gv = p_chemGraph.identifyThermoMatchedSite(f);
//        	if (gv != null) {
//        		redundance = gv.size();
//        		if (redundance > 0) {
//        			LJGroupData ga = (LJGroupData)p_library.get(f);
//        			if (ga != null) {
//        				LJGroupData temp = new LJGroupData(ga);
//        				temp.multiply(redundance);
//        				result.plus(temp);
//        				temp = null;
//        			}
//        		}
//        	}
//        }
//        p_chemGraph.getGraph().resetMatchedGC();
//        return result;
//        //#]
//    }

    /**
    Requires: the central node of p_chemGraph has been set to the thermo center atom.
    Effects: find a matched thermo functional group in the group tree for the pass-in p_chemGraph, return this functional group's thermo value.  If no leaf is found, throw  GroupNotFoundException
    Modifies:
    */
    //## operation findGAGroup(ChemGraph)
    public LJGroupData findGroup(ChemGraph p_chemGraph) throws GroupNotFoundException, MultipleGroupFoundException, InvalidCenterTypeException {
        //#[ operation findGAGroup(ChemGraph)
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



        //#]
    }



    public LJGroupData findRingGroup(ChemGraph p_chemGraph) throws GroupNotFoundException, MultipleGroupFoundException, InvalidCenterTypeException {
        //#[ operation findGAGroup(ChemGraph)
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



        //#]
    }




    //## operation read(String,String,String,String,String,String,String,String,String)
	public void read(String p_groupDictionary, String p_groupTree, String p_groupLibrary, String p_ringDictionary, String p_ringTree, String p_ringLibrary) { //,String p_solventDictionary,String p_solventLibrary) {
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

	        	// step 2: read in Ring Correction
	                // begin pey
	                readRingDictionary(p_ringDictionary);
	                readRingTree(p_ringTree);
                    readRingLibrary(p_ringLibrary);
	                // System.out.println("tree height = " + ringTree.height());
	                // end pey

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


//    //## operation readStandardCorrectionLibrary(String,HashMap)
//    protected void readStandardCorrectionLibrary(String p_fileName, HashMap p_library) throws IOException {
//        //#[ operation readStandardCorrectionLibrary(String,HashMap)
//        try {
//                FileReader in = new FileReader(p_fileName);
//                BufferedReader data = new BufferedReader(in);
//
//                String line = ChemParser.readMeaningfulLine(data);
//                 while (line != null) {
//                        // step 1: read in index and name
//                        StringTokenizer token = new StringTokenizer(line);
//                        int index = Integer.parseInt(token.nextToken());
//                        String name = token.nextToken();
//                        if (p_library == ringLibrary) {
//                                String fomula = token.nextToken();
//                                String sigma = token.nextToken();
//                        }
//
//                        // setp 2: read in LJGroupData
//                        String thermo="";
//                for (int i=0;i<12;i++) {
//                                thermo = thermo.concat(token.nextToken());
//                                thermo = thermo.concat(" ");
//                        }
//                        LJGroupData gaValue = ChemParser.parseLJGroupData(thermo);
//                        String comments = "";
//                        while (token.hasMoreTokens()) {
//                               comments = comments + " " + token.nextToken();
//                       }
//                        LJGroupData newGaValue = new LJGroupData(name,gaValue,comments);
//
//                        // step3: read in graph of the functional group
//                        Graph g = ChemParser.readFGGraph(data);
//                        if (g == null) throw new NullGraphException();
//                        FunctionalGroup fg = FunctionalGroup.make(name, g);
//
//                // step4: put in library
//                        Object previous = p_library.put(fg, newGaValue);
//                        if (previous != null) {
//                                throw new ReplaceThermoGAValueException();
//                        }
//                        line = ChemParser.readMeaningfulLine(data);
//                }
//
//                in.close();
//                return;
//        }
//        catch (IOException e) {
//                throw new IOException();
//        }







        //#]
//    }


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

                        // step 3: read in LJGroupData
                        String dTc = token.nextToken();

                        try{
                            double H = Double.parseDouble(dTc);//dummy line to check whether first entry is a number; if not, it is a reference, and will raise exception that will be caught below
                            String dPc = token.nextToken();
                            String dVc = token.nextToken();
                            String dTb = token.nextToken();
                            String comments = "";
                            while (token.hasMoreTokens()) {
                                       comments = comments + " " + token.nextToken();
                            }
                            LJGroupData gaValue = new LJGroupData(Double.parseDouble(dTc), Double.parseDouble(dPc), Double.parseDouble(dVc), Double.parseDouble(dTb), comments);
                        // if there is a set of real thermo numbers, read them in and put the thermo data into library
                       // try {
                       //         double H = Double.parseDouble(thermo);
                       //         thermo = thermo.concat(" ");
                       //         for (int i=0;i<11;i++) {
                       //                 thermo = thermo.concat(token.nextToken());
                       //                 thermo = thermo.concat(" ");
                       //         }
                       //         LJGroupData gaValue = ChemParser.parseLJGroupData(thermo);
                       //         String comments = "";
                                
                       //    }
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
                                        //throw new FunctionalGroupNotFoundException(thermo);
                                        System.out.print(index);
                                        System.out.println(": " + dTc);
                                }
                                Object previous = library.put(fg, dTc);
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

    protected static LJDatabase getINSTANCE() {
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

    protected HashMap getRingLibrary() {
        return ringLibrary;
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\LJDatabase.java
*********************************************************************/


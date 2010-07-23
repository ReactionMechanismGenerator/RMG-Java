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
import jing.rxn.*;
import jing.rxnSys.*;
import java.util.*;
import jing.chemUtil.*;
import jing.chemParser.*;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\PrimaryThermoLibrary.java
//----------------------------------------------------------------------------


//## class PrimaryThermoLibrary
public class PrimaryThermoLibrary {
protected static HashMap library;
protected static HashMap dictionary;
private static PrimaryThermoLibrary INSTANCE = new PrimaryThermoLibrary();	//## attribute INSTANCE


private PrimaryThermoLibrary(){
  library = new HashMap();
  dictionary = new HashMap();
}

// 7-Jul-2009: MRH
public PrimaryThermoLibrary(HashMap Dictionary, HashMap Library) {
	dictionary = Dictionary;
	library = Library;
}

//7-Jul-2009: MRH
public PrimaryThermoLibrary(String name, String location) {
	// Create a new PrimaryThermoLibrary.
	library = new HashMap();
	dictionary = new HashMap();
	appendPrimaryThermoLibrary(name, location);
}

//7-Jul-2009: MRH
public void appendPrimaryThermoLibrary(String name, String path) {
	String thermoLibraryDirectory = System.getProperty("jing.chem.ThermoReferenceLibrary.pathName");
	String dictionaryFile = thermoLibraryDirectory + "/" + path + "/Dictionary.txt";
	String libraryFile = thermoLibraryDirectory + "/" + path + "/Library.txt";
	try {
		read(dictionaryFile, libraryFile, name);
	}
	catch (IOException e) {
		System.err.println("RMG cannot read Primary Thermo Library: " + name + "\n" + e.getMessage());
	}
}

//7-Jul-2009: MRH
public void read(String p_dictionary, String p_library, String p_name) throws IOException, FileNotFoundException {
	String source = "Primary Thermo Library: " + p_name;
	System.out.println("Reading " + source);
	dictionary = readDictionary(p_dictionary, source);
	library = readLibrary(p_library, dictionary, source);
}

//7-Jul-2009: MRH
public HashMap readLibrary(String p_thermoFileName, HashMap p_dictionary, String source) throws IOException {
	try{
		FileReader in = new FileReader(p_thermoFileName);
		BufferedReader data = new BufferedReader(in);
		String line = ChemParser.readMeaningfulLine(data);
		while (line != null){
			StringTokenizer token = new StringTokenizer(line);
			String name = token.nextToken();
			String thermo = token.nextToken();
			try{
				double H = Double.parseDouble(thermo);
				thermo = thermo.concat(" ");
				for (int i = 0; i < 11; i++){
					thermo = thermo.concat(token.nextToken());
					thermo = thermo.concat(" ");
				}
				ThermoData thermoData = ChemParser.parseThermoFromLibrary(thermo);
				String comments = "";
				while(token.hasMoreTokens()){
					comments = comments + " " + token.nextToken();
				}
				ThermoData newThermoData = new ThermoData(name, thermoData, comments, source + " (Species ID: " + name + ")");
				Graph g = (Graph)dictionary.get(name);
				if (g != null) {
					Object old = library.get(g);
					if (old == null){
						library.put(g, newThermoData);
					}
					else {
						ThermoData oldThermoData = (ThermoData)old;
						if (!oldThermoData.equals(newThermoData)) {
				            System.out.println("Duplicate thermo data (same graph, same name) in " + source);
				            System.out.println("\tIgnoring thermo data for species: " + newThermoData.getName() + " from thermo_library " + source);
				            System.out.println("\tStill storing thermo data for species: " + oldThermoData.getName() + " from thermo_library " + oldThermoData.source);
						}
					}
				}
			}
			catch (NumberFormatException e) {
				Object o = p_dictionary.get(thermo);
				if (o == null) {
					System.out.println(name + ": "+thermo);
				}
			}
			line = ChemParser.readMeaningfulLine(data);
		}
		in.close();
	    return library;
	}
	catch (Exception e){
	   throw new IOException("Can't read thermo in primary thermo library!");
	}
}

//7-Jul-2009: MRH
public HashMap readDictionary(String p_fileName, String source) throws FileNotFoundException, IOException{
	  try{
	    FileReader in = new FileReader(p_fileName);
	    BufferedReader data = new BufferedReader(in);

	    String line = ChemParser.readMeaningfulLine(data);

	    read: while(line != null){
	      StringTokenizer st = new StringTokenizer(line);
	      String name = st.nextToken();

	      data.mark(10000);
	      line = ChemParser.readMeaningfulLine(data);
	      if (line == null) break read;
	      line = line.trim();
	      data.reset();
	      Graph graph = null;

	        graph = ChemParser.readChemGraph(data);
	        graph.addMissingHydrogen();
	        
	      Object old = dictionary.get(name);
	      if (old == null){
	          ThermoData td = (ThermoData)library.get(graph);
	          if (td == null){
	          	dictionary.put(name, graph);
	          } else {
	        	  System.out.println("Ignoring species " + name + 
	        			  " -- Graph already exists in user-defined " + td.source);
	          }
	      }
	      else{
	        Graph oldGraph = (Graph)old;
	        if (!oldGraph.equals(graph)) {
	          System.out.println("Can't replace graph in primary thermo library!");
	          System.out.println("The species name '" + name + "' is given multiple times");
	          System.out.println("in the thermo library '" + source + "', yet has");
	          System.out.println("different graphs");
	          System.exit(0);
	        }

	      }
	      line = ChemParser.readMeaningfulLine(data);
	    }
	    in.close();
	    return dictionary;
	  }
	  catch (FileNotFoundException e){
	      throw new FileNotFoundException(p_fileName);
	  }
	  catch (IOException e){
	    throw new IOException(p_fileName + ":" + e.getMessage());
	  }
}


//## operation getThermoData(ChemGraph)
public ThermoData getThermoData(Graph p_graph){
  //#[ operation getThermoData(ChemGraph)
  ThermoData td = (ThermoData)library.get(p_graph);
  if (td != null){
    return td;
  }
  Iterator iter = library.keySet().iterator();
  while (iter.hasNext()){
    Graph g = (Graph)iter.next();
    g.addMissingHydrogen();
    if (g.isEquivalent(p_graph)){
      td = (ThermoData)library.get(g);
      return td;
    }
  }
  return null;
  //#]
}

//// Added by Amrit Jalan, April 19, 2009
//// Reads in values of S,B,E and V for Platts' groups
//
////## operation getThermoData(ChemGraph)
//public AbramData getAbramData(Graph p_graph){
//  //#[ operation getThermoData(ChemGraph)
//  AbramData td = (AbramData)library.get(p_graph);
//  if (td != null){
//    return td;
//  }
//  Iterator iter = library.keySet().iterator();
//  while (iter.hasNext()){
//    Graph g = (Graph)iter.next();
//    g.addMissingHydrogen();
//    if (g.isEquivalent(p_graph)){
//      td = (AbramData)library.get(g);
//      return td;
//    }
//  }
//  return null;
//  //#]
//}



public void read(String p_dictionary, String p_library) throws IOException, FileNotFoundException {
    dictionary = readDictionary(p_dictionary);
    library = readLibrary(p_library, dictionary);
}

//## operation readDictionary(String)
public HashMap readDictionary(String p_fileName) throws FileNotFoundException, IOException{
  try{
    FileReader in = new FileReader(p_fileName);
    BufferedReader data = new BufferedReader(in);

    String line = ChemParser.readMeaningfulLine(data);

    read: while(line != null){
      StringTokenizer st = new StringTokenizer(line);
      String name = st.nextToken();

      data.mark(10000);
      line = ChemParser.readMeaningfulLine(data);
      if (line == null) break read;
      line = line.trim();
      data.reset();
      Graph graph = null;

        graph = ChemParser.readChemGraph(data);
        graph.addMissingHydrogen();
      Object old = dictionary.get(name);
      if (old == null){
        dictionary.put(name, graph);
      }
      else{
        Graph oldGraph = (Graph)old;
        if (!oldGraph.equals(graph)) {
          System.out.println("Can't replace graph in primary thermo library!");
          System.exit(0);
        }

      }
      line = ChemParser.readMeaningfulLine(data);
    }
    in.close();
    return dictionary;
  }
  catch (FileNotFoundException e){
      throw new FileNotFoundException(p_fileName);
  }
  catch (IOException e){
    throw new IOException(p_fileName + ":" + e.getMessage());
  }
  //#]
}



public HashMap readLibrary(String p_thermoFileName, HashMap p_dictionary) throws IOException {
 try{
  FileReader in = new FileReader(p_thermoFileName);
  BufferedReader data = new BufferedReader(in);
  HashMap library = new HashMap();
  String line = ChemParser.readMeaningfulLine(data);
  while (line != null){
    StringTokenizer token = new StringTokenizer(line);
    String name = token.nextToken();
    String thermo = token.nextToken();
    try{
      double H = Double.parseDouble(thermo);
      thermo = thermo.concat(" ");
      for (int i = 0; i < 11; i++){
        thermo = thermo.concat(token.nextToken());
        thermo = thermo.concat(" ");
      }
        ThermoData thermoData = ChemParser.parseThermoFromLibrary(thermo);
        String comments = "";
        while(token.hasMoreTokens()){
          comments = comments + " " + token.nextToken();
        }
        ThermoData newThermoData = new ThermoData(name, thermoData, comments);
        Graph g = (Graph)dictionary.get(name);
        Object old = library.get(g);
        if (old == null){
          library.put(g, newThermoData);
        }
        else{
          ThermoData oldThermoData = (ThermoData)old;
          if (!oldThermoData.equals(newThermoData)) {
            System.out.println("Can't replace thermoData in primary thermo library!");
            System.out.println("old thermodata:"+oldThermoData.getName());
            System.out.println("new thermodata:"+newThermoData.getName());
            System.exit(0);
          }

        }

      }

    catch (NumberFormatException e){
      Object o = p_dictionary.get(thermo);
      if (o == null){
        System.out.println(name + ": "+thermo);
      }
    }
    line = ChemParser.readMeaningfulLine(data);
  }
  in.close();
    return library;
 }
 catch (Exception e){
   throw new IOException("Can't read thermo in primary thermo library!");
 }
}

protected static PrimaryThermoLibrary getINSTANCE() {
  return INSTANCE;
}


}
/*********************************************************************
        File Path	: RMG\RMG\jing\chem\PrimaryThermoLibrary.java
*********************************************************************/


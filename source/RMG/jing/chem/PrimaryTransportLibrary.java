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

/**
 * PrimaryTransportLibrary.java
 * 	Template was PrimaryThermoLibrary.java file
 * @author MRH (mrharper@mit.edu) 17-May-2010
 *
 */
public class PrimaryTransportLibrary {
	protected static HashMap library;
	protected static HashMap dictionary;
	

	public PrimaryTransportLibrary(){
		library = new HashMap();
		dictionary = new HashMap();
	}

	public PrimaryTransportLibrary(HashMap Dictionary, HashMap Library) {
		dictionary = Dictionary;
		library = Library;
	}

	public PrimaryTransportLibrary(String name, String location) {
		library = new HashMap();
		dictionary = new HashMap();
		appendPrimaryTransportLibrary(name, location);
	}

	public void appendPrimaryTransportLibrary(String name, String path) {
		String transportLibraryDirectory = System.getProperty("jing.chem.TransportReferenceLibrary.pathName");
		String dictionaryFile = transportLibraryDirectory + "/" + path + "/Dictionary.txt";
		String libraryFile = transportLibraryDirectory + "/" + path + "/Library.txt";
		try {
			read(dictionaryFile, libraryFile, name);
		}
		catch (IOException e) {
			System.err.println("RMG cannot read Primary Transport Library: " + name + "\n" + e.getMessage());
		}
	}

	public void read(String p_dictionary, String p_library, String p_name) throws IOException, FileNotFoundException {
		String source = "Primary Transport Library: " + p_name;
		System.out.println("Reading " + source);
		dictionary = readDictionary(p_dictionary, source);
		library = readLibrary(p_library, dictionary, source);
	}

	public HashMap readLibrary(String p_transportFileName, HashMap p_dictionary, String source) throws IOException {
		try{
			FileReader in = new FileReader(p_transportFileName);
			BufferedReader data = new BufferedReader(in);
			HashMap tempLibrary = new HashMap();
			String line = ChemParser.readMeaningfulLine(data);
			while (line != null){
				StringTokenizer token = new StringTokenizer(line);
				String name = token.nextToken();
				String thermo = token.nextToken();
				try{
					TransportData transData = new TransportData(
							Integer.parseInt(thermo),
							Double.parseDouble(token.nextToken()),
							Double.parseDouble(token.nextToken()),
							Double.parseDouble(token.nextToken()),
							Double.parseDouble(token.nextToken()),
							Double.parseDouble(token.nextToken()),
							true);
					String comments = "";
					while (token.hasMoreTokens()) {
						comments = comments + " " + token.nextToken();
					}
					TransportData newTransData = new TransportData(name, transData, comments, source + " (Species ID: " + name + ")",true);
					Graph g = (Graph)dictionary.get(name);
					if (g != null) {
						Object old = tempLibrary.get(g);
						if (old == null){
							tempLibrary.put(g, newTransData);
							library.put(g, newTransData);
						}
						else {
							TransportData oldTransData = (TransportData)old;
							if (!oldTransData.equals(newTransData)) {
					            System.out.println("Duplicate transport data (same graph, different name) in " + source);
					            System.out.println("\tIgnoring thermo data for species: " + newTransData.getName());
					            System.out.println("\tStill storing thermo data for species: " + oldTransData.getName());
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
		   throw new IOException("Can't read transport in primary transport library!");
		}
	}

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
					TransportData td = (TransportData)library.get(graph);
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
						System.out.println("Can't replace graph in primary transport library!");
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


	public static TransportData getTransportData(Graph p_graph){
		if (library == null) return null;
		
		TransportData td = (TransportData)library.get(p_graph);
		if (td != null){
			return td;
		}
		Iterator iter = library.keySet().iterator();
		while (iter.hasNext()){
			Graph g = (Graph)iter.next();
			g.addMissingHydrogen();
			if (g.isEquivalent(p_graph)){
				td = (TransportData)library.get(g);
				return td;
			}
		}
		return null;
	}

}
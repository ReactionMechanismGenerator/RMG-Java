////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
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

/* Created by Amrit Jalan on December 21, 2010
 */

package jing.chem;

import java.io.*;
import java.util.*;

import jing.chemUtil.*;
import jing.chemParser.*;
import jing.rxnSys.Logger;

/**
 *
 * @author amrit
 */
public class PrimaryAbrahamLibrary {
    	protected static HashMap library;
	protected static HashMap dictionary;


	public PrimaryAbrahamLibrary(){
		library = new HashMap();
		dictionary = new HashMap();
	}

	public PrimaryAbrahamLibrary(HashMap Dictionary, HashMap Library) {
		dictionary = Dictionary;
		library = Library;
	}

	public PrimaryAbrahamLibrary(String name, String location) {
		library = new HashMap();
		dictionary = new HashMap();
		appendPrimaryAbrahamLibrary(name, location);
	}

	public void appendPrimaryAbrahamLibrary(String name, String path) {
		String abrahamLibraryDirectory = System.getProperty("jing.chem.ThermoReferenceLibrary.pathName");
		String dictionaryFile = abrahamLibraryDirectory + "/" + path + "/Dictionary.txt";
		String libraryFile = abrahamLibraryDirectory + "/" + path + "/Library.txt";
		try {
			read(dictionaryFile, libraryFile, name);
		}
		catch (IOException e) {
			System.err.println("RMG cannot read Primary Abraham Library: " + name + "\n" + e.getMessage());
		}
	}

	public void read(String p_dictionary, String p_library, String p_name) throws IOException, FileNotFoundException {
		String source = "Primary Abraham Library: " + p_name;
		System.out.println("Reading " + source);
		dictionary = readDictionary(p_dictionary, source);
		library = readLibrary(p_library, dictionary, source);
	}

	public HashMap readLibrary(String p_transportFileName, HashMap p_dictionary, String source) throws IOException {
		try{
			FileReader in = new FileReader(p_transportFileName);
			BufferedReader data = new BufferedReader(in);
			HashMap tempLibrary = new HashMap();
			String line = ChemParser.readMeaningfulLine(data, true);
			while (line != null){
				StringTokenizer token = new StringTokenizer(line);
				String name = token.nextToken();
				String abram = token.nextToken();
				try{
					AbramData abramData = new AbramData(
							Double.parseDouble(abram),
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
					AbramData newAbramData = new AbramData(name, abramData, comments, source + " (Species ID: " + name + ")",true);
					Graph g = (Graph)dictionary.get(name);
					if (g != null) {
						Object old = tempLibrary.get(g);
						if (old == null){
							tempLibrary.put(g, newAbramData);
							library.put(g, newAbramData);
						}
						else {
							AbramData oldAbramData = (AbramData)old;
							if (!oldAbramData.equals(newAbramData)) {
					            System.out.println("Duplicate Abraham data (same graph, different name) in " + source);
					            System.out.println("\tIgnoring Abraham data for species: " + newAbramData.getName());
					            System.out.println("\tStill storing Abraham data for species: " + oldAbramData.getName());
							}
						}
					}
				}
				catch (NumberFormatException e) {
					Object o = p_dictionary.get(abram);
					if (o == null) {
						System.out.println(name + ": "+abram);
					}
				}
				line = ChemParser.readMeaningfulLine(data, true);
			}
			in.close();
		    return library;
		}
		catch (Exception e){
			Logger.logStackTrace(e);
		   throw new IOException("Can't read Abraham descriptors in primary Abraham library!");
		}
	}

	public HashMap readDictionary(String p_fileName, String source) throws FileNotFoundException, IOException{
		try{
			FileReader in = new FileReader(p_fileName);
			BufferedReader data = new BufferedReader(in);

			String line = ChemParser.readMeaningfulLine(data, true);

			read: while(line != null){
				StringTokenizer st = new StringTokenizer(line);
				String name = st.nextToken();

				data.mark(10000);
				line = ChemParser.readMeaningfulLine(data, true);
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
						System.out.println("Can't replace graph in primary Abraham library!");
						System.exit(0);
					}
				}
				line = ChemParser.readMeaningfulLine(data, true);
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


	public static AbramData getAbrahamData(Graph p_graph){
		if (library == null) return null;

		AbramData ab = (AbramData)library.get(p_graph);
		if (ab != null){
			return ab;
		}
		Iterator iter = library.keySet().iterator();
		while (iter.hasNext()){
			Graph g = (Graph)iter.next();
			g.addMissingHydrogen();
			if (g.isEquivalent(p_graph)){
				ab = (AbramData)library.get(g);
				return ab;
			}
		}
		return null;
	}

}

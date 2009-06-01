import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.StringTokenizer;

import jing.chemParser.ChemParser;


public class Mol2AdjList {

	/**
	 * Mol2ChemGraph converts a .mol file to an RMG adjacency list.
	 * 	This program was originally created because the InChI executable was
	 * 	not stable on Linux systems.  Rather than go from the InChI string to
	 * 	the adjacency list in "one step" (via the poorly named InChI2ChemGraph
	 * 	function), the user now supplies the .mol file; creating a .mol file
	 * 	from the InChI on Linux is normally stable.
	 */
	public static void main(String[] args) {
        // Read in the .mol file
        FileReader in = null;
		try {
			in = new FileReader(args[0]);
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
		
		// Extract the information in the first line (Count Line) of the .mol file
		StringTokenizer st = new StringTokenizer(molFileLines[0]);
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
			// Extract radical information
			int numOfRads = Integer.parseInt(st.nextToken());
			for (int j=0; j<numOfRads; j++) {
				int atom = Integer.parseInt(st.nextToken());
				int radType = Integer.parseInt(st.nextToken());
				if (radType == 1)
					adjListRadical[atom-1] = "2S ";
				else if (radType == 2)
					adjListRadical[atom-1] = "1 ";
				else if (radType== 3)
					adjListRadical[atom-1] = "2T ";
				else
					adjListRadical[atom-1] = "3 ";
			}
		}
		
		// Construct the entire adjacency list from its individual lines
		String cgString = "";
		for (int i=0; i<numOfAtoms-1; i++) {
			cgString += (i+1) + adjListElement[i] + adjListRadical[i] + adjListConnectivity[i] + "\n";
		}
		cgString += numOfAtoms + adjListElement[numOfAtoms-1] + adjListRadical[numOfAtoms-1] + adjListConnectivity[numOfAtoms-1];
		
		System.out.println(cgString);
	}

}

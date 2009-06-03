import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.StringTokenizer;

import jing.chemParser.ChemParser;

//	InChI2AdjList.java
//	Written by: Michael Harper (mrharper@mit.edu)
//	Last modified: 26-May-2009
//	This function converts a list of InChI strings to a list
//		of adjacency lists.  The user supplies the location and
//		name of the .txt file containing the list of InChI strings.
//		The output of this function is a file named adjList_output.txt
//		which holds a list of adjacency lists, corresponding to the
//		list of InChI strings.
//	NOTE: The InChI executable must be located in the $RMG/software/InChI
//		folder AND must be version 1.01.
//	Inputs:
//		- .txt file (of any name) containing list of InChIs
//	Output:
//		- .txt file (named adjList_output) containing list of
//			adjacency lists
//
//	3-Jun-2009 (MRH)
//	Renamed from InChI2ChemGraph.java (since the output is not a ChemGraph,
//		but a string in the form of an adjacency list)

public class InChI2AdjList {
	public static void main(String[] args) {
		
		// This string will hold the InChI strings and their respective
		//	adjacency lists
		String listOfAdjLists = "";
		
		// In case a folder named InChI is not already present, make it
		File inchi = new File("InChI");
		ChemParser.deleteDir(inchi);
		inchi.mkdir();
		
		// Read in the input file
        FileReader input_file = null;
		try {
			input_file = new FileReader(args[0]);
		} catch (FileNotFoundException e) {
			String err = "Error reading input file to InChI2ChemGraph: " +
				args[0] + "\n" + e.toString();
			System.out.println(err);
		}
        
		// Read in the list of InChIs, line-by-line
		BufferedReader input_reader = new BufferedReader(input_file);
		String p_inchi = ChemParser.readMeaningfulLine(input_reader);
		
		// While more InChIs remain
		while (p_inchi != null) {
		
			String workingDirectory = System.getenv("RMG");
			File inchi_file = null;
			
			// Save InChI string in inchi.txt file
	        try {
	        	inchi_file = new File(inchi + "/inchi.txt");
	        	FileWriter fw = new FileWriter(inchi_file);
	        	fw.write(p_inchi);
	        	fw.close();
	        } catch (IOException e) {
	        	String err = "Error writing inchi.txt file for InChI-to-molFile conversion: " + e.toString();
	        	System.out.println(err);
	        }
			
	        // Call cINChI-1 executable file
	        //	The call (more specifically, the way to pass options)
	        //		depends on the operating system
	        String[] optionsArgument = new String[2];
	        if (getOs().toLowerCase().equals("windows")) {
	        	optionsArgument[0] = "/InChI2Struct";
	        	optionsArgument[1] = "/OutputSDF";
	        } else if (getOs().toLowerCase().equals("linux")) {
	        	optionsArgument[0] = "-InChI2Struct";
	        	optionsArgument[1] = "-OutputSDF";
	        } else if (getOs().toLowerCase().equals("mac")) {
	        	optionsArgument[0] = "-InChI2Struct";
	        	optionsArgument[1] = "-OutputSDF";
	        }
	
	        int exitValue = -1;
	        try {
	            String[] command = {workingDirectory + "/software/InChI/cInChI-1",
	                    "inchi.txt",
	                    "temp.txt",
	                    optionsArgument[0]};
	            File runningDir = new File("InChI");
	            Process InChI = Runtime.getRuntime().exec(command, null, runningDir);
	
	            InputStream errStream = InChI.getErrorStream();
	            InputStream inpStream = InChI.getInputStream();
	            errStream.close();
	            inpStream.close();
	
	            exitValue = InChI.waitFor();
	        }
	        catch (Exception e) {
	            String err = "Error running cINChI-1 while converting InChI string " +
	            	"to temporary file: " + e.toString();
	            System.out.println(err);
	        }

	        try {
	            String[] command = {workingDirectory + "/software/InChI/cInChI-1",
	                    "temp.txt",
	                    "temp.mol",
	                    optionsArgument[1]};
	            File runningDir = new File("InChI");
	            Process InChI = Runtime.getRuntime().exec(command, null, runningDir);
	
	            InputStream errStream = InChI.getErrorStream();
	            InputStream inpStream = InChI.getInputStream();
	            errStream.close();
	            inpStream.close();
	
	            exitValue = InChI.waitFor();
	        }
	        catch (Exception e) {
	            String err = "Error running cINChI-1 while converting temporary " +
	            	"file to .mol file: " + e.toString();
	            System.out.println(err);
	        }
		
	        // Read in the .mol file
	        FileReader in = null;
			try {
				in = new FileReader(inchi + "/temp.mol");
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
				cgString += (i+1) + adjListElement[i] + adjListRadical[i] + adjListConnectivity[i] + "\r";
			}
			cgString += numOfAtoms + adjListElement[numOfAtoms-1] + adjListRadical[numOfAtoms-1] + adjListConnectivity[numOfAtoms-1];
			
			listOfAdjLists += p_inchi + "\n" + cgString + "\n\n";
			
			p_inchi = ChemParser.readMeaningfulLine(input_reader);
		}
		
		// Write the output file adjList_output.txt
		try {
            File adjList = new File("adjList_output.txt");
            FileWriter fw = new FileWriter(adjList);
            fw.write(listOfAdjLists);
            fw.close();
            System.out.println("Program complete.  File adjList_output.txt written successfully.");
		} catch (IOException e) {
			System.out.println("Error writing adjList_output.txt file: " + e.toString());
			System.exit(0);
		}
			
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
}
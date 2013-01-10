import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;
import jing.chem.Species;
import jing.chemParser.ChemParser;

// InChI2AdjList.java
// Written by: Michael Harper (mrharper@mit.edu)
// Last modified: 26-May-2009
// This function converts a list of InChI strings to a list
// of adjacency lists. The user supplies the location and
// name of the .txt file containing the list of InChI strings.
// The output of this function is a file named adjList_output.txt
// which holds a list of adjacency lists, corresponding to the
// list of InChI strings.
// NOTE: The InChI executable must be located in the $RMG/bin
// folder AND must be version 1.01.
// Inputs:
// - .txt file (of any name) containing list of InChIs
// Output:
// - .txt file (named adjList_output) containing list of
// adjacency lists
//
// 3-Jun-2009 (MRH)
// Renamed from InChI2ChemGraph.java (since the output is not a ChemGraph,
// but a string in the form of an adjacency list)
//
// 16-Jun-2009 (MRH)
// Have module call the Species.inchi2AdjList function. This ensures there
// is only one copy of the code.
public class InChI2AdjList {
    public static void main(String[] args) {
        String inchiDirectory = System
                .getProperty("RMG.InChI_running_directory");
        // This string will hold the InChI strings and their respective
        // adjacency lists
        String listOfAdjLists = "";
        // In case a folder named InChI is not already present, make it
        File inchi = new File(inchiDirectory);
        ChemParser.deleteDir(inchi);
        inchi.mkdir();
        // Read in the input file
        FileReader input_file = null;
        try {
            input_file = new FileReader(args[0]);
        } catch (FileNotFoundException e) {
            String err = "Error reading input file to InChI2AdjList: "
                    + args[0] + "\n" + e.toString();
            System.out.println(err);
        }
        // Read in the list of InChIs, line-by-line
        BufferedReader input_reader = new BufferedReader(input_file);
        String p_inchi = ChemParser.readMeaningfulLine(input_reader, true);
        // While more InChIs remain
        while (p_inchi != null) {
            // Convert the InChI to an adjacency list
            // Write the InChI and adjacency list to the listOfAdjList string
            StringTokenizer st = new StringTokenizer(p_inchi);
            String name = st.nextToken();
            listOfAdjLists += name + "\n";
            if (st.hasMoreTokens())
                listOfAdjLists += Species.inchi2AdjList(st.nextToken())
                        + "\n\n";
            else
                listOfAdjLists += Species.inchi2AdjList(name) + "\n\n";
            p_inchi = ChemParser.readMeaningfulLine(input_reader, true);
        }
        // Write the output file adjList_output.txt
        try {
            File adjList = new File("adjList_output.txt");
            FileWriter fw = new FileWriter(adjList);
            fw.write(listOfAdjLists);
            fw.close();
            System.out
                    .println("Program complete.  File adjList_output.txt written successfully.");
        } catch (IOException e) {
            System.out.println("Error writing adjList_output.txt file: "
                    + e.toString());
            System.exit(0);
        }
    }
}
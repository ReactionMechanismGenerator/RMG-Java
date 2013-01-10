import java.io.File;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.IOException;
import jing.chem.Species;

public class Mol2AdjList {
    /**
     * Mol2ChemGraph converts a .mol file to an RMG adjacency list. This program was originally created because the
     * InChI executable was not stable on Linux systems. Rather than go from the InChI string to the adjacency list in
     * "one step" (via the poorly named InChI2ChemGraph function), the user now supplies the .mol file; creating a .mol
     * file from the InChI on Linux is normally stable.
     */
    public static void main(String[] args) {
        // String will contain the list of chemgraphs
        String allChemGraphs = "";
        String singlespeciesChemGraphs = "";
        // Assume the input is the directory containing the list of .mol files
        File folder = new File(args[0]);
        File[] listOfFiles = folder.listFiles();
        for (int numFiles = 0; numFiles < listOfFiles.length; numFiles++) {
            String filename = listOfFiles[numFiles].getAbsolutePath();
            String[] directories = filename.split("\\\\");
            String speciesName = directories[directories.length - 1];
            int speciesId = numFiles + 1;
            speciesName = speciesName.substring(0, speciesName.length() - 4);
            singlespeciesChemGraphs = speciesName + "\n"
                    + Species.mol2AdjList(speciesName, filename) + "\n";
            // allChemGraphs += Species.mol2AdjList(speciesName, filename);
            allChemGraphs += singlespeciesChemGraphs + "\n";
        }
        // String cgString = Species.mol2AdjList(args[0], args[0]);
        // System.out.println(cgString);
        System.out.println(allChemGraphs);
        // Write the output file chemgraphs_output.txt
        try {
            File ChemGraphs = new File("chemgraphs_output.txt");
            FileWriter fw = new FileWriter(ChemGraphs);
            fw.write(allChemGraphs);
            fw.close();
            System.out
                    .println("Program complete.  File chemgraphs_output.txt written successfully.");
        } catch (IOException e) {
            System.out.println("Error writing chemgraphs_output.txt file: "
                    + e.toString());
            System.exit(0);
        }
    }
}

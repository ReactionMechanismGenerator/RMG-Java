import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import jing.chem.ChemGraph;
import jing.chem.ForbiddenStructureException;
import jing.chem.Species;
import jing.chemParser.ChemParser;
import jing.chemUtil.Graph;

public class AdjList2InChI {
    public static void main(String[] args) {
// initializeSystemProperties();
        RMG.globalInitializeSystemProperties();
        String outputString = "";
        try {
            FileReader in = new FileReader(args[0]);
            BufferedReader data = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(data, true);
            while (line != null) {
                String speciesName = line;
                Graph g = ChemParser.readChemGraph(data);
                ChemGraph cg = null;
                try {
                    cg = ChemGraph.make(g);
                } catch (ForbiddenStructureException e) {
                    System.out
                            .println("Error in reading graph: Graph contains a forbidden structure.\n"
                                    + g.toString());
                    System.exit(0);
                }
                String[] inchiStrings = Species.generateInChI(cg);
                outputString += speciesName + "\t" + inchiStrings[0] + "\n";
                line = ChemParser.readMeaningfulLine(data, true);
            }
            in.close();
        } catch (FileNotFoundException e) {
            System.err.println("File was not found!\n");
        } catch (IOException e) {
            System.err.println("Something wrong with ChemParser.readChemGraph");
        }
        try {
            File inchiOutput = new File("inchi_output.txt");
            FileWriter fw = new FileWriter(inchiOutput);
            fw.write(outputString);
            fw.close();
        } catch (IOException e) {
            System.out.println("Error in writing inchi_output.txt file.");
            System.exit(0);
        }
        System.out.println("Results written in inchi_output.txt");
    };
}
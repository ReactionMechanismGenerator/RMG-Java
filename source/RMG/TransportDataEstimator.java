import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import jing.chem.ChemGraph;
import jing.chem.ForbiddenStructureException;
import jing.chem.InvalidChemGraphException;
import jing.chem.Species;
import jing.chem.TransportData;
import jing.chemParser.ChemParser;
import jing.rxnSys.ReactionModelGenerator;

public class TransportDataEstimator {
    public static void main(String[] args) {
        RMG.globalInitializeSystemProperties();
        /*
         * When a species is made (as will be done in this module), the thermo is calculated by calling GATPFit. The
         * .exe expects a working directory named GATPFit to be present in the current directory.
         */
        File gatpfit = new File(System.getProperty("RMG.GATPFitDir"));
        gatpfit.mkdir();
        String transportProperties = "";
        LinkedHashMap speciesFromInputFile = new LinkedHashMap();
        // Read in the file contained in args[0] string
        try {
            BufferedReader reader = new BufferedReader(new FileReader(new File(
                    args[0])));
            ReactionModelGenerator rmg = new ReactionModelGenerator();
            String line = ChemParser.readMeaningfulLine(reader, true);
            // Read in the Database field
            if (line.toLowerCase().startsWith("database")) {
                RMG.extractAndSetDatabasePath(line);
            } else {
                System.err.println("TransportDataEstimator: Could not"
                        + " locate the Database field");
                System.exit(0);
            }
            line = ChemParser.readMeaningfulLine(reader, true);
            // Read in the Primary Transport Library, if it exists
            if (line.toLowerCase().startsWith("primarytransportlibrary")) {
                rmg.readAndMakePTransL(reader);
            } else {
                System.err
                        .println("TransportDataEstimator: Could not locate the PrimaryTransportLibrary field."
                                + "Line read was: " + line);
                System.exit(0);
            }
            line = ChemParser.readMeaningfulLine(reader, true);
            while (line != null) {
                String speciesName = line;
                ChemGraph cg = ChemGraph.make(ChemParser.readChemGraph(reader));
                ReactionModelGenerator
                        .addChemGraphToListIfNotPresent_ElseTerminate(
                                speciesFromInputFile, cg, "");
                Species sp = Species.make(speciesName, cg);
                TransportData lj4species = sp.getChemkinTransportData();
                String whitespace = "                ";
                // Write the 6 transport properties
                if (speciesName.length() > 16) {
                    System.out
                            .println("!Species name contains more than 16 characters: "
                                    + speciesName
                                    + "\n!CHEMKIM Pre-Processor will most likely throw an error.");
                    transportProperties += speciesName + "   "
                            + lj4species.toString();
                } else
                    transportProperties += speciesName
                            + whitespace.substring(speciesName.length())
                            + "   " + lj4species.toString();
                transportProperties += " ! " + lj4species.getSource() + "\t"
                        + lj4species.getComment() + "\n";
                // (Attempt to) Read next line of input file
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } catch (InvalidChemGraphException e) {
            e.printStackTrace();
        } catch (ForbiddenStructureException e) {
            e.printStackTrace();
        } catch (IOException e) {
            System.out.println(e.toString());
        }
        try {
            File trandat = new File("tran.dat");
            FileWriter fw = new FileWriter(trandat);
            fw.write(transportProperties);
            fw.close();
            System.out.println("Results written to tran.dat");
        } catch (IOException e) {
            System.out.println("Could not write tran.dat");
            System.exit(0);
        }
    }
}
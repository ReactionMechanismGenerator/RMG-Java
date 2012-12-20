// this utility adds (modified) inchis as comments to Chemkin files
/**
 * @author gmagoon
 */
import java.util.*;
import java.io.*;
import jing.chem.*;
import jing.chemParser.*;
import jing.chemUtil.*;

// arg[0]=dictionary file
// arg[1]=original chemkin file
// arg[2]=new chemkin file
public class InchifyCHEMKINusingDictionary {
    public static void main(String[] args) {
        RMG.globalInitializeSystemProperties();
        // 1. read the Dictionary file, converting to modified inchis along the way and adding them into inchiDict
// LinkedHashMap
        LinkedHashMap inchiDict = new LinkedHashMap();
        try {
            FileReader in = new FileReader(args[0]);
            BufferedReader data = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(data, true);
            // While more species
            while (line != null) {
                System.out.println(line);// print the name of the molecule
                Graph g = ChemParser.readChemGraph(data);
                System.out.println(g);
                ChemGraph chemgraph = ChemGraph.make(g);
                // get the inchi
                String inchiString = chemgraph.getModifiedInChIAnew();
                // get the molecule name; ***this assumes that the names in the dictionary are the same as the names in
// the chemkin file; additional considerations would be needed if the CHEMKIN names are shortened names of the
// dictionary names***
                String name = line.split("\\s+")[0];// split on whitespace and take the first part of the name (in case
// the second part happens to be a pre-existing InChI
                inchiDict.put(name, inchiString);// add the chemkin name as key, and InChI as value
                line = ChemParser.readMeaningfulLine(data, true);
            }
            in.close();
        } catch (FileNotFoundException e) {
            System.err.println("File was not found!\n");
        } catch (IOException e) {
            System.err.println("Something wrong with ChemParser.readChemGraph");
        } catch (ForbiddenStructureException e) {
            System.err.println("Something wrong with ChemGraph.make");
        }
        // 2. read the CHEMKIN file, while simultaneously writing a CHEMKIN output file with InChI comments, where
// available
        try {
            FileReader chemkinOld = new FileReader(args[1]);
            FileWriter chemkinNew = new FileWriter(args[2]);
            BufferedReader chemkinReader = new BufferedReader(chemkinOld);
            BufferedWriter chemkinWriter = new BufferedWriter(chemkinNew);
            String line = chemkinReader.readLine();
            // While more InChIs remain
            while (line != null) {
                String[] spl = line.split(" "); // split on space
                int len = spl.length;
                // if the last "word" is a 1 and the first word can be found in the dictionary, write a comment line
// before copying the line; otherwise, just copy the line
                if (spl[len - 1].equals("1")) {
                    if (inchiDict.containsKey(spl[0])) {
                        chemkinWriter.write("! [_ SMILES=\""
                                + inchiDict.get(spl[0]) + "\" _]\n");
                    } else {
                        System.out
                                .println("Warning: Could not find InChI corresponding to CHEMKIN file species "
                                        + spl[0]);
                    }
                }
                chemkinWriter.write(line + "\n");
                line = chemkinReader.readLine();
            }
            chemkinOld.close();
            chemkinWriter.flush();
            chemkinWriter.close();
            chemkinNew.close();
        } catch (FileNotFoundException e) {
            System.err.println("File was not found!\n");
        } catch (IOException eIO) {
            System.err.println("IO Exception!\n");
        }
    }
}

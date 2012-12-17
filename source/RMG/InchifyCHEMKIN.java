// this utility adds inchis as comments to Chemkin files
/**
 * @author gmagoon
 */
import java.util.*;
import java.io.*;

public class InchifyCHEMKIN {
    public static void main(String[] args) {
        // 1. read the inchis into LinkedHashMap
        LinkedHashMap inchiDict = new LinkedHashMap();
        try {
            FileReader inchiDictionary = new FileReader(args[0]);
            BufferedReader inchiDictionaryBR = new BufferedReader(
                    inchiDictionary);
            String line = inchiDictionaryBR.readLine();
            // While more InChIs remain
            while (line != null) {
                String[] spl = line.split("\t"); // split on tab
                inchiDict.put(spl[0], spl[2]);// add the chemkin name as key, and InChI as value
                line = inchiDictionaryBR.readLine();
            }
            inchiDictionary.close();
        } catch (FileNotFoundException e) {
            System.err.println("File was not found!\n");
        } catch (IOException eIO) {
            System.err.println("IO Exception!\n");
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
                if (spl[len - 1].equals("1") && inchiDict.containsKey(spl[0])) {
                    chemkinWriter.write("! [_ SMILES=\""
                            + inchiDict.get(spl[0]) + "\" _]\n");
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
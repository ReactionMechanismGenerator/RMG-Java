package applications;
//this utility reads the inchiDictionary.txt file and computes thermo using QMTP without HBI/ring-only options

/**
 *
 * @author gmagoon
 */
 
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;

import jing.chem.Species;
import jing.chemParser.ChemParser;
import java.util.*;
import java.io.*;
import jing.chem.*;
import jing.chemParser.*;
import jing.param.*;
import jing.chemUtil.*;
//import bondGroups.*;
import jing.rxn.*;
import jing.rxnSys.*;
import jing.mathTool.*;

public class DictionaryReader {

    public static void main(String[] args) {
         String workingDir = System.getenv("RMG");
         System.setProperty("RMG.workingDirectory", workingDir);
         System.setProperty("jing.chem.LJDatabase.pathName", workingDir + "/databases/RMG_database/transport_groups");
         
         //block below was commented out and name was changed from inchiDictionaryReader to DictionaryReader after I realized biradicals on adjacent were causing difficulty
         
        //1. read the dictionary and create the adjacency lists (partially copied from InChI2AdjList)
//         String listOfAdjLists ="";
//         String name= "RMG_database";
//	 String workingDir = System.getenv("RMG");
//         System.setProperty("RMG.workingDirectory", workingDir);
//         // Read in the input file
//         FileReader input_file = null;
//	try {
//            input_file = new FileReader("inchiDictionary.txt");
//	} catch (FileNotFoundException e) {
//            String err = "Error reading input file to InChIDictionaryReader" + "\n" + e.toString();
//            System.out.println(err);
//	}
//        
//        // Read in the list of InChIs, line-by-line
//	BufferedReader input_reader = new BufferedReader(input_file);
//	String line = ChemParser.readMeaningfulLine(input_reader);	
//        // While more InChIs remain
//	while (line != null) {
//            //	Convert the InChI to an adjacency list
//            //	Write the InChI and adjacency list to the listOfAdjList string
//            String [] splitString = line.trim().split("\\s+");
//            String inChI = splitString[2];
//            String validInChI= inChI.replaceAll("/mult\\d+", "");
//            listOfAdjLists += inChI + "\n";
//            listOfAdjLists += Species.inchi2AdjList(validInChI) + "\n\n";
//            line = ChemParser.readMeaningfulLine(input_reader);
//        }
//        // Write the output file adjList_output.txt
//	try {
//            File adjList = new File("inchiDictionaryAdjList.txt");
//            FileWriter fw = new FileWriter(adjList);
//            fw.write(listOfAdjLists);
//            fw.close();
//            System.out.println("Program complete.  File inchiDictionaryAdjList.txt written successfully.");
//	} catch (IOException e) {
//                System.out.println("Error writing inchiDictionaryAdjList.txt file: " + e.toString());
//            System.exit(0);
//	}

        ChemGraph.useQM = true;
        ChemGraph.useQMonCyclicsOnly = false;
        Global.maxRadNumForQM = 9999;
        QMTP.qmfolder = "DictionaryQMFiles/";
        //QMTP.qmprogram = "gaussian03";
        QMTP.qmprogram = "mopac";
        QMTP.usePolar = true;

        LinkedHashMap speciesFromInputFile = new LinkedHashMap();
        
         //2. calculate the thermo using QMTP with no HBI; include non-ring species; (primaryThermoLibrary species will still not be included); for triplets, use guess=mix and use gaussian?
          try {
            //FileReader in = new FileReader("inchiDictionaryAdjList.txt");
            FileReader in = new FileReader("RMG_Dictionary.txt");
            BufferedReader data = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(data, true);	
            // While more InChIs remain
            while (line != null) {
                System.out.println(line);//print the name of the molecule
                Graph g = ChemParser.readChemGraph(data);
                System.out.println(g);
                if(!line.equals("InChI=1/H")&&!line.startsWith("HJ(")&&!line.startsWith("H(")){//{for some reason, H does not seem to work in Gaussian, even manually, without freq keyword; not sure about why MOPAC fails
                    ChemGraph chemgraph = ChemGraph.make(g);

                    ReactionModelGenerator.addChemGraphToListIfNotPresent_ElseTerminate(speciesFromInputFile,chemgraph,"");

                    Species spe = Species.make("molecule",chemgraph);
                    //System.out.println(spe.getName());
                            
                    //calculate and display Lennard-Jones estimates based on Joback correlations for critical properties
                    GATransportP LJGAPP=GATransportP.getINSTANCE();
                    TransportData ljdata = LJGAPP.generateTransportData(chemgraph);
                    System.out.println("Sums (dTc, dPc, dVc, dTb): " + ljdata.getGAData().deltaToString());
                    System.out.println("Tc (K) = " + ljdata.getGAData().calculateTc());
                    System.out.println("Pc (bar) = " + ljdata.getGAData().calculatePc());
                    System.out.println("Vc (cc/mol) = " + ljdata.getGAData().calculateVc());
                    System.out.println("Tb (K) = " + ljdata.getGAData().calculateTb());
                    System.out.println("omega = " + ljdata.getGAData().calculateOmega());
                    System.out.println("sigma (Angstroms) = " + ljdata.getSigma());
                    System.out.println("epsilon (K) = " + ljdata.getEpsilon());
                    System.out.println("\n\n");
                     /*
                      *	Following line added by MRH on 10Aug2009:
                      *		After the species is made, the chemgraph is not necessarily the same
                      *		(e.g. the species contains resonance isomers, and the adjacency list
                      *		supplied by the user is not the most stable isomer).  This was causing
                      *		a discrepancy to be displayed to the screen: the "H=" value we were
                      *		displaying corresponded to the chemgraph supplied by the user; the 
                      *		ThermData corresponded to the most stable isomer.
                      *
                      *		An example of a troublesome chemgraph:
                      *		1 C 0 {2,T}
                      *		2 C 0 {1,T} {3,S}
                      *		3 O 0 {2,S}
                      */

                    // chemgraph = spe.getChemGraph();
                    // chemgraph.getThermoData();
                }
                line = ChemParser.readMeaningfulLine(data, true);
            }
            in.close();
 }
 catch (FileNotFoundException e) {
   System.err.println("File was not found!\n");
 }
 catch(IOException e){
   System.err.println("Something wrong with ChemParser.readChemGraph");
 }
 catch(ForbiddenStructureException e){
   System.err.println("Something wrong with ChemGraph.make");
 }
        
         //after the run, unoptimized 3D-mol files should be present in 3Dmolfiles directory and Gaussian results should be in InChIDictionaryQMFiles folder
    }
    
}

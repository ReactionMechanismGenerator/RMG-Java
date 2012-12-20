// //////////////////////////////////////////////////////////////////////////////
//
// RMG - Reaction Mechanism Generator
//
// Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
// RMG Team (rmg_dev@mit.edu)
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// //////////////////////////////////////////////////////////////////////////////
import java.util.*;
import java.io.*;
import qm.QMFlags;
import jing.chem.*;
import jing.chemParser.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.rxnSys.*;

/**
 * This {@link ThermoDataEstimator} extends the functionalities of the previous <code>ThermoDataEstimator</code> by
 * allowing thermodynamic properties of chemical species to generated using quantum-mechanical methods.<BR>
 * <BR>
 * It does so by reading in a few additional lines in the input file: <LI>QM active flag <LI>QM method <LI>QM for
 * cyclics only flag <LI>max radical number for QM integer value
 */
public class ThermoDataEstimator {
    /**
     * @param args
     *            filename of input file
     */
    public static void main(String[] args) {
        RMG.globalInitializeSystemProperties();
        createFolders();
        Species.setAddID(false);// don't add IDs to species names
        /**
         * TODO program against interfaces, not implementations!
         */
        LinkedHashMap speciesFromInputFile = new LinkedHashMap();
        Map<ChemGraph, String> mappedChemGraphsToNames = null;
        QMFlags qmflags = null;
        try {
            File file = new File(args[0]);
            BufferedReader reader = new BufferedReader(new FileReader(file));
            readDatabasePath(reader);
            qmflags = readQMFlags(reader);
            Global.maxRadNumForQM = qmflags.maxRadNumForQM.intValue();
            ChemGraph.TDMETHOD = qmflags.TDSTRATEGY;
            QMTP.qmprogram = qmflags.method.toLowerCase();
            QMTP.connectivityCheck = qmflags.connectivityCheck.intValue();
            readAtomConstraints(reader);
            readPrimaryThermoLibrary(reader);
            mappedChemGraphsToNames = readChemGraphsFromFile(
                    speciesFromInputFile, reader);
        } catch (InvalidChemGraphException e) {
            e.printStackTrace();
        } catch (ForbiddenStructureException e) {
            e.printStackTrace();
        } catch (IOException e) {
            System.out.println(e.toString());
        }
        generateTDProperties(mappedChemGraphsToNames);
        Logger.info("Done!\n");
    }

    /**
     * dynamic constraints setting on the number of atoms (C,O,heavy,etc) in the species
     */
    private static void readAtomConstraints(BufferedReader reader) {
        String line;
        line = ChemParser.readMeaningfulLine(reader, true);
        while (!line.equals("END")) {
            if (line.startsWith("MaxCarbonNumber")) {
                StringTokenizer st = new StringTokenizer(line);
                String dummyString = st.nextToken();
                int maxCNum = Integer.parseInt(st.nextToken());
                ChemGraph.setMaxCarbonNumber(maxCNum);
                Logger.info("Note: Overriding default MAX_CARBON_NUM with user-defined value: "
                        + maxCNum);
            } else if (line.startsWith("MaxOxygenNumber")) {
                StringTokenizer st = new StringTokenizer(line);
                String dummyString = st.nextToken(); // This should hold "MaxOxygenNumberPerSpecies:"
                int maxONum = Integer.parseInt(st.nextToken());
                ChemGraph.setMaxOxygenNumber(maxONum);
                Logger.info("Note: Overriding default MAX_OXYGEN_NUM with user-defined value: "
                        + maxONum);
            } else if (line.startsWith("MaxRadicalNumber")) {
                StringTokenizer st = new StringTokenizer(line);
                String dummyString = st.nextToken(); // This should hold "MaxRadicalNumberPerSpecies:"
                int maxRadNum = Integer.parseInt(st.nextToken());
                ChemGraph.setMaxRadicalNumber(maxRadNum);
                Logger.info("Note: Overriding default MAX_RADICAL_NUM with user-defined value: "
                        + maxRadNum);
            } else if (line.startsWith("MaxSulfurNumber")) {
                StringTokenizer st = new StringTokenizer(line);
                String dummyString = st.nextToken(); // This should hold "MaxSulfurNumberPerSpecies:"
                int maxSNum = Integer.parseInt(st.nextToken());
                ChemGraph.setMaxSulfurNumber(maxSNum);
                Logger.info("Note: Overriding default MAX_SULFUR_NUM with user-defined value: "
                        + maxSNum);
            } else if (line.startsWith("MaxSiliconNumber")) {
                StringTokenizer st = new StringTokenizer(line);
                String dummyString = st.nextToken(); // This should hold "MaxSiliconNumberPerSpecies:"
                int maxSiNum = Integer.parseInt(st.nextToken());
                ChemGraph.setMaxSiliconNumber(maxSiNum);
                Logger.info("Note: Overriding default MAX_SILICON_NUM with user-defined value: "
                        + maxSiNum);
            } else if (line.startsWith("MaxHeavyAtom")) {
                StringTokenizer st = new StringTokenizer(line);
                String dummyString = st.nextToken(); // This should hold "MaxHeavyAtomPerSpecies:"
                int maxHANum = Integer.parseInt(st.nextToken());
                ChemGraph.setMaxHeavyAtomNumber(maxHANum);
                Logger.info("Note: Overriding default MAX_HEAVYATOM_NUM with user-defined value: "
                        + maxHANum);
            } else if (line.startsWith("MaxCycleNumber")) {
                StringTokenizer st = new StringTokenizer(line);
                String dummyString = st.nextToken(); // This should hold "MaxCycleNumberPerSpecies:"
                int maxCycleNum = Integer.parseInt(st.nextToken());
                ChemGraph.setMaxCycleNumber(maxCycleNum);
                Logger.info("Note: Overriding default MAX_CYCLE_NUM with user-defined value: "
                        + maxCycleNum);
            }
            line = ChemParser.readMeaningfulLine(reader, true);
        }
    }

    /**
     * method that iterates over all read-in ChemGraph's, generates TD properties for all of them, and writes properties
     * to the logger also, writes a chemkin format file named TDEresultsCHEMKIN.dat with the thermo for all the species
     * 
     * @param mappedChemGraphsToNames
     *            map with Chemgraphs and mapped names
     */
    private static void generateTDProperties(
            Map<ChemGraph, String> mappedChemGraphsToNames) {
        try {
            // setup the chemkin output file
            File chemkinFile = new File("TDEresultsCHEMKIN.dat");
            FileWriter fw = new FileWriter(chemkinFile);
            BufferedWriter bw = new BufferedWriter(fw);
            // iterate through all the species
            for (ChemGraph chemgraph : mappedChemGraphsToNames.keySet()) {
                Species spe = Species.make(
                        mappedChemGraphsToNames.get(chemgraph), chemgraph);
                ChemGraph stableChemGraph = spe.getChemGraph();
                writeThermoDataInfo(spe, stableChemGraph);
                fw.write(getChemkinString(spe, stableChemGraph) + "\n");// write to the Chemkin file
            }
            fw.close();// close the chemkin file
        } catch (IOException e) {
            Logger.error("Problem creating, writing, or closing Chemkin file!");
            Logger.logStackTrace(e);
        }
    }

    /**
     * Create the working folders for QMTP required folders
     */
    private static void createFolders() {
        createFolder(System.getProperty("RMG.GATPFitDir"), true);
        createFolder(System.getProperty("RMG.InChI_running_directory"), true);
        createFolder(System.getProperty("RMG.2DmolfilesDir"), true);
        createFolder(System.getProperty("RMG.3DmolfilesDir"), true);
        createFolder(System.getProperty("RMG.qmCalculationsDir"), false); // Preserving QM files between runs will speed
// things up considerably
        createFolder(System.getProperty("RMG.qmLibraryDir"), false);
    }

    /**
     * reader method for reading all QM flags
     * 
     * @param reader
     * @return
     */
    private static QMFlags readQMFlags(BufferedReader reader) {
        QMFlags qmFlags = new QMFlags();
        String line = ChemParser.readMeaningfulLine(reader, true);
        qmFlags.TDSTRATEGY = line;
        if (!qmFlags.TDSTRATEGY.equalsIgnoreCase("BensonOnly")) {
            line = ChemParser.readMeaningfulLine(reader, true);
            qmFlags.method = line;
            line = ChemParser.readMeaningfulLine(reader, true);
            qmFlags.maxRadNumForQM = Integer.parseInt(line);
            line = ChemParser.readMeaningfulLine(reader, true);
            String checkConnSetting = line.toLowerCase();
            if (checkConnSetting.equals("off")) {// no connectivity checking
                qmFlags.connectivityCheck = 0;
            } else if (checkConnSetting.equals("check")) {// print a warning if the connectivity doesn't appear to match
                qmFlags.connectivityCheck = 1;
            } else if (checkConnSetting.equals("confirm")) {// consider the run a failure if the connectivity doesn't appear
    // to match
                qmFlags.connectivityCheck = 2;
            } else {
                Logger.critical("input.txt: Inappropriate 'CheckConnectivity' value (should be 'off', 'check', or 'confirm')");
                System.exit(0);
            }
        } else {
            // Dummy values; not actually used by Benson GA method
            qmFlags.method = "";
            qmFlags.maxRadNumForQM = 0;
            qmFlags.connectivityCheck = 0;
        }
        return qmFlags;
    }

    private static void writeThermoDataInfo(Species spe,
            ChemGraph stableChemGraph) {
        Logger.info(stableChemGraph.toString());
        Logger.info("The number of resonance isomers is "
                + spe.getResonanceIsomersHashSet().size());
        Logger.info("The NASA data is \n"
                + getChemkinString(spe, stableChemGraph));
        Logger.info("ThermoData is \n"
                + stableChemGraph.getThermoData().toString());
        int symm = stableChemGraph.getSymmetryNumber();
        Logger.info(" symmetry number = " + symm);
        Temperature T = new Temperature(298.0, "K");
        String chemicalFormula = stableChemGraph.getChemicalFormula();
        Logger.info(chemicalFormula + "  H = " + stableChemGraph.calculateH(T));
        Logger.info("");
    }

    private static Map<ChemGraph, String> readChemGraphsFromFile(
            LinkedHashMap speciesFromInputFile, BufferedReader reader)
            throws IOException, ForbiddenStructureException {
        Map<ChemGraph, String> chemgraphNamesMap = new LinkedHashMap<ChemGraph, String>();
        String line = ChemParser.readMeaningfulLine(reader, true);
        while (line != null) {
            String name = line;
            Graph g = ChemParser.readChemGraph(reader);
            ChemGraph cg = ChemGraph.make(g);
            ReactionModelGenerator
                    .addChemGraphToListIfNotPresent_ElseTerminate(
                            speciesFromInputFile, cg, "");
            chemgraphNamesMap.put(cg, name);
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        return chemgraphNamesMap;
    }

    private static void readPrimaryThermoLibrary(BufferedReader reader) {
        ReactionModelGenerator rmg = new ReactionModelGenerator();
        String line;
        line = ChemParser.readMeaningfulLine(reader, true);
        if (line.toLowerCase().startsWith("primarythermolibrary")) {
            rmg.readAndMakePTL(reader);
        } else {
            Logger.error("ThermoDataEstimator: Could not locate the PrimaryThermoLibrary field."
                    + "Line read was: " + line);
            System.exit(0);
        }
    }

    private static void readDatabasePath(BufferedReader reader) {
        String line = ChemParser.readMeaningfulLine(reader, true);
        if (line.toLowerCase().startsWith("database")) {
            RMG.extractAndSetDatabasePath(line);
        } else {
            Logger.error("ThermoDataEstimator: Could not"
                    + " locate the Database field");
            System.exit(0);
        }
    }

    /**
     * Create a folder in the current (working) directory, deleting the existing folder and its contents if desired.
     * 
     * @param name
     *            The name of the folder to create
     * @param deleteExisting
     *            true to delete the existing folder (and its contents!), false to preserve it
     */
    public static void createFolder(String name, boolean deleteExisting) {
        File folder = new File(name);
        if (deleteExisting)
            ChemParser.deleteDir(folder);
        if (!folder.exists())
            folder.mkdir();
    }

    private static String getChemkinString(Species spe,
            ChemGraph stableChemGraph) {
        return "!" + stableChemGraph.getThermoComments() + "\n" + "!"
                + spe.getNasaThermoSource() + "\n" + "! [_ SMILES=\""
                + spe.getInChI() + "\" _]\n" + spe.getNasaThermoData();
    }
}

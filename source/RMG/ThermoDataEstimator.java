////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

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
 * allowing thermodynamic properties of chemical species to generated using 
 * quantum-mechanical methods.<BR><BR>
 * 
 * It does so by reading in a few additional lines in the input file:
 * <LI> QM active flag
 * <LI> QM method
 * <LI> QM for cyclics only flag
 * <LI> max radical number for QM integer value
 * 
 * 
 *
 */
public class ThermoDataEstimator {


	/**
	 * 
	 * @param args  filename of input file
	 */
	public static void main(String[] args) {

		RMG.globalInitializeSystemProperties();

		createFolders();

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
			ChemGraph.useQM = qmflags.qmActive.booleanValue();
			QMTP.qmprogram = qmflags.method.toLowerCase();
			ChemGraph.useQMonCyclicsOnly = qmflags.qmOnCyclicsOnly.booleanValue();

			readPrimaryThermoLibrary(reader);
			mappedChemGraphsToNames = readChemGraphsFromFile(speciesFromInputFile, reader);

		}
		catch (InvalidChemGraphException e) {
			e.printStackTrace();
		} catch (ForbiddenStructureException e) {
			e.printStackTrace();
		}
		catch (IOException e) {
			System.out.println(e.toString());
		}

		generateTDProperties(mappedChemGraphsToNames);

		Logger.info("Done!\n");

	}
	/**
	 * method that iterates over all read-in ChemGraph's, generates TD
	 * properties for all of them, and writes properties to the logger
	 * @param mappedChemGraphsToNames map with Chemgraphs and mapped names
	 */
	private static void generateTDProperties(
			Map<ChemGraph, String> mappedChemGraphsToNames) {
		for(ChemGraph chemgraph : mappedChemGraphsToNames.keySet()){

			Species spe = Species.make(mappedChemGraphsToNames.get(chemgraph),chemgraph);

			ChemGraph stableChemGraph = spe.getChemGraph();
			writeThermoDataInfo(spe, stableChemGraph);
		}
	}
	/**
	 * Create the working folders for QMTP required folders
	 */
	private static void createFolders() {
		createFolder("GATPFit", true);
		createFolder("InChI", true);
		createFolder("2Dmolfiles", true);   // Not sure if we should be deleting this
		createFolder("3Dmolfiles", true);   // Not sure if we should be deleting this
		createFolder("QMfiles", false);     // Preserving QM files between runs will speed things up considerably
	}

	/**
	 * reader method for reading all QM flags
	 * @param reader
	 * @return
	 */
	private static QMFlags readQMFlags(BufferedReader reader) {
		QMFlags qmFlags = new QMFlags();

		String line = ChemParser.readMeaningfulLine(reader, true);
		qmFlags.qmActive = Boolean.parseBoolean(line);

		line = ChemParser.readMeaningfulLine(reader, true);
		qmFlags.method = line;

		line = ChemParser.readMeaningfulLine(reader, true);
		qmFlags.qmOnCyclicsOnly = Boolean.parseBoolean(line);

		line = ChemParser.readMeaningfulLine(reader, true);
		qmFlags.maxRadNumForQM = Integer.parseInt(line);

		return qmFlags;

	}

	private static void writeThermoDataInfo(Species spe, ChemGraph stableChemGraph) {
		Logger.info(stableChemGraph.toString());

		Logger.info("The number of resonance isomers is " + 
				spe.getResonanceIsomersHashSet().size());

		Logger.info("The NASA data is \n!"+spe.getNasaThermoSource()+"\n"+
				"!" + stableChemGraph.getThermoComments() + "\n" +
				spe.getNasaThermoData());

		Logger.info("ThermoData is \n" + 
				stableChemGraph.getThermoData().toString());

		int symm = stableChemGraph.getSymmetryNumber();
		Logger.info(" symmetry number = "+symm);

		Temperature T = new Temperature(298.0,"K");

		String chemicalFormula = stableChemGraph.getChemicalFormula();

		Logger.info(chemicalFormula + "  H = " + stableChemGraph.calculateH(T));

		Logger.info("");
	}

	private static Map<ChemGraph, String> readChemGraphsFromFile(LinkedHashMap speciesFromInputFile,
			BufferedReader reader) throws IOException, ForbiddenStructureException {

		Map<ChemGraph,String> chemgraphNamesMap = new HashMap<ChemGraph,String>();

		String line = ChemParser.readMeaningfulLine(reader, true);

		while (line != null) {
			String name = line;
			Graph g = ChemParser.readChemGraph(reader);
			ChemGraph cg = ChemGraph.make(g);

			ReactionModelGenerator.addChemGraphToListIfNotPresent_ElseTerminate(speciesFromInputFile,cg,"");

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
		}
		else {
			Logger.error("ThermoDataEstimator: Could not locate the PrimaryThermoLibrary field." +
					"Line read was: " + line);
			System.exit(0);
		}
	}

	private static void readDatabasePath(BufferedReader reader) {
		String line = ChemParser.readMeaningfulLine(reader, true);
		if (line.toLowerCase().startsWith("database")) {
			RMG.extractAndSetDatabasePath(line);
		}
		else {
			Logger.error("ThermoDataEstimator: Could not"
					+ " locate the Database field");
			System.exit(0);
		}
	}

	/**
	 * Create a folder in the current (working) directory, deleting the
	 * existing folder and its contents if desired.
	 * 
	 * @param name The name of the folder to create
	 * @param deleteExisting true to delete the existing folder (and its contents!), false to preserve it
	 */
	public static void createFolder(String name, boolean deleteExisting) {
		File folder = new File(name);
		if (deleteExisting)
			ChemParser.deleteDir(folder);
		if (!folder.exists())
			folder.mkdir();
	}









}


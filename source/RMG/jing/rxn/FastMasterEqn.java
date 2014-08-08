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
package jing.rxn;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.StringReader;
import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.StringTokenizer;
import Jama.Matrix;
import jing.chem.Species;
import jing.chem.SpectroscopicData;
import jing.mathTool.UncertainDouble;
import jing.param.GasConstant;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxnSys.CoreEdgeReactionModel;
import jing.rxnSys.ReactionModelGenerator;
import jing.rxnSys.ReactionSystem;
import jing.rxnSys.Logger;

/**
 * Used to estimate pressure-dependent rate coefficients k(T, P) for a PDepNetwork object. The estimation is done by a
 * call to the Fortran module FAME. There are two methods for estimating k(T, P): the modified strong collision and the
 * reservoir state methods.
 * 
 * @author jwallen
 */
public class FastMasterEqn implements PDepKineticsEstimator {
    /**
     * The number of times the FAME module has been called for any network since the inception of this RMG execution.
     * Used to be used to number the FAME input and output files, but now the networks have individual IDs that are used
     * for this purpose.
     */
    private static int runCount = 0;

    /**
     * An enumeraction of pressure-dependent kinetics estimation methods:
     * <ul>
     * <li>NONE - The method could not be assessed.
     * <li>STRONGCOLLISION - The steady state/modified strong collision method of Chang, Bozzelli, and Dean.
     * <li>RESERVOIRSTATE - The steady state/reservoir state method of N. J. B. Green and Bhatti.
     * </ul>
     */
    public enum Mode {
        NONE, STRONGCOLLISION, RESERVOIRSTATE
    };

    /**
     * The mode to use for estimating the pressure-dependent kinetics. The mode represents the set of approximations
     * used to estimate the k(T, P) values.
     */
    private Mode mode;
    /**
     * A list of the temperatures at which the rate coefficient has been explicitly calculated..
     */
    private static Temperature[] temperatures;
    /**
     * A list of the pressures at which the rate coefficient has been explicitly calculated.
     */
    private static Pressure[] pressures;
    private static int numTBasisFuncs = -1;
    private static int numPBasisFuncs = -1;
    /**
     * The number of grains to use in a fame calculation (written in input file)
     */
    private static int numGrains = 251;
    /**
     * boolean, detailing whether the high-P-limit is greater than all of the fame-computed k(T,P)
     */
    private static boolean pdepRatesOK = true;
    /**
     * The number of atoms above which to skip pressure dependence. By default this is set to an arbitrarily large value
     * so that pressure dependence is always run.
     */
    private static int maxAtoms = 100000;
    private static int maxCarbonAtoms = 100000;

    // ==========================================================================
    //
    // Constructors
    //
    /**
     * Creates a new object with the desired mode. The mode represents the set of approximations used to estimate the
     * k(T, P) values.
     * 
     * @param m
     *            The mode to use for estimating the pressure-dependent kinetics.
     */
    public FastMasterEqn(Mode m) {
        setMode(m);
    }

    // ==========================================================================
    //
    // Accessors
    //
    /**
     * Returns the current pressure-dependent kinetics estimation mode.
     * 
     * @return The current pressure-dependent kinetics estimation mode.
     */
    public Mode getMode() {
        return mode;
    }

    /**
     * Sets the pressure-dependent kinetics estimation mode.
     * 
     * @param m
     *            The new pressure-dependent kinetics estimation mode.
     */
    public void setMode(Mode m) {
        mode = m;
    }

    // ==========================================================================
    //
    // Other methods
    //
    /**
     * Executes a pressure-dependent rate coefficient calculation.
     * 
     * @param pdn
     *            The pressure-dependent reaction network of interest
     * @param rxnSystem
     *            The reaction system of interest
     * @param cerm
     *            The current core/edge reaction model
     */
    public void runPDepCalculation(PDepNetwork pdn, ReactionSystem rxnSystem,
            CoreEdgeReactionModel cerm) {
        // If, for some unknown reason, the isomer list for the network is
        // empty, then generate it here using the path reactions
        LinkedList<PDepIsomer> isomerList0 = pdn.getIsomers();
        if (isomerList0.size() == 0) {
            for (int i = 0; i < pdn.getPathReactions().size(); i++) {
                PDepReaction rxn = pdn.getPathReactions().get(i);
                if (!isomerList0.contains(rxn.getReactant()))
                    isomerList0.add(rxn.getReactant());
                if (!isomerList0.contains(rxn.getProduct()))
                    isomerList0.add(rxn.getProduct());
            }
        }
        // No update needed if network is not altered
        if (pdn.getAltered() == false)
            return;
        // Don't do networks with isomers made up of only monatomic species
        boolean shouldContinue = true;
        for (ListIterator iter = pdn.getIsomers().listIterator(); iter
                .hasNext();) {
            PDepIsomer isomer = (PDepIsomer) iter.next();
            boolean allMonatomic = true;
            for (int i = 0; i < isomer.getNumSpecies(); i++) {
                if (!isomer.getSpecies(i).isMonatomic())
                    allMonatomic = false;
            }
            if (allMonatomic)
                shouldContinue = false;
        }
        // No update needed if network is only two wells and one reaction
        boolean noIncludedIsomers = true;
        for (ListIterator iter = pdn.getIsomers().listIterator(); iter
                .hasNext();) {
            PDepIsomer isomer = (PDepIsomer) iter.next();
            if (isomer.isUnimolecular() && isomer.getIncluded())
                noIncludedIsomers = false;
        }
        if (pdn.getIsomers().size() == 2 && pdn.getPathReactions().size() == 1
                && noIncludedIsomers)
            shouldContinue = false;
        if (!shouldContinue) {
            /*
             * LinkedList<PDepReaction> paths = pdn.getPathReactions(); LinkedList<PDepReaction> net =
             * pdn.getNetReactions(); net.clear(); for (int i = 0; i < paths.size(); i++) net.add(paths.get(i));
             */
            pdn.setAltered(false);
            return;
        }
        // Get working directory (to find FAME executable)
        String dir = System.getProperty("RMG.workingDirectory");
        // Determine wells and reactions; skip if no reactions in network
        LinkedList<Species> speciesList = pdn.getSpeciesList();
        LinkedList<PDepReaction> pathReactionList = pdn.getPathReactions();
        if (pathReactionList.size() == 0) {
            Logger.warning("Empty pressure-dependent network detected. Skipping.");
            return;
        }
        Logger.info("Solving PDepNetwork #" + Integer.toString(pdn.getID())
                + " (" + pdn.getSpeciesType() + ")");
        // Sort isomers such that the order is:
        // 1. Explored unimolecular isomers
        // 2. Bimolecular reactant/product channels
        // 3. Unexplored unimolecular isomers (treated as product channels)
        LinkedList<PDepIsomer> isomerList = new LinkedList<PDepIsomer>();
        int nIsom = 0, nReac = 0, nProd = 0;
        for (int i = 0; i < pdn.getIsomers().size(); i++) {
            PDepIsomer isom = pdn.getIsomers().get(i);
            if (isom.isUnimolecular() && isom.getIncluded())
                isomerList.add(isom);
        }
        nIsom = isomerList.size();
        if (nIsom == 0) {
            // We need at least one unimolecular isomer in order to perform a
            // P-dep calculation
            for (int i = 0; i < pdn.getIsomers().size(); i++) {
                PDepIsomer isom = pdn.getIsomers().get(i);
                if (isom.isUnimolecular() && !isom.getIncluded())
                    isomerList.add(isom);
            }
            nIsom = isomerList.size();
            for (int i = 0; i < pdn.getIsomers().size(); i++) {
                PDepIsomer isom = pdn.getIsomers().get(i);
                if (isom.isMultimolecular())
                    isomerList.add(isom);
            }
            nReac = isomerList.size() - nIsom;
            nProd = 0;
        } else {
            for (int i = 0; i < pdn.getIsomers().size(); i++) {
                PDepIsomer isom = pdn.getIsomers().get(i);
                if (isom.isMultimolecular())
                    isomerList.add(isom);
            }
            nReac = isomerList.size() - nIsom;
            for (int i = 0; i < pdn.getIsomers().size(); i++) {
                PDepIsomer isom = pdn.getIsomers().get(i);
                if (isom.isUnimolecular() && !isom.getIncluded())
                    isomerList.add(isom);
            }
            nProd = isomerList.size() - nIsom - nReac;
        }
        // Make sure all species have spectroscopic data
        for (ListIterator<PDepIsomer> iter = isomerList.listIterator(); iter
                .hasNext();) {
            PDepIsomer isomer = iter.next();
            for (int i = 0; i < isomer.getNumSpecies(); i++) {
                Species species = isomer.getSpecies(i);
                if (!species.hasSpectroscopicData())
                    species.generateSpectroscopicData();
            }
        }
        // Run garbage collection before we submit to FAME to prevent jobs from crashing with large P-Dep networks.
        // JDM January 22, 2010
        // if (runTime.freeMemory() < runTime.totalMemory()/3) {
        if (isomerList.size() >= 10) {
            Runtime runTime = Runtime.getRuntime();
            runTime.gc();
        }
        // Create FAME input files
        String input = writeInputString(pdn, rxnSystem, speciesList,
                isomerList, pathReactionList, nIsom, nReac, nProd);
        try {
            String fameOutputDir = System.getProperty("RMG.fameOutputDir");
            FileWriter fw = new FileWriter(new File(fameOutputDir + "/"
                    + Integer.toString(pdn.getID()) + "_input.txt"));
            fw.write(input);
            fw.close();
            }
            catch (IOException ex) {
            Logger.info("Unable to save FAME input file for pdep network.");
            
            }
        // Only used in the case of an error:
        StringBuilder output = new StringBuilder();
        int id = pdn.getID();
        // DEBUG only: Write input file
        if (false) { // Only do this if trying to debug
            try {
                String path = System.getProperty("RMG.fameOutputDir") + "/";
                if (id < 10)
                    path += "000";
                else if (id < 100)
                    path += "00";
                else if (id < 1000)
                    path += "0";
                path += Integer.toString(id);
                path += "_";
                path += Integer.toString(nIsom);
                path += "_input.txt";
                FileWriter fw = new FileWriter(new File(path));
                fw.write(input);
                fw.close();
            } catch (IOException ex) {
                Logger.warning("Unable to write FAME input file.");
                System.exit(0);
            }
        }
        // FAME system call
        try {
            String[] command = { dir + "/bin/fame.exe" };
            File runningDir = new File(System.getProperty("RMG.fameOutputDir"));
            Process fame = Runtime.getRuntime().exec(command, null, runningDir);
            if (fame == null)
                throw new PDepException("Couldn't start FAME process.");
            BufferedReader stdout = new BufferedReader(new InputStreamReader(
                    fame.getInputStream()));
            BufferedReader stderr = new BufferedReader(new InputStreamReader(
                    fame.getErrorStream()));
            PrintStream stdin = new PrintStream(new BufferedOutputStream(
                    fame.getOutputStream(), 1024), true);
            stdin.print(input);
            stdin.flush();
            if (stdin.checkError()) { // Flush the stream and check its error state.
                Logger.error("Error sending input to fame.exe");
            }
            stdin.close();
            String line;
            try {
                line = stdout.readLine().trim();
            } catch (NullPointerException e) {
                Logger.verbose("FAME Stderr:");
                String errline = null;
                while ((errline = stderr.readLine()) != null) {
                    Logger.verbose(errline);
                }
                throw new PDepException(
                        "FAME reported an error; FAME job was likely unsuccessful.");
            }
            /*
             * This was useful when FAME output started with a bunch of #IN and #DEBUG lines but now such information
             * appears in fame.log // advance to first line that's not for debug purposes while (
             * line.startsWith("#IN:") || line.contains("#DEBUG:") ) { //output.append(line).append("\n"); line =
             * stdout.readLine().trim(); }
             */
            if (!line.startsWith("#####")) { // Output looks like an error.
                // correct output begins with ######...
                // erroneous output does not begin with ######...
                // Print FAME stdout and error
                Logger.verbose("FAME Error:");
                Logger.verbose(line);
                output.append(line).append("\n");
                while (((line = stdout.readLine()) != null)
                        || ((line = stderr.readLine()) != null)) {
                    output.append(line).append("\n");
                    Logger.verbose(line);
                }
                // clean up i/o streams (we may be trying again with ModifiedStrongCollision and carrying on)
                stdout.close();
                stderr.close();
                throw new PDepException(
                        "Fame output looks like an error occurred.");
            }
            // Clean up i/o streams
            // This may be needed to release memory, which is especially
            // important for FAME since it can easily be called tens of
            // thousands of times in a single job
            while (stderr.ready() && (line = stderr.readLine()) != null) {
                Logger.error(line);
            }
            stderr.close();
            // Parse FAME output file and update accordingly
            if (parseOutputStream(stdout, pdn, rxnSystem, cerm, isomerList)) {
                // Reset altered flag
                pdn.setAltered(false);
                // Write finished indicator to console
                String formula = pdn.getSpeciesType();
                Logger.verbose("PDepNetwork #" + Integer.toString(id) + " ("
                        + formula + ") solved: " + pdn.getNetReactions().size()
                        + " included and "
                        + pdn.getNonincludedReactions().size()
                        + " nonincluded net reactions.");
            }
            stdout.close(); // close all output streams so the fame process terminates and we get past the following
// line:
            int exitValue = fame.waitFor();
        } catch (Exception e) {
            Logger.logStackTrace(e);
            Logger.error(e.getMessage());
            if (e.getCause() == null) {
                Logger.info("Could be because of insufficient memory and not actually a problem with FAME or the input files.");
                Logger.info("Try running fame.exe on its own.");
            }
            // Save bad input to file
            try {
                String fameOutputDir = System.getProperty("RMG.fameOutputDir");
                //FileWriter fw = new FileWriter(new File(fameOutputDir + "/"
                //        + Integer.toString(pdn.getID()) + "_input.txt"));
                //fw.write(input);
                //fw.close();
                // Don't need to duplicate the saving of the file
                Logger.info("Troublesome FAME input saved to fame/"
                        + Integer.toString(pdn.getID()) + "_input.txt");
                FileWriter fwo = new FileWriter(new File(fameOutputDir + "/"
                        + Integer.toString(pdn.getID()) + "_output.txt"));
                fwo.write(output.toString());
                fwo.close();
                Logger.info("Troublesome FAME result saved to fame/"
                        + Integer.toString(pdn.getID()) + "_output.txt");
            } catch (IOException ex) {
                Logger.info("Unable to save FAME input that caused the error.");
                System.exit(0);
            }
            // If using RS method, fall back to MSC
            if (mode == Mode.RESERVOIRSTATE) { // / mode is not defined if running modified strong collision
                Logger.info("Falling back to modified strong collision mode for this network.");
                mode = Mode.STRONGCOLLISION;
                runPDepCalculation(pdn, rxnSystem, cerm);
                mode = Mode.RESERVOIRSTATE;
                return;
            } else {
                Logger.critical("Error running FAME.");
                System.exit(0);
            }
        }
        /*
         * MRH 26Feb2010: Checking whether pdep rates exceed the high-P-limit Although fame converges, the computed
         * k(T,P) may exceed the high-P-limit, due to the number of grains being too small. If any of the pdep rates
         * exceed the high-P-limit by greater than a factor of 2, the "pdepRatesOK" boolean is set to false and fame
         * will be re-executed, using an increased number of grains After all pdep rates are below the high-P-limit (or
         * the number of grains exceeds 1000), we exit the while loop. The number of grains is then reset to 251 (which
         * was the standard before)
         */
        while (!pdepRatesOK) {
            // runPDepCalculation will increase numGrains if it is going to set pdepRatesOK=false.
            // and will set pdepRatesOK=true if numGrains > 1000.
            runPDepCalculation(pdn, rxnSystem, cerm);
        }
        numGrains = 251;
        runCount++;
    }

    /**
     * Creates the input file needed by FAME that represents a pressure- dependent reaction network.
     * 
     * @param pdn
     *            The reaction network of interest
     * @param rxnSystem
     *            The reaction system of interest
     * @param speciesList
     *            The set of species in the network
     * @param isomerList
     *            The set of isomers in the network
     * @param pathReactionList
     *            The set of path reactions in the network
     */
    public String writeInputString(PDepNetwork pdn, ReactionSystem rxnSystem,
            LinkedList<Species> speciesList, LinkedList<PDepIsomer> isomerList,
            LinkedList<PDepReaction> pathReactionList, int nIsom, int nReac,
            int nProd) {
        StringBuilder input = new StringBuilder(4096);
        Temperature stdTemp = new Temperature(298, "K");
        // Collect simulation parameters
        // MRH 28Feb: Commented temperature/pressure lines (variables are never read)
// Temperature temperature = rxnSystem.getPresentTemperature();
// Pressure pressure = rxnSystem.getPresentPressure();
        BathGas bathGas = new BathGas(rxnSystem);
        int numChebTempPolys, numChebPressPolys;
        if (numTBasisFuncs == -1) {
            numChebTempPolys = 4;
            numChebPressPolys = 4;
        } else {
            numChebTempPolys = numTBasisFuncs;
            numChebPressPolys = numPBasisFuncs;
        }
        // Determine reference energies
        double grainMinEnergy = getGrainMinEnergy(isomerList); // [=] kJ/mol
        double grainMaxEnergy = getGrainMaxEnergy(isomerList, 2100); // [=] kJ/mol
        double grainSize = getGrainSize(grainMinEnergy, grainMaxEnergy,
                pdn.getNumUniIsomers()); // [=] kJ/mol
        // Create the input string
        try {
            // Header
            input.append("################################################################################\n");
            input.append("#\n");
            input.append("#	FAME input file\n");
            input.append("#\n");
            input.append("################################################################################\n");
            input.append("\n");
            input.append("# All syntax in this file is case-insensitive\n");
            input.append("\n");
            // Method
            input.append("# The method to use to extract the phenomenological rate coefficients k(T, P)\n");
            input.append("# 	Options: ModifiedStrongCollision, ReservoirState\n");
            if (mode == Mode.STRONGCOLLISION)
                input.append("ModifiedStrongCollision\n");
            else if (mode == Mode.RESERVOIRSTATE)
                input.append("ReservoirState\n");
            else
                throw new Exception(
                        "Unable to determine method to use to estimate phenomenological rate coefficients.");
            input.append("\n");
            // Temperatures
            input.append("# The temperatures at which to estimate k(T, P)\n");
            input.append("# 	First item is the number of temperatures \n");
            input.append("#	Second item is the units; options are K, C, F, or R\n");
            input.append("#	Next comes the minimum and maximum temperatures in the specified units\n");
            input.append("#	Remaining items are the temperature values in the specified units\n");
            input.append(Integer.toString(temperatures.length) + " K ");
            input.append(Double.toString(PDepRateConstant.getTMin().getK())
                    + " " + Double.toString(PDepRateConstant.getTMax().getK())
                    + "\n");
            for (int i = 0; i < temperatures.length; i++)
                input.append(Double.toString(temperatures[i].getK()) + "\n");
            input.append("\n");
            // Pressures
            input.append("# The pressures at which to estimate k(T, P)\n");
            input.append("# 	First item is the number of pressures \n");
            input.append("#	Second item is the units ; options are bar, atm, Pa, or torr\n");
            input.append("#	Next comes the minimum and maximum pressures in the specified units\n");
            input.append("#	Remaining items are the temperature values in the specified units\n");
            input.append(Integer.toString(pressures.length) + " Pa ");
            input.append(Double.toString(PDepRateConstant.getPMin().getPa())
                    + " " + Double.toString(PDepRateConstant.getPMax().getPa())
                    + "\n");
            for (int i = 0; i < pressures.length; i++)
                input.append(Double.toString(pressures[i].getPa()) + "\n");
            input.append("\n");
            // Interpolation model
            input.append("# The interpolation model to use to fit k(T, P)\n");
            input.append("#	Option 1: No interpolation\n");
            input.append("#		Example: None\n");
            input.append("#	Option 2: Chebyshev polynomials\n");
            input.append("#		Option must be accompanied by two numbers, indicating the number of\n");
            input.append("#		terms in the Chebyshev polynomials for temperature and pressure, \n");
            input.append("#		respectively\n");
            input.append("#		Example: Chebyshev 4 4\n");
            input.append("#	Option 3: Pressure-dependent Arrhenius\n");
            input.append("#		Example: PDepArrhenius\n");
            if (PDepRateConstant.getDefaultMode() == PDepRateConstant.Mode.CHEBYSHEV)
                input.append("Chebyshev " + Integer.toString(numChebTempPolys)
                        + " " + Integer.toString(numChebPressPolys) + "\n");
            else if (PDepRateConstant.getDefaultMode() == PDepRateConstant.Mode.PDEPARRHENIUS)
                input.append("PDepArrhenius\n");
            else
                input.append("None\n");
            input.append("\n");
            // Number of energy grains to use (determines to an extent the accuracy and precision of the results)
            input.append("# A method for determining the number of energy grains to use\n");
            input.append("# Specify both the minimum number of grains and the maximum grain size\n");
            input.append("# Allowed units for the grain size are J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n");
            input.append("# Example:\n");
            input.append("#		NumGrains 251\n");
            input.append("#		GrainSize J/mol 4184\n");
            input.append("NumGrains " + numGrains);
            input.append("\n");
            input.append("GrainSize J/mol 4184");
            input.append("\n\n");
            // Collisional transfer probability model to use
            input.append("# Collisional transfer probability model\n");
            input.append("# 	Option 1: Single exponential down\n");
            input.append("#		Option must also be accompanied by three additional lines\n");
            input.append("#        First line is units and value of alpha parameter; allowed units are J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n");
            input.append("#        Second line is units and value of T0 parameter; allowed units are K\n");
            input.append("#		Third line is value of n parameter (without units)\n");
            input.append("SingleExpDown\n");
            input.append("J/mol "
                    + Double.toString(bathGas.getDeltaEdown().getAlpha() * 1000)
                    + "\n");
            input.append("K "
                    + Double.toString(bathGas.getDeltaEdown().getT0()) + "\n");
            input.append(Double.toString(bathGas.getDeltaEdown().getN()) + "\n");
            input.append("\n");
            // Other parameters for bath gas
            input.append("# Bath gas parameters\n");
            input.append("# 	Molecular weight ; allowed units are g/mol or u\n");
            input.append("# 	Lennard-Jones sigma parameter ; allowed units are m or A\n");
            input.append("# 	Lennard-Jones epsilon parameter ; allowed units are J or K\n");
            input.append("u " + Double.toString(bathGas.getMolecularWeight())
                    + "\n");
            input.append("m " + Double.toString(bathGas.getLJSigma()) + "\n");
            input.append("J " + Double.toString(bathGas.getLJEpsilon()) + "\n");
            input.append("\n");
            input.append("# The number of species in the network (minimum of 2)\n");
            input.append(Integer.toString(speciesList.size()) + "\n");
            input.append("\n");
            for (int i = 0; i < speciesList.size(); i++) {
                Species spec = speciesList.get(i);
                spec.calculateTransportParameters();
                input.append("# Species identifier (128 characters or less, no spaces)\n");
                input.append(spec.getFullName() + "\n");
                input.append("# Ground-state energy; allowed units are J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n");
                input.append("J/mol "
                        + Double.toString(spec.calculateH(stdTemp) * 4184)
                        + "\n");
                input.append("# Thermodynamics data:\n");
                input.append("# 	Standard enthalpy of formation; allowed units are J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n");
                input.append("# 	Standard entropy of formation; allowed units are permutations of energy (J, kJ, cal, or kcal) and temperature (K, C, F, or R)\n");
                input.append("# 	Heat capacity at 300, 400, 500, 600, 800, 1000, and 1500 K\n");
                input.append("J/mol "
                        + Double.toString(spec.calculateH(stdTemp) * 4184)
                        + "\n");
                input.append("J/mol*K "
                        + Double.toString(spec.calculateS(stdTemp) * 4.184)
                        + "\n");
                input.append("7 J/mol*K\n");
                input.append(Double.toString(spec.calculateCp(new Temperature(
                        300, "K")) * 4.184) + "\n");
                input.append(Double.toString(spec.calculateCp(new Temperature(
                        400, "K")) * 4.184) + "\n");
                input.append(Double.toString(spec.calculateCp(new Temperature(
                        500, "K")) * 4.184) + "\n");
                input.append(Double.toString(spec.calculateCp(new Temperature(
                        600, "K")) * 4.184) + "\n");
                input.append(Double.toString(spec.calculateCp(new Temperature(
                        800, "K")) * 4.184) + "\n");
                input.append(Double.toString(spec.calculateCp(new Temperature(
                        1000, "K")) * 4.184) + "\n");
                input.append(Double.toString(spec.calculateCp(new Temperature(
                        1500, "K")) * 4.184) + "\n");
                input.append("# Species gas parameters\n");
                input.append("# 	Molecular weight; allowed units are g/mol or u\n");
                input.append("# 	Lennard-Jones sigma parameter; allowed units are m or A\n");
                input.append("# 	Lennard-Jones epsilon parameter; allowed units are J or K\n");
                input.append("u " + Double.toString(spec.getMolecularWeight())
                        + "\n");
                input.append("m "
                        + Double.toString(spec.getChemkinTransportData()
                                .getSigma() * 1e-10) + "\n");
                input.append("J "
                        + Double.toString(spec.getChemkinTransportData()
                                .getEpsilon() * 1.380665e-23) + "\n");
                input.append("# Harmonic oscillators; allowed units are Hz and cm^-1\n");
                SpectroscopicData data = spec.getSpectroscopicData();
                input.append(Integer.toString(data.getVibrationCount())
                        + " cm^-1\n");
                for (int j = 0; j < data.getVibrationCount(); j++)
                    input.append(Double.toString(data.getVibration(j)) + "\n");
                input.append("# Rigid rotors; allowed units are Hz and cm^-1\n");
                input.append(Integer.toString(data.getRotationCount())
                        + " cm^-1\n");
                for (int j = 0; j < data.getRotationCount(); j++)
                    input.append(Double.toString(data.getRotation(j)));
                input.append("# Hindered rotor frequencies and barriers\n");
                input.append(Integer.toString(data.getHinderedCount())
                        + " cm^-1\n");
                for (int j = 0; j < data.getHinderedCount(); j++)
                    input.append(Double.toString(data.getHinderedFrequency(j))
                            + "\n");
                input.append(Integer.toString(data.getHinderedCount())
                        + " cm^-1\n");
                for (int j = 0; j < data.getHinderedCount(); j++)
                    input.append(Double.toString(data.getHinderedBarrier(j))
                            + "\n");
                input.append("# Symmetry number\n");
                input.append(Integer.toString(spec.getChemGraph()
                        .calculateSymmetryNumber()) + "\n");
                input.append("\n");
            }
            input.append("# The number of isomers in the network (minimum of 2)\n");
            input.append(Integer.toString(nIsom) + "\n");
            input.append("# The number of reactant channels in the network\n");
            input.append(Integer.toString(nReac) + "\n");
            input.append("# The number of product channels in the network\n");
            input.append(Integer.toString(nProd) + "\n");
            input.append("\n");
            for (int i = 0; i < isomerList.size(); i++) {
                PDepIsomer isomer = isomerList.get(i);
                input.append("# The number and identifiers of each species in the isomer\n");
                input.append(Integer.toString(isomer.getNumSpecies()));
                for (int j = 0; j < isomer.getNumSpecies(); j++) {
                    Species spec = isomer.getSpecies(j);
                    input.append(" " + spec.getFullName());
                }
                input.append("\n\n");
            }
            int numberOfPathReactions = 0;
            for (int i = 0; i < pathReactionList.size(); i++) {
                PDepReaction rxn = pathReactionList.get(i);
                numberOfPathReactions += rxn.getKinetics().length;
            }
            input.append("# The number of reactions in the network (minimum of 2)\n");
            input.append(Integer.toString(numberOfPathReactions) + "\n");
            input.append("\n");
            for (int i = 0; i < pathReactionList.size(); i++) {
                PDepReaction rxn = pathReactionList.get(i);
                // The reaction must be in the direction for which we are using the kinetics
                if (!rxn.isForward())
                    throw new PDepException(
                            "Encountered a path reaction that was not a forward reaction!");
                double A = 0.0, Ea = 0.0, n = 0.0;
                Kinetics[] k_array = rxn.getKinetics();
                for (int j = 0; j < k_array.length; j++) {
                    Kinetics kin = k_array[j];
                    A = kin.getAValue();
                    n = kin.getNValue();
                    Ea = kin.getEValue();// kin should be ArrheniusKinetics (rather than ArrheniusEPKinetics), so it
// should be correct to use getEValue here (similarly for other uses in this file)
                    if (A == 0 && n == 0 && Ea == 0)
                        throw new PDepException(
                                "Path reaction "
                                        + rxn.toString()
                                        + " has an A, n, and Ea of zero, which would cause FAME to crash.");
                    input.append("# The reaction equation, in the form A + B --> C + D\n");
                    input.append(rxn.toString() + "\n");
                    input.append("# Indices of the reactant and product isomers, starting with 1\n");
                    input.append(Integer.toString(isomerList.indexOf(rxn
                            .getReactant()) + 1) + " ");
                    input.append(Integer.toString(isomerList.indexOf(rxn
                            .getProduct()) + 1) + "\n");
                    input.append("# Ground-state energy; allowed units are J/mol, kJ/mol, cal/mol, kcal/mol, or cm^-1\n");
                    if (Ea < 0.0)
                        input.append("J/mol "
                                + Double.toString((rxn.getReactant()
                                        .calculateH(stdTemp)) * 4184) + "\n");
                    else
                        input.append("J/mol "
                                + Double.toString((Ea + rxn.getReactant()
                                        .calculateH(stdTemp)) * 4184) + "\n");
                    input.append("# High-pressure-limit kinetics model k(T):\n");
                    input.append("#	Option 1: Arrhenius\n");
                    input.append("# 	Arrhenius preexponential factor ; allowed units are combinations of volume {m^3, L, or cm^3} and time {s^-1}\n");
                    input.append("# 	Arrhenius activation energy ; allowed units are J/mol, kJ/mol, cal/mol, or kcal/mol\n");
                    input.append("# 	Arrhenius temperature exponent\n");
                    input.append("Arrhenius\n");
                    if (rxn.getReactant().isUnimolecular())
                        input.append("s^-1 " + Double.toString(A) + "\n");
                    else if (rxn.getReactant().isMultimolecular())
                        input.append("m^3/mol*s " + Double.toString(A * 1.0e-6)
                                + "\n");
                    input.append("J/mol " + Double.toString(Ea * 4184) + "\n");
                    input.append(Double.toString(n) + "\n");
                    input.append("\n");
                }
            }
        } catch (IndexOutOfBoundsException e) {
            Logger.error("Error: IndexOutOfBoundsException thrown.");
            Logger.logStackTrace(e);
        } catch (Exception e) {
            Logger.error(e.getMessage());
            Logger.logStackTrace(e);
        }
        return input.toString();
    }

    /**
     * Parses one meaningful line from a FAME output buffer. Returns null at the end of the file (as would br.readLine()
     * )
     */
    public String readMeaningfulLine(BufferedReader br) throws IOException {
        String str = "";
        boolean found = false;
        while (!found) {
            str = br.readLine();
            if (str == null)
                return null;
            str = str.trim();
            found = !(str.length() == 0 || str.startsWith("#"));
        }
        return str;
    }

    /**
     * Parses a FAME output file and updates the reaction network and system accordingly.
     * 
     * @param br
     *            The BufferedReader containing the stream that is to be parsed. THIS IS LEFT OPEN.
     * @param pdn
     *            The pressure-dependent reaction network of interest
     * @param rxnSystem
     *            The reaction system of interest
     * @param cerm
     *            The current core/edge reaction model
     * @param isomerList
     *            The set of isomers in the network
     */
    public boolean parseOutputStream(BufferedReader br, PDepNetwork pdn,
            ReactionSystem rxnSystem, CoreEdgeReactionModel cerm,
            LinkedList<PDepIsomer> isomerList) throws PDepException {
        String dir = System.getProperty("RMG.workingDirectory");
        double Tmin = 0, Tmax = 0, Pmin = 0, Pmax = 0;
        // JWA January 28 2010
        // Clear the included reactions from the last successful fame execution
        LinkedList<PDepReaction> netReactionList = pdn.getNetReactions();
        netReactionList.clear();
        // Also need to clear the nonincluded reactions - we will regenerate them
        // based on the new output from fame
        pdn.getNonincludedReactions().clear();
        // Below we are going to use netReactionList to temporarily store all
        // reactions, both included and nonincluded, and will then sort them
        // into included and nonincluded
        boolean ignoredARate = false;
        /*
         * MRH 26Feb2010: Checking whether pdep rates exceed the high-P-limit This boolean reflects whether all pdep
         * rates are below the high-P-limit. The variable is initialized to true for each call (in hope of avoiding any
         * infinite loop scenarios). If any of the pdep rates exceed the high-P-limit, the boolean is changed to false
         * (see code below)
         */
        pdepRatesOK = true;
        try {
            String str = "";
            StringTokenizer tkn;
            String method = readMeaningfulLine(br);
            str = readMeaningfulLine(br);
            tkn = new StringTokenizer(str);
            int numTemp = Integer.parseInt(tkn.nextToken());
            tkn.nextToken();
            Tmin = Double.parseDouble(tkn.nextToken());
            for (int i = 0; i < numTemp - 2; i++)
                tkn.nextToken();
            Tmax = Double.parseDouble(tkn.nextToken());
            str = readMeaningfulLine(br);
            tkn = new StringTokenizer(str);
            int numPress = Integer.parseInt(tkn.nextToken());
            tkn.nextToken();
            Pmin = Double.parseDouble(tkn.nextToken());
            for (int i = 0; i < numPress - 2; i++)
                tkn.nextToken();
            // Added by MRH on 10Feb2010 to handle the case where user only specifies one pressure
            if (numPress == 1)
                Pmax = Pmin;
            else
                Pmax = Double.parseDouble(tkn.nextToken());
            str = readMeaningfulLine(br);
            tkn = new StringTokenizer(str);
            String model = tkn.nextToken();
            int numChebT = 0, numChebP = 0;
            if (model.toLowerCase().contains("chebyshev")) {
                numChebT = Integer.parseInt(tkn.nextToken());
                numChebP = Integer.parseInt(tkn.nextToken());
            }
            str = readMeaningfulLine(br);
            int numSpecies = Integer.parseInt(str);
            str = readMeaningfulLine(br);
            int numIsomers = Integer.parseInt(str);
            str = readMeaningfulLine(br);
            int numReactants = Integer.parseInt(str);
            str = readMeaningfulLine(br);
            int numProducts = Integer.parseInt(str);
            str = readMeaningfulLine(br);
            int numReactions = Integer.parseInt(str);
            str = readMeaningfulLine(br);
            int numKinetics = Integer.parseInt(str);
            for (int i = 0; i < numKinetics; i++) {
                double[][] rates = new double[numTemp][numPress];
                double[][] chebyshev = null;
                PDepArrheniusKinetics pDepArrhenius = null;
                // Reactant and product isomers
                str = readMeaningfulLine(br);
                if (str == null) {
                    throw new PDepException(
                            String.format(
                                    "Was expecting %d sets of phenomenological rate coefficients in fame output, but output ends after %s",
                                    numKinetics, i));
                }
                tkn = new StringTokenizer(str);
                int reac = Integer.parseInt(tkn.nextToken()) - 1;
                int prod = Integer.parseInt(tkn.nextToken()) - 1;
                // Table of phenomenological rate coefficients
                boolean valid = true;
                str = readMeaningfulLine(br); // First row is list of pressures
                for (int t = 0; t < numTemp; t++) {
                    str = readMeaningfulLine(br);
                    tkn = new StringTokenizer(str);
                    tkn.nextToken();
                    for (int p = 0; p < numPress; p++) {
                        rates[t][p] = Double.parseDouble(tkn.nextToken());
                        if (rates[t][p] < 0 || Double.isNaN(rates[t][p])
                                || Double.isInfinite(rates[t][p]))
                            valid = false;
                    }
                }
                PDepRateConstant pDepRate = new PDepRateConstant(rates);
                // Chebyshev interpolation model
                if (model.toLowerCase().contains("chebyshev")) {
                    chebyshev = new double[numChebT][numChebP];
                    for (int t = 0; t < numChebT; t++) {
                        str = readMeaningfulLine(br);
                        tkn = new StringTokenizer(str);
                        for (int p = 0; p < numChebP; p++) {
                            chebyshev[t][p] = Double.parseDouble(tkn
                                    .nextToken());
                        }
                    }
                    pDepRate.setChebyshev(chebyshev);
                }
                // PDepArrhenius interpolation model
                else if (model.toLowerCase().contains("pdeparrhenius")) {
                    pDepArrhenius = new PDepArrheniusKinetics(numPress);
                    for (int p = 0; p < numPress; p++) {
                        str = readMeaningfulLine(br);
                        tkn = new StringTokenizer(str);
                        tkn.nextToken();
                        UncertainDouble A = new UncertainDouble(
                                Double.parseDouble(tkn.nextToken()), 0.0, "A");
                        UncertainDouble Ea = new UncertainDouble(
                                Double.parseDouble(tkn.nextToken()) / 4184.0,
                                0.0, "A");
                        UncertainDouble n = new UncertainDouble(
                                Double.parseDouble(tkn.nextToken()), 0.0, "A");
                        String Trange = Double.toString(Tmin) + "-"
                                + Double.toString(Tmax) + " K";
                        String Prange = Double.toString(Pmin / 1.0e5) + "-"
                                + Double.toString(Pmax / 1.0e5) + " bar";
                        ArrheniusKinetics kinetics = new ArrheniusKinetics(A,
                                n, Ea, Trange, 0, "Result of FAME calculation",
                                "Prange = " + Prange);
                        pDepArrhenius.setKinetics(p, pressures[p], kinetics);
                    }
                    pDepRate.setPDepArrheniusKinetics(pDepArrhenius);
                }
                // If the fitted rate coefficient is not valid, then don't add the net reaction
                if (!valid) {
                    ignoredARate = true;
                    continue;
                }
                // Initialize net reaction
                PDepIsomer reactant = isomerList.get(reac);
                PDepIsomer product = isomerList.get(prod);
                PDepReaction rxn = new PDepReaction(reactant, product, pDepRate);
                // Create the reverse net reaction if applicable
                if (product.getIncluded())
                    rxn.generateReverseReaction();
                else
                    rxn.setReverseReaction(null);
                /*
                 * MRH 26Feb2010: Checking whether pdep rates exceed the high-P-limit We grab all of the reactions in
                 * the pathReactionList. These are the RMG-generated reactions (not chemically-activated reactions) and
                 * thus have "natural" high-P-limit kinetics. If the current PDepReaction "rxn"'s structure matches one
                 * of the structures of the PDepReactions located in the pathReactionList, we compare the high-P-limit
                 * kinetics of the pathReactionList (these values either come from the RMG database or are "fitted"
                 * parameters, based on the reverse kinetics + equilibrium constant) with the k(T,P_max) for each T in
                 * the "temperatures" array. NOTE: Not every "rxn" will have a match in the pathReactionList. If the
                 * pdep rate is greater than 2x the high-P-limit, we consider this to be different. If the number of
                 * grains is less than 1000, we set the pdepRatesOK boolean to false, so that another fame calculation
                 * will ensue If the number of grains exceeds 1000, we continue on with the simulation, but alert the
                 * user of the discrepancy. The value of 2 was somewhat randomly chosen by MRH. For a toy case of tBuOH
                 * pyrolysis (with 1e-6 reaction time and 0.9 error tolerance): Before code addition: iC4H8/H2O final
                 * concentration = 7.830282E-10, run time < 1 min 2x : iC4H8/H2O final concentration = 6.838976E-10, run
                 * time ~ 2 min 1.5x : iC4H8/H2O final concentration = 6.555548E-10, run time ~ 4 min 1.1x : iC4H8/H2O
                 * final concentration = 6.555548E-10, run time ~ 10 min The value of 2 (for the toy model) seems to be
                 * a good balance between speed and accuracy P.S. Want to keep the name fame (fast approximate me)
                 * instead of having to change to smame (slow more accurate me). ;) JWA 01Nov2011: Note that the k(T,P)
                 * values always combine both direct and well-skipping effects. (In this sense they are not true rate
                 * coefficients, but are instead "flux" coefficients.) At low T and high P, the well-skipping effect is
                 * usually very small. However, there are many examples of isomerization reactions for which the
                 * well-skipping rate is much larger than the direct rate (e.g. due to a very high barrier for the
                 * direct reaction). For this reason, we do not apply the check to isomerization reactions, since they
                 * are not necessarily wrong if the check fails.
                 */
                LinkedList pathReactionList = pdn.getPathReactions();
                boolean foundHighPLimitRxn = false;
                Temperature stdtemp = new Temperature(298, "K");
                double Hrxn;
                for (int HighPRxNum = 0; HighPRxNum < pathReactionList.size(); HighPRxNum++) {
                    PDepReaction rxnWHighPLimit = (PDepReaction) pathReactionList
                            .get(HighPRxNum);
                    if (rxn.getStructure()
                            .equals(rxnWHighPLimit.getStructure())) {
                        if (rxn.getReactant().isUnimolecular()
                                && rxn.getProduct().isUnimolecular())
                            // Don't apply the check to isomerization reactions; see above comment
                            continue;
                        foundHighPLimitRxn = true;
                        Hrxn = rxnWHighPLimit.calculateHrxn(stdtemp);
                        double A = 0.0, Ea = 0.0, n = 0.0;
                        if (rxnWHighPLimit.isForward()) {
                            Kinetics[] k_array = rxnWHighPLimit.getKinetics();
                            Kinetics kin = computeKUsingLeastSquares(k_array,
                                    Hrxn);
                            A = kin.getAValue();
                            Ea = kin.getEValue();
                            n = kin.getNValue();
                            // While I'm here, and know which reaction was the High-P limit, set the comment in the
// P-dep reaction
                            rxn.setComments("NetReaction from PDepNetwork #"
                                    + Integer.toString(pdn.getID()) + " ("
                                    + pdn.getSpeciesType() + ")"
                                    + " High-P Limit: "
                                    + kin.getSource().toString() + " "
                                    + kin.getComment().toString());
                        } else {
                            Kinetics[] k_array = rxnWHighPLimit
                                    .getFittedReverseKinetics();
                            Kinetics kin = computeKUsingLeastSquares(k_array,
                                    -Hrxn);// gmagoon: I'm not sure, with forward/reverse reactions here whether it is
// correct to use Hrxn or -Hrxn, but in any case, getFittedReverseKinetics should return an ArrheniusKinetics (not
// ArrheniusEPKinetics) object, so it will not be used in computeKUsingLeastSquares anyway
                            A = kin.getAValue();
                            Ea = kin.getEValue();
                            n = kin.getNValue();
                            // While I'm here, and know which reaction was the High-P limit, set the comment in the
// P-dep reaction
                            Kinetics[] fwd_kin = rxnWHighPLimit.getKinetics();
                            String commentsForForwardKinetics = "";
                            if (fwd_kin.length > 1)
                                commentsForForwardKinetics += "Summation of kinetics:\n!";
                            for (int numKs = 0; numKs < fwd_kin.length; ++numKs) {
                                commentsForForwardKinetics += "High-P Limit Reverse: "
                                        + fwd_kin[numKs].getSource().toString()
                                        + fwd_kin[numKs].getComment()
                                                .toString();
                                if (numKs != fwd_kin.length - 1)
                                    commentsForForwardKinetics += "\n!";
                            }
                            rxn.setComments("NetReaction from PDepNetwork #"
                                    + Integer.toString(pdn.getID()) + " ("
                                    + pdn.getSpeciesType() + ")" + " "
                                    + commentsForForwardKinetics);
                        }
                        if (ReactionModelGenerator
                                .rerunFameWithAdditionalGrains()) {
                            double[][] all_ks = rxn.getPDepRate()
                                    .getRateConstants();
                            double T = temperatures[0].getK(); // lowest temperature will have the highest
// over_high_P_factor.
                            double k_highPlimit = A
                                    * Math.pow(T, n)
                                    * Math.exp(-Ea / GasConstant.getKcalMolK()
                                            / T);
                            double over_high_P_factor = all_ks[0][pressures.length - 1]
                                    / k_highPlimit;
                            if (over_high_P_factor > 2) {
                                Logger.info("For reaction " + rxn.toString());
                                Logger.info(String
                                        .format("Pressure-dependent rate coefficient at %.0fK %.1fBar "
                                                + "exceeds high-P-limit rate  by factor of %.1f .",
                                                T, Pmax * 1e-5,
                                                over_high_P_factor));
                                pdepRatesOK = false;
                            }
                        }
                    }
                }
                // If not found, we have a "nonIncluded" (pressure-dependent) reaction
                if (!foundHighPLimitRxn) {
                    rxn.setComments("NetReaction from PDepNetwork #"
                            + Integer.toString(pdn.getID()) + " ("
                            + pdn.getSpeciesType() + ")");
                }
                // Add net reaction to list
                netReactionList.add(rxn);
            }
            // If we ignored a rate coefficient from the FAME output for any
            // reason, then fail the FAME job
            if (ignoredARate) {
                throw new PDepException(
                        "One or more rate coefficients from the FAME output was ignored, possibly due to a NaN or Inf rate.");
            }
            // Update reaction lists (sort into included and nonincluded)
            pdn.updateReactionLists(cerm);
        } catch (PDepException e) {
            Logger.verbose(pdn.toString());
            Logger.error(e.getMessage());
            Logger.logStackTrace(e);
            throw new PDepException("Unable to parse FAME output file. "
                    + e.getMessage());
        } catch (IOException e) {
            Logger.error("Unable to read from file \"fame_output.txt\".");
            System.exit(0);
        }
        if (!pdepRatesOK) {
            if (numGrains > 1000) {
                Logger.info("Number of grains already exceeds 1000. "
                        + "Continuing with results from current fame run.");
                pdepRatesOK = true; // this function is called inside a while(!pdepRatesOK) loop.
            } else {
                numGrains = numGrains + 250;
                Logger.info(String.format("Re-running fame with %d grains.",
                        numGrains));
                return false;
            }
        }
        return true;
    }

    /**
     * Determines the maximum energy grain to use in the calculation. The maximum energy grain is chosen to be 25 * R *
     * T above the highest potential energy on the potential energy surface, which hopefully will capture the majority
     * of the equilibrium distributions even at the high temperature end of the calculation.
     * 
     * @param uniIsomers
     *            The set of unimolecular isomers in the network
     * @param multiIsomers
     *            The set of multimolecular isomers in the network
     * @param T
     *            The temperature of the calculation in K
     * @return The maximum grain energy in kJ/mol
     */
    private double getGrainMaxEnergy(LinkedList<PDepIsomer> isomerList, double T) {
        double Emax = 0.0;
        Temperature stdTemp = new Temperature(298, "K");
        for (int i = 0; i < isomerList.size(); i++) {
            PDepIsomer isomer = isomerList.get(i);
            double E = isomer.calculateH(stdTemp);
            if (E > Emax)
                Emax = E;
        }
        Emax *= 4.184;
        // Emax += 100 * 0.008314 * T;
        Emax += 50 * 0.008314 * T;
        // Round up to nearest ten kJ/mol
        return Math.ceil(Emax / 10) * 10.0; // [=] kJ/mol
    }

    /**
     * Determines the minimum energy grain to use in the calculation. The maximum energy grain is chosen to be at or
     * just below the energy of the lowest energy isomer in the system.
     * 
     * @param uniIsomers
     *            The set of unimolecular isomers in the network
     * @param multiIsomers
     *            The set of multimolecular isomers in the network
     * @return The minimum grain energy in kJ/mol
     */
    private double getGrainMinEnergy(LinkedList<PDepIsomer> isomerList) {
        double Emin = 1000000.0;
        Temperature stdTemp = new Temperature(298, "K");
        for (int i = 0; i < isomerList.size(); i++) {
            PDepIsomer isomer = isomerList.get(i);
            double E = isomer.calculateH(stdTemp);
            if (E < Emin)
                Emin = E;
        }
        Emin *= 4.184;
        // Round down to nearest ten kJ/mol
        return Math.floor(Emin / 10) * 10.0; // [=] kJ/mol
    }

    /**
     * Determines a suitable grain size for the calculation. Too small a grain size will cause the simulation to take a
     * very long time to conduct; too large a grain size will cause the simulation results to not be very precise. The
     * choice is hopefully made such that the grains get larger as the system does, so as to not get completely bogged
     * down in the larger networks.
     * 
     * @param grainMinEnergy
     *            The minimum energy grain in kJ/mol
     * @param grainMaxEnergy
     *            The maximum energy grain in kJ/mol
     * @param numUniWells
     *            The number of unimolecular isomers in the network.
     * @return
     */
    private double getGrainSize(double grainMinEnergy, double grainMaxEnergy,
            int numUniWells) {
        if (numUniWells < 5)
            return (grainMaxEnergy - grainMinEnergy) / 200;
        else if (numUniWells < 10)
            return (grainMaxEnergy - grainMinEnergy) / 100;
        else if (numUniWells < 20)
            return (grainMaxEnergy - grainMinEnergy) / 50;
        else
            return (grainMaxEnergy - grainMinEnergy) / 20;
    }

    /**
     * Creates an empty file on the hard disk for the FORTRAN executable to use to write its output. If this file is not
     * created, the FORTRAN executable will write to the file 'fort.2', which is checked by the readOutputFile()
     * function but still discouraged.
     */
    public void touchOutputFile() {
        try {
            File output = new File("fame/fame_output.txt");
            if (output.exists())
                output.delete();
            output.createNewFile();
        } catch (IOException e) {
            Logger.error("Unable to touch file \"fame_output.txt\".");
            Logger.error(e.getMessage());
            System.exit(0);
        } catch (Exception e) {
            Logger.error(e.getMessage());
        }
    }

    /**
     * Returns the mode as a string suitable for writing the FAME input file.
     * 
     * @return The mode as a string suitable for writing the FAME input file.
     */
    public String getModeString() {
        if (mode == Mode.STRONGCOLLISION)
            return "ModifiedStrongCollision";
        else if (mode == Mode.RESERVOIRSTATE)
            return "ReservoirState";
        else
            return "";
    }

    public static Temperature[] getTemperatures() {
        return temperatures;
    }

    public static void setTemperatures(Temperature[] temps) {
        temperatures = temps;
    }

    public static Pressure[] getPressures() {
        return pressures;
    }

    public static void setPressures(Pressure[] press) {
        pressures = press;
    }

    public static void setNumTBasisFuncs(int n) {
        numTBasisFuncs = n;
    }

    public static void setNumPBasisFuncs(int m) {
        numPBasisFuncs = m;
    }

    /**
     * Return the number of atoms above which pressure dependence is assumed to be negligible.
     */
    public static int getMaxAtoms() {
        return maxAtoms;
    }

    /**
     * Return the number of carbon atoms above which pressure dependence is assumed to be negligible.
     */
    public static int getMaxCarbonAtoms() {
        return maxCarbonAtoms;
    }
    
    /**
     * Set the number of atoms above which pressure dependence is assumed to be negligible. To always run pressure
     * dependence, set this to a very large number.
     */
    public static void setMaxAtoms(int atoms) {
        maxAtoms = atoms;
    }

    /**
     * Set the number of atoms above which pressure dependence is assumed to be negligible. To always run pressure
     * dependence, set this to a very large number.
     */
    public static void setMaxCarbonAtoms(int atoms) {
        maxCarbonAtoms = atoms;
    }
    
    /**
     * Return true if the given reaction should be considered as a pressure dependent reaction or false otherwise. A
     * pressure dependent reaction has the form A -> B, A -> B + C, or A + B -> C with a total number of atoms (or number of Carbon atoms) below a
     * certain threshold.
     * 
     * @param reaction
     *            The reaction to assess
     * @return true if the reaction should be considered as pressure dependent, false otherwise
     */
    public static boolean isReactionPressureDependent(Reaction reaction) {
        if (reaction.getAtomNumber() > maxAtoms || reaction.getCarbonAtomNumber() > maxCarbonAtoms)
            return false;
        else
            return (reaction.getReactantNumber() == 1 || reaction
                    .getProductNumber() == 1);
    }

    public static Kinetics computeKUsingLeastSquares(Kinetics[] k_array,
            double Hrxn) {
        /*
         * MRH 24MAR2010: If the reaction has more than one set of Arrhenius kinetics, sum all kinetics and re-fit for a
         * single set of Arrhenius kinetics
         */
        Kinetics k = null;
        if (k_array.length > 1) {
            double[] T = new double[5];
            T[0] = 300;
            T[1] = 600;
            T[2] = 900;
            T[3] = 1200;
            T[4] = 1500;
            double[] k_total = new double[5];
            for (int numKinetics = 0; numKinetics < k_array.length; ++numKinetics) {
                for (int numTemperatures = 0; numTemperatures < T.length; ++numTemperatures) {
                    double Ea = 0.0;
                    Ea = ((ArrheniusKinetics) k_array[numKinetics]).getEValue();
                    k_total[numTemperatures] += k_array[numKinetics]
                            .getAValue()
                            * Math.pow(T[numTemperatures],
                                    k_array[numKinetics].getNValue())
                            * Math.exp(-Ea / GasConstant.getKcalMolK()
                                    / T[numTemperatures]);
                }
            }
            // Construct matrix X and vector y
            double[][] y = new double[k_total.length][1];
            double[][] X = new double[k_total.length][3];
            for (int i = 0; i < 5; i++) {
                y[i][0] = Math.log(k_total[i]);
                X[i][0] = 1;
                X[i][1] = Math.log(T[i]);
                X[i][2] = 1 / T[i];
            }
            // Solve least-squares problem using inv(XT*X)*(XT*y)
            Matrix X_matrix = new Matrix(X);
            Matrix y_matrix = new Matrix(y);
            Matrix b_matrix = X_matrix.solve(y_matrix);
            UncertainDouble uA = new UncertainDouble(Math.exp(b_matrix
                    .get(0, 0)), 0.0, "Adding");
            UncertainDouble un = new UncertainDouble(b_matrix.get(1, 0), 0.0,
                    "Adding");
            UncertainDouble uE = new UncertainDouble(-GasConstant.getKcalMolK()
                    * b_matrix.get(2, 0), 0.0, "Adding");
            String commentsForFittedKinetics = "";
            for (int numKs = 0; numKs < k_array.length; ++numKs) {
                commentsForFittedKinetics += "!"
                        + k_array[numKs].getSource().toString();
                if (k_array[numKs].getComment() != null)
                    commentsForFittedKinetics += k_array[numKs].getComment()
                            .toString();
                if (numKs != k_array.length - 1)
                    commentsForFittedKinetics += "\n";
            }
            k = new ArrheniusKinetics(uA, un, uE, "300-1500K", 5,
                    "Summation of kinetics:", commentsForFittedKinetics);
        } else {
            k = k_array[0];
        }
        return k;
    }

    public static int getNumTBasisFuncs() {
        return numTBasisFuncs;
    }

    public static int getNumPBasisFuncs() {
        return numPBasisFuncs;
    }
}

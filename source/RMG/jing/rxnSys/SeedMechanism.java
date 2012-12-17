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
package jing.rxnSys;

import java.io.*;
import jing.mathTool.UncertainDouble;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxn.*;
import jing.chem.*;
import java.util.*;
import jing.chemUtil.*;
import jing.chemParser.*;

/**
 * This is a new class called SeedMechanism. SeedMechanism is the same class as the old PrimaryKineticLibrary, just with
 * a new (and more appropriate) name. RMG will automatically include every species and reaction contained in a Seed
 * Mechanism. Furthermore, the user has the option to pass multiple Seed Mechanisms to RMG. In the event of a duplicate
 * species/reaction, RMG will use the first instance it finds (i.e. the order of the Seed Mechanisms listed in the
 * condition.txt file is important). MRH 9-Jun-2009
 */
/*
 * Comments from old PrimaryKineticLibrary: This is the primary reaction set that any reaction system has to include
 * into its model. For example, in combustion system, we build a primary small molecule reaction set, and every
 * combustion/oxidation system should include such a primary reaction library. The reaction / rates are basically from
 * Leeds methane oxidation mechanism.
 */
public class SeedMechanism {
    protected String name;
    protected LinkedHashSet reactionSet = new LinkedHashSet();
    protected LinkedHashMap speciesSet = new LinkedHashMap();
    private boolean generateReactions = false;
    public LinkedList allPdepNetworks = new LinkedList();

    // Constructors
    public SeedMechanism(String p_mechName, String p_directoryPath,
            boolean p_generateReactions, boolean p_fromRestart)
            throws IOException {
        name = p_mechName;
        generateReactions = p_generateReactions;
        if (p_directoryPath == null)
            throw new NullPointerException(
                    "RMG does not recognize Seed Mechanism directory path: Value is null");
        try {
            read(p_directoryPath, p_fromRestart, p_mechName);
        } catch (IOException e) {
            throw new IOException("Error in reading Seed Mechanism: "
                    + p_mechName + '\n' + e.getMessage());
        }
    }

    public SeedMechanism() {
    }

    public void appendSeedMechanism(String new_mechName,
            String new_directoryPath, boolean p_generateReactions,
            boolean p_fromRestart) throws IOException {
        if (p_generateReactions)
            setGenerateReactions(p_generateReactions);
        setName(name + "/" + new_mechName);
        try {
            read(new_directoryPath, p_fromRestart, new_mechName);
        } catch (IOException e) {
            throw new IOException("Error in reading Seed Mechanism: "
                    + new_mechName + '\n' + e.getMessage());
        }
    }

    public LinkedHashSet getSpeciesSet() {
        return new LinkedHashSet(speciesSet.values());
    }

    public void read(String p_directoryName, boolean p_fromRestart,
            String seedMechName) throws IOException {
        Logger.info("Reading seed mechanism from directory " + p_directoryName);
        LinkedHashMap localSpecies = null;
        LinkedHashSet localReactions = null;
        try {
            if (!p_directoryName.endsWith("/"))
                p_directoryName = p_directoryName + "/";
            if (!p_fromRestart) {
                String speciesFile = p_directoryName + "species.txt";
                String reactionFile = p_directoryName + "reactions.txt";
                String pdepreactionFile = p_directoryName + "pdepreactions.txt";
                speciesSet.putAll(readSpecies(speciesFile, seedMechName,
                        "Seed Mechanism: "));
                reactionSet.addAll(readReactions(reactionFile, seedMechName,
                        speciesSet, "Seed Mechanism: "));
                reactionSet.addAll(readPdepReactions(pdepreactionFile,
                        seedMechName, speciesSet, "Seed Mechanism: "));
            } else {
                String speciesFile = p_directoryName + "coreSpecies.txt";
                String pdepreactionFile = p_directoryName + "pdepreactions.txt";
                speciesSet.putAll(readSpecies(speciesFile, seedMechName,
                        "Seed Mechanism: "));
                reactionSet.addAll(readPdepReactions(pdepreactionFile,
                        seedMechName, speciesSet, "Seed Mechanism: "));
            }
            return;
        } catch (Exception e) {
            Logger.logStackTrace(e);
            throw new IOException("RMG cannot read entire Seed Mechanism: "
                    + p_directoryName + "\n" + e.getMessage());
        }
    }

    public LinkedHashSet readReactions(String p_reactionFileName,
            String p_name, LinkedHashMap allSpecies, String source)
            throws IOException {
        LinkedHashSet localReactions = new LinkedHashSet();
        try {
            FileReader in = new FileReader(p_reactionFileName);
            BufferedReader data = new BufferedReader(in);
            double[] multipliers = parseReactionRateUnits(data);
            double A_multiplier = multipliers[0];
            double E_multiplier = multipliers[1];
            String line = ChemParser.readMeaningfulLine(data, true);
            read: while (line != null) {
                Reaction r;
                try {
                    r = ChemParser.parseArrheniusReaction(allSpecies, line,
                            A_multiplier, E_multiplier);
                    r.setKineticsSource(source + p_name, 0);
                    r.setKineticsComments(" ", 0);
                    r.setIsFromPrimaryKineticLibrary(true);
                    (r.getKinetics())[0].setFromPrimaryKineticLibrary(true);
                } catch (InvalidReactionFormatException e) {
                    throw new InvalidReactionFormatException(line + ": "
                            + e.getMessage());
                }
                if (r == null)
                    throw new InvalidReactionFormatException(line);
                line = ChemParser.readMeaningfulLine(data, true);
                localReactions = updateReactionList(r, localReactions, true,
                        line);
                if (line != null && line.toLowerCase().startsWith("dup"))
                    line = ChemParser.readMeaningfulLine(data, true);
            }
            in.close();
            // Check that all reactions labeled with "DUP" tag actually have a duplicate
            for (Iterator iter = localReactions.iterator(); iter.hasNext();) {
                Reaction r = (Reaction) iter.next();
                if (r.getExpectDuplicate()) {
                    if (r.getKinetics().length < 2) {
                        Logger.critical("The following reaction was labeled with a duplicate tag,"
                                + " but did not contain a duplicate:\n\t"
                                + r.toString()
                                + "\nPlease correct the kinetic_libraries directory reactions.txt file");
                        System.exit(0);
                    }
                }
            }
            return localReactions;
        } catch (Exception e) {
            Logger.logStackTrace(e);
            Logger.error("RMG did not read the following " + source + p_name
                    + " file: " + p_reactionFileName + " because "
                    + e.getMessage());
            return null;
        }
    }

    public LinkedHashMap readSpecies(String p_speciesFileName, String p_name,
            String source) throws IOException {
        LinkedHashMap localSpecies = new LinkedHashMap();
        try {
            FileReader in = new FileReader(p_speciesFileName);
            BufferedReader data = new BufferedReader(in);
            // step 1: read in structure
            String line = ChemParser.readMeaningfulLine(data, true);
            read: while (line != null) {
                // GJB allow unreactive species
                StringTokenizer st = new StringTokenizer(line);
                String name = st.nextToken().trim();
                boolean IsReactive = true;
                if (st.hasMoreTokens()) {
                    String reactive = st.nextToken().trim();
                    if (reactive.equalsIgnoreCase("unreactive"))
                        IsReactive = false;
                }
                Graph graph;
                try {
                    graph = ChemParser.readChemGraph(data);
                    if (graph == null)
                        throw new IOException("Graph was null");
                } catch (IOException e) {
                    throw new InvalidChemGraphException("Cannot read species '"
                            + name + "': " + e.getMessage());
                }
                ChemGraph cg = ChemGraph.make(graph);
                Species spe = Species.make(name, cg);
                // Check if species (i.e. chemgraph) already exists in localSpecies
                if (localSpecies.containsValue(spe)) {
                    for (Iterator iter2 = localSpecies.values().iterator(); iter2
                            .hasNext();) {
                        Species spcInList = (Species) iter2.next();
                        if (spcInList.equals(spe)) {
                            Logger.critical("Species '"
                                    + name
                                    + "' chemgraph already exists in library"
                                    + " with name "
                                    + spcInList.getName()
                                    + "\n\tRemove one of the instances before proceeding.");
                            System.exit(0);
                        }
                    }
                }
                // GJB: Turn off reactivity if necessary, but don't let code turn it on
                // again if was already set as unreactive from input file
                if (IsReactive == false)
                    spe.setReactivity(IsReactive);
                localSpecies.put(name, spe);
                line = ChemParser.readMeaningfulLine(data, true);
            }
            in.close();
            return localSpecies;
        } catch (Exception e) {
            Logger.logStackTrace(e);
            String message = "RMG cannot read the \"species.txt\" file in the "
                    + source + p_name + "\n" + e.getMessage();
            if (e instanceof jing.chem.ForbiddenStructureException)
                message += "\n Try adding it to a Thermo Library to make it allowed.";
            throw new IOException(message);
        }
    }

    public LinkedHashSet readPdepReactions(String pdepFileName, String p_name,
            LinkedHashMap allSpecies, String source) throws IOException {
        LinkedHashSet localReactions = new LinkedHashSet();
        LinkedList pdepNetworks = getPDepNetworks();
        try {
            FileReader in = new FileReader(pdepFileName);
            BufferedReader data = new BufferedReader(in);
            double[] multipliers = parseReactionRateUnits(data);
            double A_multiplier = multipliers[0];
            double E_multiplier = multipliers[1];
            String nextLine = ChemParser.readMeaningfulLine(data, true);
            read: while (nextLine != null) {
                Reaction r;
                try {
                    r = ChemParser.parseArrheniusReaction(allSpecies, nextLine,
                            A_multiplier, E_multiplier);
                } catch (InvalidReactionFormatException e) {
                    throw new InvalidReactionFormatException(nextLine + ": "
                            + e.getMessage());
                }
                if (r == null)
                    throw new InvalidReactionFormatException(nextLine);
                r.setIsFromPrimaryKineticLibrary(true);
                /*
                 * Read the next line and determine what to do based on the presence/absence of keywords
                 */
                nextLine = ChemParser.readMeaningfulLine(data, true);
                boolean continueToReadRxn = true;
                // Initialize all of the possible pdep variables
                LinkedHashMap thirdBodyList = new LinkedHashMap();
                UncertainDouble uA = new UncertainDouble(0.0, 0.0, "Adder");
                UncertainDouble un = new UncertainDouble(0.0, 0.0, "Adder");
                UncertainDouble uE = new UncertainDouble(0.0, 0.0, "Adder");
                ArrheniusKinetics low = new ArrheniusKinetics(uA, un, uE, "",
                        0, "", "");
                double a = 0.0;
                double T3star = 0.0;
                double Tstar = 0.0;
                double T2star = 0.0;
                boolean troe7 = false;
                int numPLOGs = 0;
                PDepArrheniusKinetics pdepkineticsPLOG = new PDepArrheniusKinetics(
                        numPLOGs);
                /*
                 * When reading in the auxillary information for the pdep reactions, let's not assume the order is fixed
                 * (i.e. third-bodies and their efficiencies, then lindemann, then troe) The order of the if statement
                 * is important as the "troe" and "low" lines will also contain a "/"; thus, the elseif contains "/"
                 * needs to be last.
                 */
                /*
                 * Additions by MRH on 26JUL2010: Allowing RMG to read-in PLOG and CHEB formatting
                 */
                while (continueToReadRxn) {
                    if (nextLine == null) {
                        continueToReadRxn = false;
                    } else if (nextLine.toLowerCase().contains("troe")) {
                        // read in troe parameters
                        StringTokenizer st = new StringTokenizer(nextLine, "/");
                        String temp = st.nextToken().trim(); // TROE
                        String troeString = st.nextToken().trim(); // List of troe parameters
                        st = new StringTokenizer(troeString);
                        int n = st.countTokens();
                        if (n != 3 && n != 4)
                            throw new InvalidKineticsFormatException(
                                    "Troe parameter number = " + n
                                            + " for reaction: " + r.toString());
                        a = Double.parseDouble(st.nextToken().trim());
                        T3star = Double.parseDouble(st.nextToken().trim());
                        Tstar = Double.parseDouble(st.nextToken().trim());
                        if (st.hasMoreTokens()) {
                            troe7 = true;
                            T2star = Double.parseDouble(st.nextToken().trim());
                        }
                        nextLine = ChemParser.readMeaningfulLine(data, true);
                    } else if (nextLine.toLowerCase().contains("low")) {
                        // read in lindemann parameters
                        StringTokenizer st = new StringTokenizer(nextLine, "/");
                        String temp = st.nextToken().trim(); // LOW
                        String lowString = st.nextToken().trim(); // Modified Arrhenius parameters
                        /*
                         * MRH 17Feb2010: The units of the k_zero (LOW) Arrhenius parameters are different from the
                         * units of k_inf Arrhenius parameters by a factor of cm3/mol, hence the getReactantNumber()+1
                         */
                        low = ChemParser.parseSimpleArrheniusKinetics(
                                lowString, A_multiplier, E_multiplier,
                                r.getReactantNumber() + 1);
                        nextLine = ChemParser.readMeaningfulLine(data, true);
                    } else if (nextLine.contains("CHEB")) {
                        /*
                         * CHEBYSHEV FORMAT
                         */
                        // Read in the Tmin/Tmax and Pmin/Pmax information
                        StringTokenizer st_cheb = new StringTokenizer(nextLine,
                                "/");
                        String nextToken = st_cheb.nextToken(); // Should be TCHEB or PCHEB
                        StringTokenizer st_minmax = new StringTokenizer(
                                st_cheb.nextToken());
                        double Tmin = 0.0;
                        double Tmax = 0.0;
                        double Pmin = 0.0;
                        double Pmax = 0.0;
                        if (nextToken.trim().equals("TCHEB")) {
                            Tmin = Double.parseDouble(st_minmax.nextToken());
                            Tmax = Double.parseDouble(st_minmax.nextToken());
                        } else {
                            Pmin = Double.parseDouble(st_minmax.nextToken());
                            Pmax = Double.parseDouble(st_minmax.nextToken());
                        }
                        nextToken = st_cheb.nextToken(); // Should be PCHEB or TCHEB
                        st_minmax = new StringTokenizer(st_cheb.nextToken());
                        if (nextToken.trim().equals("TCHEB")) {
                            Tmin = Double.parseDouble(st_minmax.nextToken());
                            Tmax = Double.parseDouble(st_minmax.nextToken());
                        } else {
                            Pmin = Double.parseDouble(st_minmax.nextToken());
                            Pmax = Double.parseDouble(st_minmax.nextToken());
                        }
                        // Read in the N/M values (number of polynomials in the Temp and Press domain)
                        nextLine = ChemParser.readMeaningfulLine(data, true);
                        st_cheb = new StringTokenizer(nextLine, "/");
                        nextToken = st_cheb.nextToken(); // Should be CHEB
                        st_minmax = new StringTokenizer(st_cheb.nextToken());
                        int numN = Integer.parseInt(st_minmax.nextToken());
                        int numM = Integer.parseInt(st_minmax.nextToken());
                        // Read in the coefficients
                        nextLine = ChemParser.readMeaningfulLine(data, true);
                        double[] unorderedChebyCoeffs = new double[numN * numM];
                        int chebyCoeffCounter = 0;
                        while (nextLine != null && nextLine.contains("CHEB")) {
                            st_cheb = new StringTokenizer(nextLine, "/");
                            nextToken = st_cheb.nextToken(); // Should be CHEB
                            st_minmax = new StringTokenizer(st_cheb.nextToken());
                            while (st_minmax.hasMoreTokens()) {
                                unorderedChebyCoeffs[chebyCoeffCounter] = Double
                                        .parseDouble(st_minmax.nextToken());
                                ++chebyCoeffCounter;
                            }
                            nextLine = ChemParser
                                    .readMeaningfulLine(data, true);
                        }
                        // Order the chebyshev coefficients
                        double[][] chebyCoeffs = new double[numN][numM];
                        for (int numRows = 0; numRows < numN; numRows++) {
                            for (int numCols = 0; numCols < numM; numCols++) {
                                chebyCoeffs[numRows][numCols] = unorderedChebyCoeffs[numM
                                        * numRows + numCols];
                            }
                        }
                        // Make the ChebyshevPolynomials, PDepReaction, and add to list of PDepReactions
                        ChebyshevPolynomials coeffs = new ChebyshevPolynomials(
                                numN, new Temperature(Tmin, "K"),
                                new Temperature(Tmax, "K"), numM, new Pressure(
                                        Pmin, "bar"),
                                new Pressure(Pmax, "bar"), chebyCoeffs);
                        PDepIsomer reactants = new PDepIsomer(r.getStructure()
                                .getReactantList());
                        PDepIsomer products = new PDepIsomer(r.getStructure()
                                .getProductList());
                        PDepRateConstant pdepRC = new PDepRateConstant(coeffs);
                        PDepReaction pdeprxn = new PDepReaction(reactants,
                                products, pdepRC);
                        pdeprxn.setHighPKinetics(r.getKinetics()[0]);
                        r = pdeprxn;
                        r.setIsFromPrimaryKineticLibrary(true);
                        allPdepNetworks.add(pdeprxn);
                        continueToReadRxn = false;
                    } else if (nextLine.contains("PLOG")) {
                        /*
                         * PLOG FORMAT
                         */
                        while (nextLine != null && nextLine.contains("PLOG")) {
                            // Increase the PLOG counter
                            ++numPLOGs;
                            // Store the previous PLOG information in temporary arrays
                            Pressure[] previousPressures = new Pressure[numPLOGs];
                            ArrheniusKinetics[] previousKinetics = new ArrheniusKinetics[numPLOGs];
                            for (int previousNumPLOG = 0; previousNumPLOG < numPLOGs - 1; ++previousNumPLOG) {
                                previousPressures[previousNumPLOG] = pdepkineticsPLOG
                                        .getPressure(previousNumPLOG);
                                previousKinetics[previousNumPLOG] = pdepkineticsPLOG
                                        .getKinetics(previousNumPLOG);
                            }
                            // Read in the new PLOG information, and add this to the temporary array
                            PDepArrheniusKinetics newpdepkinetics = parsePLOGline(
                                    nextLine, A_multiplier, E_multiplier,
                                    r.getReactantNumber());
                            previousPressures[numPLOGs - 1] = newpdepkinetics
                                    .getPressure(0);
                            previousKinetics[numPLOGs - 1] = newpdepkinetics
                                    .getKinetics(0);
                            // Re-initialize pdepkinetics and populate with stored information
                            pdepkineticsPLOG = new PDepArrheniusKinetics(
                                    numPLOGs);
                            pdepkineticsPLOG.setPressures(previousPressures);
                            pdepkineticsPLOG
                                    .setRateCoefficients(previousKinetics);
                            // Read the next line
                            nextLine = ChemParser
                                    .readMeaningfulLine(data, true);
                        }
                        // read DUPLICATE tag if it's there
                        boolean expecting_duplicate = false;
                        if (nextLine != null
                                && nextLine.toLowerCase().startsWith("dup")) {
                            expecting_duplicate = true;
                            nextLine = ChemParser
                                    .readMeaningfulLine(data, true);
                        }
                        // done with this reaction
                        continueToReadRxn = false;
                        // Make the PDepReaction
                        PDepIsomer reactants = new PDepIsomer(r.getStructure()
                                .getReactantList());
                        PDepIsomer products = new PDepIsomer(r.getStructure()
                                .getProductList());
                        PDepRateConstant pdepRC = new PDepRateConstant(
                                pdepkineticsPLOG);
                        PDepReaction pdeprxn = new PDepReaction(reactants,
                                products, pdepRC);
                        pdeprxn.setHighPKinetics(r.getKinetics()[0]);
                        r = pdeprxn;
                        r.setIsFromPrimaryKineticLibrary(true);
                        if (expecting_duplicate)
                            r.setExpectDuplicate(true);
                        // see if we've already made this reaction in this seed mechanism
                        Iterator allRxnsIter = localReactions.iterator();
                        boolean foundRxn = false;
                        while (allRxnsIter.hasNext()) {
                            Reaction old = (Reaction) allRxnsIter.next();
                            if (old.equals(r)) {
                                if (!expecting_duplicate)
                                    throw new InvalidKineticsFormatException(
                                            "Unexpected duplicate reaction found for "
                                                    + r.toString());
                                PDepReaction oldPDep = (PDepReaction) old; // we should be able to cast it
                                oldPDep.addAdditionalKinetics(
                                        r.getKinetics()[0], 1, true); // high-P limit kinetics
                                oldPDep.addPDepArrheniusKinetics(pdepkineticsPLOG); // PLOG kinetics
                                continue read; // break out of inner loops and read the next reaction
                            }
                        }
                        // Add to the list of PDepReactions
                        allPdepNetworks.add(pdeprxn);
                        // Re-initialize the pdepkinetics variable
                        numPLOGs = 0;
                        pdepkineticsPLOG = new PDepArrheniusKinetics(numPLOGs);
                    } else if (nextLine.contains("/")) {
                        /*
                         * Collision efficiencies
                         */
                        // read in third body colliders + efficiencies
                        thirdBodyList.putAll(ChemParser.parseThirdBodyList(
                                nextLine, allSpecies));
                        nextLine = ChemParser.readMeaningfulLine(data, true);
                    } else if (nextLine.toLowerCase().startsWith("dup")) {
                        // read DUPLICATE tag
                        throw new InvalidKineticsFormatException(
                                "DUPLICATE pdep kinetics only handled for PLOG form. "
                                        + r.toString());
                    } else {
                        // the nextLine is a "new" reaction, hence we need to exit the while loop
                        continueToReadRxn = false;
                    }
                }
                /*
                 * Make the pdep reaction, according to which parameters are present
                 */
                Reaction tbr;
                if (r instanceof PDepReaction) {
                    // Chebyshev or PLog
                    tbr = r;
                    // For other forms this is stored in tbr.kinetics[0].source rather than trb.comments,
                    // but the PDepReaction class is different, and the net effect is the same:
                    tbr.setComments(source + p_name);
                } else if ((a == 0.0) && (T3star == 0.0) && (Tstar == 0.0)) {
                    // Not a troe reaction
                    if (low.getAValue() == 0.0) {
                        // thirdbody reaction
                        tbr = ThirdBodyReaction.make(r, thirdBodyList);
                    } else {
                        // lindemann reaction
                        tbr = LindemannReaction.make(r, thirdBodyList, low);
                    }
                } else {
                    // troe reaction
                    tbr = TROEReaction.make(r, thirdBodyList, low, a, T3star,
                            Tstar, troe7, T2star);
                }
                tbr.setKineticsSource(source + p_name, 0);
                tbr.setKineticsComments(" ", 0);
                // RHW thinks this should always be set - it later determines if the rate is multiplied by the reaction
// path degeneracy or not.
                tbr.setIsFromPrimaryKineticLibrary(true);
                (tbr.getKinetics())[0].setFromPrimaryKineticLibrary(true);
                localReactions.add(tbr);
                Reaction reverse = tbr.getReverseReaction();
                if (reverse != null)
                    localReactions.add(reverse);
                if (!tbr.repOk())
                    throw new RuntimeException(String.format(
                            "Something wrong with reaction %s", tbr));
            }
            // = check duplicates are correctly labeled
            for (Iterator<Reaction> allRxnsIter = localReactions.iterator(); allRxnsIter
                    .hasNext();) {
                Reaction r = allRxnsIter.next();
                if (!(r instanceof PDepReaction))
                    continue; // ThirdBodyReactions can not be duplicates at this time
                PDepReaction rxn = (PDepReaction) r;
                PDepRateConstant pDepRate = rxn.getPDepRate();
                if (rxn.getExpectDuplicate()) {
                    if (pDepRate.getMode() != PDepRateConstant.Mode.PDEPARRHENIUS)
                        throw new InvalidKineticsFormatException(
                                "DUPLICATE pdep kinetics only handled for PLOG form.");
                    if (pDepRate.getPDepArrheniusKinetics().length <= 1)
                        throw new InvalidKineticsFormatException(
                                "Was expecting DULPLICATE kinetics but none found for "
                                        + rxn.toString());
                }
                if (pDepRate.getMode() == PDepRateConstant.Mode.PDEPARRHENIUS) {
                    if (pDepRate.getPDepArrheniusKinetics().length > 1) {
                        if (!rxn.getExpectDuplicate())
                            throw new InvalidKineticsFormatException(
                                    "Found unlabeled duplicate kinetics for "
                                            + rxn.toString());
                    }
                }
            }
            // =
            in.close();
            return localReactions;
        } catch (Exception e) {
            /*
             * 25Jun2009-MRH: When reading the Primary Kinetic Library, we should not require the user to supply troe
             * reactions. In the instance that no "troeReactions.txt" file exists, inform user of this but continue
             * simulation.
             */
            Logger.logStackTrace(e);
            Logger.error("RMG could not find, or read in its entirety, the pressure-dependent reactions file (pdepreactions.txt)"
                    + "\n\tin the " + source + p_name + "\n" + e.getMessage());
            return null;
        }
    }

    public int size() {
        return speciesSet.size();
    }

    public String getName() {
        return name;
    }

    public void setName(String p_name) {
        name = p_name;
    }

    public LinkedHashSet getReactionSet() {
        return reactionSet;
    }

    public boolean shouldGenerateReactions() {
        return generateReactions;
    }

    public void setGenerateReactions(boolean generateReactions) {
        this.generateReactions = generateReactions;
    }

    public double[] parseReactionRateUnits(BufferedReader data) {
        double[] multipliers = new double[2];
        String line = ChemParser.readMeaningfulLine(data, true);
        if (line.startsWith("Unit")) {
            line = ChemParser.readMeaningfulLine(data, true);
            unit: while (!(line.startsWith("Reaction"))) {
                if (line.startsWith("A")) {
                    StringTokenizer st = new StringTokenizer(line);
                    String temp = st.nextToken();
                    String unit = st.nextToken().trim();
                    if (unit.compareToIgnoreCase("mol/cm3/s") == 0) {
                        multipliers[0] = 1;
                    } else if (unit.compareToIgnoreCase("mol/liter/s") == 0) {
                        multipliers[0] = 1e3;
                    } else if (unit.compareToIgnoreCase("molecule/cm3/s") == 0) {
                        multipliers[0] = 6.022e23;
                    }
                } else if (line.startsWith("E")) {
                    StringTokenizer st = new StringTokenizer(line);
                    String temp = st.nextToken();
                    String unit = st.nextToken().trim();
                    if (unit.compareToIgnoreCase("kcal/mol") == 0) {
                        multipliers[1] = 1;
                    } else if (unit.compareToIgnoreCase("cal/mol") == 0) {
                        multipliers[1] = 1e-3;
                    } else if (unit.compareToIgnoreCase("kJ/mol") == 0) {
                        multipliers[1] = 1 / 4.186;
                    } else if (unit.compareToIgnoreCase("J/mol") == 0) {
                        multipliers[1] = 1 / 4186;
                    } else if (unit.compareToIgnoreCase("Kelvin") == 0) {
                        multipliers[1] = 1.987e-3;
                    }
                }
                line = ChemParser.readMeaningfulLine(data, true);
            }
        }
        return multipliers;
    }

    public LinkedHashSet updateReactionList(Reaction r,
            LinkedHashSet listOfRxns, boolean generateReverse, String nextLine) {
        Iterator allRxnsIter = listOfRxns.iterator();
        boolean foundRxn = false;
        while (allRxnsIter.hasNext()) {
            Reaction old = (Reaction) allRxnsIter.next();
            if (old.equals(r)) {
                if (old.getExpectDuplicate()
                        && nextLine.toLowerCase().startsWith("dup")) {
                    old.addAdditionalKinetics(r.getKinetics()[0], 1, true);
                    foundRxn = true;
                    break;
                } else {
                    Logger.critical("Unspecified duplicate reaction found in library:\n\t"
                            + r.toString()
                            + "\nPlease add 'DUP' tag after each instance of desired duplicate reaction");
                    System.exit(0);
                }
            }
        }
        if (!foundRxn) {
            if (nextLine != null && nextLine.toLowerCase().startsWith("dup")) {
                r.setExpectDuplicate(true);
            }
            listOfRxns.add(r);
            if (generateReverse) {
                Reaction reverse = r.getReverseReaction();
                if (reverse != null) {
                    listOfRxns.add(reverse);
                }
            }
        }
        return listOfRxns;
    }

    public PDepArrheniusKinetics parsePLOGline(String line,
            double A_multiplier, double E_multiplier, int numReac) {
        PDepRateConstant pdepk = new PDepRateConstant();
        // Delimit the PLOG line by "/"
        StringTokenizer st_slash = new StringTokenizer(line, "/");
        String dummyString = st_slash.nextToken();
        // Delimit the data of the PLOG line by whitespace
        StringTokenizer st = new StringTokenizer(st_slash.nextToken());
        // If reading a chemkin plog line, pressure MUST be in atmospheres
        Pressure p = new Pressure(Double.parseDouble(st.nextToken().trim()),
                "atm");
        // ***note: PLOG uses the same units for Ea and A as Arrhenius expressions; this has been a persistent source of
// confusion; see
// https://github.com/GreenGroup/RMG-Java/commit/2947e7b8d5b1e3e19543f2489990fa42e43ecad2#commitcomment-844009
        double A = Double.parseDouble(st.nextToken().trim());
        // convert the units (cf. similar code in ChemParser)
        if (numReac == 1) {
            // do nothing, no conversion needed
        } else if (numReac == 2) {
            A = A * A_multiplier;
        } else if (numReac == 3) {
            A = A * A_multiplier * A_multiplier;
        } else {
            Logger.error("Unsupported number of reactants (" + numReac + "): "
                    + line);
            System.exit(0);
        }
        UncertainDouble dA = new UncertainDouble(A, 0.0, "A");
        UncertainDouble dn = new UncertainDouble(Double.parseDouble(st
                .nextToken().trim()), 0.0, "A");
        double Ea = Double.parseDouble(st.nextToken().trim());
        UncertainDouble dE = new UncertainDouble(Ea * E_multiplier, 0.0, "A");
        ArrheniusKinetics k = new ArrheniusKinetics(dA, dn, dE, "", 1, "", "");
        PDepArrheniusKinetics pdepAK = new PDepArrheniusKinetics(1);
        pdepAK.setKinetics(0, p, k);
        return pdepAK;
    }

    public LinkedList getPDepNetworks() {
        return allPdepNetworks;
    }
}

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
import jing.rxnSys.ReactionSystem;
import jing.rxn.*;
import jing.chem.*;
import java.util.*;
import jing.mathTool.UncertainDouble;
import jing.param.*;
import jing.chemUtil.*;
import jing.chemParser.*;

// ## package jing::rxnSys
// ----------------------------------------------------------------------------
// jing\rxnSys\ReactionModelGenerator.java
// ----------------------------------------------------------------------------
// ## class ReactionModelGenerator
public class ReactionModelGenerator {
    protected LinkedList timeStep;
    protected ReactionModel reactionModel;
    protected String workingDirectory;
    protected LinkedList reactionSystemList;
    protected int paraInfor;
    protected boolean error;
    protected boolean sensitivity;
    protected LinkedList species;
    protected LinkedList initialStatusList;
    protected double rtol;
    protected static double atol;
    protected PrimaryKineticLibrary primaryKineticLibrary;
    protected ReactionLibrary ReactionLibrary;
    protected ReactionModelEnlarger reactionModelEnlarger;
    protected LinkedHashSet speciesSeed;
    protected ReactionGenerator reactionGenerator;
    protected LibraryReactionGenerator lrg;
    protected LinkedList tempList;
    protected LinkedList presList;
    protected LinkedList validList;
    protected LinkedList initList = new LinkedList();
    protected LinkedList beginList = new LinkedList();
    protected LinkedList endList = new LinkedList();
    protected LinkedList lastTList = new LinkedList();
    protected LinkedList currentTList = new LinkedList();
    protected LinkedList lastPList = new LinkedList();
    protected LinkedList currentPList = new LinkedList();
    protected LinkedList conditionChangedList = new LinkedList();
    protected LinkedList reactionChangedList = new LinkedList();
    protected int numConversions;
    protected String equationOfState;
    // 24Jun2009 MRH: variable stores the first temperature encountered in the condition.txt file
    // This temperature is used to select the "best" kinetics from the rxn library
    protected static Temperature temp4BestKinetics;
    protected static boolean useDiffusion;
    protected static boolean useSolvation;
    protected SeedMechanism seedMechanism = null;
    protected PrimaryThermoLibrary primaryThermoLibrary;
    protected PrimaryTransportLibrary primaryTransportLibrary;
    protected PrimaryAbrahamLibrary primaryAbrahamLibrary;
    protected static SolventData solvent;
    protected SolventLibrary solventLibrary;
    protected static double viscosity;
    protected boolean readrestart = false;
    protected boolean writerestart = false;
    protected LinkedHashSet restartCoreSpcs = new LinkedHashSet();
    protected LinkedHashSet restartEdgeSpcs = new LinkedHashSet();
    protected LinkedHashSet restartCoreRxns = new LinkedHashSet();
    protected LinkedHashSet restartEdgeRxns = new LinkedHashSet();
    // Constructors
    private LinkedHashSet specs = new LinkedHashSet();
    // public static native long getCpuTime();
    // static {System.loadLibrary("cpuTime");}
    public static boolean rerunFame = false;
    public static boolean keepQMfiles = false;
    protected static double tolerance;// can be interpreted as "coreTol" (vs. edgeTol)
    protected static double termTol;
    protected static double edgeTol;
    protected static int minSpeciesForPruning;
    protected static int maxEdgeSpeciesAfterPruning;
    public int limitingReactantID = 1;
    public int numberOfEquivalenceRatios = 0;

    // ## operation ReactionModelGenerator()
    public ReactionModelGenerator() {
        workingDirectory = System.getProperty("RMG.workingDirectory");
    }

    // ## operation initializeReactionSystem()
    // 10/24/07 gmagoon: changed name to initializeReactionSystems
    public void initializeReactionSystems() throws InvalidSymbolException,
            IOException {
        // #[ operation initializeReactionSystem()
        try {
            String initialConditionFile = System
                    .getProperty("jing.rxnSys.ReactionModelGenerator.conditionFile");
            if (initialConditionFile == null) {
                Logger.critical("undefined system property: jing.rxnSys.ReactionModelGenerator.conditionFile");
                System.exit(0);
            }
            // double sandeep = getCpuTime();
            // System.out.println(getCpuTime()/1e9/60);
            FileReader in = new FileReader(initialConditionFile);
            BufferedReader reader = new BufferedReader(in);
            // TemperatureModel temperatureModel = null;//10/27/07 gmagoon: commented out
            // PressureModel pressureModel = null;//10/27/07 gmagoon: commented out
            // ReactionModelEnlarger reactionModelEnlarger = null;//10/9/07 gmagoon: commented out: unneeded now and
// causes scope problems
            FinishController finishController = null;
            // DynamicSimulator dynamicSimulator = null;//10/27/07 gmagoon: commented out and replaced with following
// line
            LinkedList dynamicSimulatorList = new LinkedList();
            setPrimaryKineticLibrary(null);// 10/14/07 gmagoon: changed to use setPrimaryReactionLibrary
            double[] conversionSet = new double[50];
            String line = ChemParser.readMeaningfulLine(reader, true);
            /*
             * if (line.startsWith("Restart")){ StringTokenizer st = new StringTokenizer(line); String token =
             * st.nextToken(); token = st.nextToken(); if (token.equalsIgnoreCase("true")) {
             * //Runtime.getRuntime().exec("cp Restart/allSpecies.txt Restart/allSpecies1.txt");
             * //Runtime.getRuntime().exec("echo  >> allSpecies.txt"); restart = true; } else if
             * (token.equalsIgnoreCase("false")) { Runtime.getRuntime().exec("rm Restart/allSpecies.txt"); restart =
             * false; } else throw new InvalidSymbolException("UnIdentified Symbol "+token+" after Restart:"); } else
             * throw new InvalidSymbolException("Can't find Restart!");
             */
            // line = ChemParser.readMeaningfulLine(reader);
            if (line.startsWith("Database")) {// svp
                line = ChemParser.readMeaningfulLine(reader, true);
            } else
                throw new InvalidSymbolException("Can't find database!");
            // if (line.startsWith("PrimaryThermoLibrary")){//svp
            // line = ChemParser.readMeaningfulLine(reader);
            // }
            // else throw new InvalidSymbolException("Can't find primary thermo library!");
            /*
             * Added by MRH on 15-Jun-2009 Give user the option to change the maximum carbon, oxygen, and/or radical
             * number for all species. These lines will be optional in the condition.txt file. Values are hard- coded
             * into RMG (in ChemGraph.java), but any user- defined input will override these values.
             */
            /*
             * Moved from before InitialStatus to before PrimaryThermoLibary by MRH on 27-Oct-2009 Overriding default
             * values of maximum number of "X" per chemgraph should come before RMG attempts to make any chemgraph. The
             * first instance RMG will attempt to make a chemgraph is in reading the primary thermo library.
             */
            line = readMaxAtomTypes(line, reader);
// if (line.startsWith("MaxCarbonNumber")) {
// StringTokenizer st = new StringTokenizer(line);
// String dummyString = st.nextToken(); // This should hold "MaxCarbonNumberPerSpecies:"
// int maxCNum = Integer.parseInt(st.nextToken());
// ChemGraph.setMaxCarbonNumber(maxCNum);
// System.out.println("Note: Overriding RMG-defined MAX_CARBON_NUM with user-defined value: " + maxCNum);
// line = ChemParser.readMeaningfulLine(reader);
// }
// if (line.startsWith("MaxOxygenNumber")) {
// StringTokenizer st = new StringTokenizer(line);
// String dummyString = st.nextToken(); // This should hold "MaxOxygenNumberPerSpecies:"
// int maxONum = Integer.parseInt(st.nextToken());
// ChemGraph.setMaxOxygenNumber(maxONum);
// System.out.println("Note: Overriding RMG-defined MAX_OXYGEN_NUM with user-defined value: " + maxONum);
// line = ChemParser.readMeaningfulLine(reader);
// }
// if (line.startsWith("MaxRadicalNumber")) {
// StringTokenizer st = new StringTokenizer(line);
// String dummyString = st.nextToken(); // This should hold "MaxRadicalNumberPerSpecies:"
// int maxRadNum = Integer.parseInt(st.nextToken());
// ChemGraph.setMaxRadicalNumber(maxRadNum);
// System.out.println("Note: Overriding RMG-defined MAX_RADICAL_NUM with user-defined value: " + maxRadNum);
// line = ChemParser.readMeaningfulLine(reader);
// }
// if (line.startsWith("MaxSulfurNumber")) {
// StringTokenizer st = new StringTokenizer(line);
// String dummyString = st.nextToken(); // This should hold "MaxSulfurNumberPerSpecies:"
// int maxSNum = Integer.parseInt(st.nextToken());
// ChemGraph.setMaxSulfurNumber(maxSNum);
// System.out.println("Note: Overriding RMG-defined MAX_SULFUR_NUM with user-defined value: " + maxSNum);
// line = ChemParser.readMeaningfulLine(reader);
// }
// if (line.startsWith("MaxSiliconNumber")) {
// StringTokenizer st = new StringTokenizer(line);
// String dummyString = st.nextToken(); // This should hold "MaxSiliconNumberPerSpecies:"
// int maxSiNum = Integer.parseInt(st.nextToken());
// ChemGraph.setMaxSiliconNumber(maxSiNum);
// System.out.println("Note: Overriding RMG-defined MAX_SILICON_NUM with user-defined value: " + maxSiNum);
// line = ChemParser.readMeaningfulLine(reader);
// }
// if (line.startsWith("MaxHeavyAtom")) {
// StringTokenizer st = new StringTokenizer(line);
// String dummyString = st.nextToken(); // This should hold "MaxHeavyAtomPerSpecies:"
// int maxHANum = Integer.parseInt(st.nextToken());
// ChemGraph.setMaxHeavyAtomNumber(maxHANum);
// System.out.println("Note: Overriding RMG-defined MAX_HEAVYATOM_NUM with user-defined value: " + maxHANum);
// line = ChemParser.readMeaningfulLine(reader);
// }
            /*
             * Read in the Primary Thermo Library MRH 7-Jul-2009
             */
            if (line.startsWith("PrimaryThermoLibrary:")) {
                /*
                 * MRH 27Feb2010: Changing the "read in Primary Thermo Library information" code into it's own method.
                 * Other modules (e.g. PopulateReactions) will be utilizing the exact code. Rather than copying and
                 * pasting code into other modules, just have everything call this new method: readAndMakePTL
                 */
                readAndMakePTL(reader);
            } else
                throw new InvalidSymbolException(
                        "Error reading condition.txt file: "
                                + "Could not locate PrimaryThermoLibrary field");
            line = ChemParser.readMeaningfulLine(reader, true);
            /*
             * MRH 17-May-2010: Added primary transport library field
             */
            if (line.toLowerCase().startsWith("primarytransportlibrary")) {
                readAndMakePTransL(reader);
            } else
                throw new InvalidSymbolException(
                        "Error reading condition.txt file: "
                                + "Could not locate PrimaryTransportLibrary field.");
            line = ChemParser.readMeaningfulLine(reader, true);
            // Extra forbidden structures may be specified after the Primary Thermo Library
            if (line.startsWith("ForbiddenStructures:")) {
                readExtraForbiddenStructures(reader);
                line = ChemParser.readMeaningfulLine(reader, true);
            }
            if (line.toLowerCase().startsWith("readrestart")) {
                StringTokenizer st = new StringTokenizer(line);
                String tempString = st.nextToken(); // "ReadRestart:"
                tempString = st.nextToken();
                if (tempString.toLowerCase().equals("yes")) {
                    readrestart = true;
                    readRestartSpecies();
                    readRestartReactions();
                } else
                    readrestart = false;
                line = ChemParser.readMeaningfulLine(reader, true);
            } else
                throw new InvalidSymbolException(
                        "Cannot locate ReadRestart field");
            if (line.toLowerCase().startsWith("writerestart")) {
                StringTokenizer st = new StringTokenizer(line);
                String tempString = st.nextToken(); // "WriteRestart:"
                tempString = st.nextToken();
                if (tempString.toLowerCase().equals("yes"))
                    writerestart = true;
                else
                    writerestart = false;
                line = ChemParser.readMeaningfulLine(reader, true);
            } else
                throw new InvalidSymbolException(
                        "Cannot locate WriteRestart field");
            // read temperature model
            // gmagoon 10/23/07: modified to handle multiple temperatures; note that this requires different formatting
// of units in condition.txt
            if (line.startsWith("TemperatureModel:")) {
                createTModel(line);
// StringTokenizer st = new StringTokenizer(line);
// String name = st.nextToken();
// String modelType = st.nextToken();
// //String t = st.nextToken();
// String unit = st.nextToken();
// unit = ChemParser.removeBrace(unit);
// if (modelType.equals("Constant")) {
// tempList = new LinkedList();
// //read first temperature
// double t = Double.parseDouble(st.nextToken());
// tempList.add(new ConstantTM(t, unit));
// Temperature temp = new Temperature(t, unit);//10/29/07 gmagoon: added this line and next two lines to set
// Global.lowTemperature and Global.highTemperature
// Global.lowTemperature = (Temperature)temp.clone();
// Global.highTemperature = (Temperature)temp.clone();
// //read remaining temperatures
// while (st.hasMoreTokens()) {
// t = Double.parseDouble(st.nextToken());
// tempList.add(new ConstantTM(t, unit));
// temp = new Temperature(t,unit);//10/29/07 gmagoon: added this line and next two "if" statements to set
// Global.lowTemperature and Global.highTemperature
// if(temp.getK() < Global.lowTemperature.getK())
// Global.lowTemperature = (Temperature)temp.clone();
// if(temp.getK() > Global.highTemperature.getK())
// Global.highTemperature = (Temperature)temp.clone();
// }
// // Global.temperature = new Temperature(t,unit);
// }
                // 10/23/07 gmagoon: commenting out; further updates needed to get this to work
                // else if (modelType.equals("Curved")) {
                // String t = st.nextToken();
                // // add reading curved temperature function here
                // temperatureModel = new CurvedTM(new LinkedList());
                // }
// else {
// throw new InvalidSymbolException("condition.txt: Unknown TemperatureModel = " + modelType);
// }
            } else
                throw new InvalidSymbolException(
                        "condition.txt: can't find TemperatureModel!");
            // read in pressure model
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("PressureModel:")) {
                createPModel(line);
// StringTokenizer st = new StringTokenizer(line);
// String name = st.nextToken();
// String modelType = st.nextToken();
// //String p = st.nextToken();
// String unit = st.nextToken();
// unit = ChemParser.removeBrace(unit);
// if (modelType.equals("Constant")) {
// presList = new LinkedList();
// //read first pressure
// double p = Double.parseDouble(st.nextToken());
// Pressure pres = new Pressure(p, unit);
// Global.lowPressure = (Pressure)pres.clone();
// Global.highPressure = (Pressure)pres.clone();
// presList.add(new ConstantPM(p, unit));
// //read remaining temperatures
// while (st.hasMoreTokens()) {
// p = Double.parseDouble(st.nextToken());
// presList.add(new ConstantPM(p, unit));
// pres = new Pressure(p, unit);
// if(pres.getBar() < Global.lowPressure.getBar())
// Global.lowPressure = (Pressure)pres.clone();
// if(pres.getBar() > Global.lowPressure.getBar())
// Global.highPressure = (Pressure)pres.clone();
// }
// //Global.pressure = new Pressure(p, unit);
// }
// //10/23/07 gmagoon: commenting out; further updates needed to get this to work
// //else if (modelType.equals("Curved")) {
// // // add reading curved pressure function here
// // pressureModel = new CurvedPM(new LinkedList());
// //}
// else {
// throw new InvalidSymbolException("condition.txt: Unknown PressureModel = " + modelType);
// }
            } else
                throw new InvalidSymbolException(
                        "condition.txt: can't find PressureModel!");
            // after PressureModel comes an optional line EquationOfState
            // if "EquationOfState: Liquid" is found then initial concentrations are assumed to be correct
            // if it is ommited, then initial concentrations are normalised to ensure PV=NRT (ideal gas law)
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("EquationOfState")) {
                StringTokenizer st = new StringTokenizer(line);
                String name = st.nextToken();
                String eosType = st.nextToken().toLowerCase();
                if (eosType.equals("liquid")) {
                    equationOfState = "Liquid";
                    Logger.info("Equation of state: Liquid. Relying on concentrations in input file to get density correct; not checking PV=NRT");
                }
                line = ChemParser.readMeaningfulLine(reader, true);
            }
            // Read in InChI generation
            if (line.startsWith("InChIGeneration:")) {
                StringTokenizer st = new StringTokenizer(line);
                String name = st.nextToken();
                String inchiOnOff = st.nextToken().toLowerCase();
                if (inchiOnOff.equals("on")) {
                    Species.useInChI = true;
                } else if (inchiOnOff.equals("off")) {
                    Species.useInChI = false;
                } else
                    throw new InvalidSymbolException(
                            "condition.txt: Unknown InChIGeneration flag: "
                                    + inchiOnOff);
                line = ChemParser.readMeaningfulLine(reader, true);
            }
            // Read in Solvation effects
            if (line.startsWith("Solvation:")) {
                StringTokenizer st = new StringTokenizer(line);
                String name = st.nextToken();
                String solvationOnOff = st.nextToken().toLowerCase();
                if (solvationOnOff.equals("on")) {
                    setUseSolvation(true);
                    Species.useSolvation = true;
                    readAndMakePAL();
                    String solventname = st.nextToken().toLowerCase();
                    readAndMakeSL(solventname);
                    Logger.info(String
                            .format("Using solvation corrections to thermochemsitry with solvent properties of %s",
                                    solventname));
                } else if (solvationOnOff.startsWith("off")) {
                    setUseSolvation(false);
                    Species.useSolvation = false;
                } else
                    throw new InvalidSymbolException(
                            "condition.txt: Unknown solvation flag: "
                                    + solvationOnOff);
                line = ChemParser.readMeaningfulLine(reader, true);
            }
            // Read in Diffusion effects
            // If 'Diffusion' is 'on' then override the settings made by the solvation flag and sets solvation 'on'
            if (line.startsWith("Diffusion:")) {
                StringTokenizer st = new StringTokenizer(line);
                String name = st.nextToken();
                String diffusionOnOff = st.nextToken().toLowerCase();
                if (diffusionOnOff.equals("on")) {
                    String viscosity_str = st.nextToken();
                    viscosity = Double.parseDouble(viscosity_str);
                    setUseDiffusion(true);
                    Logger.info(String
                            .format("Using diffusion corrections to kinetics with solvent viscosity of %.3g Pa.s.",
                                    viscosity));
                } else if (diffusionOnOff.equals("off")) {
                    setUseDiffusion(false);
                } else
                    throw new InvalidSymbolException(
                            "condition.txt: Unknown diffusion flag: "
                                    + diffusionOnOff);
                line = ChemParser.readMeaningfulLine(reader, true);// read in reactants or thermo line
            }
            /*
             * AJ 12JULY2010: Right now we do not want RMG to throw an exception if it cannot find a diffusion flag
             */
            // else throw new InvalidSymbolException("condition.txt: Cannot find diffusion flag.");
            // Should have already read in reactants or thermo line into variable 'line'
            // Read in optional QM thermo generation
            if (line.startsWith("ThermoMethod:")) {
                StringTokenizer st = new StringTokenizer(line);
                String name = st.nextToken();
                String thermoMethod = st.nextToken().toLowerCase();
                ChemGraph.TDMETHOD = thermoMethod;
                if (thermoMethod.contains("qm")
                        || thermoMethod.contains("hybrid")) {
                    if (st.hasMoreTokens()) {// override the default qmprogram ("both") if there are more; current
// options: "gaussian03" and "mopac" and of course, "both"
                        QMTP.qmprogram = st.nextToken().toLowerCase();
                    }
                    // line=ChemParser.readMeaningfulLine(reader, true);
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("MaxRadNumForQM:")) {
                        StringTokenizer st3 = new StringTokenizer(line);
                        String nameRadNum = st3.nextToken();
                        Global.maxRadNumForQM = Integer.parseInt(st3
                                .nextToken());
                    } else {
                        Logger.critical("condition.txt: Can't find 'MaxRadNumForQM:' field");
                        System.exit(0);
                    }
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("CheckConnectivity:")) {
                        StringTokenizer st4 = new StringTokenizer(line);
                        String nameCheckConnectivity = st4.nextToken();
                        String checkConnSetting = st4.nextToken().toLowerCase();
                        if (checkConnSetting.equals("off")) {// no connectivity checking
                            QMTP.connectivityCheck = 0;
                        } else if (checkConnSetting.equals("check")) {// print a warning if the connectivity doesn't
// appear to match
                            QMTP.connectivityCheck = 1;
                        } else if (checkConnSetting.equals("confirm")) {// consider the run a failure if the
// connectivity doesn't appear to match
                            QMTP.connectivityCheck = 2;
                        } else {
                            Logger.critical("condition.txt: Inappropriate 'CheckConnectivity' value (should be 'off', 'check', or 'confirm')");
                            System.exit(0);
                        }
                    } else {
                        Logger.critical("condition.txt: Can't find 'CheckConnectivity:' field (should be 'off', 'check', or 'confirm')");
                        System.exit(0);
                    }
                    line = ChemParser.readMeaningfulLine(reader, true); // read in either QM 'KeepQMFiles:' option or
// 'InitialStatus' line.
                    if (line.startsWith("KeepQMFiles:")) {
                        StringTokenizer st5 = new StringTokenizer(line);
                        String nameQmVerbose = st5.nextToken(); // String Verbose
                        String checkQmVerbose = st5.nextToken().toLowerCase();
                        if (checkQmVerbose.equals("no")) {
                            QMTP.keepQMfiles = false;
                        } else if (!checkQmVerbose.equals("yes")) {
                            Logger.critical("condition.txt: QMTP 'KeepQMFiles' field should be 'yes' or 'no'");
                            System.exit(0);
                        }
                        // Read another line
                        line = ChemParser.readMeaningfulLine(reader, true);
                    } else {
                        Logger.critical("Can't find QMTP 'KeepQMFiles:' field. Defaulting to 'yes'.");
                        QMTP.keepQMfiles = true;
                    }
                }// otherwise, the flag useQM will remain false by default and the traditional group additivity approach
// will be used
            }
            // // Read in Solvation effects
            // if (line.startsWith("Solvation:")) {
            // StringTokenizer st = new StringTokenizer(line);
            // String name = st.nextToken();
            // String solvationOnOff = st.nextToken().toLowerCase();
            // if (solvationOnOff.equals("on")) {
            // Species.useSolvation = true;
            // } else if (solvationOnOff.equals("off")) {
            // Species.useSolvation = false;
            // }
            // else throw new InvalidSymbolException("condition.txt: Unknown solvation flag: " + solvationOnOff);
            // }
            // else throw new InvalidSymbolException("condition.txt: Cannot find solvation flag.");
            // read in reactants
            //
            // 10/4/07 gmagoon: moved to initializeCoreEdgeReactionModel
            // LinkedHashSet p_speciesSeed = new LinkedHashSet();//gmagoon 10/4/07: changed to p_speciesSeed
            // setSpeciesSeed(p_speciesSeed);//gmagoon 10/4/07: added
            LinkedHashMap speciesSet = new LinkedHashMap();
            /*
             * 7/Apr/2010: MRH Neither of these variables are utilized
             */
// LinkedHashMap speciesStatus = new LinkedHashMap();
// int speciesnum = 1;
            // System.out.println(line);
            if (line.startsWith("InitialStatus")) {
                speciesSet = populateInitialStatusListWithReactiveSpecies(reader);
// line = ChemParser.readMeaningfulLine(reader);
// while (!line.equals("END")) {
// StringTokenizer st = new StringTokenizer(line);
// String index = st.nextToken();
// String name = null;
// if (!index.startsWith("(")) name = index;
// else name = st.nextToken();
// //if (restart) name += "("+speciesnum+")";
// // 24Jun2009: MRH
// // Check if the species name begins with a number.
// // If so, terminate the program and inform the user to choose
// // a different name. This is implemented so that the chem.inp
// // file generated will be valid when run in Chemkin
// try {
// int doesNameBeginWithNumber = Integer.parseInt(name.substring(0,1));
// System.out.println("\nA species name should not begin with a number." +
// " Please rename species: " + name + "\n");
// System.exit(0);
// } catch (NumberFormatException e) {
// // We're good
// }
// speciesnum ++;
// if (!(st.hasMoreTokens())) throw new InvalidSymbolException("Couldn't find concentration of species: "+name);
// String conc = st.nextToken();
// double concentration = Double.parseDouble(conc);
// String unit = st.nextToken();
// unit = ChemParser.removeBrace(unit);
// if (unit.equals("mole/l") || unit.equals("mol/l") || unit.equals("mole/liter") || unit.equals("mol/liter")) {
// concentration /= 1000;
// unit = "mol/cm3";
// }
// else if (unit.equals("mole/m3") || unit.equals("mol/m3")) {
// concentration /= 1000000;
// unit = "mol/cm3";
// }
// else if (unit.equals("molecule/cm3") || unit.equals("molecules/cm3")) {
// concentration /= 6.022e23;
// }
// else if (!unit.equals("mole/cm3") && !unit.equals("mol/cm3")) {
// throw new InvalidUnitException("Species Concentration in condition.txt!");
// }
//
// //GJB to allow "unreactive" species that only follow user-defined library reactions.
// // They will not react according to RMG reaction families
// boolean IsReactive = true;
// boolean IsConstantConcentration = false;
// while (st.hasMoreTokens()) {
// String reactive = st.nextToken().trim();
// if (reactive.equalsIgnoreCase("unreactive"))
// IsReactive = false;
// if (reactive.equalsIgnoreCase("constantconcentration"))
// IsConstantConcentration=true;
// }
//
// Graph g = ChemParser.readChemGraph(reader);
// ChemGraph cg = null;
// try {
// cg = ChemGraph.make(g);
// }
// catch (ForbiddenStructureException e) {
// System.out.println("Forbidden Structure:\n" + e.getMessage());
// throw new InvalidSymbolException("A species in the input file has a forbidden structure.");
// }
// //System.out.println(name);
// Species species = Species.make(name,cg);
// species.setReactivity(IsReactive); // GJB
// species.setConstantConcentration(IsConstantConcentration);
// speciesSet.put(name, species);
// getSpeciesSeed().add(species);
// double flux = 0;
// int species_type = 1; // reacted species
// SpeciesStatus ss = new SpeciesStatus(species,species_type,concentration,flux);
// speciesStatus.put(species, ss);
// line = ChemParser.readMeaningfulLine(reader);
// }
// ReactionTime initial = new ReactionTime(0,"S");
// //10/23/07 gmagoon: modified for handling multiple temperature, pressure conditions; note: concentration within
// speciesStatus (and list of conversion values) should not need to be modified for each T,P since this is done within
// isTPCconsistent in ReactionSystem
// initialStatusList = new LinkedList();
// for (Iterator iter = tempList.iterator(); iter.hasNext(); ) {
// TemperatureModel tm = (TemperatureModel)iter.next();
// for (Iterator iter2 = presList.iterator(); iter2.hasNext(); ){
// PressureModel pm = (PressureModel)iter2.next();
// // LinkedHashMap speStat = (LinkedHashMap)speciesStatus.clone();//10/31/07 gmagoon: trying creating multiple
// instances of speciesStatus to address issues with concentration normalization (last normalization seems to apply to
// all)
// Set ks = speciesStatus.keySet();
// LinkedHashMap speStat = new LinkedHashMap();
// for (Iterator iter3 = ks.iterator(); iter3.hasNext();){//11/1/07 gmagoon: perform deep copy; (is there an easier or
// more elegant way to do this?)
// SpeciesStatus ssCopy = (SpeciesStatus)speciesStatus.get(iter3.next());
// speStat.put(ssCopy.getSpecies(),new
// SpeciesStatus(ssCopy.getSpecies(),ssCopy.getSpeciesType(),ssCopy.getConcentration(),ssCopy.getFlux()));
// }
// initialStatusList.add(new InitialStatus(speStat,tm.getTemperature(initial),pm.getPressure(initial)));
// }
// }
            } else
                throw new InvalidSymbolException(
                        "condition.txt: can't find InitialStatus!");
            // read in inert gas concentration
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("InertGas:")) {
                populateInitialStatusListWithInertSpecies(reader);
// line = ChemParser.readMeaningfulLine(reader);
// while (!line.equals("END")) {
// StringTokenizer st = new StringTokenizer(line);
// String name = st.nextToken().trim();
// String conc = st.nextToken();
// double inertConc = Double.parseDouble(conc);
// String unit = st.nextToken();
// unit = ChemParser.removeBrace(unit);
// if (unit.equals("mole/l") || unit.equals("mol/l") || unit.equals("mole/liter") || unit.equals("mol/liter")) {
// inertConc /= 1000;
// unit = "mol/cm3";
// }
// else if (unit.equals("mole/m3") || unit.equals("mol/m3")) {
// inertConc /= 1000000;
// unit = "mol/cm3";
// }
// else if (unit.equals("molecule/cm3") || unit.equals("molecules/cm3")) {
// inertConc /= 6.022e23;
// unit = "mol/cm3";
// }
// else if (!unit.equals("mole/cm3") && !unit.equals("mol/cm3")) {
// throw new InvalidUnitException("Inert Gas Concentration not recognized: " + unit);
// }
//
// //SystemSnapshot.putInertGas(name,inertConc);
// for(Iterator iter=initialStatusList.iterator();iter.hasNext(); ){//6/23/09 gmagoon: needed to change this to
// accommodate non-static inertConc
// ((InitialStatus)iter.next()).putInertGas(name,inertConc);
// }
// line = ChemParser.readMeaningfulLine(reader);
// }
            } else
                throw new InvalidSymbolException(
                        "condition.txt: can't find Inert gas concentration!");
            // read in spectroscopic data estimator
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("SpectroscopicDataEstimator:")) {
                setSpectroscopicDataMode(line);
// StringTokenizer st = new StringTokenizer(line);
// String name = st.nextToken();
// String sdeType = st.nextToken().toLowerCase();
// if (sdeType.equals("frequencygroups") || sdeType.equals("default")) {
// SpectroscopicData.mode = SpectroscopicData.Mode.FREQUENCYGROUPS;
// }
// else if (sdeType.equals("therfit") || sdeType.equals("threefrequencymodel")) {
// SpectroscopicData.mode = SpectroscopicData.Mode.THREEFREQUENCY;
// }
// else if (sdeType.equals("off") || sdeType.equals("none")) {
// SpectroscopicData.mode = SpectroscopicData.Mode.OFF;
// }
// else throw new InvalidSymbolException("condition.txt: Unknown SpectroscopicDataEstimator = " + sdeType);
            } else
                throw new InvalidSymbolException(
                        "condition.txt: can't find SpectroscopicDataEstimator!");
            // pressure dependence and related flags
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.toLowerCase().startsWith("pressuredependence:"))
                line = setPressureDependenceOptions(line, reader);
            else
                throw new InvalidSymbolException(
                        "condition.txt: can't find PressureDependence flag!");
            if (readrestart)
                if (PDepNetwork.generateNetworks)
                    readPDepNetworks();
            // include species (optional)
            /*
             * MRH 3-APR-2010: This if statement is no longer necessary and was causing an error when the
             * PressureDependence field was set to "off"
             */
// if (!PDepRateConstant.getMode().name().equals("CHEBYSHEV") &&
// !PDepRateConstant.getMode().name().equals("PDEPARRHENIUS"))
// line = ChemParser.readMeaningfulLine(reader);
            // read in finish controller
            if (line.startsWith("FinishController")) {
                line = ChemParser.readMeaningfulLine(reader, true);
                StringTokenizer st = new StringTokenizer(line);
                String index = st.nextToken();
                String goal = st.nextToken();
                String type = st.nextToken();
                TerminationTester tt;
                if (type.startsWith("Conversion")) {
                    LinkedList spc = new LinkedList();
                    while (st.hasMoreTokens()) {
                        String name = st.nextToken();
                        Species spe = (Species) speciesSet.get(name);
                        if (spe == null)
                            throw new InvalidConversionException(
                                    "Unknown reactant in 'Goal Conversion' field of input file : "
                                            + name);
                        setLimitingReactantID(spe.getID());
                        String conv = st.nextToken();
                        double conversion;
                        try {
                            if (conv.endsWith("%")) {
                                conversion = Double.parseDouble(conv.substring(
                                        0, conv.length() - 1)) / 100;
                            } else {
                                conversion = Double.parseDouble(conv);
                            }
                            conversionSet[49] = conversion;
                        } catch (NumberFormatException e) {
                            throw new NumberFormatException(
                                    "wrong number format for conversion in initial condition file!");
                        }
                        SpeciesConversion sc = new SpeciesConversion(spe,
                                conversion);
                        spc.add(sc);
                    }
                    tt = new ConversionTT(spc);
                } else if (type.startsWith("ReactionTime")) {
                    double time = Double.parseDouble(st.nextToken());
                    String unit = ChemParser.removeBrace(st.nextToken());
                    ReactionTime rt = new ReactionTime(time, unit);
                    tt = new ReactionTimeTT(rt);
                } else {
                    throw new InvalidSymbolException(
                            "condition.txt: Unknown FinishController = " + type);
                }
                line = ChemParser.readMeaningfulLine(reader, true);
                st = new StringTokenizer(line, ":");
                String temp = st.nextToken();
                String tol = st.nextToken();
                try {
                    if (tol.endsWith("%")) {
                        tolerance = Double.parseDouble(tol.substring(0,
                                tol.length() - 1)) / 100;
                    } else {
                        tolerance = Double.parseDouble(tol);
                    }
                } catch (NumberFormatException e) {
                    throw new NumberFormatException(
                            "wrong number format for conversion in initial condition file!");
                }
                ValidityTester vt = null;
                if (reactionModelEnlarger instanceof RateBasedRME)
                    vt = new RateBasedVT(tolerance);
                else if (reactionModelEnlarger instanceof RateBasedPDepRME)
                    vt = new RateBasedPDepVT(tolerance);
                else
                    throw new InvalidReactionModelEnlargerException();
                finishController = new FinishController(tt, vt);
            } else
                throw new InvalidSymbolException(
                        "condition.txt: can't find FinishController!");
            // read in dynamic simulator
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("DynamicSimulator")) {
                StringTokenizer st = new StringTokenizer(line, ":");
                String temp = st.nextToken();
                String simulator = st.nextToken().trim();
                // read in non-negative option if it exists: syntax would be something like this:
// "DynamicSimulator: DASSL: non-negative"
                if (st.hasMoreTokens()) {
                    if (st.nextToken().trim().toLowerCase()
                            .equals("non-negative")) {
                        if (simulator.toLowerCase().equals("dassl"))
                            JDAS.nonnegative = true;
                        else {
                            Logger.critical("Non-negative option is currently only supported for DASSL. Switch to DASSL solver or remove non-negative option.");
                            System.exit(0);
                        }
                    }
                }
                numConversions = 0;// 5/6/08 gmagoon: moved declaration from initializeReactionSystem() to be an
// attribute so it can be accessed by modelGenerator()
                // int numConversions = 0;
                boolean autoflag = false;// 5/2/08 gmagoon: updating the following if/else-if block to consider input
// where we want to check model validity within the ODE solver at each time step; this will be indicated by the use of a
// string beginning with "AUTO" after the "TimeStep" or "Conversions" line
                // read in time step
                line = ChemParser.readMeaningfulLine(reader, true);
                if (line.startsWith("TimeStep:")) {
                    if (!(finishController.terminationTester instanceof ReactionTimeTT))
                        throw new InvalidSymbolException(
                                "'TimeStep:' specified but finish controller goal is not reaction time.");
                    st = new StringTokenizer(line);
                    temp = st.nextToken();
                    while (st.hasMoreTokens()) {
                        temp = st.nextToken();
                        if (temp.startsWith("AUTO")) {// note potential opportunity for making case insensitive by
// using: temp.toUpperCase().startsWith("AUTO")
                            autoflag = true;
                        } else if (!autoflag) {// use "else if" to make sure additional numbers are not read in case
// numbers are erroneously used following AUTO; note that there could still be a problem if numbers come before "AUTO"
                            double tStep = Double.parseDouble(temp);
                            String unit = "sec";
                            setTimeStep(new ReactionTime(tStep, unit));
                        }
                    }
                    ((ReactionTimeTT) finishController.terminationTester)
                            .setTimeSteps(timeStep);
                } else if (line.startsWith("Conversions:")) {
                    if (!(finishController.terminationTester instanceof ConversionTT))
                        throw new InvalidSymbolException(
                                "'Conversions:' specified but finish controller goal is not conversion.");
                    st = new StringTokenizer(line);
                    temp = st.nextToken();
                    int i = 0;
                    SpeciesConversion sc = (SpeciesConversion) ((ConversionTT) finishController.terminationTester).speciesGoalConversionSet
                            .get(0);
                    Species convSpecies = sc.species;
                    Iterator iter = ((InitialStatus) (initialStatusList.get(0)))
                            .getSpeciesStatus();// 10/23/07 gmagoon: changed to use first element of initialStatusList,
// as subsequent operations should not be affected by which one is chosen
                    double initialConc = 0;
                    while (iter.hasNext()) {
                        SpeciesStatus sps = (SpeciesStatus) iter.next();
                        if (sps.species.equals(convSpecies))
                            initialConc = sps.concentration;
                    }
                    while (st.hasMoreTokens()) {
                        temp = st.nextToken();
                        if (temp.startsWith("AUTO")) {
                            autoflag = true;
                        } else if (!autoflag) {
                            double conv = Double.parseDouble(temp);
                            conversionSet[i] = (1 - conv) * initialConc;
                            i++;
                        }
                    }
                    conversionSet[i] = (1 - conversionSet[49]) * initialConc;
                    numConversions = i + 1;
                } else
                    throw new InvalidSymbolException(
                            "In condition file can't find 'TimeStep:' or 'Conversions:' for dynamic simulator.");
                //
                if (temp.startsWith("AUTOPRUNE")) {// for the AUTOPRUNE case, read in additional lines for termTol and
// edgeTol
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("TerminationTolerance:")) {
                        st = new StringTokenizer(line);
                        temp = st.nextToken();
                        termTol = Double.parseDouble(st.nextToken());
                    } else {
                        Logger.critical("Cannot find TerminationTolerance in condition.txt");
                        System.exit(0);
                    }
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("PruningTolerance:")) {
                        st = new StringTokenizer(line);
                        temp = st.nextToken();
                        edgeTol = Double.parseDouble(st.nextToken());
                    } else {
                        Logger.critical("Cannot find PruningTolerance in condition.txt");
                        System.exit(0);
                    }
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("MinSpeciesForPruning:")) {
                        st = new StringTokenizer(line);
                        temp = st.nextToken();
                        minSpeciesForPruning = Integer.parseInt(st.nextToken());
                    } else {
                        Logger.critical("Cannot find MinSpeciesForPruning in condition.txt");
                        System.exit(0);
                    }
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("MaxEdgeSpeciesAfterPruning:")) {
                        st = new StringTokenizer(line);
                        temp = st.nextToken();
                        maxEdgeSpeciesAfterPruning = Integer.parseInt(st
                                .nextToken());
                    } else {
                        Logger.critical("Cannot find MaxEdgeSpeciesAfterPruning in condition.txt");
                        System.exit(0);
                    }
                    // print header for pruning log (based on restart format)
                    BufferedWriter bw = null;
                    try {
                        File f = new File(System.getProperty("RMG.PruningDir"),
                                "edgeReactions.txt");
                        bw = new BufferedWriter(new FileWriter(f, true));
                        String EaUnits = ArrheniusKinetics.getEaUnits();
                        bw.write("UnitsOfEa: " + EaUnits);
                        bw.newLine();
                    } catch (IOException ex) {
                        Logger.logStackTrace(ex);
                    } finally {
                        try {
                            if (bw != null) {
                                bw.flush();
                                bw.close();
                            }
                        } catch (IOException ex) {
                            Logger.logStackTrace(ex);
                        }
                    }
                } else if (temp.startsWith("AUTO")) {// in the non-autoprune case (i.e. original AUTO functionality), we
// set the new parameters to values that should reproduce original functionality
                    termTol = tolerance;
                    edgeTol = 0;
                    minSpeciesForPruning = 999999;// arbitrary high number (actually, the value here should not matter,
// since pruning should not be done)
                    maxEdgeSpeciesAfterPruning = 999999;
                }
                // read in atol
                line = ChemParser.readMeaningfulLine(reader, true);
                if (line.startsWith("Atol:")) {
                    st = new StringTokenizer(line);
                    temp = st.nextToken();
                    atol = Double.parseDouble(st.nextToken());
                } else
                    throw new InvalidSymbolException(
                            "condition.txt: can't find Atol for dynamic simulator!");
                // read in rtol
                line = ChemParser.readMeaningfulLine(reader, true);
                if (line.startsWith("Rtol:")) {
                    st = new StringTokenizer(line);
                    temp = st.nextToken();
                    String rel_tol = st.nextToken();
                    if (rel_tol.endsWith("%"))
                        rtol = Double.parseDouble(rel_tol.substring(0,
                                rel_tol.length() - 1));
                    else
                        rtol = Double.parseDouble(rel_tol);
                } else
                    throw new InvalidSymbolException(
                            "condition.txt: can't find Rtol for dynamic simulator!");
                if (simulator.equals("DASPK")) {
                    paraInfor = 0;// svp
                    // read in SA
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("Error bars")) {// svp
                        st = new StringTokenizer(line, ":");
                        temp = st.nextToken();
                        String sa = st.nextToken().trim();
                        if (sa.compareToIgnoreCase("on") == 0) {
                            paraInfor = 1;
                            error = true;
                        } else if (sa.compareToIgnoreCase("off") == 0) {
                            paraInfor = 0;
                            error = false;
                        } else
                            throw new InvalidSymbolException(
                                    "condition.txt: can't find error on/off information!");
                    } else
                        throw new InvalidSymbolException(
                                "condition.txt: can't find SA information!");
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("Display sensitivity coefficients")) {// svp
                        st = new StringTokenizer(line, ":");
                        temp = st.nextToken();
                        String sa = st.nextToken().trim();
                        if (sa.compareToIgnoreCase("on") == 0) {
                            paraInfor = 1;
                            sensitivity = true;
                        } else if (sa.compareToIgnoreCase("off") == 0) {
                            if (paraInfor != 1) {
                                paraInfor = 0;
                            }
                            sensitivity = false;
                        } else
                            throw new InvalidSymbolException(
                                    "condition.txt: can't find SA on/off information!");
                        // 10/23/07 gmagoon: changed below from dynamicSimulator to dynamicSimulatorList
                        // 6/25/08 gmagoon: changed loop to use i index, and updated DASPK constructor to pass i
// (mirroring changes to DASSL
                        // 6/25/08 gmagoon: updated to pass autoflag and validity tester; this requires FinishController
// block of input file to be present before DynamicSimulator block, but this requirement may have already existed
// anyway, particularly in construction of conversion/time step lists; *perhaps we should formalize this requirement by
// checking to make sure validityTester is not null?
                        for (int i = 0; i < initialStatusList.size(); i++) {
                            dynamicSimulatorList.add(new JDASPK(rtol, atol, 0,
                                    (InitialStatus) initialStatusList.get(i),
                                    i, finishController.getValidityTester(),
                                    autoflag, termTol, tolerance));
                        }
                    }
                    species = new LinkedList();
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("Display sensitivity information")) {
                        line = ChemParser.readMeaningfulLine(reader, true);
                        Logger.info(line);
                        while (!line.equals("END")) {
                            st = new StringTokenizer(line);
                            String name = st.nextToken();
                            if (name.toUpperCase().equals("ALL"))
                                ReactionSystem.printAllSens = true; // gmagoon 12/22/09: if the line contains the word
// "all", turn on the flag to print out sensitivity information for everything
                            species.add(name);
                            line = ChemParser.readMeaningfulLine(reader, true);
                        }
                    }
                } else if (simulator.equals("DASSL")) {
                    // 10/23/07 gmagoon: changed below from dynamicSimulator to dynamicSimulatorList
                    // for (Iterator iter = initialStatusList.iterator(); iter.hasNext(); ) {
                    // dynamicSimulatorList.add(new JDASSL(rtol, atol, 0, (InitialStatus)iter.next()));
                    // }
                    // 11/1/07 gmagoon: changed loop to use i index, and updated DASSL constructor to pass i
                    // 5/5/08 gmagoon: updated to pass autoflag and validity tester; this requires FinishController
// block of input file to be present before DynamicSimulator block, but this requirement may have already existed
// anyway, particularly in construction of conversion/time step lists; *perhaps we should formalize this requirement by
// checking to make sure validityTester is not null?
                    for (int i = 0; i < initialStatusList.size(); i++) {
                        dynamicSimulatorList.add(new JDASSL(rtol, atol, 0,
                                (InitialStatus) initialStatusList.get(i), i,
                                finishController.getValidityTester(), autoflag,
                                termTol, tolerance));
                    }
                } else if (simulator.equals("Chemkin")) {
                    line = ChemParser.readMeaningfulLine(reader, true);
                    if (line.startsWith("ReactorType")) {
                        st = new StringTokenizer(line, ":");
                        temp = st.nextToken();
                        String reactorType = st.nextToken().trim();
                        // 10/23/07 gmagoon: changed below from dynamicSimulator to dynamicSimulatorList
                        for (Iterator iter = initialStatusList.iterator(); iter
                                .hasNext();) {
                            // dynamicSimulatorList.add(new JDASPK(rtol, atol, 0, (InitialStatus)iter.next()));
                            dynamicSimulatorList.add(new Chemkin(rtol, atol,
                                    reactorType));// 11/4/07 gmagoon: fixing apparent cut/paste error
                        }
                    }
                } else
                    throw new InvalidSymbolException(
                            "condition.txt: Unknown DynamicSimulator = "
                                    + simulator);
                // 10/23/07 gmagoon: changed below from dynamicSimulator to dynamicSimulatorList; note: although
// conversionSet should actually be different for each T,P condition, it will be modified in isTPCconsistent within
// ReactionSystem
                for (Iterator iter = dynamicSimulatorList.iterator(); iter
                        .hasNext();) {
                    double[] cs = conversionSet.clone();// 11/1/07 gmagoon: trying to make sure multiple instances of
// conversionSet are used
                    ((DynamicSimulator) (iter.next())).addConversion(cs,
                            numConversions);
                }
            } else
                throw new InvalidSymbolException(
                        "condition.txt: can't find DynamicSimulator!");
            // read in reaction model enlarger
            /*
             * Read in the Primary Kinetic Library The user can specify as many PKLs, including none, as they like.
             */
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("PrimaryKineticLibrary:")) {
                readAndMakePKL(reader);
            } else
                throw new InvalidSymbolException(
                        "condition.txt: can't find PrimaryKineticLibrary");
            // Reaction Library
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("ReactionLibrary:")) {
                readAndMakeReactionLibrary(reader);
            } else
                throw new InvalidSymbolException(
                        "condition.txt: can't find ReactionLibrary");
            /*
             * Added by MRH 12-Jun-2009 The SeedMechanism acts almost exactly as the old PrimaryKineticLibrary did.
             * Whatever is in the SeedMechanism will be placed in the core at the beginning of the simulation. The user
             * can specify as many seed mechanisms as they like, with the priority (in the case of duplicates) given to
             * the first instance. There is no on/off flag.
             */
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("SeedMechanism:")) {
                line = ChemParser.readMeaningfulLine(reader, true);
                while (!line.equals("END")) {
                    String name = extractLibraryName(line);
                    line = ChemParser.readMeaningfulLine(reader, true);
                    String[] tempString = line.split("Location: ");
                    String location = tempString[tempString.length - 1].trim();
                    line = ChemParser.readMeaningfulLine(reader, true);
                    tempString = line.split("GenerateReactions: ");
                    String generateStr = tempString[tempString.length - 1]
                            .trim();
                    boolean generate = true;
                    if (generateStr.equalsIgnoreCase("yes")
                            || generateStr.equalsIgnoreCase("on")
                            || generateStr.equalsIgnoreCase("true")) {
                        generate = true;
                        Logger.info("Will generate cross-reactions between species in seed mechanism "
                                + name);
                    } else if (generateStr.equalsIgnoreCase("no")
                            || generateStr.equalsIgnoreCase("off")
                            || generateStr.equalsIgnoreCase("false")) {
                        generate = false;
                        Logger.info("Will NOT initially generate cross-reactions between species in seed mechanism "
                                + name);
                        Logger.info("This may have unintended consequences");
                    } else {
                        Logger.critical("Input file invalid");
                        Logger.critical("Please include a 'GenerateReactions: yes/no' line for seed mechanism "
                                + name);
                        System.exit(0);
                    }
                    String path = System
                            .getProperty("jing.rxn.ReactionLibrary.pathName");
                    path += "/" + location;
                    if (getSeedMechanism() == null)
                        setSeedMechanism(new SeedMechanism(name, path,
                                generate, false));
                    else
                        getSeedMechanism().appendSeedMechanism(name, path,
                                generate, false);
                    line = ChemParser.readMeaningfulLine(reader, true);
                }
                if (getSeedMechanism() != null)
                    Logger.info("Seed Mechanisms in use: "
                            + getSeedMechanism().getName());
                else
                    setSeedMechanism(null);
            } else
                throw new InvalidSymbolException(
                        "Error reading condition.txt file: "
                                + "Could not locate SeedMechanism field");
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("ChemkinUnits")) {
                line = ChemParser.readMeaningfulLine(reader, true);
                if (line.startsWith("Verbose:")) {
                    StringTokenizer st = new StringTokenizer(line);
                    String dummyString = st.nextToken();
                    String OnOff = st.nextToken().toLowerCase();
                    if (OnOff.equals("off")) {
                        ArrheniusKinetics.setVerbose(false);
                    } else if (OnOff.equals("on")) {
                        ArrheniusKinetics.setVerbose(true);
                    }
                    line = ChemParser.readMeaningfulLine(reader, true);
                }
                /*
                 * MRH 3MAR2010: Adding user option regarding chemkin file New field: If user would like the empty
                 * SMILES string printed with each species in the thermochemistry portion of the generated chem.inp file
                 */
                if (line.toUpperCase().startsWith("SMILES")) {
                    StringTokenizer st = new StringTokenizer(line);
                    String dummyString = st.nextToken(); // Should be "SMILES:"
                    String OnOff = st.nextToken().toLowerCase();
                    if (OnOff.equals("off")) {
                        Chemkin.setSMILES(false);
                    } else if (OnOff.equals("on")) {
                        Chemkin.setSMILES(true);
                        /*
                         * MRH 9MAR2010: MRH decided not to generate an InChI for every new species during an RMG
                         * simulation (especially since it is not used for anything). Instead, they will only be
                         * generated in the post-processing, if the user asked for InChIs.
                         */
                        // Species.useInChI = true;
                    }
                    line = ChemParser.readMeaningfulLine(reader, true);
                }
                if (line.startsWith("A")) {
                    StringTokenizer st = new StringTokenizer(line);
                    String dummyString = st.nextToken(); // Should be "A:"
                    String units = st.nextToken();
                    if (units.equals("moles") || units.equals("molecules"))
                        ArrheniusKinetics.setAUnits(units);
                    else {
                        Logger.critical("Units for A were not recognized: "
                                + units);
                        System.exit(0);
                    }
                } else
                    throw new InvalidSymbolException(
                            "Error reading condition.txt file: "
                                    + "Could not locate Chemkin units A field.");
                line = ChemParser.readMeaningfulLine(reader, true);
                if (line.startsWith("Ea")) {
                    StringTokenizer st = new StringTokenizer(line);
                    String dummyString = st.nextToken(); // Should be "Ea:"
                    String units = st.nextToken();
                    if (units.equals("kcal/mol") || units.equals("cal/mol")
                            || units.equals("kJ/mol") || units.equals("J/mol")
                            || units.equals("Kelvins"))
                        ArrheniusKinetics.setEaUnits(units);
                    else {
                        Logger.critical("Units for Ea were not recognized: "
                                + units);
                        System.exit(0);
                    }
                } else
                    throw new InvalidSymbolException(
                            "Error reading condition.txt file: "
                                    + "Could not locate Chemkin units Ea field.");
            } else
                throw new InvalidSymbolException(
                        "Error reading condition.txt file: "
                                + "Could not locate ChemkinUnits field.");
            in.close();
            // 11/6/07 gmagoon: initializing temperatureArray and pressureArray before libraryReactionGenerator is
// initialized (initialization calls PDepNetwork and performs initializekLeak); UPDATE: moved after initialStatusList
// initialization (in case primaryKineticLibrary calls the similar pdep functions
            // LinkedList temperatureArray = new LinkedList();
            // LinkedList pressureArray = new LinkedList();
            // Iterator iterIS = initialStatusList.iterator();
            // for (Iterator iter = tempList.iterator(); iter.hasNext(); ) {
            // TemperatureModel tm = (TemperatureModel)iter.next();
            // for (Iterator iter2 = presList.iterator(); iter2.hasNext(); ){
            // PressureModel pm = (PressureModel)iter2.next();
            // InitialStatus is = (InitialStatus)iterIS.next();
            // temperatureArray.add(tm.getTemperature(is.getTime()));
            // pressureArray.add(pm.getPressure(is.getTime()));
            // }
            // }
            // PDepNetwork.setTemperatureArray(temperatureArray);
            // PDepNetwork.setPressureArray(pressureArray);
            // 10/4/07 gmagoon: moved to modelGeneration()
            // ReactionGenerator p_reactionGenerator = new TemplateReactionGenerator();//10/4/07 gmagoon: changed to
// p_reactionGenerator from reactionGenerator
            // setReactionGenerator(p_reactionGenerator);//10/4/07 gmagoon: added
            /*
             * MRH 12-Jun-2009 A TemplateReactionGenerator now requires a Temperature be passed to it. This allows RMG
             * to determine the "best" kinetic parameters to use in the mechanism generation. For now, I choose to pass
             * the first temperature in the list of temperatures. RMG only outputs one mechanism, even for multiple
             * temperature/pressure systems, so we can only have one set of kinetics.
             */
            Temperature t = new Temperature(300, "K");
            for (Iterator iter = tempList.iterator(); iter.hasNext();) {
                TemperatureModel tm = (TemperatureModel) iter.next();
                t = tm.getTemperature(new ReactionTime(0, "sec"));
                setTemp4BestKinetics(t);
                break;
            }
            setReactionGenerator(new TemplateReactionGenerator()); // 11/4/07 gmagoon: moved from modelGeneration;
// mysteriously, moving this later moves "Father" lines up in output at runtime, immediately after condition file (as in
// original code); previously, these Father lines were just before "Can't read primary kinetic library files!"
            lrg = new LibraryReactionGenerator(ReactionLibrary);// 10/10/07 gmagoon: moved from modelGeneration
// (sequence lrg increases species id, and the different sequence was causing problems as main species id was 6 instead
// of 1); //10/31/07 gmagoon: restored this line from 10/10/07 backup: somehow it got lost along the way; 11/5/07
// gmagoon: changed to use "lrg =" instead of setLibraryReactionGenerator
            // 10/24/07 gmagoon: updated to use multiple reactionSystem variables
            reactionSystemList = new LinkedList();
            // LinkedList temperatureArray = new LinkedList();//10/30/07 gmagoon: added temperatureArray variable for
// passing to PDepNetwork; 11/6/07 gmagoon: moved before initialization of lrg;
            // LinkedList pressureArray = new LinkedList();//10/30/07 gmagoon: same for pressure;//UPDATE: commenting
// out: not needed if updateKLeak is done for one temperature/pressure at a time; 11/1-2/07 restored;11/6/07 gmagoon:
// moved before initialization of lrg;
            Iterator iter3 = initialStatusList.iterator();
            Iterator iter4 = dynamicSimulatorList.iterator();
            int i = 0;// 10/30/07 gmagoon: added
            for (Iterator iter = tempList.iterator(); iter.hasNext();) {
                TemperatureModel tm = (TemperatureModel) iter.next();
                // InitialStatus is = (InitialStatus)iter3.next();//10/31/07 gmagoon: fixing apparent bug by moving
// these inside inner "for loop"
                // DynamicSimulator ds = (DynamicSimulator)iter4.next();
                for (Iterator iter2 = presList.iterator(); iter2.hasNext();) {
                    PressureModel pm = (PressureModel) iter2.next();
                    for (int numConcList = 0; numConcList < initialStatusList
                            .size() / tempList.size() / presList.size(); ++numConcList) {
// InitialStatus is = (InitialStatus)iter3.next();//10/31/07 gmagoon: moved from outer "for loop""
                        InitialStatus is = (InitialStatus) initialStatusList
                                .get(i);
                        DynamicSimulator ds = (DynamicSimulator) iter4.next();
                        // temperatureArray.add(tm.getTemperature(is.getTime()));//10/30/07 gmagoon: added; //10/31/07
// added .getTemperature(is.getTime()); 11/6/07 gmagoon: moved before initialization of lrg;
                        // pressureArray.add(pm.getPressure(is.getTime()));//10/30/07 gmagoon: added//UPDATE: commenting
// out: not needed if updateKLeak is done for one temperature/pressure at a time;11/1-2/07 restored with
// .getTemperature(is.getTime()) added;11/6/07 gmagoon: moved before initialization of lrg;
                        // 11/1/07 gmagoon: trying to make a deep copy of terminationTester when it is instance of
// ConversionTT
                        // UPDATE: actually, I don't think this deep copy was necessary; original case with
// FinishController fc = new FinishController(finishController.getTerminationTester(),
// finishController.getValidityTester()) is probably OK; (in any case, this didn't do completetly deep copy (references
// to speciesConversion element in LinkedList were the same);
                        // TerminationTester termTestCopy;
                        // if (finishController.getTerminationTester() instanceof ConversionTT){
                        // ConversionTT termTest = (ConversionTT)finishController.getTerminationTester();
                        // LinkedList spcCopy = (LinkedList)(termTest.getSpeciesGoalConversionSetList().clone());
                        // termTestCopy = new ConversionTT(spcCopy);
                        // }
                        // else{
                        // termTestCopy = finishController.getTerminationTester();
                        // }
                        FinishController fc = new FinishController(
                                finishController.getTerminationTester(),
                                finishController.getValidityTester());// 10/31/07 gmagoon: changed to create new
// finishController instance in each case (apparently, the finish controller becomes associated with reactionSystem in
// setFinishController within ReactionSystem); alteratively, could use clone, but might need to change FinishController
// to be "cloneable"
                        // FinishController fc = new FinishController(termTestCopy,
// finishController.getValidityTester());
                        i++;// 10/30/07 gmagoon: added
                        Logger.info("Creating reaction system " + i);
                        reactionSystemList.add(new ReactionSystem(tm, pm,
                                reactionModelEnlarger, fc, ds,
                                getPrimaryKineticLibrary(),
                                getReactionGenerator(), getSpeciesSeed(), is,
                                getReactionModel(), lrg, i, equationOfState));
                        Logger.info((initialStatusList.get(i - 1)).toString()
                                + "\n");
                    }
                }
            }
            // PDepNetwork.setTemperatureArray(temperatureArray);//10/30/07 gmagoon: passing temperatureArray to
// PDepNetwork; 11/6/07 gmagoon: moved before initialization of lrg;
            // PDepNetwork.setPressureArray(pressureArray);//10/30/07 gmagoon: same for pressure;//UPDATE: commenting
// out: not needed if updateKLeak is done for one temperature/pressure at a time; 11/1-2/07 restored; 11/6/07 gmagoon:
// moved before initialization of lrg;
        } catch (IOException e) {
            Logger.error("Error reading reaction system initialization file.");
            throw new IOException("Input file error: " + e.getMessage());
        }
        Logger.info("");
    }

    public void setReactionModel(ReactionModel p_ReactionModel) {
        reactionModel = p_ReactionModel;
    }

    public void modelGeneration() {
        // long begin_t = System.currentTimeMillis();
        try {
            ChemGraph.readForbiddenStructure();
            setSpeciesSeed(new LinkedHashSet());// 10/4/07 gmagoon moved from initializeCoreEdgeReactionModel
            // setReactionGenerator(new TemplateReactionGenerator());//10/4/07 gmagoon: moved inside
// initializeReactionSystem; 11/3-4/07 gmagoon: probably reverted on or before 10/10/07 (although I have not
// investigated this change in detail); //11/4/07 gmagoon: moved inside initializeReactionSystems
            // setLibraryReactionGenerator(new LibraryReactionGenerator());//10/10/07 gmagoon: moved after
// initializeReactionSystem
            // initializeCoreEdgeReactionModel();//10/4/07 gmagoon moved from below to run
// initializeCoreEdgeReactionModel before initializeReactionSystem; 11/3-4/07 gmagoon: probably reverted on or before
// 10/10/07
            initializeReactionSystems();
        } catch (IOException e) {
            Logger.logStackTrace(e);
            Logger.critical(e.getMessage());
            System.exit(0);
        } catch (InvalidSymbolException e) {
            Logger.logStackTrace(e);
            Logger.critical(e.getMessage());
            System.exit(0);
        }
        // 10/31/07 gmagoon: initialize validList (to false) before initializeCoreEdgeReactionModel is called
        validList = new LinkedList();
        for (Integer i = 0; i < reactionSystemList.size(); i++) {
            validList.add(false);
        }
        initializeCoreEdgeReactionModel();// 10/4/07 gmagoon: moved before initializeReactionSystem; 11/3-4/07 gmagoon:
// probably reverted on or before 10/10/07
        // 10/24/07 gmagoon: changed to use reactionSystemList
        // LinkedList initList = new LinkedList();//10/25/07 gmagoon: moved these variables to apply to entire class
        // LinkedList beginList = new LinkedList();
        // LinkedList endList = new LinkedList();
        // LinkedList lastTList = new LinkedList();
        // LinkedList currentTList = new LinkedList();
        // LinkedList lastPList = new LinkedList();
        // LinkedList currentPList = new LinkedList();
        // LinkedList conditionChangedList = new LinkedList();
        // LinkedList reactionChangedList = new LinkedList();
        // 5/6/08 gmagoon: determine whether there are intermediate time/conversion steps, type of termination tester is
// based on characteristics of 1st reaction system (it is assumed that they are all identical in terms of type of
// termination tester)
        boolean intermediateSteps = true;
        ReactionSystem rs0 = (ReactionSystem) reactionSystemList.get(0);
        if (rs0.finishController.terminationTester instanceof ReactionTimeTT) {
            if (timeStep == null) {
                intermediateSteps = false;
            }
        } else if (numConversions == 1) { // if we get to this block, we presumably have a conversion terminationTester;
// this required moving numConversions to be attribute...alternative to using numConversions is to access one of the
// DynamicSimulators and determine conversion length
            intermediateSteps = false;
        }
        // 10/24/07 gmagoon: note: each element of for loop could be done in parallel if desired; some modifications
// would be needed
        for (Iterator iter = reactionSystemList.iterator(); iter.hasNext();) {
            ReactionSystem rs = (ReactionSystem) iter.next();
            if ((reactionModelEnlarger instanceof RateBasedPDepRME)) {// 1/2/09 gmagoon and rwest: only call
// initializePDepNetwork for P-dep cases
                if (!readrestart)
                    rs.initializePDepNetwork();
            }
            ReactionTime init = rs.getInitialReactionTime();
            initList.add(init);
            ReactionTime begin = init;
            beginList.add(begin);
            ReactionTime end;
            if (rs.finishController.terminationTester instanceof ReactionTimeTT) {
                // 5/5/08 gmagoon: added below if statement to avoid null pointer exception in cases where there are no
// intermediate time steps specified
                if (!(timeStep == null)) {
                    end = (ReactionTime) timeStep.get(0);
                } else {
                    end = ((ReactionTimeTT) rs.finishController.terminationTester).finalTime;
                }
                // end = (ReactionTime)timeStep.get(0);
                endList.add(end);
            } else {
                end = new ReactionTime(1e6, "sec");
                endList.add(end);
            }
            // int iterationNumber = 1;
            lastTList.add(rs.getTemperature(init));
            currentTList.add(rs.getTemperature(init));
            lastPList.add(rs.getPressure(init));
            currentPList.add(rs.getPressure(init));
            conditionChangedList.add(false);
            reactionChangedList.add(false);// 10/31/07 gmagoon: added
            // Chemkin.writeChemkinInputFile(reactionSystem.getReactionModel(),reactionSystem.getPresentStatus());
        }
        int iterationNumber = 1;
        LinkedList terminatedList = new LinkedList();// 10/24/07 gmagoon: this may not be necessary, as if one
// reactionSystem is terminated, I think all should be terminated
        // validList = new LinkedList();//10/31/07 gmagoon: moved before initializeCoreEdgeReactionModel
        // 10/24/07 gmagoon: initialize allTerminated and allValid to true; these variables keep track of whether all
// the reactionSystem variables satisfy termination and validity, respectively
        boolean allTerminated = true;
        boolean allValid = true;
        // IF RESTART IS TURNED ON
        // Update the systemSnapshot for each ReactionSystem in the reactionSystemList
        if (readrestart) {
            for (Integer i = 0; i < reactionSystemList.size(); i++) {
                ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
                InitialStatus is = rs.getInitialStatus();
                putRestartSpeciesInInitialStatus(is, i);
                rs.appendUnreactedSpeciesStatus(
                        (InitialStatus) initialStatusList.get(i),
                        rs.getPresentTemperature());
            }
        }
        printModelSize();
        // Write Chemkin input file only for the LAST reaction system (preserving old behaviour from when it used to be
// overwritten N times).
        Chemkin.writeChemkinInputFile((ReactionSystem) reactionSystemList
                .getLast());
        Logger.info(String.format(
                "Running time: %.3f min",
                +(System.currentTimeMillis() - Global.tAtInitialization) / 1000. / 60.));
        printMemoryUsed();
        Logger.flush();
        // 10/24/07 gmagoon: note: each element of for loop could be done in parallel if desired; some modifications
// would be needed
        for (Integer i = 0; i < reactionSystemList.size(); i++) {
            ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
            ReactionTime begin = (ReactionTime) beginList.get(i);
            ReactionTime end = (ReactionTime) endList.get(i);
            endList.set(i, rs.solveReactionSystem(begin, end, true, true, true,
                    iterationNumber - 1));
            boolean terminated = rs.isReactionTerminated();
            terminatedList.add(terminated);
            if (!terminated)
                allTerminated = false;
            boolean valid = rs.isModelValid();
            // validList.add(valid);
            validList.set(i, valid);// 10/31/07 gmagoon: validList initialization moved before
// initializeCoreEdgeReactionModel
            if (!valid)
                allValid = false;
            reactionChangedList.set(i, false);
        }
        // 9/1/09 gmagoon: if we are using QM, output a file with the CHEMKIN name, the RMG name, the (modified) InChI,
// and the (modified) InChIKey
        if (!ChemGraph.TDMETHOD.contains("benson")) {
            writeInChIs(getReactionModel());
        }
        writeDictionary(getReactionModel());
        // System.exit(0);
        StringBuilder print_info = Global.diagnosticInfo;
        print_info
                .append("\nMolecule \t Flux\t\tTime\t \t\t \t Core \t \t Edge \t \t memory\n");
        print_info
                .append(" \t moleular \t characteristic \t findspecies \t moveUnreactedToReacted \t enlarger \t restart1 \t totalEnlarger \t resetSystem  \t readSolverFile\t writeSolverFile \t justSolver \t SolverIterations \t solverSpeciesStatus \t Totalsolver \t gc  \t restart+diagnosis \t chemkin thermo \t chemkin reactions \t validitytester \t Species \t Reactions\t Species\t Reactions \t memory used  \t allSpecies \t TotalTime \t findRateConstant\t identifyReactedSites \t reactChemGraph \t makespecies\t CheckReverseReaction \t makeTemplateReaction \t getReactionfromStruc \t genReverseFromReac");
        print_info.append("\t\t\t\t\t\t\t"
                + ((CoreEdgeReactionModel) getReactionModel())
                        .getReactedSpeciesSet().size()
                + "\t"
                + ((CoreEdgeReactionModel) getReactionModel())
                        .getReactedReactionSet().size()
                + "\t"
                + ((CoreEdgeReactionModel) getReactionModel())
                        .getUnreactedSpeciesSet().size()
                + "\t"
                + ((CoreEdgeReactionModel) getReactionModel())
                        .getUnreactedReactionSetIncludingReverseSize() + "\t"
                + Global.makeSpecies + "\n");
        double solverMin = 0;
        double vTester = 0;
        /*
         * if (!restart){ writeRestartFile(); writeCoreReactions(); writeAllReactions(); }
         */
        // System.exit(0);
        SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
        Logger.debug("Species dictionary size: " + dictionary.size());
        // boolean reactionChanged = false;//10/24/07 gmagoon: I don't know if this is even required, but I will change
// to use reactionChangedList (I put analogous line of code for list in above for loop); update: yes, it is required; I
// had been thinking of conditionChangedList
        // 10/24/07: changed to use allTerminated and allValid
        // step 2: iteratively grow reaction system
        while (!allTerminated || !allValid) {
            while (!allValid) {
                // writeCoreSpecies();
                double pt = System.currentTimeMillis();
                // Grab all species from primary kinetics / reaction libraries
                // WE CANNOT PRUNE THESE SPECIES
                LinkedHashMap unprunableSpecies = new LinkedHashMap();
                if (getPrimaryKineticLibrary() != null) {
                    unprunableSpecies
                            .putAll(getPrimaryKineticLibrary().speciesSet);
                }
                if (getReactionLibrary() != null) {
                    unprunableSpecies.putAll(getReactionLibrary()
                            .getDictionary());
                }
                // prune the reaction model (this will only do something in the AUTO case)
                pruneReactionModel(unprunableSpecies);
                garbageCollect();
                // System.out.println("After pruning:");
                // printModelSize();
                // ENLARGE THE MODEL!!! (this is where the good stuff happens)
                enlargeReactionModel();
                double totalEnlarger = (System.currentTimeMillis() - pt) / 1000 / 60;
                // PDepNetwork.completeNetwork(reactionSystem.reactionModel.getSpeciesSet());
                // 10/24/07 gmagoon: changed to use reactionSystemList
                if ((reactionModelEnlarger instanceof RateBasedPDepRME)) {// 1/2/09 gmagoon and rwest: only call
// initializePDepNetwork for P-dep cases
                    for (Iterator iter = reactionSystemList.iterator(); iter
                            .hasNext();) {
                        ReactionSystem rs = (ReactionSystem) iter.next();
                        rs.initializePDepNetwork();
                    }
                    // reactionSystem.initializePDepNetwork();
                }
                printModelSize();
                // Write Chemkin input file only for the LAST reaction system (preserving old behaviour from when it
// used to be overwritten N times).
                Chemkin.writeChemkinInputFile((ReactionSystem) reactionSystemList
                        .getLast());
                Logger.info(String.format(
                        "Running time: %.3f min",
                        +(System.currentTimeMillis() - Global.tAtInitialization) / 1000. / 60.));
                printMemoryUsed();
                Logger.flush();
                pt = System.currentTimeMillis();
                // 10/24/07 gmagoon: changed to use reactionSystemList
                for (Iterator iter = reactionSystemList.iterator(); iter
                        .hasNext();) {
                    ReactionSystem rs = (ReactionSystem) iter.next();
                    rs.resetSystemSnapshot();
                }
                // reactionSystem.resetSystemSnapshot();
                double resetSystem = (System.currentTimeMillis() - pt) / 1000 / 60;
                // 10/24/07 gmagoon: changed to use reactionSystemList
                for (Integer i = 0; i < reactionSystemList.size(); i++) {
                    // reactionChanged = true;
                    ReactionSystem rs = (ReactionSystem) reactionSystemList
                            .get(i);
                    reactionChangedList.set(i, true);
                    // begin = init;
                    beginList.set(i, (ReactionTime) initList.get(i));
                    if (rs.finishController.terminationTester instanceof ReactionTimeTT) {
                        // 5/5/08 gmagoon: added below if statement to avoid null pointer exception in cases where there
// are no intermediate time steps specified
                        if (!(timeStep == null)) {
                            endList.set(i, (ReactionTime) timeStep.get(0));
                        } else {
                            endList.set(
                                    i,
                                    ((ReactionTimeTT) rs.finishController.terminationTester).finalTime);
                        }
                        // endList.set(i, (ReactionTime)timeStep.get(0));
                        // end = (ReactionTime)timeStep.get(0);
                    } else
                        endList.set(i, new ReactionTime(1e6, "sec"));
                    // end = new ReactionTime(1e6,"sec");
                    // iterationNumber = 1;//10/24/07 gmagoon: moved outside of loop
                    currentTList.set(i,
                            rs.getTemperature((ReactionTime) beginList.get(i)));
                    currentPList.set(i,
                            rs.getPressure((ReactionTime) beginList.get(i)));
                    conditionChangedList
                            .set(i,
                                    !(((Temperature) currentTList.get(i))
                                            .equals((Temperature) lastTList
                                                    .get(i)))
                                            || !(((Pressure) currentPList
                                                    .get(i))
                                                    .equals((Pressure) lastPList
                                                            .get(i))));
                    // currentT = reactionSystem.getTemperature(begin);
                    // currentP = reactionSystem.getPressure(begin);
                    // conditionChanged = (!currentT.equals(lastT) || !currentP.equals(lastP));
                }
                iterationNumber = 1;
                double startTime = System.currentTimeMillis();
                // 10/24/07 gmagoon: changed to use reactionSystemList
                for (Integer i = 0; i < reactionSystemList.size(); i++) {
                    ReactionSystem rs = (ReactionSystem) reactionSystemList
                            .get(i);
                    boolean reactionChanged = (Boolean) reactionChangedList
                            .get(i);
                    boolean conditionChanged = (Boolean) conditionChangedList
                            .get(i);
                    ReactionTime begin = (ReactionTime) beginList.get(i);
                    ReactionTime end = (ReactionTime) endList.get(i);
                    endList.set(i, rs.solveReactionSystem(begin, end, false,
                            reactionChanged, conditionChanged,
                            iterationNumber - 1));
                    // end = reactionSystem.solveReactionSystem(begin, end, false, reactionChanged, conditionChanged,
// iterationNumber-1);
                }
                solverMin = solverMin
                        + (System.currentTimeMillis() - startTime) / 1000 / 60;
                startTime = System.currentTimeMillis();
                // 9/1/09 gmagoon: if we are using QM, output a file with the CHEMKIN name, the RMG name, the (modified)
// InChI, and the (modified) InChIKey
                if (!ChemGraph.TDMETHOD.contains("benson")) {
                    writeInChIs(getReactionModel());
                }
                writeDictionary(getReactionModel());
                double chemkint = (System.currentTimeMillis() - startTime) / 1000 / 60;
                // this is cheap to do, so do it even if not saving other restart files
                writeRestartConditionFile();
                if (writerestart) {
                    /*
                     * Rename current restart files: In the event RMG fails while writing the restart files, user won't
                     * lose any information
                     */
                    String[] restartFiles = { "coreReactions.txt",
                            "coreSpecies.txt", "edgeReactions.txt",
                            "edgeSpecies.txt", "pdepnetworks.txt",
                            "pdepreactions.txt" };
                    writeBackupRestartFiles(restartFiles);
                    writeCoreSpecies();
                    writeCoreReactions();
                    writeEdgeSpecies();
                    writeEdgeReactions();
                    if (PDepNetwork.generateNetworks == true)
                        writePDepNetworks();
                    /*
                     * Remove backup restart files from Restart folder
                     */
                    removeBackupRestartFiles(restartFiles);
                }
                // 10/24/07 gmagoon: changed to use reactionSystemList
                Logger.info("");
                for (Integer i = 0; i < reactionSystemList.size(); i++) {
                    ReactionSystem rs = (ReactionSystem) reactionSystemList
                            .get(i);
                    Logger.info(String.format(
                            "For reaction system: %d out of %d", i + 1,
                            reactionSystemList.size()));
                    Logger.info(String.format("At this time: %10.4e s",
                            ((ReactionTime) endList.get(i)).getTime()));
                    Species spe = SpeciesDictionary
                            .getSpeciesFromID(getLimitingReactantID());
                    double conv = rs.getPresentConversion(spe);
                    Logger.info(String.format("Conversion of %s is: %-10.4g",
                            spe.getFullName(), conv));
                }
                Logger.info("");
                startTime = System.currentTimeMillis();
                double mU = Runtime.getRuntime().totalMemory()
                        - Runtime.getRuntime().freeMemory();
                Logger.info("");
                double gc = (System.currentTimeMillis() - startTime) / 1000. / 60.;
                startTime = System.currentTimeMillis();
                // 10/24/07 gmagoon: updating to use reactionSystemList
                allValid = true;
                for (Integer i = 0; i < reactionSystemList.size(); i++) {
                    ReactionSystem rs = (ReactionSystem) reactionSystemList
                            .get(i);
                    boolean valid = rs.isModelValid();
                    if (!valid)
                        allValid = false;
                    validList.set(i, valid);
                    // valid = reactionSystem.isModelValid();
                }
                vTester = vTester + (System.currentTimeMillis() - startTime)
                        / 1000 / 60;
                startTime = System.currentTimeMillis();
                writeDiagnosticInfo();
                writeEnlargerInfo();
                double restart2 = (System.currentTimeMillis() - startTime) / 1000 / 60;
                int allSpecies, allReactions;
                allSpecies = SpeciesDictionary.getInstance().size();
                print_info
                        .append(totalEnlarger
                                + "\t"
                                + resetSystem
                                + "\t"
                                + Global.readSolverFile
                                + "\t"
                                + Global.writeSolverFile
                                + "\t"
                                + Global.solvertime
                                + "\t"
                                + Global.solverIterations
                                + "\t"
                                + Global.speciesStatusGenerator
                                + "\t"
                                + solverMin
                                + "\t"
                                + gc
                                + "\t"
                                + restart2
                                + "\t"
                                + Global.chemkinThermo
                                + '\t'
                                + Global.chemkinReaction
                                + "\t"
                                + vTester
                                + "\t"
                                + ((CoreEdgeReactionModel) getReactionModel())
                                        .getReactedSpeciesSet().size()
                                + "\t"
                                + ((CoreEdgeReactionModel) getReactionModel())
                                        .getReactedReactionSet().size()
                                + "\t"
                                + ((CoreEdgeReactionModel) getReactionModel())
                                        .getUnreactedSpeciesSet().size()
                                + "\t"
                                + ((CoreEdgeReactionModel) getReactionModel())
                                        .getUnreactedReactionSetIncludingReverseSize()
                                + "\t"
                                + mU
                                + "\t"
                                + allSpecies
                                + "\t"
                                + (System.currentTimeMillis() - Global.tAtInitialization)
                                / 1000 / 60 + "\t"
                                + String.valueOf(Global.RT_findRateConstant)
                                + "\t" + Global.RT_identifyReactedSites + "\t"
                                + Global.RT_reactChemGraph + "\t"
                                + Global.makeSpecies + "\t"
                                + Global.checkReactionReverse + "\t"
                                + Global.makeTR + "\t"
                                + Global.getReacFromStruc + "\t"
                                + Global.generateReverse + "\n");
            }
            // 5/6/08 gmagoon: in order to handle cases where no intermediate time/conversion steps are used, only
// evaluate the next block of code when there are intermediate time/conversion steps
            double startTime = System.currentTimeMillis();
            if (intermediateSteps) {
                for (Integer i = 0; i < reactionSystemList.size(); i++) {
                    ReactionSystem rs = (ReactionSystem) reactionSystemList
                            .get(i);
                    reactionChangedList.set(i, false);
                    // reactionChanged = false;
                    Temperature currentT = (Temperature) currentTList.get(i);
                    Pressure currentP = (Pressure) currentPList.get(i);
                    lastTList.set(i, (Temperature) currentT.clone());
                    lastPList.set(i, (Pressure) currentP.clone());
                    // lastT = (Temperature)currentT.clone();
                    // lastP = (Pressure)currentP.clone();
                    currentTList.set(i,
                            rs.getTemperature((ReactionTime) beginList.get(i)));// 10/24/07 gmagoon: ****I think this
// should actually be at end? (it shouldn't matter for isothermal/isobaric case)
                    currentPList.set(i,
                            rs.getPressure((ReactionTime) beginList.get(i)));
                    conditionChangedList
                            .set(i,
                                    !(((Temperature) currentTList.get(i))
                                            .equals((Temperature) lastTList
                                                    .get(i)))
                                            || !(((Pressure) currentPList
                                                    .get(i))
                                                    .equals((Pressure) lastPList
                                                            .get(i))));
                    // currentT = reactionSystem.getTemperature(begin);//10/24/07 gmagoon: ****I think this should
// actually be at end? (it shouldn't matter for isothermal/isobaric case)
                    // currentP = reactionSystem.getPressure(begin);
                    // conditionChanged = (!currentT.equals(lastT) || !currentP.equals(lastP));
                    beginList
                            .set(i, ((SystemSnapshot) (rs
                                    .getSystemSnapshotEnd().next())).time);
                    // begin=((SystemSnapshot)(reactionSystem.getSystemSnapshotEnd().next())).time;
                    if (rs.finishController.terminationTester instanceof ReactionTimeTT) {
                        if (iterationNumber < timeStep.size()) {
                            endList.set(i, (ReactionTime) timeStep
                                    .get(iterationNumber));
                            // end = (ReactionTime)timeStep.get(iterationNumber);
                        } else
                            endList.set(
                                    i,
                                    ((ReactionTimeTT) rs.finishController.terminationTester).finalTime);
                        // end = ((ReactionTimeTT)reactionSystem.finishController.terminationTester).finalTime;
                    } else
                        endList.set(i, new ReactionTime(1e6, "sec"));
                    // end = new ReactionTime(1e6,"sec");
                }
                iterationNumber++;
                startTime = System.currentTimeMillis();// 5/6/08 gmagoon: moved declaration outside of if statement so
// it can be accessed in subsequent vTester line; previous steps are probably so fast that I could eliminate this line
// without much effect on normal operation with intermediate steps
                // double startTime = System.currentTimeMillis();
                // 10/24/07 gmagoon: changed to use reactionSystemList
                for (Integer i = 0; i < reactionSystemList.size(); i++) {
                    ReactionSystem rs = (ReactionSystem) reactionSystemList
                            .get(i);
                    boolean reactionChanged = (Boolean) reactionChangedList
                            .get(i);
                    boolean conditionChanged = (Boolean) conditionChangedList
                            .get(i);
                    ReactionTime begin = (ReactionTime) beginList.get(i);
                    ReactionTime end = (ReactionTime) endList.get(i);
                    endList.set(i, rs.solveReactionSystem(begin, end, false,
                            reactionChanged, false, iterationNumber - 1));
                    // end = reactionSystem.solveReactionSystem(begin, end, false, reactionChanged, false,
// iterationNumber-1);
                }
                solverMin = solverMin
                        + (System.currentTimeMillis() - startTime) / 1000 / 60;
                startTime = System.currentTimeMillis();
                // 5/6/08 gmagoon: changed to separate validity and termination testing, and termination testing is done
// last...termination testing should be done even if there are no intermediate conversions; however, validity is
// guaranteed if there are no intermediate conversions based on previous conditional if statement
                allValid = true;
                for (Integer i = 0; i < reactionSystemList.size(); i++) {
                    ReactionSystem rs = (ReactionSystem) reactionSystemList
                            .get(i);
                    boolean valid = rs.isModelValid();
                    validList.set(i, valid);
                    if (!valid)
                        allValid = false;
                }
            }// 5/6/08 gmagoon: end of block for intermediateSteps
            allTerminated = true;
            for (Integer i = 0; i < reactionSystemList.size(); i++) {
                ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
                boolean terminated = rs.isReactionTerminated();
                terminatedList.set(i, terminated);
                if (!terminated) {
                    allTerminated = false;
                    Logger.error("Reaction System " + (i + 1)
                            + " has not reached its termination criterion");
                    if (rs.isModelValid()
                            && runKillableToPreventInfiniteLoop(
                                    intermediateSteps, iterationNumber)) {
                        Logger.info("although it seems to be valid (complete), so it was not interrupted for being invalid.");
                        Logger.info("This probably means there was an error with the ODE solver, and we risk entering an endless loop.");
                        Logger.info("Stopping.");
                        throw new Error();
                    }
                }
            }
            // 10/24/07 gmagoon: changed to use reactionSystemList, allValid
            if (allValid) {
                Logger.info("Model generation completed!");
                printModelSize();
            }
            vTester = vTester + (System.currentTimeMillis() - startTime) / 1000
                    / 60;// 5/6/08 gmagoon: for case where intermediateSteps = false, this will use startTime declared
// just before intermediateSteps loop, and will only include termination testing, but no validity testing
        }
        // always write this, even if not saving other restart files
        writeRestartConditionFile();
        // System.out.println("Performing model reduction");
        if (writerestart) {
            /*
             * Rename current restart files: In the event RMG fails while writing the restart files, user won't lose any
             * information
             */
            String[] restartFiles = { "coreReactions.txt", "coreSpecies.txt",
                    "edgeReactions.txt", "edgeSpecies.txt", "pdepnetworks.txt",
                    "pdepreactions.txt" };
            writeBackupRestartFiles(restartFiles);
            writeCoreSpecies();
            writeCoreReactions();
            writeEdgeSpecies();
            writeEdgeReactions();
            if (PDepNetwork.generateNetworks == true)
                writePDepNetworks();
            /*
             * Remove backup restart files from Restart folder
             */
            removeBackupRestartFiles(restartFiles);
        }
        if (paraInfor != 0) {
            Logger.info("Model Generation performed. Now generating sensitivity data.");
            // 10/24/07 gmagoon: updated to use reactionSystemList
            LinkedList dynamicSimulator2List = new LinkedList();
            for (Integer i = 0; i < reactionSystemList.size(); i++) {
                ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
                // 6/25/08 gmagoon: updated to pass index i
                // 6/25/08 gmagoon: updated to pass (dummy) finishController and autoflag (set to false here);
                dynamicSimulator2List.add(new JDASPK(rtol, atol, paraInfor,
                        (InitialStatus) initialStatusList.get(i), i));
                // DynamicSimulator dynamicSimulator2 = new JDASPK(rtol, atol, paraInfor, initialStatus);
                ((DynamicSimulator) dynamicSimulator2List.get(i))
                        .addConversion(
                                ((JDASPK) rs.dynamicSimulator).conversionSet,
                                ((JDASPK) rs.dynamicSimulator).conversionSet.length);
                // dynamicSimulator2.addConversion(((JDASPK)reactionSystem.dynamicSimulator).conversionSet,
// ((JDASPK)reactionSystem.dynamicSimulator).conversionSet.length);
                rs.setDynamicSimulator((DynamicSimulator) dynamicSimulator2List
                        .get(i));
                // reactionSystem.setDynamicSimulator(dynamicSimulator2);
                int numSteps = rs.systemSnapshot.size() - 1;
                rs.resetSystemSnapshot();
                beginList.set(i, (ReactionTime) initList.get(i));
                // begin = init;
                if (rs.finishController.terminationTester instanceof ReactionTimeTT) {
                    endList.set(
                            i,
                            ((ReactionTimeTT) rs.finishController.terminationTester).finalTime);
                    // end = ((ReactionTimeTT)reactionSystem.finishController.terminationTester).finalTime;
                } else {
                    ReactionTime end = (ReactionTime) endList.get(i);
                    endList.set(i, end.add(end));
                    // end = end.add(end);
                }
                terminatedList.set(i, false);
                // terminated = false;
                ReactionTime begin = (ReactionTime) beginList.get(i);
                ReactionTime end = (ReactionTime) endList.get(i);
                rs.solveReactionSystemwithSEN(begin, end, true, false, false);
                // reactionSystem.solveReactionSystemwithSEN(begin, end, true, false, false);
            }
        }
        // All of the reaction systems are the same, so just write the chemkin
        // file for the first reaction system
        ReactionSystem rs = (ReactionSystem) reactionSystemList.get(0);
        Chemkin.writeChemkinInputFile(getReactionModel(), rs.getPresentStatus());
        // 9/1/09 gmagoon: if we are using QM, output a file with the CHEMKIN name, the RMG name, the (modified) InChI,
// and the (modified) InChIKey
        if (!ChemGraph.TDMETHOD.contains("benson")) {
            writeInChIs(getReactionModel());
        }
        writeDictionary(getReactionModel());
    }

    /*
     * Write a condition file suitable for restarting the job. Basically, the same condition file that you started with,
     * but with all the core species added.
     */
    private void writeRestartConditionFile() {
        StringBuilder restartConditionFile = new StringBuilder(
                "// This is based on the original input condition file, but with core species added.");
        // Read in the condition file
        try {
            String initialConditionFile = System
                    .getProperty("jing.rxnSys.ReactionModelGenerator.conditionFile");
            FileReader condition_file = new FileReader(initialConditionFile);
            BufferedReader condition_file_br = new BufferedReader(
                    condition_file);
            String line = ChemParser.readMeaningfulAndEmptyLine(
                    condition_file_br, true);
            // Read lines from condition file and store them in the string till
            // you reach the species definition
            while (!line.equals("InitialStatus:")) {
                restartConditionFile.append(line + "\n");
                line = ChemParser.readMeaningfulAndEmptyLine(condition_file_br,
                        true);
            }
            restartConditionFile.append(line + "\n");
            // Read in the next line
            line = ChemParser.readMeaningfulLine(condition_file_br, true);
            // Read in the species in the initial condition file
            // and figure out the number of concentrations specified and units
            String unit = "";
            int numSpeciesInitialConditionFile = 0;
            int numConcentrations = 0;
            // This blocks reads in the initial species in condition file
            while (line != null && !line.equals("END")) {
                StringTokenizer st = new StringTokenizer(line);
                String index = st.nextToken();
                String name = null;
                if (!index.startsWith("("))
                    name = index;
                else
                    name = st.nextToken();
                // The next token will be the concentration units
                unit = st.nextToken();
                // The remaining tokens are the either token unreactive, constant concentrations
                // or concentrations for the species
                // The number of concentrations
                numConcentrations = 0;
                while (st.hasMoreTokens()) {
                    String reactive = st.nextToken().trim();
                    // Dont count the unreactive and constantconcentration tokens
                    if (!reactive.equalsIgnoreCase("unreactive")
                            && !reactive
                                    .equalsIgnoreCase("constantconcentration"))
                        ++numConcentrations;
                }
                restartConditionFile.append(line + "\n");
                Graph g;
                try {
                    g = ChemParser.readChemGraph(condition_file_br);
                    ChemGraph cg = null;
                    try {
                        cg = ChemGraph.make(g);
                        restartConditionFile.append(cg.toStringWithoutH(0)
                                + "\n");
                    } catch (ForbiddenStructureException e) {
                        Logger.error("Forbidden Structure in restart condition file! Will not be written.");
                    }
                } catch (IOException e1) {
                    Logger.error("Something wrong with writing the graph to the restartConditionFile");
                    System.exit(0);
                }
                // increment species counter
                ++numSpeciesInitialConditionFile;
                line = ChemParser.readMeaningfulLine(condition_file_br, true);
            }
            // Creates an array initialized with default value of double 0.0, with size of numConcentrations
            double[] ConcArray = new double[numConcentrations];
            String ConcString = ArraytoString(ConcArray);
            // Then print the core and append the units and 0.0 for number of concentrations
            // *** Strictly assuming that species in condition file are the initial species
            // in the RMG_Dictionary ***
            int speciesCounter = 1;
            Iterator iter = getReactionModel().getSpecies();
            while (iter.hasNext()) {
                Species spe = (Species) iter.next();
                if (speciesCounter > numSpeciesInitialConditionFile) {
                    // Append the name + units + concentrations as 0.0
                    restartConditionFile.append(spe.getChemkinName() + " "
                            + unit + " " + ConcString + "\n");
                    // Append the Chemgraph
                    restartConditionFile.append(spe.getChemGraph().toString(0)
                            + "\n");
                }
                ++speciesCounter;
            }
            // Read lines till end of condition file
            while (line != null) {
                restartConditionFile.append(line + "\n");
                line = ChemParser.readMeaningfulAndEmptyLine(condition_file_br,
                        true);
            }
            // Write the restartConditionFile.txt output file
            try {
                File rxns = new File(System.getProperty("RMG.jobOutputDir"),
                        "RestartConditionFile.txt");
                FileWriter fw_rxns = new FileWriter(rxns);
                fw_rxns.write(restartConditionFile.toString());
                fw_rxns.close();
            } catch (IOException e) {
                System.out
                        .println("Could not write RestartConditionFile.txt file");
                System.exit(0);
            }
        } catch (FileNotFoundException e) {
            Logger.error("Condition file not found");
        }
    }

    private String ArraytoString(double[] concArray) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < concArray.length; i++) {
            if (i != 0)
                sb.append(" ");
            sb.append(concArray[i]);
        }
        return sb.toString();
    }

    // 9/1/09 gmagoon: this function writes a "dictionary" with Chemkin name, RMG name, (modified) InChI, and InChIKey
    // this is based off of writeChemkinFile in ChemkinInputFile.java
    private void writeInChIs(ReactionModel p_reactionModel) {
        StringBuilder result = new StringBuilder();
        for (Iterator iter = ((CoreEdgeReactionModel) p_reactionModel).core
                .getSpecies(); iter.hasNext();) {
            Species species = (Species) iter.next();
            result.append(species.getChemkinName() + "\t" + species.getName()
                    + "\t" + species.getChemGraph().getModifiedInChIAnew()
                    + "\t" + species.getChemGraph().getModifiedInChIKeyAnew()
                    + "\n");
        }
        String file = new File(System.getProperty("RMG.jobOutputDir"),
                "inchiDictionary.txt").getPath();
        try {
            FileWriter fw = new FileWriter(file);
            fw.write(result.toString());
            fw.close();
        } catch (Exception e) {
            Logger.logStackTrace(e);
            Logger.critical("Error in writing InChI file inchiDictionary.txt!");
            Logger.critical(e.getMessage());
            System.exit(0);
        }
    }

    // 9/14/09 gmagoon: function to write dictionary, based on code copied from RMG.java
    private void writeDictionary(ReactionModel rm) {
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) rm;
        // Write core species to RMG_Dictionary.txt
        StringBuilder coreSpecies = new StringBuilder(1024);
        Iterator iter = cerm.getSpecies();
        while (iter.hasNext()) {
            int i = 1;
            Species spe = (Species) iter.next();
            coreSpecies.append(spe.getChemkinName());
            if (Species.useInChI)
                coreSpecies.append(" "
                        + spe.getChemGraph().getModifiedInChIAnew());
            coreSpecies.append("\n");
            coreSpecies.append(spe.getChemGraph().toString(i) + "\n\n");
        }
        try {
            File rmgDictionary = new File(
                    System.getProperty("RMG.jobOutputDir"),
                    "RMG_Dictionary.txt");
            FileWriter fw = new FileWriter(rmgDictionary);
            fw.write(coreSpecies.toString());
            fw.close();
        } catch (IOException e) {
            Logger.critical("Could not write RMG_Dictionary.txt");
            System.exit(0);
        }
        // If we have solvation on, then every time we write the dictionary, also write the solvation properties
        if (Species.useSolvation) {
            writeSolvationProperties(rm);
        }
    }

    private void writeSolvationProperties(ReactionModel rm) {
        // Write core species to RMG_Solvation_Properties.txt
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) rm;
        StringBuilder result = new StringBuilder();
        result.append("ChemkinName\tChemicalFormula\tMolecularWeight\tRadius\tDiffusivity\tAbrahamS\tAbrahamB\tAbrahamE\tAbrahamL\tAbrahamA\tAbrahamV\tChemkinName\n\n");
        Iterator iter = cerm.getSpecies();
        while (iter.hasNext()) {
            Species spe = (Species) iter.next();
            result.append(spe.getChemkinName() + "\t");
            result.append(spe.getChemGraph().getChemicalFormula() + "\t");
            result.append(spe.getMolecularWeight() + "\t");
            result.append(spe.getChemGraph().getRadius() + "\t");
            result.append(spe.getChemGraph().getDiffusivity() + "\t");
            result.append(spe.getChemGraph().getAbramData().toString() + "\t");
            result.append(spe.getChemkinName() + "\n");
        }
        try {
            File rmgSolvationProperties = new File(
                    System.getProperty("RMG.jobOutputDir"),
                    "RMG_Solvation_Properties.txt");
            FileWriter fw = new FileWriter(rmgSolvationProperties);
            fw.write(result.toString());
            fw.close();
        } catch (IOException e) {
            Logger.critical("Could not write RMG_Solvation_Properties.txt");
            System.exit(0);
        }
    }

    /*
     * MRH 23MAR2010: Commenting out deprecated parseRestartFiles method
     */
// private void parseRestartFiles() {
// parseAllSpecies();
// parseCoreSpecies();
// parseEdgeSpecies();
// parseAllReactions();
// parseCoreReactions();
//
// }
    /*
     * MRH 23MAR2010: Commenting out deprecated parseEdgeReactions method
     */
// private void parseEdgeReactions() {
// SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
// //HasMap speciesMap = dictionary.dictionary;
// try{
// File coreReactions = new File(System.getProperty("RMG.RestartDir"),"edgeReactions.txt");
// FileReader fr = new FileReader(coreReactions);
// BufferedReader reader = new BufferedReader(fr);
// String line = ChemParser.readMeaningfulLine(reader);
// boolean found = false;
// LinkedHashSet reactionSet = new LinkedHashSet();
// while (line != null){
// Reaction reaction = ChemParser.parseEdgeArrheniusReaction(dictionary,line,1,1);
// boolean added = reactionSet.add(reaction);
// if (!added){
// if (reaction.hasResonanceIsomerAsReactant()){
// //Structure reactionStructure = reaction.getStructure();
// found = getResonanceStructure(reaction,"reactants", reactionSet);
// }
// if (reaction.hasResonanceIsomerAsProduct() && !found){
// //Structure reactionStructure = reaction.getStructure();
// found = getResonanceStructure(reaction,"products", reactionSet);
// }
//
// if (!found){
// System.out.println("Cannot add reaction "+line+" to the Reaction Edge. All resonance isomers have already been added");
// System.exit(0);
// }
// else found = false;
// }
// //Reaction reverse = reaction.getReverseReaction();
// //if (reverse != null) reactionSet.add(reverse);
// line = ChemParser.readMeaningfulLine(reader);
// }
// ((CoreEdgeReactionModel)getReactionModel()).addReactionSet(reactionSet);
// }
// catch (IOException e){
// System.out.println("Could not read the corespecies restart file");
// System.exit(0);
// }
//
// }
    /*
     * MRH 23MAR2010: Commenting out deprecated parseAllSpecies method
     */
// public void parseCoreReactions() {
// SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
// int i=1;
// //HasMap speciesMap = dictionary.dictionary;
// try{
// File coreReactions = new File(System.getProperty("RMG.RestartDir"),"coreReactions.txt");
// FileReader fr = new FileReader(coreReactions);
// BufferedReader reader = new BufferedReader(fr);
// String line = ChemParser.readMeaningfulLine(reader);
// boolean found = false;
// LinkedHashSet reactionSet = new LinkedHashSet();
// while (line != null){
// Reaction reaction =
// ChemParser.parseCoreArrheniusReaction(dictionary,line,1,1);//,((CoreEdgeReactionModel)reactionSystem.reactionModel));
// boolean added = reactionSet.add(reaction);
// if (!added){
// if (reaction.hasResonanceIsomerAsReactant()){
// //Structure reactionStructure = reaction.getStructure();
// found = getResonanceStructure(reaction,"reactants", reactionSet);
// }
// if (reaction.hasResonanceIsomerAsProduct() && !found){
// //Structure reactionStructure = reaction.getStructure();
// found = getResonanceStructure(reaction,"products", reactionSet);
// }
//
// if (!found){
// System.out.println("Cannot add reaction "+line+" to the Reaction Core. All resonance isomers have already been added");
// //System.exit(0);
// }
// else found = false;
// }
//
// Reaction reverse = reaction.getReverseReaction();
// if (reverse != null) {
// reactionSet.add(reverse);
// //System.out.println(2 + "\t " + line);
// }
//
// //else System.out.println(1 + "\t" + line);
// line = ChemParser.readMeaningfulLine(reader);
// i=i+1;
// }
// ((CoreEdgeReactionModel)getReactionModel()).addReactedReactionSet(reactionSet);
// }
// catch (IOException e){
// System.out.println("Could not read the coreReactions restart file");
// System.exit(0);
// }
//
// }
    /*
     * MRH 23MAR2010: Commenting out deprecated parseAllSpecies method
     */
// private void parseAllReactions() {
// SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
// int i=1;
// //HasMap speciesMap = dictionary.dictionary;
// try{
// File allReactions = new File(System.getProperty("RMG.RestartDir"),"allReactions.txt");
// FileReader fr = new FileReader(allReactions);
// BufferedReader reader = new BufferedReader(fr);
// String line = ChemParser.readMeaningfulLine(reader);
// boolean found = false;
// LinkedHashSet reactionSet = new LinkedHashSet();
// OuterLoop:
// while (line != null){
// Reaction reaction =
// ChemParser.parseArrheniusReaction(dictionary,line,1,1,((CoreEdgeReactionModel)getReactionModel()));
//
// if (((CoreEdgeReactionModel)getReactionModel()).categorizeReaction(reaction)==-1){
// boolean added = reactionSet.add(reaction);
// if (!added){
// found = false;
// if (reaction.hasResonanceIsomerAsReactant()){
// //Structure reactionStructure = reaction.getStructure();
// found = getResonanceStructure(reaction,"reactants", reactionSet);
// }
// if (reaction.hasResonanceIsomerAsProduct() && !found){
// //Structure reactionStructure = reaction.getStructure();
// found = getResonanceStructure(reaction,"products", reactionSet);
// }
//
// if (!found){
// Iterator iter = reactionSet.iterator();
// while (iter.hasNext()){
// Reaction reacTemp = (Reaction)iter.next();
// if (reacTemp.equals(reaction)){
// reactionSet.remove(reacTemp);
// reactionSet.add(reaction);
// break;
// }
// }
//
// //System.out.println("Cannot add reaction "+line+" to the Reaction Core. All resonance isomers have already been added");
// //System.exit(0);
// }
//
// //else found = false;
// }
// }
//
//
// /*Reaction reverse = reaction.getReverseReaction();
// if (reverse != null && ((CoreEdgeReactionModel)reactionSystem.reactionModel).isReactedReaction(reaction)) {
// reactionSet.add(reverse);
// //System.out.println(2 + "\t " + line);
// }*/
// //else System.out.println(1 + "\t" + line);
//
// i=i+1;
//
// line = ChemParser.readMeaningfulLine(reader);
// }
// ((CoreEdgeReactionModel)getReactionModel()).addReactionSet(reactionSet);
// }
// catch (IOException e){
// System.out.println("Could not read the corespecies restart file");
// System.exit(0);
// }
//
// }
    private boolean getResonanceStructure(Reaction p_Reaction, String rOrP,
            LinkedHashSet reactionSet) {
        Structure reactionStructure = p_Reaction.getStructure();
        // Structure tempreactionStructure = new
// Structure(reactionStructure.getReactantList(),reactionStructure.getProductList());
        boolean found = false;
        if (rOrP.equals("reactants")) {
            Iterator originalreactants = reactionStructure.getReactants();
            LinkedHashSet tempHashSet = new LinkedHashSet();
            while (originalreactants.hasNext()) {
                tempHashSet.add(originalreactants.next());
            }
            Iterator reactants = tempHashSet.iterator();
            while (reactants.hasNext() && !found) {
                ChemGraph reactant = (ChemGraph) reactants.next();
                if (reactant.getSpecies().hasResonanceIsomers()) {
                    Iterator chemGraphIterator = reactant.getSpecies()
                            .getResonanceIsomers();
                    ChemGraph newChemGraph;// = (ChemGraph)chemGraphIterator.next();
                    while (chemGraphIterator.hasNext() && !found) {
                        newChemGraph = (ChemGraph) chemGraphIterator.next();
                        reactionStructure.removeReactants(reactant);
                        reactionStructure.addReactants(newChemGraph);
                        reactant = newChemGraph;
                        if (reactionSet.add(p_Reaction)) {
                            found = true;
                        }
                    }
                }
            }
        } else {
            Iterator originalproducts = reactionStructure.getProducts();
            LinkedHashSet tempHashSet = new LinkedHashSet();
            while (originalproducts.hasNext()) {
                tempHashSet.add(originalproducts.next());
            }
            Iterator products = tempHashSet.iterator();
            while (products.hasNext() && !found) {
                ChemGraph product = (ChemGraph) products.next();
                if (product.getSpecies().hasResonanceIsomers()) {
                    Iterator chemGraphIterator = product.getSpecies()
                            .getResonanceIsomers();
                    ChemGraph newChemGraph;// = (ChemGraph)chemGraphIterator.next();
                    while (chemGraphIterator.hasNext() && !found) {
                        newChemGraph = (ChemGraph) chemGraphIterator.next();
                        reactionStructure.removeProducts(product);
                        reactionStructure.addProducts(newChemGraph);
                        product = newChemGraph;
                        if (reactionSet.add(p_Reaction)) {
                            found = true;
                        }
                    }
                }
            }
        }
        return found;
    }

    public void parseCoreSpecies() {
        // String restartFileContent ="";
        // int speciesCount = 0;
        // boolean added;
        SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
        try {
            File coreSpecies = new File(System.getProperty("RMG.RestartDir"),
                    "coreSpecies.txt");
            FileReader fr = new FileReader(coreSpecies);
            BufferedReader reader = new BufferedReader(fr);
            String line = ChemParser.readMeaningfulLine(reader, true);
            // LinkedHashSet speciesSet = new LinkedHashSet();
            // if (reactionSystem == null){//10/24/07 gmagoon: commenting out since contents of if was already commented
// out anyway
            // //ReactionSystem reactionSystem = new ReactionSystem();
            // }
            setReactionModel(new CoreEdgeReactionModel());// 10/4/07 gmagoon:changed to setReactionModel
            while (line != null) {
                StringTokenizer st = new StringTokenizer(line);
                String index = st.nextToken();
                int ID = Integer.parseInt(index);
                Species spe = dictionary.getSpeciesFromID(ID);
                if (spe == null)
                    Logger.warning("There was no species with ID " + ID
                            + " in the species dictionary");
                ((CoreEdgeReactionModel) getReactionModel())
                        .addReactedSpecies(spe);
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } catch (IOException e) {
            Logger.critical("Could not read the corespecies restart file");
            System.exit(0);
        }
    }

    public static void garbageCollect() {
        System.gc();
    }

    public static void printMemoryUsed() {
        garbageCollect();
        Runtime rT = Runtime.getRuntime();
        double uM, tM, fM;
        tM = rT.totalMemory() / 1.0e6;
        fM = rT.freeMemory() / 1.0e6;
        uM = tM - fM;
        Logger.debug("After garbage collection:");
        Logger.info(String.format("Memory used: %.2f MB / %.2f MB (%.2f%%)",
                uM, tM, uM / tM * 100.));
    }

    /*
     * MRH 23MAR2010: Commenting out deprecated parseAllSpecies method
     */
// public LinkedHashSet parseAllSpecies() {
// // String restartFileContent ="";
// int speciesCount = 0;
// LinkedHashSet speciesSet = new LinkedHashSet();
//
// boolean added;
// try{
// long initialTime = System.currentTimeMillis();
//
// File coreSpecies = new File ("allSpecies.txt");
// BufferedReader reader = new BufferedReader(new FileReader(coreSpecies));
// String line = ChemParser.readMeaningfulLine(reader);
// int i=0;
// while (line!=null) {
// i++;
// StringTokenizer st = new StringTokenizer(line);
// String index = st.nextToken();
// String name = null;
// if (!index.startsWith("(")) name = index;
// else name = st.nextToken().trim();
// int ID = getID(name);
//
// name = getName(name);
// Graph g = ChemParser.readChemGraph(reader);
// ChemGraph cg = null;
// try {
// cg = ChemGraph.make(g);
// }
// catch (ForbiddenStructureException e) {
// System.out.println("Forbidden Structure:\n" + e.getMessage());
// System.exit(0);
// }
// Species species;
// if (ID == 0)
// species = Species.make(name,cg);
// else
// species = Species.make(name,cg,ID);
// speciesSet.add(species);
// double flux = 0;
// int species_type = 1;
// line = ChemParser.readMeaningfulLine(reader);
// System.out.println(line);
// }
// }
// catch (IOException e){
// System.out.println("Could not read the allSpecies restart file");
// System.exit(0);
// }
// return speciesSet;
// }
    private String getName(String name) {
        // int id;
        String number = "";
        int index = 0;
        if (!name.endsWith(")"))
            return name;
        else {
            char[] nameChars = name.toCharArray();
            String temp = String.copyValueOf(nameChars);
            int i = name.length() - 2;
            // char test = "(";
            while (i > 0) {
                if (name.charAt(i) == '(') {
                    index = i;
                    i = 0;
                } else
                    i = i - 1;
            }
        }
        number = name.substring(0, index);
        return number;
    }

    private int getID(String name) {
        int id;
        String number = "";
        if (!name.endsWith(")"))
            return 0;
        else {
            char[] nameChars = name.toCharArray();
            int i = name.length() - 2;
            // char test = "(";
            while (i > 0) {
                if (name.charAt(i) == '(')
                    i = 0;
                else {
                    number = name.charAt(i) + number;
                    i = i - 1;
                }
            }
        }
        id = Integer.parseInt(number);
        return id;
    }

    /*
     * MRH 23MAR2010: Commenting out deprecated parseAllSpecies method
     */
// private void parseEdgeSpecies() {
// // String restartFileContent ="";
// SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
// try{
// File edgeSpecies = new File (System.getProperty("RMG.RestartDir"),"edgeSpecies.txt");
// FileReader fr = new FileReader(edgeSpecies);
// BufferedReader reader = new BufferedReader(fr);
// String line = ChemParser.readMeaningfulLine(reader);
// //LinkedHashSet speciesSet = new LinkedHashSet();
//
// while (line!=null) {
//
// StringTokenizer st = new StringTokenizer(line);
// String index = st.nextToken();
// int ID = Integer.parseInt(index);
// Species spe = dictionary.getSpeciesFromID(ID);
// if (spe == null)
// System.out.println("There was no species with ID "+ID +" in the species dictionary");
// //reactionSystem.reactionModel = new CoreEdgeReactionModel();
// ((CoreEdgeReactionModel)getReactionModel()).addUnreactedSpecies(spe);
// line = ChemParser.readMeaningfulLine(reader);
// }
//
// }
// catch (IOException e){
// System.out.println("Could not read the edgepecies restart file");
// System.exit(0);
// }
//
// }
    /*
     * private int calculateAllReactionsinReactionTemplate() { int totalnum = 0; TemplateReactionGenerator trg =
     * (TemplateReactionGenerator)reactionSystem.reactionGenerator; Iterator iter =
     * trg.getReactionTemplateLibrary().getReactionTemplate(); while (iter.hasNext()){ ReactionTemplate rt =
     * (ReactionTemplate)iter.next(); totalnum += rt.getNumberOfReactions(); } return totalnum; }
     */
    private void writeEnlargerInfo() {
        try {
            File diagnosis = new File(System.getProperty("RMG.jobOutputDir"),
                    "enlarger.xls");
            FileWriter fw = new FileWriter(diagnosis);
            fw.write(Global.enlargerInfo.toString());
            fw.close();
        } catch (IOException e) {
            Logger.critical("Cannot write enlarger file");
            System.exit(0);
        }
    }

    private void writeDiagnosticInfo() {
        try {
            File diagnosis = new File(System.getProperty("RMG.jobOutputDir"),
                    "diagnosis.xls");
            FileWriter fw = new FileWriter(diagnosis);
            fw.write(Global.diagnosticInfo.toString());
            fw.close();
        } catch (IOException e) {
            Logger.critical("Cannot write diagnosis file");
            System.exit(0);
        }
    }

    // 10/25/07 gmagoon: I don't think this is used, but I will update to use reactionSystem and reactionTime as
// parameter to access temperature; commented-out usage of writeRestartFile will need to be modified
    // Is still incomplete.
    public void writeRestartFile(ReactionSystem p_rs, ReactionTime p_time) {
        // writeCoreSpecies(p_rs);
        // writeCoreReactions(p_rs, p_time);
        // writeEdgeSpecies();
        // writeAllReactions(p_rs, p_time);
        // writeEdgeReactions(p_rs, p_time);
        // String restartFileName;
        // String restartFileContent="";
        /*
         * File restartFile; try { restartFileName="Restart.txt"; restartFile = new File(restartFileName); FileWriter fw
         * = new FileWriter(restartFile); restartFileContent = restartFileContent+ "TemperatureModel: Constant " +
         * reactionSystem.temperatureModel.getTemperature(reactionSystem.getInitialReactionTime()).getStandard()+
         * " (K)\n"; restartFileContent = restartFileContent+ "PressureModel: Constant " +
         * reactionSystem.pressureModel.getPressure(reactionSystem.getInitialReactionTime()).getAtm() + " (atm)\n\n";
         * restartFileContent = restartFileContent+ "InitialStatus: \n"; int speciesCount = 1; for(Iterator
         * iter=reactionSystem.reactionModel.getSpecies();iter.hasNext();){ Species species = (Species) iter.next();
         * restartFileContent = restartFileContent + "("+ speciesCount + ") "+species.getChemkinName() + "  " +
         * reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n"; restartFileContent = restartFileContent +
         * species.toStringWithoutH(1) + "\n"; speciesCount = speciesCount + 1; } restartFileContent =
         * restartFileContent + "\n\n END \n\nInertGas:\n"; for (Iterator
         * iter=reactionSystem.getPresentStatus().getInertGas(); iter.hasNext();){ String inertName= (String)
         * iter.next(); restartFileContent = restartFileContent + inertName + " " +
         * reactionSystem.getPresentStatus().getInertGas(inertName) + " mol/cm3 \n"; } restartFileContent =
         * restartFileContent + "END\n"; double tolerance; if (reactionSystem.reactionModelEnlarger instanceof
         * RateBasedPDepRME){ restartFileContent = restartFileContent +
         * "ReactionModelEnlarger: RateBasedPDepModelEnlarger\n"; restartFileContent = restartFileContent +
         * "FinishController: RateBasedPDepFinishController\n"; } else { restartFileContent = restartFileContent +
         * "ReactionModelEnlarger: RateBasedModelEnlarger\n"; restartFileContent = restartFileContent +
         * "FinishController: RateBasedFinishController\n"; } if (reactionSystem.finishController.terminationTester
         * instanceof ConversionTT){ restartFileContent = restartFileContent + "(1) Goal Conversion: "; for (Iterator
         * iter =
         * ((ConversionTT)reactionSystem.finishController.terminationTester).getSpeciesGoalConversionSet();iter.hasNext
         * ();){ SpeciesConversion sc = (SpeciesConversion) iter.next(); restartFileContent = restartFileContent +
         * sc.getSpecies().getName()+" "+sc.conversion + "\n"; } } else { restartFileContent = restartFileContent +
         * "(1) Goal ReactionTime: "
         * +((ReactionTimeTT)reactionSystem.finishController.terminationTester).finalTime.toString(); }
         * restartFileContent = restartFileContent + "Error Tolerance: " +((RateBasedVT)
         * reactionSystem.finishController.validityTester).tolerance + "\n"; fw.write(restartFileContent); fw.close(); }
         * catch (IOException e) { System.out.println("Could not write the restart file"); System.exit(0); }
         */
    }

    /*
     * MRH 25MAR2010 This method is no longer used
     */
    /*
     * Only write the forward reactions in the model core. The reverse reactions are generated from the forward
     * reactions.
     */
    // 10/25/07 gmagoon: added reaction system and reaction time as parameters and eliminated use of Global.temperature
// private void writeEdgeReactions(ReactionSystem p_rs, ReactionTime p_time) {
// StringBuilder restartFileContent =new StringBuilder();
// int reactionCount = 1;
// try{
// File coreSpecies = new File (System.getProperty("RMG.RestartDir"),"edgeReactions.txt");
// FileWriter fw = new FileWriter(coreSpecies);
// for(Iterator iter=((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSet().iterator();iter.hasNext();){
//
// Reaction reaction = (Reaction) iter.next();
// //if (reaction.getDirection()==1){
// //restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " +
// reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
// restartFileContent = restartFileContent.append(reaction.toRestartString(p_rs.getTemperature(p_time)) + "\n");
// reactionCount = reactionCount + 1;
// //}
// }
// //restartFileContent += "\nEND";
// fw.write(restartFileContent.toString());
// fw.close();
// }
// catch (IOException e){
// System.out.println("Could not write the restart edgereactions file");
// System.exit(0);
// }
//
// }
    /*
     * MRH 25MAR2010: This method is no longer used
     */
    // 10/25/07 gmagoon: added reaction system and reaction time as parameters and eliminated use of Global.temperature
// private void writeAllReactions(ReactionSystem p_rs, ReactionTime p_time) {
// StringBuilder restartFileContent = new StringBuilder();
// int reactionCount = 1;
// try{
// File allReactions = new File (System.getProperty("RMG.RestartDir"),"allReactions.txt");
// FileWriter fw = new FileWriter(allReactions);
// for(Iterator iter=getReactionModel().getReaction();iter.hasNext();){
//
// Reaction reaction = (Reaction) iter.next();
//
// //restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " +
// reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
// restartFileContent = restartFileContent.append(reaction.toRestartString(p_rs.getTemperature(p_time)) + "\n");
//
// }
//
// for(Iterator iter=((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSet().iterator();iter.hasNext();){
//
// Reaction reaction = (Reaction) iter.next();
// //if (reaction.getDirection()==1){
// //restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " +
// reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
// restartFileContent = restartFileContent.append(reaction.toRestartString(p_rs.getTemperature(p_time)) + "\n");
//
// }
// //restartFileContent += "\nEND";
// fw.write(restartFileContent.toString());
// fw.close();
// }
// catch (IOException e){
// System.out.println("Could not write the restart edgereactions file");
// System.exit(0);
// }
//
// }
    private void writeEdgeSpecies() {
        BufferedWriter bw = null;
        Logger.info("Writing Restart Edge Species");
        try {
            bw = new BufferedWriter(new FileWriter(new File(
                    System.getProperty("RMG.RestartDir"), "edgeSpecies.txt")));
            for (Iterator iter = ((CoreEdgeReactionModel) getReactionModel())
                    .getUnreactedSpeciesSet().iterator(); iter.hasNext();) {
                Species species = (Species) iter.next();
                bw.write(species.getFullName());
                bw.newLine();
                int dummyInt = 0;
                bw.write(species.getChemGraph().toString(dummyInt));
                bw.newLine();
            }
        } catch (IOException ex) {
            Logger.logStackTrace(ex);
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                Logger.logStackTrace(ex);
            }
        }
    }

    private void writePrunedEdgeSpecies(Species species) {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new FileWriter(new File(
                    System.getProperty("RMG.PruningDir"), "edgeSpecies.txt"),
                    true));
            bw.write(species.getChemkinName());
            bw.newLine();
            int dummyInt = 0;
            bw.write(species.getChemGraph().toString(dummyInt));
            bw.newLine();
        } catch (IOException ex) {
            Logger.logStackTrace(ex);
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                Logger.logStackTrace(ex);
            }
        }
    }

    /*
     * MRH 25MAR2010: This method is no longer used
     */
    // 10/25/07 gmagoon: added reaction system and reaction time as parameters and eliminated use of Global.temperature
// private void writeCoreReactions(ReactionSystem p_rs, ReactionTime p_time) {
// StringBuilder restartFileContent = new StringBuilder();
// int reactionCount = 0;
// try{
// File coreSpecies = new File (System.getProperty("RMG.RestartDir"),"coreReactions.txt");
// FileWriter fw = new FileWriter(coreSpecies);
// for(Iterator iter=getReactionModel().getReaction();iter.hasNext();){
//
// Reaction reaction = (Reaction) iter.next();
// if (reaction.getDirection()==1){
// //restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " +
// reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
// restartFileContent = restartFileContent.append(reaction.toRestartString(p_rs.getTemperature(p_time)) + "\n");
// reactionCount = reactionCount + 1;
// }
// }
// //restartFileContent += "\nEND";
// fw.write(restartFileContent.toString());
// fw.close();
// }
// catch (IOException e){
// System.out.println("Could not write the restart corereactions file");
// System.exit(0);
// }
// }
    private void writeCoreSpecies() {
        BufferedWriter bw = null;
        Logger.info("Writing Restart Core Species");
        try {
            bw = new BufferedWriter(new FileWriter(new File(
                    System.getProperty("RMG.RestartDir"), "coreSpecies.txt")));
            for (Iterator iter = getReactionModel().getSpecies(); iter
                    .hasNext();) {
                Species species = (Species) iter.next();
                bw.write(species.getFullName());
                bw.newLine();
                int dummyInt = 0;
                bw.write(species.getChemGraph().toString(dummyInt));
                bw.newLine();
            }
        } catch (IOException ex) {
            Logger.logStackTrace(ex);
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                Logger.logStackTrace(ex);
            }
        }
    }

    private void writeCoreReactions() {
        BufferedWriter bw_rxns = null;
        BufferedWriter bw_pdeprxns = null;
        Logger.info("Writing Restart Core Reactions");
        try {
            bw_rxns = new BufferedWriter(new FileWriter(new File(
                    System.getProperty("RMG.RestartDir"), "coreReactions.txt")));
            bw_pdeprxns = new BufferedWriter(new FileWriter(new File(
                    System.getProperty("RMG.RestartDir"), "pdepreactions.txt")));
            String EaUnits = ArrheniusKinetics.getEaUnits();
            String AUnits = ArrheniusKinetics.getAUnits();
            bw_rxns.write("UnitsOfEa: " + EaUnits);
            bw_rxns.newLine();
            bw_pdeprxns.write("Unit:\nA: mol/cm3/s\nE: " + EaUnits
                    + "\n\nReactions:");
            bw_pdeprxns.newLine();
            CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) getReactionModel();
            LinkedHashSet allcoreRxns = cerm.core.reaction;
            for (Iterator iter = allcoreRxns.iterator(); iter.hasNext();) {
                Reaction reaction = (Reaction) iter.next();
                if (reaction.isForward()) {
                    if (reaction instanceof TROEReaction) {
                        TROEReaction troeRxn = (TROEReaction) reaction;
                        bw_pdeprxns.write(troeRxn
                                .toRestartString(new Temperature(298, "K")));
                        bw_pdeprxns.newLine();
                    } else if (reaction instanceof LindemannReaction) {
                        LindemannReaction lindeRxn = (LindemannReaction) reaction;
                        bw_pdeprxns.write(lindeRxn
                                .toRestartString(new Temperature(298, "K")));
                        bw_pdeprxns.newLine();
                    } else if (reaction instanceof ThirdBodyReaction) {
                        ThirdBodyReaction tbRxn = (ThirdBodyReaction) reaction;
                        bw_pdeprxns.write(tbRxn
                                .toRestartString(new Temperature(298, "K")));
                        bw_pdeprxns.newLine();
                    } else {
                        // bw.write(reaction.toChemkinString(new Temperature(298,"K")));
                        bw_rxns.write(reaction.toRestartString(new Temperature(
                                298, "K"), false));
                        bw_rxns.newLine();
                    }
                }
            }
        } catch (IOException ex) {
            Logger.logStackTrace(ex);
        } finally {
            try {
                if (bw_rxns != null) {
                    bw_rxns.flush();
                    bw_rxns.close();
                }
                if (bw_pdeprxns != null) {
                    bw_pdeprxns.flush();
                    bw_pdeprxns.close();
                }
            } catch (IOException ex) {
                Logger.logStackTrace(ex);
            }
        }
    }

    private void writeEdgeReactions() {
        BufferedWriter bw = null;
        Logger.info("Writing Restart Edge Reactions");
        try {
            bw = new BufferedWriter(new FileWriter(new File(
                    System.getProperty("RMG.RestartDir"), "edgeReactions.txt")));
            String EaUnits = ArrheniusKinetics.getEaUnits();
            bw.write("UnitsOfEa: " + EaUnits);
            bw.newLine();
            CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) getReactionModel();
            LinkedHashSet alledgeRxns = cerm.edge.reaction;
            for (Iterator iter = alledgeRxns.iterator(); iter.hasNext();) {
                Reaction reaction = (Reaction) iter.next();
                if (reaction.isForward()) {
                    // bw.write(reaction.toChemkinString(new Temperature(298,"K")));
                    bw.write(reaction.toRestartString(
                            new Temperature(298, "K"), false));
                    bw.newLine();
                } else if (reaction.getReverseReaction().isForward()) {
                    // bw.write(reaction.getReverseReaction().toChemkinString(new Temperature(298,"K")));
                    bw.write(reaction.getReverseReaction().toRestartString(
                            new Temperature(298, "K"), false));
                    bw.newLine();
                } else
                    Logger.warning("Could not determine forward direction for following rxn: "
                            + reaction.toString());
            }
        } catch (IOException ex) {
            Logger.logStackTrace(ex);
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                Logger.logStackTrace(ex);
            }
        }
    }

    // gmagoon 4/5/10: based on Mike's writeEdgeReactions
    private void writePrunedEdgeReaction(Reaction reaction) {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new FileWriter(new File(
                    System.getProperty("RMG.PruningDir"), "edgeReactions.txt"),
                    true));
            if (reaction.isForward()) {
                bw.write(reaction.toChemkinString(new Temperature(298, "K")));
                // bw.write(reaction.toRestartString(new Temperature(298,"K")));
                bw.newLine();
            } else if (reaction.getReverseReaction().isForward()) {
                bw.write(reaction.getReverseReaction().toChemkinString(
                        new Temperature(298, "K")));
                // bw.write(reaction.getReverseReaction().toRestartString(new Temperature(298,"K")));
                bw.newLine();
            } else
                Logger.warning("Could not determine forward direction for following rxn: "
                        + reaction.toString());
        } catch (IOException ex) {
            Logger.logStackTrace(ex);
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                Logger.logStackTrace(ex);
            }
        }
    }

    private void writePDepNetworks() {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new FileWriter(new File(
                    System.getProperty("RMG.RestartDir"), "pdepnetworks.txt")));
            int numFameTemps = PDepRateConstant.getTemperatures().length;
            int numFamePress = PDepRateConstant.getPressures().length;
            int numChebyTemps = ChebyshevPolynomials.getDefaultNT();
            int numChebyPress = ChebyshevPolynomials.getDefaultNP();
            // int numPlog = PDepArrheniusKinetics.getNumPressures();
            int numPlog = numFamePress; // probably often the case for FAME-generated PLOG rates (but not for seed or
// library reactions)
            String EaUnits = ArrheniusKinetics.getEaUnits();
            bw.write("UnitsOfEa: " + EaUnits);
            bw.newLine();
            bw.write("NumberOfFameTemps: " + numFameTemps);
            bw.newLine();
            bw.write("NumberOfFamePress: " + numFamePress);
            bw.newLine();
            bw.write("NumberOfChebyTemps: " + numChebyTemps);
            bw.newLine();
            bw.write("NumberOfChebyPress: " + numChebyPress);
            bw.newLine();
            bw.write("NumberOfPLogs: Could be different for seed and library reactions, but default for FAME-generated rates is probably "
                    + numPlog);
            bw.newLine();
            bw.newLine();
            LinkedList allNets = PDepNetwork.getNetworks();
            for (Iterator iter = allNets.iterator(); iter.hasNext();) {
                PDepNetwork pdepnet = (PDepNetwork) iter.next();
                bw.write("PDepNetwork #" + pdepnet.getID());
                bw.newLine();
                // Write netReactionList
                LinkedList netRxns = pdepnet.getNetReactions();
                bw.write("netReactionList:");
                bw.newLine();
                for (Iterator iter2 = netRxns.iterator(); iter2.hasNext();) {
                    PDepReaction currentPDepRxn = (PDepReaction) iter2.next();
                    bw.write(currentPDepRxn.toString());
                    bw.newLine();
                    bw.write(writeRatesAndParameters(currentPDepRxn,
                            numFameTemps, numFamePress, numChebyTemps,
                            numChebyPress, numPlog));
                }
                // Write nonincludedReactionList
                LinkedList nonIncludeRxns = pdepnet.getNonincludedReactions();
                bw.write("nonIncludedReactionList:");
                bw.newLine();
                for (Iterator iter2 = nonIncludeRxns.iterator(); iter2
                        .hasNext();) {
                    PDepReaction currentPDepRxn = (PDepReaction) iter2.next();
                    bw.write(currentPDepRxn.toString());
                    bw.newLine();
                    bw.write(writeRatesAndParameters(currentPDepRxn,
                            numFameTemps, numFamePress, numChebyTemps,
                            numChebyPress, numPlog));
                    PDepReaction currentPDepReverseRxn = currentPDepRxn
                            .getReverseReaction();
                    // Not all nonIncludedReactions are reversible
                    // (MRH FEB-16-2011) ... and not all reverse reactions
                    // have pressure-dependent rate coefficients (apparently)
                    if (currentPDepReverseRxn != null
                            && currentPDepReverseRxn.getPDepRate() != null) {
                        bw.write(currentPDepReverseRxn.toString());
                        bw.newLine();
                        bw.write(writeRatesAndParameters(currentPDepReverseRxn,
                                numFameTemps, numFamePress, numChebyTemps,
                                numChebyPress, numPlog));
                    }
                }
                // Write pathReactionList
                LinkedList pathRxns = pdepnet.getPathReactions();
                bw.write("pathReactionList:");
                bw.newLine();
                for (Iterator iter2 = pathRxns.iterator(); iter2.hasNext();) {
                    PDepReaction currentPDepRxn = (PDepReaction) iter2.next();
                    bw.write(currentPDepRxn.getDirection()
                            + "\t"
                            + currentPDepRxn.toRestartString(new Temperature(
                                    298, "K")));
                    bw.newLine();
                }
                bw.newLine();
                bw.newLine();
            }
        } catch (IOException ex) {
            Logger.logStackTrace(ex);
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                Logger.logStackTrace(ex);
            }
        }
    }

    public String writeRatesAndParameters(PDepReaction pdeprxn,
            int numFameTemps, int numFamePress, int numChebyTemps,
            int numChebyPress, int numPlog) {
        StringBuilder sb = new StringBuilder();
        // Write the rate coefficients
        double[][] rateConstants = pdeprxn.getPDepRate().getRateConstants();
        for (int i = 0; i < numFameTemps; i++) {
            for (int j = 0; j < numFamePress; j++) {
                sb.append(rateConstants[i][j] + "\t");
            }
            sb.append("\n");
        }
        sb.append("\n");
        // If chebyshev polynomials are present, write them
        if (numChebyTemps != 0) {
            ChebyshevPolynomials chebyPolys = pdeprxn.getPDepRate()
                    .getChebyshev();
            for (int i = 0; i < numChebyTemps; i++) {
                for (int j = 0; j < numChebyPress; j++) {
                    sb.append(chebyPolys.getAlpha(i, j) + "\t");
                }
                sb.append("\n");
            }
            sb.append("\n");
        }
        // If plog parameters are present, write them
        PDepArrheniusKinetics[] kinetics = pdeprxn.getPDepRate()
                .getPDepArrheniusKinetics();
        if (kinetics != null) {
            for (int i = 0; i < kinetics.length; i++) {
                sb.append(kinetics[i].toChemkinString(pdeprxn
                        .getReactantNumber()));
                sb.append("\n");
            }
        }
        /*
         * else if (numPlog != 0) { PDepArrheniusKinetics kinetics = pdeprxn.getPDepRate().getPDepArrheniusKinetics();
         * for (int i=0; i<numPlog; i++) { double Hrxn = pdeprxn.calculateHrxn(new Temperature(298,"K"));
         * sb.append(kinetics.pressures[i].getPa() + "\t" + kinetics.getKinetics(i).toChemkinString(Hrxn,new
         * Temperature(298,"K"),false) + "\n"); } sb.append("\n"); }
         */
        return sb.toString();
    }

    public LinkedList getTimeStep() {
        return timeStep;
    }

    public static boolean getUseDiffusion() {
        return useDiffusion;
    }

    public static boolean getUseSolvation() {
        return useSolvation;
    }

    public void setUseDiffusion(Boolean p_boolean) {
        useDiffusion = p_boolean;
    }

    public void setUseSolvation(Boolean p_boolean) {
        useSolvation = p_boolean;
    }

    public void setTimeStep(ReactionTime p_timeStep) {
        if (timeStep == null)
            timeStep = new LinkedList();
        timeStep.add(p_timeStep);
    }

    public String getWorkingDirectory() {
        return workingDirectory;
    }

    public void setWorkingDirectory(String p_workingDirectory) {
        workingDirectory = p_workingDirectory;
    }

    // svp
    public boolean getError() {
        return error;
    }

    // svp
    public boolean getSensitivity() {
        return sensitivity;
    }

    public LinkedList getSpeciesList() {
        return species;
    }

    // gmagoon 10/24/07: commented out getReactionSystem and setReactionSystem
    // public ReactionSystem getReactionSystem() {
    // return reactionSystem;
    // }
    // 11/2/07 gmagoon: adding accessor method for reactionSystemList
    public LinkedList getReactionSystemList() {
        return reactionSystemList;
    }

    // added by gmagoon 9/24/07
    // public void setReactionSystem(ReactionSystem p_ReactionSystem) {
    // reactionSystem = p_ReactionSystem;
    // }
    // copied from ReactionSystem.java by gmagoon 9/24/07
    public ReactionModel getReactionModel() {
        return reactionModel;
    }

    public void readRestartSpecies() {
        Logger.info("Reading in species from Restart folder");
        // Read in core species -- NOTE code is almost duplicated in Read in edge species (second part of procedure)
        try {
            FileReader in = new FileReader(new File(
                    System.getProperty("RMG.RestartDir"), "coreSpecies.txt"));
            BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader, true);
            while (line != null) {
                // The first line of a new species is the user-defined name
                String totalSpeciesName = line;
                String[] splitString1 = totalSpeciesName.split("[(]");
                String[] splitString2 = splitString1[splitString1.length - 1]
                        .split("[)]");
                // The remaining lines are the graph
                Graph g = ChemParser.readChemGraph(reader);
                // Make the ChemGraph, assuming it does not contain a forbidden structure
                ChemGraph cg = null;
                try {
                    cg = ChemGraph.make(g);
                } catch (ForbiddenStructureException e) {
                    Logger.critical("Error reading graph: Graph contains a forbidden structure.\n"
                            + g.toString());
                    System.exit(0);
                }
                // Make the species
                int intLocation = totalSpeciesName.indexOf("("
                        + splitString2[0] + ")");
                Species species = Species.make(
                        totalSpeciesName.substring(0, intLocation), cg,
                        Integer.parseInt(splitString2[0]));
                // Add the new species to the set of species
                restartCoreSpcs.add(species);
                /*
                 * int species_type = 1; // reacted species for (int i=0; i<numRxnSystems; i++) { SpeciesStatus ss = new
                 * SpeciesStatus(species,species_type,y[i],yprime[i]); speciesStatus[i].put(species, ss); }
                 */
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } catch (IOException e) {
            Logger.logStackTrace(e);
        }
        // Read in edge species
        try {
            FileReader in = new FileReader(new File(
                    System.getProperty("RMG.RestartDir"), "edgeSpecies.txt"));
            BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader, true);
            while (line != null) {
                // The first line of a new species is the user-defined name
                String totalSpeciesName = line;
                String[] splitString1 = totalSpeciesName.split("[(]");
                String[] splitString2 = splitString1[splitString1.length - 1]
                        .split("[)]"); // Change JDM to reflect MRH 2-11-2010
                // The remaining lines are the graph
                Graph g = ChemParser.readChemGraph(reader);
                // Make the ChemGraph, assuming it does not contain a forbidden structure
                ChemGraph cg = null;
                try {
                    cg = ChemGraph.make(g);
                } catch (ForbiddenStructureException e) {
                    Logger.critical("Error reading graph: Graph contains a forbidden structure.\n"
                            + g.toString());
                    System.exit(0);
                }
                // Rewrite the species name ... with the exception of the (#)
                String speciesName = splitString1[0];
                for (int numTokens = 1; numTokens < splitString1.length - 1; ++numTokens) {
                    speciesName += "(" + splitString1[numTokens];
                }
                // Make the species
                Species species = Species.make(speciesName, cg,
                        Integer.parseInt(splitString2[0]));
                // Add the new species to the set of species
                restartEdgeSpcs.add(species);
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } catch (IOException e) {
            Logger.logStackTrace(e);
        }
    }

    public void readRestartReactions() {
        // Grab the IDs from the core species
        int[] coreSpcsIds = new int[restartCoreSpcs.size()];
        int i = 0;
        for (Iterator iter = restartCoreSpcs.iterator(); iter.hasNext();) {
            Species spcs = (Species) iter.next();
            coreSpcsIds[i] = spcs.getID();
            ++i;
        }
        Logger.info("Reading reactions from Restart folder");
        // Read in core reactions
        try {
            FileReader in = new FileReader(new File(
                    System.getProperty("RMG.RestartDir"), "coreReactions.txt"));
            BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader, true);
            // Determine units of Ea
            StringTokenizer st = new StringTokenizer(line);
            String tempString = st.nextToken();
            String EaUnits = st.nextToken();
            line = ChemParser.readMeaningfulLine(reader, true);
            while (line != null) {
                if (!line.trim().equals("DUP")) {
                    Reaction r = ChemParser.parseRestartReaction(line,
                            coreSpcsIds, "core", EaUnits);
                    Iterator rxnIter = restartCoreRxns.iterator();
                    boolean foundRxn = false;
                    while (rxnIter.hasNext()) {
                        Reaction old = (Reaction) rxnIter.next();
                        if (old.equals(r)) {
                            old.addAdditionalKinetics(r.getKinetics()[0], 1,
                                    true);
                            foundRxn = true;
                            break;
                        }
                    }
                    if (!foundRxn) {
                        if (r.hasReverseReaction())
                            r.generateReverseReaction();
                        restartCoreRxns.add(r);
                    }
                }
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } catch (IOException e) {
            Logger.logStackTrace(e);
        }
        /*
         * Read in the pdepreactions.txt file: This file contains third-body, lindemann, and troe reactions A RMG
         * mechanism would only have these reactions if a user specified them in a Seed Mechanism, meaning they are core
         * species & reactions. Place these reactions in a new Seed Mechanism, using the coreSpecies.txt file as the
         * species.txt file.
         */
        SeedMechanism restart_seed_mechanism = null;
        try {
            String path = System.getProperty("RMG.RestartDir");
            restart_seed_mechanism = new SeedMechanism("Restart", path, false,
                    true);
        } catch (IOException e1) {
            Logger.logStackTrace(e1);
        }
        restartCoreRxns.addAll(restart_seed_mechanism.getReactionSet());
        // Read in edge reactions
        try {
            FileReader in = new FileReader(new File(
                    System.getProperty("RMG.RestartDir"), "edgeReactions.txt"));
            BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader, true);
            // Determine units of Ea
            StringTokenizer st = new StringTokenizer(line);
            String tempString = st.nextToken();
            String EaUnits = st.nextToken();
            line = ChemParser.readMeaningfulLine(reader, true);
            while (line != null) {
                if (!line.trim().equals("DUP")) {
                    Reaction r = ChemParser.parseRestartReaction(line,
                            coreSpcsIds, "edge", EaUnits);
                    Iterator rxnIter = restartEdgeRxns.iterator();
                    boolean foundRxn = false;
                    while (rxnIter.hasNext()) {
                        Reaction old = (Reaction) rxnIter.next();
                        if (old.equals(r)) {
                            old.addAdditionalKinetics(r.getKinetics()[0], 1,
                                    true);
                            foundRxn = true;
                            break;
                        }
                    }
                    if (!foundRxn) {
                        r.generateReverseReaction();
                        restartEdgeRxns.add(r);
                    }
                }
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } catch (IOException e) {
            Logger.logStackTrace(e);
        }
    }

    public LinkedHashMap getRestartSpeciesStatus(int i) {
        LinkedHashMap speciesStatus = new LinkedHashMap();
        try {
            FileReader in = new FileReader(new File(
                    System.getProperty("RMG.RestartDir"), "coreSpecies.txt"));
            BufferedReader reader = new BufferedReader(in);
            Integer numRxnSystems = Integer.parseInt(ChemParser
                    .readMeaningfulLine(reader, true));
            String line = ChemParser.readMeaningfulLine(reader, true);
            while (line != null) {
                // The first line of a new species is the user-defined name
                String totalSpeciesName = line;
                String[] splitString1 = totalSpeciesName.split("[(]");
                String[] splitString2 = splitString1[1].split("[)]");
                double y = 0.0;
                double yprime = 0.0;
                for (int j = 0; j < numRxnSystems; j++) {
                    StringTokenizer st = new StringTokenizer(
                            ChemParser.readMeaningfulLine(reader, true));
                    if (j == i) {
                        y = Double.parseDouble(st.nextToken());
                        yprime = Double.parseDouble(st.nextToken());
                    }
                }
                // The remaining lines are the graph
                Graph g = ChemParser.readChemGraph(reader);
                // Make the ChemGraph, assuming it does not contain a forbidden structure
                ChemGraph cg = null;
                try {
                    cg = ChemGraph.make(g);
                } catch (ForbiddenStructureException e) {
                    Logger.critical("Error reading graph: Graph contains a forbidden structure.\n"
                            + g.toString());
                    System.exit(0);
                }
                // Make the species
                Species species = Species.make(splitString1[0], cg);
                // Add the new species to the set of species
                // restartCoreSpcs.add(species);
                int species_type = 1; // reacted species
                SpeciesStatus ss = new SpeciesStatus(species, species_type, y,
                        yprime);
                speciesStatus.put(species, ss);
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } catch (IOException e) {
            Logger.logStackTrace(e);
        }
        return speciesStatus;
    }

    public void putRestartSpeciesInInitialStatus(InitialStatus is, int i) {
        try {
            FileReader in = new FileReader(new File(
                    System.getProperty("RMG.RestartDir"), "coreSpecies.txt"));
            BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader, true);
            while (line != null) {
                // The first line of a new species is the user-defined name
                String totalSpeciesName = line;
                String[] splitString1 = totalSpeciesName.split("[(]");
                String[] splitString2 = splitString1[1].split("[)]");
                // The remaining lines are the graph
                Graph g = ChemParser.readChemGraph(reader);
                // Make the ChemGraph, assuming it does not contain a forbidden structure
                ChemGraph cg = null;
                try {
                    cg = ChemGraph.make(g);
                } catch (ForbiddenStructureException e) {
                    Logger.critical("Error reading graph: Graph contains a forbidden structure.\n"
                            + g.toString());
                    System.exit(0);
                }
                // Make the species
                Species species = Species.make(splitString1[0], cg);
                // Add the new species to the set of species
                // restartCoreSpcs.add(species);
                if (is.getSpeciesStatus(species) == null) {
                    SpeciesStatus ss = new SpeciesStatus(species, 1, 0.0, 0.0);
                    is.putSpeciesStatus(ss);
                }
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } catch (IOException e) {
            Logger.logStackTrace(e);
        }
    }

    public void readPDepNetworks() {
        LinkedList allNetworks = PDepNetwork.getNetworks();
        try {
            FileReader in = new FileReader(new File(
                    System.getProperty("RMG.RestartDir"), "pdepnetworks.txt"));
            BufferedReader reader = new BufferedReader(in);
            StringTokenizer st = new StringTokenizer(
                    ChemParser.readMeaningfulLine(reader, true));
            String tempString = st.nextToken();
            String EaUnits = st.nextToken();
            st = new StringTokenizer(
                    ChemParser.readMeaningfulLine(reader, true));
            tempString = st.nextToken();
            int numFameTs = Integer.parseInt(st.nextToken());
            st = new StringTokenizer(
                    ChemParser.readMeaningfulLine(reader, true));
            tempString = st.nextToken();
            int numFamePs = Integer.parseInt(st.nextToken());
            st = new StringTokenizer(
                    ChemParser.readMeaningfulLine(reader, true));
            tempString = st.nextToken();
            int numChebyTs = Integer.parseInt(st.nextToken());
            st = new StringTokenizer(
                    ChemParser.readMeaningfulLine(reader, true));
            tempString = st.nextToken();
            int numChebyPs = Integer.parseInt(st.nextToken());
            st = new StringTokenizer(
                    ChemParser.readMeaningfulLine(reader, true));
            tempString = st.nextToken();
            int numPlogs = Integer.parseInt(st.nextToken());
            double[][] rateCoefficients = new double[numFameTs][numFamePs];
            double[][] chebyPolys = new double[numChebyTs][numChebyPs];
            Kinetics[] plogKinetics = new Kinetics[numPlogs];
            String line = ChemParser.readMeaningfulLine(reader, true); // line should be "PDepNetwork #"
            while (line != null) {
                line = ChemParser.readMeaningfulLine(reader, true); // line should now be "netReactionList:"
                PDepNetwork newNetwork = new PDepNetwork();
                LinkedList netRxns = newNetwork.getNetReactions();
                LinkedList nonincludeRxns = newNetwork
                        .getNonincludedReactions();
                line = ChemParser.readMeaningfulLine(reader, true); // line is either data or "nonIncludedReactionList"
                // If line is "nonincludedreactionlist", we need to skip over this while loop
                if (!line.toLowerCase().startsWith("nonincludedreactionlist")) {
                    while (!line.toLowerCase().startsWith(
                            "nonincludedreactionlist")) {
                        // Read in the forward rxn
                        String[] reactsANDprods = line.split("\\-->");
                        /*
                         * Determine if netReaction is reversible or irreversible
                         */
                        boolean reactionIsReversible = true;
                        if (reactsANDprods.length == 2)
                            reactionIsReversible = false;
                        else
                            reactsANDprods = line.split("\\<=>");
                        PDepIsomer Reactants = parseIsomerFromRestartFile(reactsANDprods[0]
                                .trim());
                        PDepIsomer Products = parseIsomerFromRestartFile(reactsANDprods[1]
                                .trim());
                        newNetwork.addIsomer(Reactants);
                        newNetwork.addIsomer(Products);
                        rateCoefficients = parseRateCoeffsFromRestartFile(
                                numFameTs, numFamePs, reader);
                        PDepRateConstant pdepk = parsePDepRateConstantFromRestartFile(
                                reader, numChebyTs, numChebyPs,
                                rateCoefficients, numPlogs, EaUnits);
                        PDepReaction forward = new PDepReaction(Reactants,
                                Products, pdepk);
                        // Read in the reverse reaction
                        if (reactionIsReversible) {
                            PDepReaction reverse = new PDepReaction(Products,
                                    Reactants, new PDepRateConstant());
                            reverse.setPDepRateConstant(null);
                            reverse.setReverseReaction(forward);
                            forward.setReverseReaction(reverse);
                        } else {
                            PDepReaction reverse = null;
                            forward.setReverseReaction(reverse);
                        }
                        netRxns.add(forward);
                        line = ChemParser.readMeaningfulLine(reader, true);
                    }
                }
                // This loop ends once line == "nonIncludedReactionList"
                line = ChemParser.readMeaningfulLine(reader, true); // line is either data or "pathReactionList"
                if (!line.toLowerCase().startsWith("pathreactionList")) {
                    while (!line.toLowerCase().startsWith("pathreactionlist")) {
                        // Read in the forward rxn
                        String[] reactsANDprods = line.split("\\-->");
                        /*
                         * Determine if nonIncludedReaction is reversible or irreversible
                         */
                        boolean reactionIsReversible = true;
                        if (reactsANDprods.length == 2)
                            reactionIsReversible = false;
                        else
                            reactsANDprods = line.split("\\<=>");
                        PDepIsomer Reactants = parseIsomerFromRestartFile(reactsANDprods[0]
                                .trim());
                        PDepIsomer Products = parseIsomerFromRestartFile(reactsANDprods[1]
                                .trim());
                        newNetwork.addIsomer(Reactants);
                        newNetwork.addIsomer(Products);
                        rateCoefficients = parseRateCoeffsFromRestartFile(
                                numFameTs, numFamePs, reader);
                        PDepRateConstant pdepk = parsePDepRateConstantFromRestartFile(
                                reader, numChebyTs, numChebyPs,
                                rateCoefficients, numPlogs, EaUnits);
                        PDepReaction forward = new PDepReaction(Reactants,
                                Products, pdepk);
                        // Read in the reverse reaction
                        if (reactionIsReversible) {
                            PDepReaction reverse = new PDepReaction(Products,
                                    Reactants, new PDepRateConstant());
                            reverse.setPDepRateConstant(null);
                            reverse.setReverseReaction(forward);
                            forward.setReverseReaction(reverse);
                        } else {
                            PDepReaction reverse = null;
                            forward.setReverseReaction(reverse);
                        }
                        nonincludeRxns.add(forward);
                        line = ChemParser.readMeaningfulLine(reader, true);
                    }
                }
                // This loop ends once line == "pathReactionList"
                line = ChemParser.readMeaningfulLine(reader, true); // line is either data or "PDepNetwork #_" or null
// (end of file)
                while (line != null
                        && !line.toLowerCase().startsWith("pdepnetwork")) {
                    st = new StringTokenizer(line);
                    int direction = Integer.parseInt(st.nextToken());
                    // First token is the rxn structure: A+B=C+D
                    // Note: Up to 3 reactants/products allowed
                    // : Either "=" or "=>" will separate reactants and products
                    String structure = st.nextToken();
                    // Separate the reactants from the products
                    boolean generateReverse = false;
                    String[] reactsANDprods = structure.split("\\=>");
                    if (reactsANDprods.length == 1) {
                        reactsANDprods = structure.split("[=]");
                        generateReverse = true;
                    }
                    SpeciesDictionary sd = SpeciesDictionary.getInstance();
                    LinkedList r = ChemParser.parseReactionSpecies(sd,
                            reactsANDprods[0]);
                    LinkedList p = ChemParser.parseReactionSpecies(sd,
                            reactsANDprods[1]);
                    Structure s = new Structure(r, p);
                    s.setDirection(direction);
                    // Next three tokens are the modified Arrhenius parameters
                    double rxn_A = Double.parseDouble(st.nextToken());
                    double rxn_n = Double.parseDouble(st.nextToken());
                    double rxn_E = Double.parseDouble(st.nextToken());
                    if (EaUnits.equals("cal/mol"))
                        rxn_E = rxn_E / 1000;
                    else if (EaUnits.equals("J/mol"))
                        rxn_E = rxn_E / 4.184 / 1000;
                    else if (EaUnits.equals("kJ/mol"))
                        rxn_E = rxn_E / 4.184;
                    else if (EaUnits.equals("Kelvins"))
                        rxn_E = rxn_E * 1.987;
                    UncertainDouble uA = new UncertainDouble(rxn_A, 0.0, "A");
                    UncertainDouble un = new UncertainDouble(rxn_n, 0.0, "A");
                    UncertainDouble uE = new UncertainDouble(rxn_E, 0.0, "A");
                    // The remaining tokens are comments
                    String comments = "";
                    if (st.hasMoreTokens()) {
                        String beginningOfComments = st.nextToken();
                        int startIndex = line.indexOf(beginningOfComments);
                        comments = line.substring(startIndex);
                    }
                    if (comments.startsWith("!"))
                        comments = comments.substring(1);
// while (st.hasMoreTokens()) {
// comments += st.nextToken();
// }
                    ArrheniusKinetics[] k = new ArrheniusKinetics[1];
                    k[0] = new ArrheniusKinetics(uA, un, uE, "", 1, "",
                            comments);
                    Reaction pathRxn = new Reaction();
                    // if (direction == 1)
                    // pathRxn = Reaction.makeReaction(s,k,generateReverse);
                    // else
                    // pathRxn = Reaction.makeReaction(s.generateReverseStructure(),k,generateReverse);
                    pathRxn = Reaction.makeReaction(s, k, generateReverse);
                    PDepIsomer Reactants = new PDepIsomer(r);
                    PDepIsomer Products = new PDepIsomer(p);
                    PDepReaction pdeppathrxn = new PDepReaction(Reactants,
                            Products, pathRxn);
                    newNetwork.addReaction(pdeppathrxn, true);
                    newNetwork.addIsomer(Reactants);
                    newNetwork.addIsomer(Products);
                    line = ChemParser.readMeaningfulLine(reader, true);
                }
                newNetwork.setAltered(false);
                PDepNetwork.getNetworks().add(newNetwork);
            }
        } catch (IOException e) {
            Logger.logStackTrace(e);
        }
    }

    public PDepIsomer parseIsomerFromRestartFile(String p_string) {
        SpeciesDictionary sd = SpeciesDictionary.getInstance();
        PDepIsomer isomer = null;
        if (p_string.contains("+")) {
            String[] indivReacts = p_string.split("[+]");
            String name = indivReacts[0].trim();
            Species spc1 = sd.getSpeciesFromNameID(name);
            if (spc1 == null) {
                spc1 = getSpeciesBySPCName(name, sd);
            }
            name = indivReacts[1].trim();
            String[] nameANDincluded = name.split("\\(included =");
            Species spc2 = sd.getSpeciesFromNameID(nameANDincluded[0].trim());
            if (spc2 == null) {
                spc2 = getSpeciesBySPCName(name, sd);
            }
            boolean isIncluded = Boolean.parseBoolean(nameANDincluded[1]
                    .substring(0, nameANDincluded[1].length() - 1));
            isomer = new PDepIsomer(spc1, spc2, isIncluded);
        } else {
            String name = p_string.trim();
            /*
             * Separate the (included =boolean) portion of the string from the name of the Isomer
             */
            String[] nameANDincluded = name.split("\\(included =");
            Species spc = sd.getSpeciesFromNameID(nameANDincluded[0].trim());
            if (spc == null) {
                spc = getSpeciesBySPCName(name, sd);
            }
            boolean isIncluded = Boolean.parseBoolean(nameANDincluded[1]
                    .substring(0, nameANDincluded[1].length() - 1));
            isomer = new PDepIsomer(spc, isIncluded);
        }
        return isomer;
    }

    public double[][] parseRateCoeffsFromRestartFile(int numFameTs,
            int numFamePs, BufferedReader reader) {
        double[][] rateCoefficients = new double[numFameTs][numFamePs];
        for (int i = 0; i < numFameTs; i++) {
            StringTokenizer st = new StringTokenizer(
                    ChemParser.readMeaningfulLine(reader, true));
            for (int j = 0; j < numFamePs; j++) {
                rateCoefficients[i][j] = Double.parseDouble(st.nextToken());
            }
        }
        return rateCoefficients;
    }

    public PDepRateConstant parsePDepRateConstantFromRestartFile(
            BufferedReader reader, int numChebyTs, int numChebyPs,
            double[][] rateCoefficients, int numPlogs, String EaUnits) {
        PDepRateConstant pdepk = null;
        if (numChebyTs > 0) {
            double chebyPolys[][] = new double[numChebyTs][numChebyPs];
            for (int i = 0; i < numChebyTs; i++) {
                StringTokenizer st = new StringTokenizer(
                        ChemParser.readMeaningfulLine(reader, true));
                for (int j = 0; j < numChebyPs; j++) {
                    chebyPolys[i][j] = Double.parseDouble(st.nextToken());
                }
            }
            ChebyshevPolynomials chebyshev = new ChebyshevPolynomials(
                    numChebyTs, ChebyshevPolynomials.getDefaultTlow(),
                    ChebyshevPolynomials.getDefaultTup(), numChebyPs,
                    ChebyshevPolynomials.getDefaultPlow(),
                    ChebyshevPolynomials.getDefaultPup(), chebyPolys);
            pdepk = new PDepRateConstant(rateCoefficients, chebyshev);
        } else if (numPlogs > 0) {
            // ***note: PLOG uses the same units for Ea and A as Arrhenius expressions; see
// https://github.com/GreenGroup/RMG-Java/commit/2947e7b8d5b1e3e19543f2489990fa42e43ecad2#commitcomment-844009
            // gmagoon 1/6/12: restart feature doesn't currently seem to support alternative A units, so I haven't made
// any updates to fix this for PLOG in my current round of changes
            PDepArrheniusKinetics pdepAK = new PDepArrheniusKinetics(numPlogs);
            for (int i = 0; i < numPlogs; i++) {
                StringTokenizer st = new StringTokenizer(
                        ChemParser.readMeaningfulLine(reader, true));
                Pressure p = new Pressure(Double.parseDouble(st.nextToken()),
                        "Pa");
                double A = Double.parseDouble(st.nextToken());
                UncertainDouble dA = new UncertainDouble(A, 0.0, "A");
                UncertainDouble dn = new UncertainDouble(Double.parseDouble(st
                        .nextToken()), 0.0, "A");
                double Ea = Double.parseDouble(st.nextToken());
                if (EaUnits.equals("cal/mol"))
                    Ea = Ea / 1000;
                else if (EaUnits.equals("J/mol"))
                    Ea = Ea / 4.184 / 1000;
                else if (EaUnits.equals("kJ/mol"))
                    Ea = Ea / 4.184;
                else if (EaUnits.equals("Kelvins"))
                    Ea = Ea * 1.987;
                UncertainDouble dE = new UncertainDouble(Ea, 0.0, "A");
                ArrheniusKinetics k = new ArrheniusKinetics(dA, dn, dE, "", 1,
                        "", "");
                pdepAK.setKinetics(i, p, k);
                pdepk = new PDepRateConstant(rateCoefficients, pdepAK);
            }
        }
        return pdepk;
    }

    /**
     * MRH 14Jan2010 getSpeciesBySPCName Input: String name - Name of species, normally chemical formula followed by
     * "J"s for radicals, and then (#) SpeciesDictionary sd This method was originally written as a complement to the
     * method readPDepNetworks. jdmo found a bug with the readrestart option. The bug was that the method was attempting
     * to add a null species to the Isomer list. The null species resulted from searching the SpeciesDictionary by
     * chemkinName (e.g. C4H8OJJ(48)), when the chemkinName present in the dictionary was SPC(48).
     */
    public Species getSpeciesBySPCName(String name, SpeciesDictionary sd) {
        String[] nameFromNumber = name.split("\\(");
        String newName = "SPC(" + nameFromNumber[1];
        return sd.getSpeciesFromChemkinName(newName);
    }

    /**
     * MRH 12-Jun-2009 Function initializes the model's core and edge. The initial core species always consists of the
     * species contained in the condition.txt file. If seed mechanisms exist, those species (and the reactions given in
     * the seed mechanism) are also added to the core. The initial edge species/reactions are determined by reacting the
     * core species by one full iteration.
     */
    public void initializeCoreEdgeModel() {
        LinkedHashSet allInitialCoreSpecies = new LinkedHashSet();
        LinkedHashSet allInitialCoreRxns = new LinkedHashSet();
        if (readrestart) {
            allInitialCoreSpecies.addAll(restartCoreSpcs);
            allInitialCoreRxns.addAll(restartCoreRxns);
        }
        // Add the species from the condition.txt (input) file
        allInitialCoreSpecies.addAll(getSpeciesSeed());
        // Add the species from the seed mechanisms, if they exist
        if (hasSeedMechanisms()) {
            allInitialCoreSpecies.addAll(getSeedMechanism().getSpeciesSet());
            allInitialCoreRxns.addAll(getSeedMechanism().getReactionSet());
        }
        CoreEdgeReactionModel cerm = new CoreEdgeReactionModel(
                allInitialCoreSpecies, allInitialCoreRxns);
        // Store the seed mechanism in the CERM so that we can access it from there when writing chemkin files:
        if (hasSeedMechanisms()) {
            cerm.setSeedMechanism(getSeedMechanism());
        }
        if (readrestart) {
            cerm.addUnreactedSpeciesSet(restartEdgeSpcs);
            cerm.addUnreactedReactionSet(restartEdgeRxns);
        }
        setReactionModel(cerm);
        PDepNetwork.reactionModel = getReactionModel();
        PDepNetwork.reactionSystem = (ReactionSystem) getReactionSystemList()
                .get(0);
        // Determine initial set of reactions and edge species using only the
        // species enumerated in the input file and the seed mechanisms as the core
        if (!readrestart) {
            LinkedHashSet reactionSet_withdup;
            LinkedHashSet reactionSet;
            // If Seed Mechanism is present and Generate Reaction is set on
            if (hasSeedMechanisms()
                    && getSeedMechanism().shouldGenerateReactions()) {
                reactionSet_withdup = getLibraryReactionGenerator().react(
                        allInitialCoreSpecies);
                reactionSet_withdup.addAll(getReactionGenerator().react(
                        allInitialCoreSpecies));
                // Removing Duplicates instances of reaction if present
                reactionSet = getLibraryReactionGenerator()
                        .RemoveDuplicateReac(reactionSet_withdup);
            } else {
                reactionSet_withdup = new LinkedHashSet();
                LinkedHashSet tempnewReactionSet = getLibraryReactionGenerator()
                        .react(allInitialCoreSpecies);
                if (tempnewReactionSet.isEmpty()) {
                    Logger.info("No reactions of initial core and seed species found in Reaction Library");
                } else {
                    Logger.info("Reactions of initial core and seed species found in Reaction Library:");
                    Logger.info(tempnewReactionSet.toString());
                }
                // Adds Reactions Found in Library Reaction Generator to Reaction Set
                reactionSet_withdup.addAll(tempnewReactionSet);
                // Generates Reaction from the Reaction Generator and adds them to Reaction Set
                Logger.info("Generating reactions of initial core species using reaction families.");
                Logger.info("(Condition file species will react with seed species, but seed species will not \n"
                        + "react with each other because GenerateReactions is off)");
                // nb. speciesSeed contains the condition file species (not the seed mechanism species)
                for (Iterator iter = speciesSeed.iterator(); iter.hasNext();) {
                    Species spec = (Species) iter.next();
                    reactionSet_withdup.addAll(getReactionGenerator().react(
                            allInitialCoreSpecies, spec, "All"));
                }
                reactionSet = getLibraryReactionGenerator()
                        .RemoveDuplicateReac(reactionSet_withdup);
            }
            // Set initial core-edge reaction model based on above results
            if (reactionModelEnlarger instanceof RateBasedRME) {
                Iterator iter = reactionSet.iterator();
                while (iter.hasNext()) {
                    Reaction r = (Reaction) iter.next();
                    cerm.addReaction(r);
                }
            } else {
                // Only keep the reactions involving bimolecular reactants and bimolecular products
                Iterator iter = reactionSet.iterator();
                while (iter.hasNext()) {
                    Reaction r = (Reaction) iter.next();
                    if (FastMasterEqn.isReactionPressureDependent(r)) {
                        cerm.categorizeReaction(r.getStructure());
                        PDepNetwork.addReactionToNetworks(r);
                    } else {
                        cerm.addReaction(r);
                    }
                }
            }
        }
        /*
         * 22-SEPT-2010 ELSE: If reading in restart files, at the very least, we should react all species present in the
         * input file (speciesSeed) with all of the coreSpecies. Before, when the following else statement was not
         * present, we would completely skip this step. Thus, RMG would come to the first ODE solve, integrate to a very
         * large time, and conclude that the model was both valid and terminated, thereby not adding any new reactions
         * to the core regardless of the conditions stated in the input file. EXAMPLE: A user runs RMG for iso-octane
         * with Restart turned on. The simulation converges and now the user would like to add a small amount of
         * 1-butanol to the input file, while reading in from the Restart files. What should happen, at the very least,
         * is 1-butanol reacts with the other species present in the input file and with the already-known coreSpecies.
         * This will, at a minimum, add these reactions to the core. Whether the model remains validated and terminated
         * depends on the conditions stated in the input file. MRH (mrharper@mit.edu)
         */
        else {
            LinkedHashSet reactionSet_withdup;
            LinkedHashSet reactionSet;
            /*
             * If the user has specified a Seed Mechanism, and that the cross reactions should be generated, generate
             * those here NOTE: Since the coreSpecies from the Restart files are treated as a Seed Mechanism, MRH is
             * inclined to comment out the following lines. Depending on how large the Seed Mechanism and/or Restart
             * files are, RMG could get "stuck" cross-reacting hundreds of species against each other.
             */
// if (hasSeedMechanisms() && getSeedMechanism().shouldGenerateReactions()) {
// reactionSet_withdup = getLibraryReactionGenerator().react(getSeedMechanism().getSpeciesSet());
// reactionSet_withdup.addAll(getReactionGenerator().react(getSeedMechanism().getSpeciesSet()));
// reactionSet = getLibraryReactionGenerator().RemoveDuplicateReac(reactionSet_withdup);
// }
            /*
             * If not, react the species present in the input file against any reaction libraries, and then against all
             * RMG-defined reaction families
             */
// else {
            reactionSet_withdup = new LinkedHashSet();
            LinkedHashSet tempnewReactionSet = getLibraryReactionGenerator()
                    .react(speciesSeed);
            if (!tempnewReactionSet.isEmpty()) {
                Logger.info("Reaction Set Found from Reaction Library "
                        + tempnewReactionSet);
            }
            // Adds Reactions Found in Library Reaction Generator to Reaction Set
            reactionSet_withdup.addAll(tempnewReactionSet);
            // Generates Reaction from the Reaction Generator and adds them to Reaction Set
            for (Iterator iter = speciesSeed.iterator(); iter.hasNext();) {
                Species spec = (Species) iter.next();
                reactionSet_withdup.addAll(getReactionGenerator().react(
                        allInitialCoreSpecies, spec, "All"));
            }
            reactionSet = getLibraryReactionGenerator().RemoveDuplicateReac(
                    reactionSet_withdup);
// }
            // Set initial core-edge reaction model based on above results
            if (reactionModelEnlarger instanceof RateBasedRME) {
                Iterator iter = reactionSet.iterator();
                while (iter.hasNext()) {
                    Reaction r = (Reaction) iter.next();
                    cerm.addReaction(r);
                }
            } else {
                // Only keep the reactions involving bimolecular reactants and bimolecular products
                Iterator iter = reactionSet.iterator();
                while (iter.hasNext()) {
                    Reaction r = (Reaction) iter.next();
                    if (FastMasterEqn.isReactionPressureDependent(r)) {
                        cerm.categorizeReaction(r.getStructure());
                        PDepNetwork.addReactionToNetworks(r);
                    } else {
                        cerm.addReaction(r);
                    }
                }
            }
        }
        for (Integer i = 0; i < reactionSystemList.size(); i++) {
            ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
            rs.setReactionModel(getReactionModel());
        }
        // We cannot return a system with no core reactions, so if this is a case we must add to the core
        while (getReactionModel().isEmpty()
                && !PDepNetwork
                        .hasCoreReactions((CoreEdgeReactionModel) getReactionModel())) {
            for (Integer i = 0; i < reactionSystemList.size(); i++) {
                ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
                if (reactionModelEnlarger instanceof RateBasedPDepRME)
                    rs.initializePDepNetwork();
                rs.appendUnreactedSpeciesStatus(
                        (InitialStatus) initialStatusList.get(i),
                        rs.getPresentTemperature());
            }
            enlargeReactionModel();
        }
        for (Integer i = 0; i < reactionSystemList.size(); i++) {
            ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
            rs.setReactionModel(getReactionModel());
        }
        return;
    }

    // 9/24/07 gmagoon: moved from ReactionSystem.java
    public void initializeCoreEdgeModelWithPKL() {
        initializeCoreEdgeModelWithoutPKL();
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) getReactionModel();
        LinkedHashSet primarySpeciesSet = getPrimaryKineticLibrary()
                .getSpeciesSet(); // 10/14/07 gmagoon: changed to use getPrimaryReactionLibrary
        LinkedHashSet primaryKineticSet = getPrimaryKineticLibrary()
                .getReactionSet();
        cerm.addReactedSpeciesSet(primarySpeciesSet);
        cerm.addPrimaryKineticSet(primaryKineticSet);
        LinkedHashSet newReactions = getReactionGenerator().react(
                cerm.getReactedSpeciesSet());
        if (reactionModelEnlarger instanceof RateBasedRME)
            cerm.addReactionSet(newReactions);
        else {
            Iterator iter = newReactions.iterator();
            while (iter.hasNext()) {
                Reaction r = (Reaction) iter.next();
                if (r.getReactantNumber() == 2 && r.getProductNumber() == 2) {
                    cerm.addReaction(r);
                }
            }
        }
        return;
        // #]
    }

    // 9/24/07 gmagoon: moved from ReactionSystem.java
    protected void initializeCoreEdgeModelWithoutPKL() {
        CoreEdgeReactionModel cerm = new CoreEdgeReactionModel(
                new LinkedHashSet(getSpeciesSeed()));
        setReactionModel(cerm);
        PDepNetwork.reactionModel = getReactionModel();
        PDepNetwork.reactionSystem = (ReactionSystem) getReactionSystemList()
                .get(0);
        // Determine initial set of reactions and edge species using only the
        // species enumerated in the input file as the core
        LinkedHashSet reactionSet = getReactionGenerator().react(
                getSpeciesSeed());
        reactionSet.addAll(getLibraryReactionGenerator()
                .react(getSpeciesSeed()));
        // Set initial core-edge reaction model based on above results
        if (reactionModelEnlarger instanceof RateBasedRME) {
            // Only keep the reactions involving bimolecular reactants and bimolecular products
            Iterator iter = reactionSet.iterator();
            while (iter.hasNext()) {
                Reaction r = (Reaction) iter.next();
                cerm.addReaction(r);
            }
        } else {
            // Only keep the reactions involving bimolecular reactants and bimolecular products
            Iterator iter = reactionSet.iterator();
            while (iter.hasNext()) {
                Reaction r = (Reaction) iter.next();
                if (FastMasterEqn.isReactionPressureDependent(r)) {
                    cerm.categorizeReaction(r.getStructure());
                    PDepNetwork.addReactionToNetworks(r);
                } else {
                    cerm.addReaction(r);
                }
            }
        }
        // 10/9/07 gmagoon: copy reactionModel to reactionSystem; there may still be scope problems, particularly in
// above elseif statement
        // 10/24/07 gmagoon: want to copy same reaction model to all reactionSystem variables; should probably also make
// similar modifications elsewhere; may or may not need to copy in ...WithPRL function
        for (Integer i = 0; i < reactionSystemList.size(); i++) {
            ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
            rs.setReactionModel(getReactionModel());
        }
        // reactionSystem.setReactionModel(getReactionModel());
        // We cannot return a system with no core reactions, so if this is a case we must add to the core
        while (getReactionModel().isEmpty()
                && !PDepNetwork
                        .hasCoreReactions((CoreEdgeReactionModel) getReactionModel())) {
            for (Integer i = 0; i < reactionSystemList.size(); i++) {
                ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
                if (reactionModelEnlarger instanceof RateBasedPDepRME)
                    rs.initializePDepNetwork();
                rs.appendUnreactedSpeciesStatus(
                        (InitialStatus) initialStatusList.get(i),
                        rs.getPresentTemperature());
            }
            enlargeReactionModel();
        }
        for (Integer i = 0; i < reactionSystemList.size(); i++) {
            ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
            rs.setReactionModel(getReactionModel());
        }
        return;
        // #]
    }

    // ## operation initializeCoreEdgeReactionModel()
    // 9/24/07 gmagoon: moved from ReactionSystem.java
    public void initializeCoreEdgeReactionModel() {
        Logger.info("Initializing core-edge reaction model");
        initializeCoreEdgeModel();
    }

    // 9/24/07 gmagoon: copied from ReactionSystem.java
    public ReactionGenerator getReactionGenerator() {
        return reactionGenerator;
    }

    // 10/4/07 gmagoon: moved from ReactionSystem.java
    public void setReactionGenerator(ReactionGenerator p_ReactionGenerator) {
        reactionGenerator = p_ReactionGenerator;
    }

    // 9/25/07 gmagoon: moved from ReactionSystem.java
    // 10/24/07 gmagoon: changed to use reactionSystemList
    // ## operation enlargeReactionModel()
    public void enlargeReactionModel() {
        // #[ operation enlargeReactionModel()
        if (reactionModelEnlarger == null)
            throw new NullPointerException("ReactionModelEnlarger");
        Logger.info("");
        Logger.info("Enlarging reaction model");
        reactionModelEnlarger.enlargeReactionModel(reactionSystemList,
                reactionModel, validList);
        return;
        // #]
    }

    public void pruneReactionModel(LinkedHashMap unprunableSpecies) {
        Runtime runtime = Runtime.getRuntime();
        LinkedHashMap prunableSpeciesMap = new LinkedHashMap();
        // check whether all the reaction systems reached target conversion/time
        boolean allReachedTarget = true;
        for (Integer i = 0; i < reactionSystemList.size(); i++) {
            JDAS ds = (JDAS) ((ReactionSystem) reactionSystemList.get(i))
                    .getDynamicSimulator();
            if (!ds.targetReached)
                allReachedTarget = false;
        }
        JDAS ds0 = (JDAS) ((ReactionSystem) reactionSystemList.get(0))
                .getDynamicSimulator(); // get the first reactionSystem dynamic simulator
        // prune the reaction model if AUTO is being used, and all reaction systems have reached target time/conversion,
// and edgeTol is non-zero (and positive, obviously), and if there are a sufficient number of species in the reaction
// model (edge + core)
        if (JDAS.autoflag
                && allReachedTarget
                && edgeTol > 0
                && (((CoreEdgeReactionModel) reactionModel).getEdge()
                        .getSpeciesNumber() + reactionModel.getSpeciesNumber()) >= minSpeciesForPruning) {
            int numberToBePruned = ((CoreEdgeReactionModel) reactionModel)
                    .getEdge().getSpeciesNumber() - maxEdgeSpeciesAfterPruning;
            // System.out.println("PDep Pruning DEBUG:\nThe number of species in the model's edge, before pruning: " +
// ((CoreEdgeReactionModel)reactionModel).getEdge().getSpeciesNumber());
            // System.out.println("PDep Pruning DEBUG:\nRMG thinks the following number of species" +
            // " needs to be pruned: " + numberToBePruned);
            Iterator iter = JDAS.edgeID.keySet().iterator();// determine the maximum edge flux ratio for each edge
// species
            while (iter.hasNext()) {
                Species spe = (Species) iter.next();
                Integer id = (Integer) JDAS.edgeID.get(spe);
                double maxmaxRatio = ds0.maxEdgeFluxRatio[id - 1];
                boolean prunable = ds0.prunableSpecies[id - 1];
                // go through the rest of the reaction systems to see if there are higher max flux ratios
                for (Integer i = 1; i < reactionSystemList.size(); i++) {
                    JDAS ds = (JDAS) ((ReactionSystem) reactionSystemList
                            .get(i)).getDynamicSimulator();
                    if (ds.maxEdgeFluxRatio[id - 1] > maxmaxRatio)
                        maxmaxRatio = ds.maxEdgeFluxRatio[id - 1];
                    if (!ds.prunableSpecies[id - 1])
                        prunable = false;// probably redundant: if the conc. is zero in one system, it should be zero in
// all systems, but it is included for completeness
                }
                // if the species is "prunable" (i.e. it doesn't have any reactions producing it with zero flux), add it
// to the prunableSpeciesMap
                if (prunable) {
                    prunableSpeciesMap.put(spe, maxmaxRatio);
                }
            }
            // repeat with the edgeLeakID; if a species appears in both lists, it will be prunable only if it is
// prunable in both cases, and the sum of maximum edgeFlux + maximum edgeLeakFlux (for each reaction system) will be
// considered; this will be a conservative overestimate of maximum (edgeFlux+edgeLeakFlux)
            iter = JDAS.edgeLeakID.keySet().iterator();
            while (iter.hasNext()) {
                Species spe = (Species) iter.next();
                Integer id = (Integer) JDAS.edgeLeakID.get(spe);
                // check whether the same species is in edgeID
                if (JDAS.edgeID.containsKey(spe)) {// the species exists in edgeID
                    if (prunableSpeciesMap.containsKey(spe)) {// the species was determined to be "prunable" based on
// edgeID
                        Integer idEdge = (Integer) JDAS.edgeID.get(spe);
                        double maxmaxRatio = ds0.maxEdgeFluxRatio[id - 1]
                                + ds0.maxEdgeFluxRatio[idEdge - 1];
                        boolean prunable = ds0.prunableSpecies[id - 1];
                        // go through the rest of the reaction systems to see if there are higher max flux ratios
                        for (Integer i = 1; i < reactionSystemList.size(); i++) {
                            JDAS ds = (JDAS) ((ReactionSystem) reactionSystemList
                                    .get(i)).getDynamicSimulator();
                            if (ds.maxEdgeFluxRatio[id - 1]
                                    + ds.maxEdgeFluxRatio[idEdge - 1] > maxmaxRatio)
                                maxmaxRatio = ds.maxEdgeFluxRatio[id - 1]
                                        + ds.maxEdgeFluxRatio[idEdge - 1];
                            if (!ds.prunableSpecies[id - 1])
                                prunable = false;// probably redundant: if the conc. is zero in one system, it should be
// zero in all systems, but it is included for completeness
                        }
                        if (prunable) {// if the species is "prunable" (i.e. it doesn't have any reactions producing it
// with zero flux), replace with the newly determined maxmaxRatio
                            prunableSpeciesMap.remove(spe);
                            prunableSpeciesMap.put(spe, maxmaxRatio);
                        } else {// otherwise, the species is not prunable in both edgeID and edgeLeakID and should be
// removed from the prunable species map
                            prunableSpeciesMap.remove(spe);
                        }
                    }
                } else {// the species is new
                    double maxmaxRatio = ds0.maxEdgeFluxRatio[id - 1];
                    boolean prunable = ds0.prunableSpecies[id - 1];
                    // go through the rest of the reaction systems to see if there are higher max flux ratios
                    for (Integer i = 1; i < reactionSystemList.size(); i++) {
                        JDAS ds = (JDAS) ((ReactionSystem) reactionSystemList
                                .get(i)).getDynamicSimulator();
                        if (ds.maxEdgeFluxRatio[id - 1] > maxmaxRatio)
                            maxmaxRatio = ds.maxEdgeFluxRatio[id - 1];
                        if (!ds.prunableSpecies[id - 1])
                            prunable = false;// probably redundant: if the conc. is zero in one system, it should be
// zero in all systems, but it is included for completeness
                    }
                    // if the species is "prunable" (i.e. it doesn't have any reactions producing it with zero flux),
// add it to the prunableSpeciesMap
                    if (prunable) {
                        prunableSpeciesMap.put(spe, maxmaxRatio);
                    }
                }
            }
            // at this point prunableSpeciesMap includes ALL prunable species, no matter how large their flux
            // System.out.println("PDep Pruning DEBUG:\nRMG has marked the following number of species" +
            // " as prunable, before checking against explored (included) species: " + prunableSpeciesMap.size());
            // Pressure dependence only: Species that are included in any
            // PDepNetwork are not eligible for pruning, so they must be removed
            // from the map of prunable species
            if (reactionModelEnlarger instanceof RateBasedPDepRME) {
                LinkedList speciesToRemove = new LinkedList();
                for (iter = prunableSpeciesMap.keySet().iterator(); iter
                        .hasNext();) {
                    Species spec = (Species) iter.next();
                    if (PDepNetwork.isSpeciesIncludedInAnyNetwork(spec))
                        speciesToRemove.add(spec);
                }
                for (iter = speciesToRemove.iterator(); iter.hasNext();) {
                    prunableSpeciesMap.remove(iter.next());
                }
            }
            // System.out.println("PDep Pruning DEBUG:\nRMG now reduced the number of prunable species," +
            // " after checking against explored (included) species, to: " + prunableSpeciesMap.size());
            // sort the prunableSpecies by maxmaxRatio
            // i.e. sort the map by values
            List prunableSpeciesList = new LinkedList(
                    prunableSpeciesMap.entrySet());
            Collections.sort(prunableSpeciesList, new Comparator() {
                public int compare(Object o1, Object o2) {
                    return ((Comparable) ((Map.Entry) (o1)).getValue())
                            .compareTo(((Map.Entry) (o2)).getValue());
                }
            });
            List speciesToPrune = new LinkedList();
            int belowThreshold = 0;
            int lowMaxFlux = 0;
            for (Iterator it = prunableSpeciesList.iterator(); it.hasNext();) {
                Map.Entry entry = (Map.Entry) it.next();
                Species spe = (Species) entry.getKey();
                double maxmaxRatio = (Double) entry.getValue();
                if (maxmaxRatio < edgeTol) {
                    Logger.info(String
                            .format("Edge species %s has a maximum flux ratio (%.2g) lower than edge inclusion threshhold (%.2g) and will be pruned.",
                                    spe.getChemkinName(), maxmaxRatio, edgeTol));
                    speciesToPrune.add(spe);
                    ++belowThreshold;
                } else if (numberToBePruned - speciesToPrune.size() > 0) {
                    if (maxmaxRatio > tolerance) {
                        Logger.warning(String
                                .format("To reach the requested maximum edge size after pruning (%d) would require pruning species with a maximum flux ratio above the requested overall tolerance (%.2g) which is inconsistent.",
                                        maxEdgeSpeciesAfterPruning, tolerance));
                        Logger.warning(String
                                .format("No more species will be pruned this iteration, leaving an edge size of %d. Please increase your overall Error Tolerance if you'd like more species to be pruned.",
                                        maxEdgeSpeciesAfterPruning
                                                + numberToBePruned
                                                - speciesToPrune.size()));
                        break;
                    }
                    Logger.info(String
                            .format("Edge species %s has a low maximum flux ratio (%.2g) and will be pruned to reduce the edge size to the maximum (%d).",
                                    spe.getChemkinName(), maxmaxRatio,
                                    maxEdgeSpeciesAfterPruning));
                    speciesToPrune.add(spe);
                    ++lowMaxFlux;
                } else
                    break; // no more to be pruned
            }
            // System.out.println("PDep Pruning DEBUG:\nRMG has marked the following number of species" +
            // " to be pruned due to max flux ratio lower than threshold: " + belowThreshold);
            // System.out.println("PDep Pruning DEBUG:\nRMG has marked the following number of species" +
            // " to be pruned due to low max flux ratio : " + lowMaxFlux);
            runtime.gc();
            double memoryUsedBeforePruning = (runtime.totalMemory() - runtime
                    .freeMemory()) / 1.0e6;
            // now, speciesToPrune has been filled with species that should be pruned from the edge
            Logger.info("Pruning...");
            // prune species from the edge
            // remove species from the edge and from the species dictionary and from edgeID
            iter = speciesToPrune.iterator();
            while (iter.hasNext()) {
                Species spe = (Species) iter.next();
                writePrunedEdgeSpecies(spe);
                ((CoreEdgeReactionModel) getReactionModel())
                        .getUnreactedSpeciesSet().remove(spe);
                // SpeciesDictionary.getInstance().getSpeciesSet().remove(spe);
                if (!unprunableSpecies.containsValue(spe))
                    SpeciesDictionary.getInstance().remove(spe);
                else
                    Logger.info("Pruning Message: Not removing the following species "
                            + "from the SpeciesDictionary\nas it is present in a Primary Kinetic / Reaction"
                            + " Library\nThe species will still be removed from the Edge of the "
                            + "Reaction Mechanism\n" + spe.toString());
                JDAS.edgeID.remove(spe);
                JDAS.edgeLeakID.remove(spe); // this would get cleaned up in another iteration when edgeLeakID is
// rebuilt, but debugging memory leaks is simpler if we explicitly clear it here.
            }
            // remove reactions from the edge involving pruned species
            ReactionTemplateLibrary rtl = ReactionTemplateLibrary.getINSTANCE();
            iter = ((CoreEdgeReactionModel) getReactionModel())
                    .getUnreactedReactionSet().iterator();
            LinkedHashSet toRemove = new LinkedHashSet();
            while (iter.hasNext()) {
                Reaction reaction = (Reaction) iter.next();
                if (reactionPrunableQ(reaction, speciesToPrune))
                    toRemove.add(reaction);
            }
            iter = toRemove.iterator();
            while (iter.hasNext()) {
                Reaction reaction = (Reaction) iter.next();
                writePrunedEdgeReaction(reaction);
                Reaction reverse = reaction.getReverseReaction();
                ((CoreEdgeReactionModel) reactionModel)
                        .removeFromUnreactedReactionSet(reaction);
                ((CoreEdgeReactionModel) reactionModel)
                        .removeFromUnreactedReactionSet(reverse);
                int number_removed;
                number_removed = rtl
                        .removeFromAllReactionDictionariesByStructure(reaction
                                .getStructure());
                if (number_removed != 1)
                    Logger.info(String.format(
                            "Removed edge forward %s from %s dictionaries",
                            reaction.getStructure(), number_removed));
                if (reverse != null) {
                    number_removed = rtl
                            .removeFromAllReactionDictionariesByStructure(reverse
                                    .getStructure());
                    if (number_removed != 1)
                        Logger.info(String.format(
                                "Removed edge reverse %s from %s dictionaries",
                                reverse.getStructure(), number_removed));
                }
                reaction.prune();
                if (reverse != null)
                    reverse.prune();
            }
            // remove reactions from PDepNetworks in PDep cases
            if (reactionModelEnlarger instanceof RateBasedPDepRME) {
                iter = PDepNetwork.getNetworks().iterator();
                LinkedHashSet pdnToRemove = new LinkedHashSet();
                LinkedHashSet toRemovePath;
                LinkedHashSet toRemoveNet;
                LinkedHashSet toRemoveNonincluded;
                LinkedHashSet toRemoveIsomer;
                while (iter.hasNext()) {
                    PDepNetwork pdn = (PDepNetwork) iter.next();
                    // identify path reactions to remove
                    Iterator rIter = pdn.getPathReactions().iterator();
                    toRemovePath = new LinkedHashSet();
                    while (rIter.hasNext()) {
                        Reaction reaction = (Reaction) rIter.next();
                        try {
                            if (reactionPrunableQ(reaction, speciesToPrune))
                                toRemovePath.add(reaction);
                        } catch (NullPointerException e) {
                            Logger.error("NullPointerException when inspecting Path Reaction");
                            Logger.logStackTrace(e);
                            Logger.error("Path reaction will not be pruned. Here is the network:");
                            Logger.error("List of species to be pruned:");
                            Iterator errIter = speciesToPrune.iterator();
                            while (errIter.hasNext()) {
                                Species errSpe = (Species) errIter.next();
                                Logger.error(errSpe.getName() + "\t"
                                        + errSpe.getFullName() + "\t"
                                        + errSpe.getChemkinName());
                            }
                            try {
                                Logger.error(pdn.toString());
                            } catch (NullPointerException e2) {
                                Logger.error("NullPointerException trying to print PDEpNetwork");
                                Logger.logStackTrace(e2);
                            }
                            Logger.flush();
                        }
                    }
                    // identify net reactions to remove
                    rIter = pdn.getNetReactions().iterator();
                    toRemoveNet = new LinkedHashSet();
                    while (rIter.hasNext()) {
                        Reaction reaction = (Reaction) rIter.next();
                        if (reactionPrunableQ(reaction, speciesToPrune))
                            toRemoveNet.add(reaction);
                    }
                    // identify nonincluded reactions to remove
                    rIter = pdn.getNonincludedReactions().iterator();
                    toRemoveNonincluded = new LinkedHashSet();
                    while (rIter.hasNext()) {
                        Reaction reaction = (Reaction) rIter.next();
                        if (reactionPrunableQ(reaction, speciesToPrune))
                            toRemoveNonincluded.add(reaction);
                    }
                    // identify isomers to remove
                    Iterator iIter = pdn.getIsomers().iterator();
                    toRemoveIsomer = new LinkedHashSet();
                    while (iIter.hasNext()) {
                        PDepIsomer pdi = (PDepIsomer) iIter.next();
                        Iterator isIter = pdi.getSpeciesListIterator();
                        while (isIter.hasNext()) {
                            Species spe = (Species) isIter.next();
                            if (speciesToPrune.contains(spe)
                                    && !toRemove.contains(pdi))
                                toRemoveIsomer.add(pdi);
                        }
                        if (pdi.getSpeciesList().size() == 0
                                && !toRemove.contains(pdi))
                            toRemoveIsomer.add(pdi);// if the pdi doesn't contain any species, schedule it for removal
                    }
                    // remove path reactions
                    Iterator iterRem = toRemovePath.iterator();
                    while (iterRem.hasNext()) {
                        Reaction reaction = (Reaction) iterRem.next();
                        Reaction reverse = reaction.getReverseReaction();
                        pdn.removeFromPathReactionList((PDepReaction) reaction);
                        pdn.removeFromPathReactionList((PDepReaction) reverse);
                        Structure fwd_structure = reaction.getStructure();
                        Structure rev_structure = null;
                        if (reverse != null) {
                            rev_structure = reverse.getStructure();
                        } else {
                            rev_structure = fwd_structure
                                    .generateReverseStructure();
                        }
                        int number_removed;
                        number_removed = rtl
                                .removeFromAllReactionDictionariesByStructure(fwd_structure);
                        if (number_removed != 1)
                            Logger.info(String
                                    .format("Removed path forward %s from %s dictionaries",
                                            fwd_structure, number_removed));
                        number_removed = rtl
                                .removeFromAllReactionDictionariesByStructure(rev_structure);
                        if (number_removed != 1)
                            Logger.info(String
                                    .format("Removed path reverse %s from %s dictionaries",
                                            rev_structure, number_removed));
                        reaction.prune();
                        if (reverse != null)
                            reverse.prune();
                    }
                    // remove net reactions
                    iterRem = toRemoveNet.iterator();
                    while (iterRem.hasNext()) {
                        Reaction reaction = (Reaction) iterRem.next();
                        Reaction reverse = reaction.getReverseReaction();
                        pdn.removeFromNetReactionList((PDepReaction) reaction);
                        pdn.removeFromNetReactionList((PDepReaction) reverse);
                        reaction.prune();
                        if (reverse != null)
                            reverse.prune();
                    }
                    // remove nonincluded reactions
                    iterRem = toRemoveNonincluded.iterator();
                    while (iterRem.hasNext()) {
                        Reaction reaction = (Reaction) iterRem.next();
                        Reaction reverse = reaction.getReverseReaction();
                        pdn.removeFromNonincludedReactionList((PDepReaction) reaction);
                        pdn.removeFromNonincludedReactionList((PDepReaction) reverse);
                        reaction.prune();
                        if (reverse != null)
                            reverse.prune();
                    }
                    // remove isomers
                    iterRem = toRemoveIsomer.iterator();
                    while (iterRem.hasNext()) {
                        PDepIsomer pdi = (PDepIsomer) iterRem.next();
                        pdn.removeFromIsomerList(pdi);
                    }
                    // remove the entire network if the network has no path or net reactions
                    if (pdn.getPathReactions().size() == 0
                            && pdn.getNetReactions().size() == 0)
                        pdnToRemove.add(pdn);
                }
                iter = pdnToRemove.iterator();
                while (iter.hasNext()) {
                    PDepNetwork pdn = (PDepNetwork) iter.next();
                    PDepNetwork.getNetworks().remove(pdn);
                }
            }
            runtime.gc();
            double memoryUsedAfterPruning = (runtime.totalMemory() - runtime
                    .freeMemory()) / 1.0e6;
            Logger.info(String.format("Number of species pruned:    %d",
                    speciesToPrune.size()));
            Logger.info(String.format("Memory used before pruning:  %10.2f MB",
                    memoryUsedBeforePruning));
            Logger.info(String.format("Memory used after pruning:   %10.2f MB",
                    memoryUsedAfterPruning));
            if (memoryUsedAfterPruning < memoryUsedBeforePruning)
                Logger.info(String.format(
                        "Memory recovered by pruning: %10.2f MB",
                        memoryUsedBeforePruning - memoryUsedAfterPruning));
            else if (speciesToPrune.size() > 100)
                // There were a significant number of species pruned, but we didn't recover any memory
                Logger.warning("No memory recovered due to pruning!");
        }
        // System.out.println("PDep Pruning DEBUG:\nThe number of species in the model's edge, after pruning: " +
// ((CoreEdgeReactionModel)reactionModel).getEdge().getSpeciesNumber());
        return;
    }

    // determines whether a reaction can be removed; returns true ; cf. categorizeReaction() in CoreEdgeReactionModel
    // returns true if the reaction involves reactants or products that are in p_prunableSpecies; otherwise returns
// false
    public boolean reactionPrunableQ(Reaction p_reaction,
            Collection p_prunableSpecies) {
        Iterator iter = p_reaction.getReactants();
        while (iter.hasNext()) {
            Species spe = (Species) iter.next();
            if (p_prunableSpecies.contains(spe))
                return true;
        }
        iter = p_reaction.getProducts();
        while (iter.hasNext()) {
            Species spe = (Species) iter.next();
            if (p_prunableSpecies.contains(spe))
                return true;
        }
        return false;
    }

    public boolean hasPrimaryKineticLibrary() {
        if (primaryKineticLibrary == null)
            return false;
        return (primaryKineticLibrary.size() > 0);
    }

    public boolean hasSeedMechanisms() {
        if (getSeedMechanism() == null)
            return false;
        return (seedMechanism.size() > 0);
    }

    // 9/25/07 gmagoon: moved from ReactionSystem.java
    public PrimaryKineticLibrary getPrimaryKineticLibrary() {
        return primaryKineticLibrary;
    }

    // 9/25/07 gmagoon: moved from ReactionSystem.java
    public void setPrimaryKineticLibrary(
            PrimaryKineticLibrary p_PrimaryKineticLibrary) {
        primaryKineticLibrary = p_PrimaryKineticLibrary;
    }

    public ReactionLibrary getReactionLibrary() {
        return ReactionLibrary;
    }

    public void setReactionLibrary(ReactionLibrary p_ReactionLibrary) {
        ReactionLibrary = p_ReactionLibrary;
    }

    // 10/4/07 gmagoon: added
    public LinkedHashSet getSpeciesSeed() {
        return speciesSeed;
    }

    // 10/4/07 gmagoon: added
    public void setSpeciesSeed(LinkedHashSet p_speciesSeed) {
        speciesSeed = p_speciesSeed;
    }

    // 10/4/07 gmagoon: added
    public LibraryReactionGenerator getLibraryReactionGenerator() {
        return lrg;
    }

    // 10/4/07 gmagoon: added
    public void setLibraryReactionGenerator(LibraryReactionGenerator p_lrg) {
        lrg = p_lrg;
    }

    public static Temperature getTemp4BestKinetics() {
        return temp4BestKinetics;
    }

    public static void setTemp4BestKinetics(Temperature firstSysTemp) {
        temp4BestKinetics = firstSysTemp;
    }

    public SeedMechanism getSeedMechanism() {
        return seedMechanism;
    }

    public void setSeedMechanism(SeedMechanism p_seedMechanism) {
        seedMechanism = p_seedMechanism;
    }

    public PrimaryThermoLibrary getPrimaryThermoLibrary() {
        return primaryThermoLibrary;
    }

    public void setPrimaryThermoLibrary(
            PrimaryThermoLibrary p_primaryThermoLibrary) {
        primaryThermoLibrary = p_primaryThermoLibrary;
    }

    public static double getAtol() {
        return atol;
    }

    public boolean runKillableToPreventInfiniteLoop(boolean intermediateSteps,
            int iterationNumber) {
        ReactionSystem rs0 = (ReactionSystem) reactionSystemList.get(0);
        if (!intermediateSteps)// if there are no intermediate steps (for example when using AUTO method), return true;
            return true;
        // if there are intermediate steps, the run is killable if the iteration number exceeds the number of time steps
// / conversions
        else if (rs0.finishController.terminationTester instanceof ReactionTimeTT) {
            if (iterationNumber - 1 > timeStep.size()) { // -1 correction needed since when this is called, iteration
// number has been incremented
                return true;
            }
        } else // the case where intermediate conversions are specified
        if (iterationNumber - 1 > numConversions) { // see above; it is possible there is an off-by-one error here, so
// further testing will be needed
            return true;
        }
        return false; // return false if none of the above criteria are met
    }

    public void readAndMakePKL(BufferedReader reader) throws IOException {
        int Ilib = 0;
        String line = ChemParser.readMeaningfulLine(reader, true);
        while (!line.equals("END")) {
            String name = extractLibraryName(line);
            line = ChemParser.readMeaningfulLine(reader, true);
            String[] tempString = line.split("Location: ");
            String location = tempString[tempString.length - 1].trim();
            String path = System
                    .getProperty("jing.rxn.ReactionLibrary.pathName");
            path += "/" + location;
            if (Ilib == 0) {
                setPrimaryKineticLibrary(new PrimaryKineticLibrary(name, path));
                Ilib++;
            } else {
                getPrimaryKineticLibrary().appendPrimaryKineticLibrary(name,
                        path);
                Ilib++;
            }
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (Ilib == 0) {
            setPrimaryKineticLibrary(null);
        } else
            Logger.info("Primary Kinetic Libraries in use: "
                    + getPrimaryKineticLibrary().getName());
    }

    public void readAndMakeReactionLibrary(BufferedReader reader)
            throws IOException {
        int Ilib = 0;
        String line = ChemParser.readMeaningfulLine(reader, true);
        while (!line.equals("END")) {
            String name = extractLibraryName(line);
            line = ChemParser.readMeaningfulLine(reader, true);
            String[] tempString = line.split("Location: ");
            String location = tempString[tempString.length - 1].trim();
            String path = System
                    .getProperty("jing.rxn.ReactionLibrary.pathName");
            path += "/" + location;
            if (Ilib == 0) {
                setReactionLibrary(new ReactionLibrary(name, path));
                Ilib++;
            } else {
                getReactionLibrary().appendReactionLibrary(name, path);
                Ilib++;
            }
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (Ilib == 0) {
            setReactionLibrary(null);
        } else
            Logger.info("Reaction Libraries in use: "
                    + getReactionLibrary().getName());
    }

    public void readAndMakePTL(BufferedReader reader) {
        int numPTLs = 0;
        String line = ChemParser.readMeaningfulLine(reader, true);
        while (!line.equals("END")) {
            String name = extractLibraryName(line);
            line = ChemParser.readMeaningfulLine(reader, true);
            String[] tempString = line.split("Location: ");
            String path = tempString[tempString.length - 1].trim();
            if (numPTLs == 0) {
                setPrimaryThermoLibrary(new PrimaryThermoLibrary(name, path));
                ++numPTLs;
            } else {
                getPrimaryThermoLibrary()
                        .appendPrimaryThermoLibrary(name, path);
                ++numPTLs;
            }
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (numPTLs == 0)
            setPrimaryThermoLibrary(null);
    }

    public void readExtraForbiddenStructures(BufferedReader reader)
            throws IOException {
        Logger.info("Reading extra forbidden structures from input file.");
        String line = ChemParser.readMeaningfulLine(reader, true);
        while (!line.equals("END")) {
            StringTokenizer token = new StringTokenizer(line);
            String fgname = token.nextToken();
            Graph fgGraph = null;
            try {
                fgGraph = ChemParser.readFGGraph(reader);
            } catch (InvalidGraphFormatException e) {
                Logger.error("Invalid functional group in " + fgname);
                throw new InvalidFunctionalGroupException(fgname + ": "
                        + e.getMessage());
            }
            if (fgGraph == null)
                throw new InvalidFunctionalGroupException(fgname);
            FunctionalGroup fg = FunctionalGroup.makeForbiddenStructureFG(
                    fgname, fgGraph);
            ChemGraph.addForbiddenStructure(fg);
            line = ChemParser.readMeaningfulLine(reader, true);
            Logger.debug(" Forbidden structure: " + fgname);
        }
    }

    public void setSpectroscopicDataMode(String line) {
        StringTokenizer st = new StringTokenizer(line);
        String name = st.nextToken();
        String sdeType = st.nextToken().toLowerCase();
        if (sdeType.equals("frequencygroups") || sdeType.equals("default")) {
            SpectroscopicData.mode = SpectroscopicData.Mode.FREQUENCYGROUPS;
        } else if (sdeType.equals("therfit")
                || sdeType.equals("threefrequencymodel")) {
            SpectroscopicData.mode = SpectroscopicData.Mode.THREEFREQUENCY;
        } else if (sdeType.equals("off") || sdeType.equals("none")) {
            SpectroscopicData.mode = SpectroscopicData.Mode.OFF;
        } else
            throw new InvalidSymbolException(
                    "condition.txt: Unknown SpectroscopicDataEstimator = "
                            + sdeType);
    }

    /**
     * Sets the pressure dependence options to on or off. If on, checks for more options and sets them as well.
     * 
     * @param line
     *            The current line in the condition file; should start with "PressureDependence:"
     * @param reader
     *            The reader currently being used to parse the condition file
     */
    public String setPressureDependenceOptions(String line,
            BufferedReader reader) throws InvalidSymbolException {
        // Determine pressure dependence mode
        StringTokenizer st = new StringTokenizer(line);
        String name = st.nextToken(); // Should be "PressureDependence:"
        String pDepType = st.nextToken();
        if (pDepType.toLowerCase().equals("off")) {
            // No pressure dependence
            reactionModelEnlarger = new RateBasedRME();
            PDepNetwork.generateNetworks = false;
            /*
             * If the Spectroscopic Data Estimator field is set to "Frequency Groups," terminate the RMG job and inform
             * the user to either: a) Set the Spectroscopic Data Estimator field to "off," OR b) Select a
             * pressure-dependent model Before, RMG would read in "Frequency Groups" with no pressure-dependence and
             * carry on. However, the calculated frequencies would not be stored / reported (plus increase the runtime),
             * so no point in calculating them.
             */
            if (SpectroscopicData.mode != SpectroscopicData.mode.OFF) {
                Logger.critical("Terminating RMG simulation: User requested frequency estimation, "
                        + "yet no pressure-dependence.\nSUGGESTION: Set the "
                        + "SpectroscopicDataEstimator field in the input file to 'off'.");
                System.exit(0);
            }
            line = ChemParser.readMeaningfulLine(reader, true);
        } else if (pDepType.toLowerCase().equals("modifiedstrongcollision")
                || pDepType.toLowerCase().equals("reservoirstate")
                || pDepType.toLowerCase().equals("chemdis")) {
            reactionModelEnlarger = new RateBasedPDepRME();
            PDepNetwork.generateNetworks = true;
            // Set pressure dependence method
            if (pDepType.toLowerCase().equals("reservoirstate"))
                ((RateBasedPDepRME) reactionModelEnlarger)
                        .setPDepKineticsEstimator(new FastMasterEqn(
                                FastMasterEqn.Mode.RESERVOIRSTATE));
            else if (pDepType.toLowerCase().equals("modifiedstrongcollision"))
                ((RateBasedPDepRME) reactionModelEnlarger)
                        .setPDepKineticsEstimator(new FastMasterEqn(
                                FastMasterEqn.Mode.STRONGCOLLISION));
            // else if (pDepType.toLowerCase().equals("chemdis"))
            // ((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new Chemdis());
            else
                throw new InvalidSymbolException(
                        "condition.txt: Unknown PressureDependence mode = "
                                + pDepType);
            RateBasedPDepRME pdepModelEnlarger = (RateBasedPDepRME) reactionModelEnlarger;
            // Turn on spectroscopic data estimation if not already on
            if (pdepModelEnlarger.getPDepKineticsEstimator() instanceof FastMasterEqn
                    && SpectroscopicData.mode == SpectroscopicData.Mode.OFF) {
                Logger.warning("Spectroscopic data needed for pressure dependence; switching SpectroscopicDataEstimator to FrequencyGroups.");
                SpectroscopicData.mode = SpectroscopicData.Mode.FREQUENCYGROUPS;
            } else if (pdepModelEnlarger.getPDepKineticsEstimator() instanceof Chemdis
                    && SpectroscopicData.mode != SpectroscopicData.Mode.THREEFREQUENCY) {
                Logger.warning("Switching SpectroscopicDataEstimator to three-frequency model.");
                SpectroscopicData.mode = SpectroscopicData.Mode.THREEFREQUENCY;
            }
            // Optional: MaxAtomsForPressureDependence
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.toLowerCase().startsWith("maxatomsforpressuredependence:")) {
                st = new StringTokenizer(line);
                name = st.nextToken();
                int atoms = Integer.parseInt(st.nextToken());
                FastMasterEqn.setMaxAtoms(atoms);
                line = ChemParser.readMeaningfulLine(reader, true);
            }
            
            // Optional: MaxCarbonAtomsForPressureDependence
            if (line.toLowerCase().startsWith("maxcarbonatomsforpressuredependence:")) {
                st = new StringTokenizer(line);
                name = st.nextToken();
                int atoms = Integer.parseInt(st.nextToken());
                FastMasterEqn.setMaxCarbonAtoms(atoms);
                line = ChemParser.readMeaningfulLine(reader, true);
            }
            
            // Next line must be PDepKineticsModel
            if (line.toLowerCase().startsWith("pdepkineticsmodel:")) {
                st = new StringTokenizer(line);
                name = st.nextToken();
                String pDepKinType = st.nextToken();
                if (pDepKinType.toLowerCase().equals("chebyshev")) {
                    PDepRateConstant
                            .setDefaultMode(PDepRateConstant.Mode.CHEBYSHEV);
                    // Default is to cubic order for basis functions
                    FastMasterEqn.setNumTBasisFuncs(4);
                    FastMasterEqn.setNumPBasisFuncs(4);
                } else if (pDepKinType.toLowerCase().equals("pdeparrhenius"))
                    PDepRateConstant
                            .setDefaultMode(PDepRateConstant.Mode.PDEPARRHENIUS);
                else if (pDepKinType.toLowerCase().equals("rate"))
                    PDepRateConstant.setDefaultMode(PDepRateConstant.Mode.RATE);
                else
                    throw new InvalidSymbolException(
                            "condition.txt: Unknown PDepKineticsModel = "
                                    + pDepKinType);
                // For Chebyshev polynomials, optionally specify the number of
                // temperature and pressure basis functions
                // Such a line would read, e.g.: "PDepKineticsModel: Chebyshev 4 4"
                if (st.hasMoreTokens()
                        && PDepRateConstant.getDefaultMode() == PDepRateConstant.Mode.CHEBYSHEV) {
                    try {
                        int numTBasisFuncs = Integer.parseInt(st.nextToken());
                        int numPBasisFuncs = Integer.parseInt(st.nextToken());
                        FastMasterEqn.setNumTBasisFuncs(numTBasisFuncs);
                        FastMasterEqn.setNumPBasisFuncs(numPBasisFuncs);
                    } catch (NoSuchElementException e) {
                        throw new InvalidSymbolException(
                                "condition.txt: Missing number of pressure basis functions for Chebyshev polynomials.");
                    }
                }
            } else
                throw new InvalidSymbolException(
                        "condition.txt: Missing PDepKineticsModel after PressureDependence line.");
            // Determine temperatures and pressures to use
            // These can be specified automatically using TRange and PRange or
            // manually using Temperatures and Pressures
            Temperature[] temperatures = null;
            Pressure[] pressures = null;
            String Tunits = "K";
            Temperature Tmin = new Temperature(300.0, "K");
            Temperature Tmax = new Temperature(2000.0, "K");
            int Tnumber = 8;
            String Punits = "bar";
            Pressure Pmin = new Pressure(0.01, "bar");
            Pressure Pmax = new Pressure(100.0, "bar");
            int Pnumber = 5;
            // Read next line of input
            line = ChemParser.readMeaningfulLine(reader, true);
            boolean done = !(line.toLowerCase().startsWith("trange:")
                    || line.toLowerCase().startsWith("prange:")
                    || line.toLowerCase().startsWith("temperatures:") || line
                    .toLowerCase().startsWith("pressures:"));
            // Parse lines containing pressure dependence options
            // Possible options are "TRange:", "PRange:", "Temperatures:", and "Pressures:"
            // You must specify either TRange or Temperatures and either PRange or Pressures
            // The order does not matter
            while (!done) {
                st = new StringTokenizer(line);
                name = st.nextToken();
                if (line.toLowerCase().startsWith("trange:")) {
                    Tunits = ChemParser.removeBrace(st.nextToken());
                    Tmin = new Temperature(Double.parseDouble(st.nextToken()),
                            Tunits);
                    Tmax = new Temperature(Double.parseDouble(st.nextToken()),
                            Tunits);
                    Tnumber = Integer.parseInt(st.nextToken());
                } else if (line.toLowerCase().startsWith("prange:")) {
                    Punits = ChemParser.removeBrace(st.nextToken());
                    Pmin = new Pressure(Double.parseDouble(st.nextToken()),
                            Punits);
                    Pmax = new Pressure(Double.parseDouble(st.nextToken()),
                            Punits);
                    Pnumber = Integer.parseInt(st.nextToken());
                } else if (line.toLowerCase().startsWith("temperatures:")) {
                    Tnumber = Integer.parseInt(st.nextToken());
                    Tunits = ChemParser.removeBrace(st.nextToken());
                    temperatures = new Temperature[Tnumber];
                    for (int i = 0; i < Tnumber; i++) {
                        temperatures[i] = new Temperature(Double.parseDouble(st
                                .nextToken()), Tunits);
                    }
                    Tmin = temperatures[0];
                    Tmax = temperatures[Tnumber - 1];
                } else if (line.toLowerCase().startsWith("pressures:")) {
                    Pnumber = Integer.parseInt(st.nextToken());
                    Punits = ChemParser.removeBrace(st.nextToken());
                    pressures = new Pressure[Pnumber];
                    for (int i = 0; i < Pnumber; i++) {
                        pressures[i] = new Pressure(Double.parseDouble(st
                                .nextToken()), Punits);
                    }
                    Pmin = pressures[0];
                    Pmax = pressures[Pnumber - 1];
                }
                // Read next line of input
                line = ChemParser.readMeaningfulLine(reader, true);
                done = !(line.toLowerCase().startsWith("trange:")
                        || line.toLowerCase().startsWith("prange:")
                        || line.toLowerCase().startsWith("temperatures:") || line
                        .toLowerCase().startsWith("pressures:"));
            }
            // Set temperatures and pressures (if not already set manually)
            if (temperatures == null) {
                temperatures = new Temperature[Tnumber];
                if (PDepRateConstant.getDefaultMode() == PDepRateConstant.Mode.CHEBYSHEV) {
                    // Use the Gauss-Chebyshev points
                    // The formula for the Gauss-Chebyshev points was taken from
                    // the Chemkin theory manual
                    for (int i = 1; i <= Tnumber; i++) {
                        double T = -Math.cos((2 * i - 1) * Math.PI
                                / (2 * Tnumber));
                        T = 2.0 / ((1.0 / Tmax.getK() - 1.0 / Tmin.getK()) * T
                                + 1.0 / Tmax.getK() + 1.0 / Tmin.getK());
                        temperatures[i - 1] = new Temperature(T, "K");
                    }
                } else {
                    // Distribute equally on a 1/T basis
                    double slope = (1.0 / Tmax.getK() - 1.0 / Tmin.getK())
                            / (Tnumber - 1);
                    for (int i = 0; i < Tnumber; i++) {
                        double T = 1.0 / (slope * i + 1.0 / Tmin.getK());
                        temperatures[i] = new Temperature(T, "K");
                    }
                }
            }
            if (pressures == null) {
                pressures = new Pressure[Pnumber];
                if (PDepRateConstant.getDefaultMode() == PDepRateConstant.Mode.CHEBYSHEV) {
                    // Use the Gauss-Chebyshev points
                    // The formula for the Gauss-Chebyshev points was taken from
                    // the Chemkin theory manual
                    for (int i = 1; i <= Pnumber; i++) {
                        double P = -Math.cos((2 * i - 1) * Math.PI
                                / (2 * Pnumber));
                        P = Math.pow(
                                10,
                                0.5 * ((Math.log10(Pmax.getBar()) - Math
                                        .log10(Pmin.getBar()))
                                        * P
                                        + Math.log10(Pmax.getBar()) + Math
                                        .log10(Pmin.getBar())));
                        pressures[i - 1] = new Pressure(P, "bar");
                    }
                } else {
                    // Distribute equally on a log P basis
                    double slope = (Math.log10(Pmax.getBar()) - Math.log10(Pmin
                            .getBar())) / (Pnumber - 1);
                    for (int i = 0; i < Pnumber; i++) {
                        double P = Math.pow(10,
                                slope * i + Math.log10(Pmin.getBar()));
                        pressures[i] = new Pressure(P, "bar");
                    }
                }
            }
            FastMasterEqn.setTemperatures(temperatures);
            PDepRateConstant.setTemperatures(temperatures);
            PDepRateConstant.setTMin(Tmin);
            PDepRateConstant.setTMax(Tmax);
            ChebyshevPolynomials.setDefaultTlow(Tmin);
            ChebyshevPolynomials.setDefaultTup(Tmax);
            FastMasterEqn.setPressures(pressures);
            PDepRateConstant.setPressures(pressures);
            PDepRateConstant.setPMin(Pmin);
            PDepRateConstant.setPMax(Pmax);
            ChebyshevPolynomials.setDefaultPlow(Pmin);
            ChebyshevPolynomials.setDefaultPup(Pmax);
            /*
             * New option for input file: DecreaseGrainSize User now has the option to re-run fame with additional
             * grains (smaller grain size) when the p-dep rate exceeds the high-P-limit rate. Default value: off
             */
            if (line.toLowerCase().startsWith("decreasegrainsize")) {
                st = new StringTokenizer(line);
                String tempString = st.nextToken(); // "DecreaseGrainSize:"
                tempString = st.nextToken().trim().toLowerCase();
                if (tempString.equals("on") || tempString.equals("yes")
                        || tempString.equals("true")) {
                    rerunFame = true;
                } else
                    rerunFame = false;
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } else {
            throw new InvalidSymbolException(
                    "condition.txt: Unknown PressureDependence = " + pDepType);
        }
        return line;
    }

    public void createTModel(String line) {
        StringTokenizer st = new StringTokenizer(line);
        String name = st.nextToken();
        String modelType = st.nextToken();
        String unit = st.nextToken();
        unit = ChemParser.removeBrace(unit);
        if (modelType.equals("Constant")) {
            tempList = new LinkedList();
            // read first temperature
            double t = Double.parseDouble(st.nextToken());
            tempList.add(new ConstantTM(t, unit));
            Temperature temp = new Temperature(t, unit);
            Global.lowTemperature = (Temperature) temp.clone();
            Global.highTemperature = (Temperature) temp.clone();
            // read remaining temperatures
            while (st.hasMoreTokens()) {
                t = Double.parseDouble(st.nextToken());
                tempList.add(new ConstantTM(t, unit));
                temp = new Temperature(t, unit);
                if (temp.getK() < Global.lowTemperature.getK())
                    Global.lowTemperature = (Temperature) temp.clone();
                if (temp.getK() > Global.highTemperature.getK())
                    Global.highTemperature = (Temperature) temp.clone();
            }
        } else {
            throw new InvalidSymbolException(
                    "condition.txt: Unknown TemperatureModel = " + modelType);
        }
    }

    public void createPModel(String line) {
        StringTokenizer st = new StringTokenizer(line);
        String name = st.nextToken();
        String modelType = st.nextToken();
        String unit = st.nextToken();
        unit = ChemParser.removeBrace(unit);
        if (modelType.equals("Constant")) {
            presList = new LinkedList();
            // read first pressure
            double p = Double.parseDouble(st.nextToken());
            Pressure pres = new Pressure(p, unit);
            Global.lowPressure = (Pressure) pres.clone();
            Global.highPressure = (Pressure) pres.clone();
            presList.add(new ConstantPM(p, unit));
            // read remaining temperatures
            while (st.hasMoreTokens()) {
                p = Double.parseDouble(st.nextToken());
                presList.add(new ConstantPM(p, unit));
                pres = new Pressure(p, unit);
                if (pres.getBar() < Global.lowPressure.getBar())
                    Global.lowPressure = (Pressure) pres.clone();
                if (pres.getBar() > Global.lowPressure.getBar())
                    Global.highPressure = (Pressure) pres.clone();
            }
        } else {
            throw new InvalidSymbolException(
                    "condition.txt: Unknown PressureModel = " + modelType);
        }
    }

    public LinkedHashMap populateInitialStatusListWithReactiveSpecies(
            BufferedReader reader) throws IOException {
        LinkedHashMap speciesSet = new LinkedHashMap();
        LinkedHashMap speciesStatus = new LinkedHashMap();
        LinkedHashMap speciesFromInputFileSet = new LinkedHashMap();
        int numSpeciesStatus = 0;
        String line = ChemParser.readMeaningfulLine(reader, true);
        while (line != null && !line.equals("END")) {
            StringTokenizer st = new StringTokenizer(line);
            String index = st.nextToken();
            String name = null;
            if (!index.startsWith("("))
                name = index;
            else
                name = st.nextToken();
            // if (restart) name += "("+speciesnum+")";
            // 24Jun2009: MRH
            // Check if the species name begins with a number.
            // If so, terminate the program and inform the user to choose
            // a different name. This is implemented so that the chem.inp
            // file generated will be valid when run in Chemkin
            try {
                int doesNameBeginWithNumber = Integer.parseInt(name.substring(
                        0, 1));
                Logger.critical("\nA species name should not begin with a number."
                        + " Please rename species: " + name + "\n");
                System.exit(0);
            } catch (NumberFormatException e) {
                // We're good
            }
            if (!(st.hasMoreTokens()))
                throw new InvalidSymbolException(
                        "Couldn't find concentration of species: " + name);
            // The next token will be the concentration units
            String unit = st.nextToken();
            unit = ChemParser.removeBrace(unit);
            // Read in the graph (and make the species) before we parse all of the concentrations
            // Will store each SpeciesStatus (NXM where N is the # of species and M is the # of
            // concentrations) in a LinkedHashMap, with a (int) counter as the unique key
            Graph g = ChemParser.readChemGraph(reader);
            ChemGraph cg = null;
            try {
                cg = ChemGraph.make(g);
            } catch (ForbiddenStructureException e) {
                Logger.error("Forbidden Structure:\n" + e.getMessage());
                throw new InvalidSymbolException(
                        "A species in the input file has a forbidden structure.");
            }
            // Make sure to kekulize species in the input file to convert adjlists
            // containing benzene bonds to single/double bonds
            cg.kekulize();
            // System.out.println(name);
            // Check to see if chemgraph already appears in the input file
            addChemGraphToListIfNotPresent_ElseTerminate(
                    speciesFromInputFileSet, cg, name);
            Species species = Species.make(name, cg);
            int numConcentrations = 0; // The number of concentrations read-in for each species
            // The remaining tokens are either:
            // The desired concentrations
            // The flag "unreactive"
            // The flag "constantconcentration"
            // GJB to allow "unreactive" species that only follow user-defined library reactions.
            // They will not react according to RMG reaction families
            boolean IsReactive = true;
            boolean IsConstantConcentration = false;
            double concentration = 0.0;
            while (st.hasMoreTokens()) {
                String reactive = st.nextToken().trim();
                if (reactive.equalsIgnoreCase("unreactive"))
                    IsReactive = false;
                else if (reactive.equalsIgnoreCase("constantconcentration"))
                    IsConstantConcentration = true;
                else {
                    try {
                        concentration = Double.parseDouble(reactive);
                    } catch (NumberFormatException e) {
                        Logger.error(String
                                .format("Unable to read concentration value '%s'. Check syntax of input line '%s'.",
                                        reactive, line));
                        throw e;
                    }
                    if (unit.equals("mole/l") || unit.equals("mol/l")
                            || unit.equals("mole/liter")
                            || unit.equals("mol/liter")) {
                        concentration /= 1000;
                        unit = "mol/cm3";
                    } else if (unit.equals("mole/m3") || unit.equals("mol/m3")) {
                        concentration /= 1000000;
                        unit = "mol/cm3";
                    } else if (unit.equals("molecule/cm3")
                            || unit.equals("molecules/cm3")) {
                        concentration /= 6.022e23;
                    } else if (!unit.equals("mole/cm3")
                            && !unit.equals("mol/cm3")) {
                        Logger.error(String
                                .format("Unable to read concentration units '%s'. Check syntax of input line '%s'.",
                                        unit, line));
                        throw new InvalidUnitException(
                                "Species Concentration in condition.txt!");
                    }
                    SpeciesStatus ss = new SpeciesStatus(species, 1,
                            concentration, 0.0);
                    speciesStatus.put(numSpeciesStatus, ss);
                    ++numSpeciesStatus;
                    ++numConcentrations;
                }
            }
            // Check if the number of concentrations read in is consistent with all previous
            // concentration counts. The first time this function is called, the variable
            // numberOfEquivalenceRatios will be initialized.
            boolean goodToGo = areTheNumberOfConcentrationsConsistent(numConcentrations);
            if (!goodToGo) {
                Logger.critical("\n\nThe number of concentrations ("
                        + numConcentrations
                        + ") supplied for species "
                        + species.getName()
                        + "\nis not consistent with the number of concentrations ("
                        + numberOfEquivalenceRatios + ") "
                        + "supplied for all previously read-in species \n\n"
                        + "Terminating RMG simulation.");
                System.exit(0);
            }
            // Make a SpeciesStatus and store it in the LinkedHashMap
// double flux = 0;
// int species_type = 1; // reacted species
            species.setReactivity(IsReactive); // GJB
            species.setConstantConcentration(IsConstantConcentration);
            speciesSet.put(name, species);
            getSpeciesSeed().add(species);
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (line == null)
            throw new IOException(
                    "Unexpected end of input file before END of reactive SpeciesStatus section");
        ReactionTime initial = new ReactionTime(0, "S");
        // 10/23/07 gmagoon: modified for handling multiple temperature, pressure conditions; note: concentration within
// speciesStatus (and list of conversion values) should not need to be modified for each T,P since this is done within
// isTPCconsistent in ReactionSystem
        initialStatusList = new LinkedList();
        for (Iterator iter = tempList.iterator(); iter.hasNext();) {
            TemperatureModel tm = (TemperatureModel) iter.next();
            for (Iterator iter2 = presList.iterator(); iter2.hasNext();) {
                PressureModel pm = (PressureModel) iter2.next();
                // LinkedHashMap speStat = (LinkedHashMap)speciesStatus.clone();//10/31/07 gmagoon: trying creating
// multiple instances of speciesStatus to address issues with concentration normalization (last normalization seems to
// apply to all)
// Set ks = speciesStatus.keySet();
// for (Iterator iter3 = ks.iterator(); iter3.hasNext();){//11/1/07 gmagoon: perform deep copy; (is there an easier or
// more elegant way to do this?)
                int reactionSystemCounter = 0;
                for (int totalConcs = 0; totalConcs < numSpeciesStatus
                        / speciesSet.size(); ++totalConcs) {
                    LinkedHashMap speStat = new LinkedHashMap();
                    for (int totalSpecs = 0; totalSpecs < speciesSet.size(); ++totalSpecs) {
                        SpeciesStatus ssCopy = (SpeciesStatus) speciesStatus
                                .get(totalConcs + totalSpecs * numSpeciesStatus
                                        / speciesSet.size());
                        String blah = ssCopy.getSpecies().getName();
                        speStat.put(ssCopy.getSpecies(), new SpeciesStatus(
                                ssCopy.getSpecies(), ssCopy.getSpeciesType(),
                                ssCopy.getConcentration(), ssCopy.getFlux()));
                    }
                    initialStatusList.add(new InitialStatus(speStat, tm
                            .getTemperature(initial), pm.getPressure(initial)));
                }
            }
        }
        return speciesSet;
    }

    public void populateInitialStatusListWithInertSpecies(BufferedReader reader) {
        String line = ChemParser.readMeaningfulLine(reader, true);
        while (!line.equals("END")) {
            StringTokenizer st = new StringTokenizer(line);
            String name = st.nextToken().trim();
            // The next token is the units of concentration
            String unit = st.nextToken();
            unit = ChemParser.removeBrace(unit);
            // The remaining tokens are concentrations
            double inertConc = 0.0;
            int counter = 0;
            int numberOfConcentrations = 0;
            while (st.hasMoreTokens()) {
                String conc = st.nextToken();
                inertConc = Double.parseDouble(conc);
                if (unit.equals("mole/l") || unit.equals("mol/l")
                        || unit.equals("mole/liter")
                        || unit.equals("mol/liter")) {
                    inertConc /= 1000;
                    unit = "mol/cm3";
                } else if (unit.equals("mole/m3") || unit.equals("mol/m3")) {
                    inertConc /= 1000000;
                    unit = "mol/cm3";
                } else if (unit.equals("molecule/cm3")
                        || unit.equals("molecules/cm3")) {
                    inertConc /= 6.022e23;
                    unit = "mol/cm3";
                } else if (!unit.equals("mole/cm3") && !unit.equals("mol/cm3")) {
                    throw new InvalidUnitException(
                            "Inert Gas Concentration not recognized: " + unit);
                }
                // SystemSnapshot.putInertGas(name,inertConc);
// for(Iterator iter=initialStatusList.iterator();iter.hasNext(); ){//6/23/09 gmagoon: needed to change this to
// accommodate non-static inertConc
// ((InitialStatus)iter.next()).putInertGas(name,inertConc);
// }
                int numberOfDiffTP = tempList.size() * presList.size();
                for (int isIndex = counter; isIndex < initialStatusList.size(); isIndex += initialStatusList
                        .size() / numberOfDiffTP) {
                    InitialStatus is = ((InitialStatus) initialStatusList
                            .get(isIndex));
                    ((InitialStatus) initialStatusList.get(isIndex))
                            .putInertGas(name, inertConc);
                }
                ++counter;
                ++numberOfConcentrations;
            }
            // Check if the number of concentrations read in is consistent with all previous
            // concentration counts. The first time this function is called, the variable
            // numberOfEquivalenceRatios will be initialized.
            boolean goodToGo = areTheNumberOfConcentrationsConsistent(numberOfConcentrations);
            if (!goodToGo) {
                Logger.critical("\n\nThe number of concentrations ("
                        + numberOfConcentrations
                        + ") supplied for species "
                        + name
                        + "\nis not consistent with the number of concentrations ("
                        + numberOfEquivalenceRatios + ") "
                        + "supplied for all previously read-in species \n\n"
                        + "Terminating RMG simulation.");
                System.exit(0);
            }
            line = ChemParser.readMeaningfulLine(reader, true);
        }
    }

    public String readMaxAtomTypes(String line, BufferedReader reader) {
        if (line.startsWith("MaxCarbonNumber")) {
            StringTokenizer st = new StringTokenizer(line);
            String dummyString = st.nextToken(); // This should hold "MaxCarbonNumberPerSpecies:"
            int maxCNum = Integer.parseInt(st.nextToken());
            ChemGraph.setMaxCarbonNumber(maxCNum);
            Logger.info("Note: Overriding default MAX_CARBON_NUM with user-defined value: "
                    + maxCNum);
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (line.startsWith("MaxOxygenNumber")) {
            StringTokenizer st = new StringTokenizer(line);
            String dummyString = st.nextToken(); // This should hold "MaxOxygenNumberPerSpecies:"
            int maxONum = Integer.parseInt(st.nextToken());
            ChemGraph.setMaxOxygenNumber(maxONum);
            Logger.info("Note: Overriding default MAX_OXYGEN_NUM with user-defined value: "
                    + maxONum);
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (line.startsWith("MaxRadicalNumber")) {
            StringTokenizer st = new StringTokenizer(line);
            String dummyString = st.nextToken(); // This should hold "MaxRadicalNumberPerSpecies:"
            int maxRadNum = Integer.parseInt(st.nextToken());
            ChemGraph.setMaxRadicalNumber(maxRadNum);
            Logger.info("Note: Overriding default MAX_RADICAL_NUM with user-defined value: "
                    + maxRadNum);
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (line.startsWith("MaxSulfurNumber")) {
            StringTokenizer st = new StringTokenizer(line);
            String dummyString = st.nextToken(); // This should hold "MaxSulfurNumberPerSpecies:"
            int maxSNum = Integer.parseInt(st.nextToken());
            ChemGraph.setMaxSulfurNumber(maxSNum);
            Logger.info("Note: Overriding default MAX_SULFUR_NUM with user-defined value: "
                    + maxSNum);
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (line.startsWith("MaxSiliconNumber")) {
            StringTokenizer st = new StringTokenizer(line);
            String dummyString = st.nextToken(); // This should hold "MaxSiliconNumberPerSpecies:"
            int maxSiNum = Integer.parseInt(st.nextToken());
            ChemGraph.setMaxSiliconNumber(maxSiNum);
            Logger.info("Note: Overriding default MAX_SILICON_NUM with user-defined value: "
                    + maxSiNum);
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (line.startsWith("MaxHeavyAtom")) {
            StringTokenizer st = new StringTokenizer(line);
            String dummyString = st.nextToken(); // This should hold "MaxHeavyAtomPerSpecies:"
            int maxHANum = Integer.parseInt(st.nextToken());
            ChemGraph.setMaxHeavyAtomNumber(maxHANum);
            Logger.info("Note: Overriding default MAX_HEAVYATOM_NUM with user-defined value: "
                    + maxHANum);
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (line.startsWith("MaxCycleNumber")) {
            StringTokenizer st = new StringTokenizer(line);
            String dummyString = st.nextToken(); // This should hold "MaxCycleNumberPerSpecies:"
            int maxCycleNum = Integer.parseInt(st.nextToken());
            ChemGraph.setMaxCycleNumber(maxCycleNum);
            Logger.info("Note: Overriding default MAX_CYCLE_NUM with user-defined value: "
                    + maxCycleNum);
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        return line;
    }

    public ReactionModelEnlarger getReactionModelEnlarger() {
        return reactionModelEnlarger;
    }

    public LinkedList getTempList() {
        return tempList;
    }

    public LinkedList getPressList() {
        return presList;
    }

    public LinkedList getInitialStatusList() {
        return initialStatusList;
    }

    public void writeBackupRestartFiles(String[] listOfFiles) {
        for (int i = 0; i < listOfFiles.length; i++) {
            // delete the backup if it exists (it shouldn't)
            File temporaryRestartFile = new File(
                    System.getProperty("RMG.RestartDir"), listOfFiles[i] + "~");
            if (temporaryRestartFile.exists())
                temporaryRestartFile.delete();
            // rename the original to the backup
            temporaryRestartFile = new File(
                    System.getProperty("RMG.RestartDir"), listOfFiles[i]);
            if (temporaryRestartFile.exists())
                temporaryRestartFile.renameTo(new File(System
                        .getProperty("RMG.RestartDir"), listOfFiles[i] + "~"));
        }
    }

    public void removeBackupRestartFiles(String[] listOfFiles) {
        for (int i = 0; i < listOfFiles.length; i++) {
            // delete the backup
            File temporaryRestartFile = new File(
                    System.getProperty("RMG.RestartDir"), listOfFiles[i] + "~");
            temporaryRestartFile.delete();
        }
    }

    public static boolean rerunFameWithAdditionalGrains() {
        return rerunFame;
    }

    public void setLimitingReactantID(int id) {
        limitingReactantID = id;
    }

    public int getLimitingReactantID() {
        return limitingReactantID;
    }

    public void readAndMakePTransL(BufferedReader reader) {
        int numPTLs = 0;
        String line = ChemParser.readMeaningfulLine(reader, true);
        while (!line.equals("END")) {
            String name = extractLibraryName(line);
            line = ChemParser.readMeaningfulLine(reader, true);
            String[] tempString = line.split("Location: ");
            String path = tempString[tempString.length - 1].trim();
            if (numPTLs == 0) {
                setPrimaryTransportLibrary(new PrimaryTransportLibrary(name,
                        path));
                ++numPTLs;
            } else {
                getPrimaryTransportLibrary().appendPrimaryTransportLibrary(
                        name, path);
                ++numPTLs;
            }
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        if (numPTLs == 0)
            setPrimaryTransportLibrary(null);
    }

    // Added by Amrit Jalan on December 21, 2010
    public void readAndMakePAL() {
        String name = "primaryAbrahamLibrary";
        String path = "primaryAbrahamLibrary";
        setPrimaryAbrahamLibrary(new PrimaryAbrahamLibrary(name, path));
        getPrimaryAbrahamLibrary().appendPrimaryAbrahamLibrary(name, path);
    }

    public void readAndMakeSL(String solventname) {
        String name = "SolventLibrary";
        String path = "SolventLibrary";
        setSolventLibrary(new SolventLibrary(name, path)); // the constructor with (name,path) reads in the library at
// construction time.
        SolventData solvent = getSolventLibrary().getSolventData(solventname);
        setSolvent(solvent);
    }

    public PrimaryTransportLibrary getPrimaryTransportLibrary() {
        return primaryTransportLibrary;
    }

    public void setPrimaryTransportLibrary(
            PrimaryTransportLibrary p_primaryTransportLibrary) {
        primaryTransportLibrary = p_primaryTransportLibrary;
    }

    public PrimaryAbrahamLibrary getPrimaryAbrahamLibrary() {
        return primaryAbrahamLibrary;
    }

    public static SolventData getSolvent() {
        return solvent;
    }

    public static double getViscosity() {
        return viscosity;
    }

    public SolventLibrary getSolventLibrary() {
        return solventLibrary;
    }

    public void setPrimaryAbrahamLibrary(
            PrimaryAbrahamLibrary p_primaryAbrahamLibrary) {
        primaryAbrahamLibrary = p_primaryAbrahamLibrary;
    }

    public void setSolventLibrary(SolventLibrary p_solventLibrary) {
        solventLibrary = p_solventLibrary;
    }

    public void setSolvent(SolventData p_solvent) {
        solvent = p_solvent;
    }

    /**
     * Print the current numbers of core and edge species and reactions to the console.
     */
    public void printModelSize() {
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) getReactionModel();
        int numberOfCoreSpecies = cerm.getReactedSpeciesSet().size();
        int numberOfEdgeSpecies = cerm.getUnreactedSpeciesSet().size();
        int numberOfCoreReactions = 0;
        int numberOfEdgeReactions = 0;
        double count = 0.0;
        for (Iterator iter = cerm.getReactedReactionSet().iterator(); iter
                .hasNext();) {
            Reaction rxn = (Reaction) iter.next();
            // The model core stores reactions in both directions
            // To avoid double-counting we must count each of these as 1/2
            if (rxn.hasReverseReaction())
                count += 0.5;
            else
                count += 1;
        }
        numberOfCoreReactions = (int) Math.round(count);
        count = 0.0;
        for (Iterator iter = cerm.getUnreactedReactionSet().iterator(); iter
                .hasNext();) {
            Reaction rxn = (Reaction) iter.next();
            // The model edge stores reactions in only one direction, so each
            // edge reaction counts as 1
            count += 1;
        }
        numberOfEdgeReactions = (int) Math.round(count);
        if (reactionModelEnlarger instanceof RateBasedPDepRME) {
            numberOfCoreReactions += PDepNetwork.getNumCoreReactions(cerm);
            numberOfEdgeReactions += PDepNetwork.getNumEdgeReactions(cerm);
        }
        Logger.info("");
        Logger.info("The model core has "
                + Integer.toString(numberOfCoreReactions) + " reactions and "
                + Integer.toString(numberOfCoreSpecies) + " species.");
        Logger.info("The model edge has "
                + Integer.toString(numberOfEdgeReactions) + " reactions and "
                + Integer.toString(numberOfEdgeSpecies) + " species.");
        // If pressure dependence is on, print some information about the networks
        if (reactionModelEnlarger instanceof RateBasedPDepRME) {
            int numberOfNetworks = PDepNetwork.getNetworks().size();
            int numberOfPathReactions = PDepNetwork.getNumPathReactions(cerm);
            int numberOfNetReactions = PDepNetwork.getNumNetReactions(cerm);
            Logger.info("There are " + Integer.toString(numberOfNetworks)
                    + " partial pressure-dependent networks containing "
                    + Integer.toString(numberOfPathReactions) + " path and "
                    + Integer.toString(numberOfNetReactions)
                    + " net reactions.");
        }
    }

    public boolean areTheNumberOfConcentrationsConsistent(int number) {
        if (numberOfEquivalenceRatios == 0)
            numberOfEquivalenceRatios = number;
        else {
            if (number == numberOfEquivalenceRatios)
                return true;
            else
                return false;
        }
        return true;
    }

    public static void addChemGraphToListIfNotPresent_ElseTerminate(
            LinkedHashMap speciesMap, ChemGraph cg, String name) {
        if (speciesMap.containsKey(cg)) {
            String message = "The same ChemGraph appears multiple times in the user-specified input file\n"
                    + String.format("Species %s has the same ChemGraph as %s",
                            name, speciesMap.get(cg));
            Logger.error(message);
            throw new RuntimeException(message);
        } else
            speciesMap.put(cg, name);
    }

    public String extractLibraryName(String line) {
        String[] tempString = line.split("Name: ");
        String name = tempString[tempString.length - 1].trim();
        if (name.contains("//")) {
            tempString = line.split("//");
            name = tempString[0].trim();
        }
        return name;
    }
}
/*********************************************************************
 * File Path : RMG\RMG\jing\rxnSys\ReactionModelGenerator.java
 *********************************************************************/

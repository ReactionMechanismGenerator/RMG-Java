////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
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


//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\ReactionModelGenerator.java
//----------------------------------------------------------------------------

//## class ReactionModelGenerator
public class ReactionModelGenerator {
	
    protected LinkedList timeStep;		//## attribute timeStep
    protected ReactionModel reactionModel;      //gmagoon 9/24/07
    protected String workingDirectory;		//## attribute workingDirectory
	
	// protected ReactionSystem reactionSystem;
    protected LinkedList reactionSystemList; //10/24/07 gmagoon: changed from reactionSystem to reactionSystemList
    
    protected int paraInfor;//svp
    protected boolean error;//svp
    protected boolean sensitivity;//svp
    protected LinkedList species;//svp
	//  protected InitialStatus initialStatus;//svp
    protected LinkedList initialStatusList; //10/23/07 gmagoon: changed from initialStatus to initialStatusList
    protected double rtol;//svp
    protected static double atol;
    protected PrimaryReactionLibrary primaryReactionLibrary;//9/24/07 gmagoon
    protected ReactionModelEnlarger reactionModelEnlarger;//9/24/07 gmagoon
    protected LinkedHashSet speciesSeed;//9/24/07 gmagoon;
    protected ReactionGenerator reactionGenerator;//9/24/07 gmagoon
    protected LibraryReactionGenerator lrg;// = new LibraryReactionGenerator();//9/24/07 gmagoon: moved from ReactionSystem.java;10/4/07 gmagoon: postponed initialization of lrg til later
    //10/23/07 gmagoon: added additional variables
    protected LinkedList tempList;
    protected LinkedList presList;
    protected LinkedList validList;//10/24/07 gmagoon: added
    //10/25/07 gmagoon: moved variables from modelGeneration()
    protected LinkedList initList = new LinkedList();
    protected LinkedList beginList = new LinkedList();
    protected LinkedList endList = new LinkedList();
    protected LinkedList lastTList = new LinkedList();
    protected LinkedList currentTList = new LinkedList();
    protected LinkedList lastPList = new LinkedList();
    protected LinkedList currentPList = new LinkedList();
    protected LinkedList conditionChangedList = new LinkedList();
    protected LinkedList reactionChangedList = new LinkedList();
    protected int numConversions;//5/6/08 gmagoon: moved from initializeReactionSystem() to be an attribute so it can be accessed by modelGenerator()
    protected String equationOfState;
    // 24Jun2009 MRH: variable stores the first temperature encountered in the condition.txt file
    //	This temperature is used to select the "best" kinetics from the rxn library
    protected static Temperature temp4BestKinetics; 
    // This is the new "PrimaryReactionLibrary"
    protected SeedMechanism seedMechanism;
    protected PrimaryThermoLibrary primaryThermoLibrary;
	
	protected boolean restart = false;
	protected boolean readrestart = false;
	protected boolean writerestart = false;
	
	protected LinkedHashSet restartCoreSpcs = new LinkedHashSet();
	protected LinkedHashSet restartEdgeSpcs = new LinkedHashSet();
	protected LinkedHashSet restartCoreRxns = new LinkedHashSet();
	protected LinkedHashSet restartEdgeRxns = new LinkedHashSet();
    // Constructors
	
	private HashSet specs = new HashSet();
	//public static native long getCpuTime();
	//static {System.loadLibrary("cpuTime");}
	
	//## operation ReactionModelGenerator()
    public  ReactionModelGenerator() {
        //#[ operation ReactionModelGenerator()
        workingDirectory = System.getProperty("RMG.workingDirectory");
		
		
        //#]
    }
	
    //## operation initializeReactionSystem()
    //10/24/07 gmagoon: changed name to initializeReactionSystems
    public void initializeReactionSystems() throws InvalidSymbolException, IOException {
        //#[ operation initializeReactionSystem()
        try {
        	String initialConditionFile = System.getProperty("jing.rxnSys.ReactionModelGenerator.conditionFile");
        	if (initialConditionFile == null) {
        		System.out.println("undefined system property: jing.rxnSys.ReactionModelGenerator.conditionFile");
        		System.exit(0);
        	}
			//double sandeep = getCpuTime();
			//System.out.println(getCpuTime()/1e9/60);
        	FileReader in = new FileReader(initialConditionFile);
        	BufferedReader reader = new BufferedReader(in);
			
        	//TemperatureModel temperatureModel = null;//10/27/07 gmagoon: commented out
        	//PressureModel pressureModel = null;//10/27/07 gmagoon: commented out
			//      	ReactionModelEnlarger reactionModelEnlarger = null;//10/9/07 gmagoon: commented out: unneeded now and causes scope problems
        	
        	FinishController finishController = null;
        	//DynamicSimulator dynamicSimulator = null;//10/27/07 gmagoon: commented out and replaced with following line
			LinkedList dynamicSimulatorList = new LinkedList();
        	//PrimaryReactionLibrary primaryReactionLibrary = null;//10/14/07 gmagoon: see below
			setPrimaryReactionLibrary(null);//10/14/07 gmagoon: changed to use setPrimaryReactionLibrary
        	double [] conversionSet = new double[50];
			String line = ChemParser.readMeaningfulLine(reader);
        	/*if (line.startsWith("Restart")){
			 StringTokenizer st = new StringTokenizer(line);
			 String token = st.nextToken();
			 token = st.nextToken();
			 if (token.equalsIgnoreCase("true")) {
			 //Runtime.getRuntime().exec("cp Restart/allSpecies.txt Restart/allSpecies1.txt");
			 //Runtime.getRuntime().exec("echo  >> allSpecies.txt");
			 restart = true;
			 }
			 else if (token.equalsIgnoreCase("false")) {
			 Runtime.getRuntime().exec("rm Restart/allSpecies.txt");
			 restart = false;
			 }
			 else throw new InvalidSymbolException("UnIdentified Symbol "+token+" after Restart:");
			 }
			 else throw new InvalidSymbolException("Can't find Restart!");*/
			
			//line = ChemParser.readMeaningfulLine(reader);
			
			if (line.startsWith("Database")){//svp
                line = ChemParser.readMeaningfulLine(reader);
			}
			else throw new InvalidSymbolException("Can't find database!");
			//              if (line.startsWith("PrimaryThermoLibrary")){//svp
			//                line = ChemParser.readMeaningfulLine(reader);
			//              }
			//              else throw new InvalidSymbolException("Can't find primary thermo library!");
			
			/*
             * Added by MRH on 15-Jun-2009
             * 	Give user the option to change the maximum carbon, oxygen,
             *		and/or radical number for all species.  These lines will be
             * 		optional in the condition.txt file.  Values are hard-
             * 		coded into RMG (in ChemGraph.java), but any user-
             * 		defined input will override these values.
             */
            /*
             * Moved from before InitialStatus to before PrimaryThermoLibary
             * by MRH on 27-Oct-2009
             * 	Overriding default values of maximum number of "X" per
             * 	chemgraph should come before RMG attempts to make any 
             *     chemgraph.  The first instance RMG will attempt to make a 
             *     chemgraph is in reading the primary thermo library.
             */
            if (line.startsWith("MaxCarbonNumber")) {
            	StringTokenizer st = new StringTokenizer(line);
            	String dummyString = st.nextToken();	// This should hold "MaxCarbonNumberPerSpecies:"
            	int maxCNum = Integer.parseInt(st.nextToken());
            	ChemGraph.setMaxCarbonNumber(maxCNum);
            	System.out.println("Note: Overriding RMG-defined MAX_CARBON_NUM with user-defined value: " + maxCNum);
            	line = ChemParser.readMeaningfulLine(reader);
            }
            if (line.startsWith("MaxOxygenNumber")) {
            	StringTokenizer st = new StringTokenizer(line);
            	String dummyString = st.nextToken();	// This should hold "MaxOxygenNumberPerSpecies:"
            	int maxONum = Integer.parseInt(st.nextToken());
            	ChemGraph.setMaxOxygenNumber(maxONum);
            	System.out.println("Note: Overriding RMG-defined MAX_OXYGEN_NUM with user-defined value: " + maxONum);
            	line = ChemParser.readMeaningfulLine(reader);
            }
            if (line.startsWith("MaxRadicalNumber")) {
            	StringTokenizer st = new StringTokenizer(line);
            	String dummyString = st.nextToken();	// This should hold "MaxRadicalNumberPerSpecies:"
            	int maxRadNum = Integer.parseInt(st.nextToken());
            	ChemGraph.setMaxRadicalNumber(maxRadNum);
            	System.out.println("Note: Overriding RMG-defined MAX_RADICAL_NUM with user-defined value: " + maxRadNum);
            	line = ChemParser.readMeaningfulLine(reader);
            }
            if (line.startsWith("MaxSulfurNumber")) {
            	StringTokenizer st = new StringTokenizer(line);
            	String dummyString = st.nextToken();	// This should hold "MaxSulfurNumberPerSpecies:"
            	int maxSNum = Integer.parseInt(st.nextToken());
            	ChemGraph.setMaxSulfurNumber(maxSNum);
            	System.out.println("Note: Overriding RMG-defined MAX_SULFUR_NUM with user-defined value: " + maxSNum);
            	line = ChemParser.readMeaningfulLine(reader);
            }
            if (line.startsWith("MaxSiliconNumber")) {
            	StringTokenizer st = new StringTokenizer(line);
            	String dummyString = st.nextToken();	// This should hold "MaxSiliconNumberPerSpecies:"
            	int maxSiNum = Integer.parseInt(st.nextToken());
            	ChemGraph.setMaxSiliconNumber(maxSiNum);
            	System.out.println("Note: Overriding RMG-defined MAX_SILICON_NUM with user-defined value: " + maxSiNum);
            	line = ChemParser.readMeaningfulLine(reader);
            }
            if (line.startsWith("MaxHeavyAtom")) {
            	StringTokenizer st = new StringTokenizer(line);
            	String dummyString = st.nextToken();	// This should hold "MaxHeavyAtomPerSpecies:"
            	int maxHANum = Integer.parseInt(st.nextToken());
            	ChemGraph.setMaxHeavyAtomNumber(maxHANum);
            	System.out.println("Note: Overriding RMG-defined MAX_HEAVYATOM_NUM with user-defined value: " + maxHANum);
            	line = ChemParser.readMeaningfulLine(reader);
            }
			
         	/*
         	 * Read in the Primary Thermo Library
         	 * MRH 7-Jul-2009
         	 */
			if (line.startsWith("PrimaryThermoLibrary:")) {
				/*
				 * MRH 27Feb2010:
				 * Changing the "read in Primary Thermo Library information" code
				 * 	into it's own method.
				 * 
				 * Other modules (e.g. PopulateReactions) will be utilizing the exact code.
				 * 	Rather than copying and pasting code into other modules, just have
				 * 	everything call this new method: readAndMakePTL
				 */
				readAndMakePTL(reader);
//             	int numPTLs = 0;
//             	line = ChemParser.readMeaningfulLine(reader);
//             	while (!line.equals("END")) {
//             		String[] tempString = line.split("Name: ");
//             		String name = tempString[tempString.length-1].trim();
//					line = ChemParser.readMeaningfulLine(reader);
//					tempString = line.split("Location: ");
//					String path = tempString[tempString.length-1].trim();
//					if (numPTLs==0) {
//                     	setPrimaryThermoLibrary(new PrimaryThermoLibrary(name,path));
//                     	++numPTLs; 	
//					}
//					else {
//                     	getPrimaryThermoLibrary().appendPrimaryThermoLibrary(name,path);
//                     	++numPTLs;
//					}
//					line = ChemParser.readMeaningfulLine(reader);
//             	}
//             	if (numPTLs == 0) setPrimaryThermoLibrary(null);
			} else throw new InvalidSymbolException("Error reading condition.txt file: "
													+ "Could not locate PrimaryThermoLibrary field");
			line = ChemParser.readMeaningfulLine(reader);
			
 			if (line.toLowerCase().startsWith("readrestart")) {
				StringTokenizer st = new StringTokenizer(line);
				String tempString = st.nextToken();	// "ReadRestart:"
				tempString = st.nextToken();
				if (tempString.toLowerCase().equals("yes")) {
					readrestart = true;
					readRestartSpecies();
				} else readrestart = false;				
				line = ChemParser.readMeaningfulLine(reader);
			} else throw new InvalidSymbolException("Cannot locate ReadRestart field");
			
			if (line.toLowerCase().startsWith("writerestart")) {
				StringTokenizer st = new StringTokenizer(line);
				String tempString = st.nextToken();	// "WriteRestart:"
				tempString = st.nextToken();
				if (tempString.toLowerCase().equals("yes"))
					writerestart = true;
				else writerestart = false;
				line = ChemParser.readMeaningfulLine(reader);
			} else throw new InvalidSymbolException("Cannot locate WriteRestart field");
			
        	// read temperature model
			//gmagoon 10/23/07: modified to handle multiple temperatures; note that this requires different formatting of units in condition.txt
        	if (line.startsWith("TemperatureModel:")) {
        		StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		String modelType = st.nextToken();
        		//String t = st.nextToken();
        		String unit = st.nextToken();
				unit = ChemParser.removeBrace(unit);
        		if (modelType.equals("Constant")) {
					tempList = new LinkedList();
					//read first temperature
					double t = Double.parseDouble(st.nextToken());
					tempList.add(new ConstantTM(t, unit));
					Temperature temp = new Temperature(t, unit);//10/29/07 gmagoon: added this line and next two lines to set Global.lowTemperature and Global.highTemperature
					Global.lowTemperature = (Temperature)temp.clone();
					Global.highTemperature = (Temperature)temp.clone();
					//read remaining temperatures
        			while (st.hasMoreTokens()) {
						t = Double.parseDouble(st.nextToken());
						tempList.add(new ConstantTM(t, unit));
						temp = new Temperature(t,unit);//10/29/07 gmagoon: added this line and next two "if" statements to set Global.lowTemperature and Global.highTemperature
						if(temp.getK() < Global.lowTemperature.getK())
							Global.lowTemperature = (Temperature)temp.clone();
						if(temp.getK() > Global.highTemperature.getK())
							Global.highTemperature = (Temperature)temp.clone();
					}
					// Global.temperature = new Temperature(t,unit);
        		}
				//10/23/07 gmagoon: commenting out; further updates needed to get this to work
        		//else if (modelType.equals("Curved")) {
				//        String t = st.nextToken();
        		//	// add reading curved temperature function here
        		//	temperatureModel = new CurvedTM(new LinkedList());
        		//}
        		else {
        			throw new InvalidSymbolException("condition.txt: Unknown TemperatureModel = " + modelType);
        		}
        	}
        	else throw new InvalidSymbolException("condition.txt: can't find TemperatureModel!");
			
        	// read in pressure model
        	line = ChemParser.readMeaningfulLine(reader);
        	if (line.startsWith("PressureModel:")) {
        		StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		String modelType = st.nextToken();
        		//String p = st.nextToken();
        		String unit = st.nextToken();
				unit = ChemParser.removeBrace(unit);
        		if (modelType.equals("Constant")) {
					presList = new LinkedList();
					//read first pressure
					double p = Double.parseDouble(st.nextToken());
					Pressure pres = new Pressure(p, unit);
					Global.lowPressure = (Pressure)pres.clone();
					Global.highPressure = (Pressure)pres.clone();
					presList.add(new ConstantPM(p, unit));
					//read remaining temperatures
        			while (st.hasMoreTokens()) {
						p = Double.parseDouble(st.nextToken());
						presList.add(new ConstantPM(p, unit));
						pres = new Pressure(p, unit);
						if(pres.getBar() < Global.lowPressure.getBar())
							Global.lowPressure = (Pressure)pres.clone();
						if(pres.getBar() > Global.lowPressure.getBar())
							Global.highPressure = (Pressure)pres.clone();
					}	
        			//Global.pressure = new Pressure(p, unit);
        		}
				//10/23/07 gmagoon: commenting out; further updates needed to get this to work
        		//else if (modelType.equals("Curved")) {
        		//	// add reading curved pressure function here
        		//	pressureModel = new CurvedPM(new LinkedList());
        		//}
        		else {
        			throw new InvalidSymbolException("condition.txt: Unknown PressureModel = " + modelType);
        		}
        	}
        	else throw new InvalidSymbolException("condition.txt: can't find PressureModel!");
            
            // after PressureModel comes an optional line EquationOfState
            // if "EquationOfState: Liquid" is found then initial concentrations are assumed to be correct
            // if it is ommited, then initial concentrations are normalised to ensure PV=NRT (ideal gas law)
            line = ChemParser.readMeaningfulLine(reader);
            if (line.startsWith("EquationOfState")) {
                StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		String eosType = st.nextToken();
                if (eosType.equals("Liquid")) {
                    equationOfState="Liquid";
                    System.out.println("Equation of state: Liquid. Relying on concentrations in input file to get density correct; not checking PV=NRT");
                }
                line = ChemParser.readMeaningfulLine(reader);
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
        		}
        		else throw new InvalidSymbolException("condition.txt: Unknown InChIGeneration flag: " + inchiOnOff);
        		line = ChemParser.readMeaningfulLine(reader);
        	}
			
            // Read in Solvation effects
            if (line.startsWith("Solvation:")) {
        		StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		String solvationOnOff = st.nextToken().toLowerCase();
        		if (solvationOnOff.equals("on")) {
        			Species.useSolvation = true;
        		} else if (solvationOnOff.equals("off")) {
        			Species.useSolvation = false;
        		}
        		else throw new InvalidSymbolException("condition.txt: Unknown solvation flag: " + solvationOnOff);
        		line = ChemParser.readMeaningfulLine(reader);
        	}
			
			//line = ChemParser.readMeaningfulLine(reader);//read in reactants or thermo line
			// Read in optional QM thermo  generation
        	if (line.startsWith("ThermoMethod:")) {
        		StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		String thermoMethod = st.nextToken().toLowerCase();
        		if (thermoMethod.equals("qm")) {
        			ChemGraph.useQM = true;
					line=ChemParser.readMeaningfulLine(reader);
					if(line.startsWith("QMForCyclicsOnly:")){
						StringTokenizer st2 = new StringTokenizer(line);
						String nameCyc = st2.nextToken();
						String option = st2.nextToken().toLowerCase();
						if (option.equals("on")) {
							ChemGraph.useQMonCyclicsOnly = true;
						}
					}
					else{
						System.out.println("condition.txt: Can't find 'QMForCyclicsOnly:' field");
						System.exit(0);
					}
					line=ChemParser.readMeaningfulLine(reader);
					if(line.startsWith("MaxRadNumForQM:")){
						StringTokenizer st3 = new StringTokenizer(line);
						String nameRadNum = st3.nextToken();
						Global.maxRadNumForQM = Integer.parseInt(st3.nextToken());
						
					}
					else{
						System.out.println("condition.txt: Can't find 'MaxRadNumForQM:' field");
						System.exit(0);
					}
        		}//otherwise, the flag useQM will remain false by default and the traditional group additivity approach will be used
				line = ChemParser.readMeaningfulLine(reader);//read in reactants
			}
            
			//            // Read in Solvation effects
			//            if (line.startsWith("Solvation:")) {
			//        		StringTokenizer st = new StringTokenizer(line);
			//        		String name = st.nextToken();
			//        		String solvationOnOff = st.nextToken().toLowerCase();
			//        		if (solvationOnOff.equals("on")) {
			//        			Species.useSolvation = true;
			//        		} else if (solvationOnOff.equals("off")) {
			//        			Species.useSolvation = false;
			//        		}
			//        		else throw new InvalidSymbolException("condition.txt: Unknown solvation flag: " + solvationOnOff);
			//        	}
			//        	else throw new InvalidSymbolException("condition.txt: Cannot find solvation flag.");
			
			
			// read in reactants
        	//
			
			//10/4/07 gmagoon: moved to initializeCoreEdgeReactionModel
        	//LinkedHashSet p_speciesSeed = new LinkedHashSet();//gmagoon 10/4/07: changed to p_speciesSeed
			//setSpeciesSeed(p_speciesSeed);//gmagoon 10/4/07: added
        	LinkedHashMap speciesSet = new LinkedHashMap();
        	LinkedHashMap speciesStatus = new LinkedHashMap();
			int speciesnum = 1;
			//System.out.println(line);
        	if (line.startsWith("InitialStatus")) {
        		line = ChemParser.readMeaningfulLine(reader);
        		while (!line.equals("END")) {
        			StringTokenizer st = new StringTokenizer(line);
        			String index = st.nextToken();
        			String name = null;
        			if (!index.startsWith("(")) name = index;
        			else name = st.nextToken();
					//if (restart) name += "("+speciesnum+")";
        			// 24Jun2009: MRH
        			//	Check if the species name begins with a number.
        			//	If so, terminate the program and inform the user to choose
        			//		a different name.  This is implemented so that the chem.inp
        			//		file generated will be valid when run in Chemkin
        			try {
        				int doesNameBeginWithNumber = Integer.parseInt(name.substring(0,1));
        				System.out.println("\nA species name should not begin with a number." +
										   " Please rename species: " + name + "\n");
        				System.exit(0);
        			} catch (NumberFormatException e) {
        				// We're good
        			}
					speciesnum ++;
					if (!(st.hasMoreTokens())) throw new InvalidSymbolException("Couldn't find concentration of species: "+name);
        			String conc = st.nextToken();
        			double concentration = Double.parseDouble(conc);
        			String unit = st.nextToken();
        			unit = ChemParser.removeBrace(unit);
        			if (unit.equals("mole/l") || unit.equals("mol/l") || unit.equals("mole/liter") || unit.equals("mol/liter")) {
        				concentration /= 1000;
        				unit = "mol/cm3";
        			}
        			else if (unit.equals("mole/m3") || unit.equals("mol/m3")) {
        				concentration /= 1000000;
        				unit = "mol/cm3";
        			}
        			else if (unit.equals("molecule/cm3") || unit.equals("molecules/cm3")) {
        				concentration /= 6.022e23;
        			}
        			else if (!unit.equals("mole/cm3") && !unit.equals("mol/cm3")) {
        				throw new InvalidUnitException("Species Concentration in condition.txt!");
        			}
					
        			//GJB to allow "unreactive" species that only follow user-defined library reactions.  
        			// They will not react according to RMG reaction families 
					boolean IsReactive = true;
                    boolean IsConstantConcentration = false;
					while (st.hasMoreTokens()) {
						String reactive = st.nextToken().trim();
						if (reactive.equalsIgnoreCase("unreactive"))
							IsReactive = false;
                        if (reactive.equalsIgnoreCase("constantconcentration"))
                            IsConstantConcentration=true;
					}
        			
        			Graph g = ChemParser.readChemGraph(reader);
        			ChemGraph cg = null;
        			try {
        				cg = ChemGraph.make(g);
        			}
        			catch (ForbiddenStructureException e) {
        				System.out.println("Forbidden Structure:\n" + e.getMessage());
        				System.exit(0);
        			}
					//System.out.println(name);
        			Species species = Species.make(name,cg);
        			species.setReactivity(IsReactive); // GJB
                    species.setConstantConcentration(IsConstantConcentration);
           			speciesSet.put(name, species);
        			getSpeciesSeed().add(species);
        			double flux = 0;
        			int species_type = 1; // reacted species
        			SpeciesStatus ss = new SpeciesStatus(species,species_type,concentration,flux);
        			speciesStatus.put(species, ss);
        			line = ChemParser.readMeaningfulLine(reader);
        		}
				ReactionTime initial = new ReactionTime(0,"S");
				//10/23/07 gmagoon: modified for handling multiple temperature, pressure conditions; note: concentration within speciesStatus (and list of conversion values) should not need to be modified for each T,P since this is done within isTPCconsistent in ReactionSystem
				initialStatusList = new LinkedList();
				for (Iterator iter = tempList.iterator(); iter.hasNext(); ) {
					TemperatureModel tm = (TemperatureModel)iter.next();
					for (Iterator iter2 = presList.iterator(); iter2.hasNext(); ){
						PressureModel pm = (PressureModel)iter2.next();
						//   LinkedHashMap speStat = (LinkedHashMap)speciesStatus.clone();//10/31/07 gmagoon: trying creating multiple instances of speciesStatus to address issues with concentration normalization (last normalization seems to apply to all)
						Set ks = speciesStatus.keySet();
						LinkedHashMap speStat = new LinkedHashMap();
						for (Iterator iter3 = ks.iterator(); iter3.hasNext();){//11/1/07 gmagoon: perform deep copy; (is there an easier or more elegant way to do this?)
							SpeciesStatus ssCopy = (SpeciesStatus)speciesStatus.get(iter3.next());
							speStat.put(ssCopy.getSpecies(),new SpeciesStatus(ssCopy.getSpecies(),ssCopy.getSpeciesType(),ssCopy.getConcentration(),ssCopy.getFlux()));
						}
						initialStatusList.add(new InitialStatus(speStat,tm.getTemperature(initial),pm.getPressure(initial)));
					}
				}
        	}
        	else throw new InvalidSymbolException("condition.txt: can't find InitialStatus!");
			
        	// read in inert gas concentration
        	line = ChemParser.readMeaningfulLine(reader);
            if (line.startsWith("InertGas:")) {
           		line = ChemParser.readMeaningfulLine(reader);
           		while (!line.equals("END")) {
        	    	StringTokenizer st = new StringTokenizer(line);
        	    	String name = st.nextToken().trim();
        			String conc = st.nextToken();
        			double inertConc = Double.parseDouble(conc);
        			String unit = st.nextToken();
        			unit = ChemParser.removeBrace(unit);
        			if (unit.equals("mole/l") || unit.equals("mol/l") || unit.equals("mole/liter") || unit.equals("mol/liter")) {
        				inertConc /= 1000;
        				unit = "mol/cm3";
        			}
        			else if (unit.equals("mole/m3") || unit.equals("mol/m3")) {
        				inertConc /= 1000000;
        				unit = "mol/cm3";
        			}
        			else if (unit.equals("molecule/cm3") || unit.equals("molecules/cm3")) {
        				inertConc /= 6.022e23;
        				unit = "mol/cm3";
        			}
        			else if (!unit.equals("mole/cm3") && !unit.equals("mol/cm3")) {
        				throw new InvalidUnitException("Inert Gas Concentration not recognized: " + unit);
        			}
					
        			//SystemSnapshot.putInertGas(name,inertConc);
					for(Iterator iter=initialStatusList.iterator();iter.hasNext(); ){//6/23/09 gmagoon: needed to change this to accommodate non-static inertConc
						((InitialStatus)iter.next()).putInertGas(name,inertConc);
					}
        	   		line = ChemParser.readMeaningfulLine(reader);
        		}
           	}
        	else throw new InvalidSymbolException("condition.txt: can't find Inert gas concentration!");
			
        	// read in spectroscopic data estimator
			line = ChemParser.readMeaningfulLine(reader);
        	if (line.startsWith("SpectroscopicDataEstimator:")) {
        		setSpectroscopicDataMode(line);
//        		StringTokenizer st = new StringTokenizer(line);
//        		String name = st.nextToken();
//        		String sdeType = st.nextToken().toLowerCase();
//        		if (sdeType.equals("frequencygroups") || sdeType.equals("default")) {
//        			SpectroscopicData.mode = SpectroscopicData.Mode.FREQUENCYGROUPS;
//        		}
//				else if (sdeType.equals("therfit") || sdeType.equals("threefrequencymodel")) {
//        			SpectroscopicData.mode = SpectroscopicData.Mode.THREEFREQUENCY;
//        		}
//				else if (sdeType.equals("off") || sdeType.equals("none")) {
//        			SpectroscopicData.mode = SpectroscopicData.Mode.OFF;
//        		}
//				else throw new InvalidSymbolException("condition.txt: Unknown SpectroscopicDataEstimator = " + sdeType);
        	}
        	else throw new InvalidSymbolException("condition.txt: can't find SpectroscopicDataEstimator!");
			
        	// pressure dependence flag
			line = ChemParser.readMeaningfulLine(reader);
        	if (line.toLowerCase().startsWith("pressuredependence:")) {
        		setPressureDependenceType(line);
//				StringTokenizer st = new StringTokenizer(line);
//        		String name = st.nextToken();
//        		String pDepType = st.nextToken();
//        		if (pDepType.toLowerCase().equals("modifiedstrongcollision") ||
//					pDepType.toLowerCase().equals("reservoirstate") ||
//					pDepType.toLowerCase().equals("chemdis")) {
//					
//					reactionModelEnlarger = new RateBasedPDepRME();
//        			PDepNetwork.generateNetworks = true;
//					
//					if (pDepType.toLowerCase().equals("reservoirstate")) {
//						((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new FastMasterEqn(FastMasterEqn.Mode.RESERVOIRSTATE));
//						if (SpectroscopicData.mode == SpectroscopicData.Mode.OFF) {
//							System.out.println("Warning: Spectroscopic data needed for pressure dependence; switching SpectroscopicDataEstimator to FrequencyGroups.");
//							SpectroscopicData.mode = SpectroscopicData.Mode.FREQUENCYGROUPS;
//						}
//					}
//					else if (pDepType.toLowerCase().equals("modifiedstrongcollision")) {
//						((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new FastMasterEqn(FastMasterEqn.Mode.STRONGCOLLISION));
//						if (SpectroscopicData.mode == SpectroscopicData.Mode.OFF) {
//							System.out.println("Warning: Spectroscopic data needed for pressure dependence; switching SpectroscopicDataEstimator to FrequencyGroups.");
//							SpectroscopicData.mode = SpectroscopicData.Mode.FREQUENCYGROUPS;
//						}
//					}
//					else if (pDepType.toLowerCase().equals("chemdis")) {
//						((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new Chemdis());
//						if (SpectroscopicData.mode != SpectroscopicData.Mode.THREEFREQUENCY) {
//							System.out.println("Warning: Switching SpectroscopicDataEstimator to three-frequency model.");
//							SpectroscopicData.mode = SpectroscopicData.Mode.THREEFREQUENCY;
//						}
//					}
//					else {
//						throw new InvalidSymbolException("condition.txt: Unknown PDepKineticsEstimator = " + pDepType);
//					}
//					
//					// Set temperatures and pressures to use in PDep kinetics estimation
//					Temperature[] temperatures = new Temperature[8];
//					temperatures[0] = new Temperature(300, "K");
//					temperatures[1] = new Temperature(400, "K");
//					temperatures[2] = new Temperature(600, "K");
//					temperatures[3] = new Temperature(900, "K");
//					temperatures[4] = new Temperature(1200, "K");
//					temperatures[5] = new Temperature(1500, "K");
//					temperatures[6] = new Temperature(1800, "K");
//					temperatures[7] = new Temperature(2100, "K");
//					FastMasterEqn.setTemperatures(temperatures);
//					PDepRateConstant.setTemperatures(temperatures);
//					ChebyshevPolynomials.setTlow(temperatures[0]);
//					ChebyshevPolynomials.setTup(temperatures[7]);
//					
//					Pressure[] pressures = new Pressure[5];
//					pressures[0] = new Pressure(0.01, "bar");
//					pressures[1] = new Pressure(0.1, "bar");
//					pressures[2] = new Pressure(1, "bar");
//					pressures[3] = new Pressure(10, "bar");
//					pressures[4] = new Pressure(100, "bar");
//					FastMasterEqn.setPressures(pressures);
//					PDepRateConstant.setPressures(pressures);
//					ChebyshevPolynomials.setPlow(pressures[0]);
//					ChebyshevPolynomials.setPup(pressures[4]);
//					
//				}
//				else if (pDepType.toLowerCase().equals("off")) {
//					// No pressure dependence
//					reactionModelEnlarger = new RateBasedRME();
//        			PDepNetwork.generateNetworks = false;
//				}
//				else {
//					throw new InvalidSymbolException("condition.txt: Unknown PressureDependence = " + pDepType);
//				}
			}
			else throw new InvalidSymbolException("condition.txt: can't find PressureDependence flag!");
			
			// pressure dependence flag
			if (reactionModelEnlarger instanceof RateBasedPDepRME) {
				line = ChemParser.readMeaningfulLine(reader);
				if (line.toLowerCase().startsWith("pdepkineticsmodel:")) {
					line = setPDepKineticsModel(line,reader);
//					StringTokenizer st = new StringTokenizer(line);
//					String name = st.nextToken();
//					String pDepKinType = st.nextToken();
//					if (pDepKinType.toLowerCase().equals("chebyshev") ||
//						pDepKinType.toLowerCase().equals("pdeparrhenius") ||
//						pDepKinType.toLowerCase().equals("rate")) {
//						
//						if (pDepKinType.toLowerCase().equals("chebyshev")) {
//							PDepRateConstant.setMode(PDepRateConstant.Mode.CHEBYSHEV);
//							line = ChemParser.readMeaningfulLine(reader);
//							if (line.startsWith("TRange")) {
//								st = new StringTokenizer(line);
//								String temp = st.nextToken(); // Should be "TRange:"
//								String TUNITS = ChemParser.removeBrace(st.nextToken());
//								double tLow = Double.parseDouble(st.nextToken());
//								Temperature TMIN = new Temperature(tLow,TUNITS);
//								ChebyshevPolynomials.setTlow(TMIN);
//								double tHigh = Double.parseDouble(st.nextToken());
//								if (tHigh <= tLow) {
//									System.err.println("Tmax is less than or equal to Tmin");
//									System.exit(0);
//								}
//								Temperature TMAX = new Temperature(tHigh,TUNITS);
//								ChebyshevPolynomials.setTup(TMAX);
//								int tResolution = Integer.parseInt(st.nextToken());
//								int tbasisFuncs = Integer.parseInt(st.nextToken());
//								if (tbasisFuncs > tResolution) {
//									System.err.println("The number of basis functions cannot exceed the number of grid points");
//									System.exit(0);
//								}
//								FastMasterEqn.setNumTBasisFuncs(tbasisFuncs);
//								line = ChemParser.readMeaningfulLine(reader);
//								String PUNITS = "";
//								Pressure PMIN = new Pressure();
//								Pressure PMAX = new Pressure();
//								int pResolution = 0;
//								int pbasisFuncs = 0;
//								if (line.startsWith("PRange")) {
//									st = new StringTokenizer(line);
//									temp = st.nextToken(); // Should be "PRange:"
//									PUNITS = ChemParser.removeBrace(st.nextToken());
//									double pLow = Double.parseDouble(st.nextToken());
//									PMIN = new Pressure(pLow,PUNITS);
//									ChebyshevPolynomials.setPlow(PMIN);
//									double pHigh = Double.parseDouble(st.nextToken());
//									if (pHigh <= pLow) {
//										System.err.println("Pmax is less than or equal to Pmin");
//										System.exit(0);
//									}
//									PMAX = new Pressure(pHigh,PUNITS);
//									ChebyshevPolynomials.setPup(PMAX);
//									pResolution = Integer.parseInt(st.nextToken());
//									pbasisFuncs = Integer.parseInt(st.nextToken());
//									if (pbasisFuncs > pResolution) {
//										System.err.println("The number of basis functions cannot exceed the number of grid points");
//										System.exit(0);
//									}
//									FastMasterEqn.setNumPBasisFuncs(pbasisFuncs);
//								}
//								else {
//									System.err.println("RMG cannot locate PRange field for Chebyshev polynomials.");
//									System.exit(0);
//								}
//								
//								// Set temperatures and pressures to use in PDep kinetics estimation
//								Temperature[] temperatures = new Temperature[tResolution];
//								for (int i=0; i<temperatures.length; i++) {
//									double tempValueTilda = Math.cos((2*(i+1)-1)*Math.PI/(2*temperatures.length));
//									double tempValue = 2 / (tempValueTilda * ((1/TMAX.getK()) - (1/TMIN.getK())) + (1/TMIN.getK()) + (1/TMAX.getK()));
//									temperatures[temperatures.length-i-1] = new Temperature(tempValue,TUNITS);
//								}
//								FastMasterEqn.setTemperatures(temperatures);
//								PDepRateConstant.setTemperatures(temperatures);
//								
//								Pressure[] pressures = new Pressure[pResolution];
//								for (int j=0; j<pressures.length; j++) {
//									double pressValueTilda = Math.cos((2*(j+1)-1)*Math.PI/(2*pressures.length));
//									double pressValue = Math.pow(10,(pressValueTilda*(Math.log10(PMAX.getBar())-Math.log10(PMIN.getBar()))+Math.log10(PMIN.getBar())+Math.log10(PMAX.getBar()))/2);
//									pressures[pressures.length-j-1] = new Pressure(pressValue,PUNITS);
//								}
//								FastMasterEqn.setPressures(pressures);
//								PDepRateConstant.setPressures(pressures);
//								line = ChemParser.readMeaningfulLine(reader);
//							}
//						}
//						else if (pDepKinType.toLowerCase().equals("pdeparrhenius")) {
//							//PDepRateConstant.setMode(PDepRateConstant.Mode.PDEPARRHENIUS);
//							/*
//							 *  Updated by MRH on 10Feb2010
//							 *  	Allow user to specify # of T's/P's solved for in fame &
//							 *  	# of PLOG's to report
//							 */
//							PDepRateConstant.setMode(PDepRateConstant.Mode.PDEPARRHENIUS);
//							line = ChemParser.readMeaningfulLine(reader);
//							if (line.startsWith("TRange")) {
//								st = new StringTokenizer(line);
//								String temp = st.nextToken(); // Should be "TRange:"
//								String TUNITS = ChemParser.removeBrace(st.nextToken());
//								int numT = Integer.parseInt(st.nextToken());
//								Temperature[] listOfTs = new Temperature[numT];
//								int counter = 0;
//								while (st.hasMoreTokens()) {
//									listOfTs[counter] = new Temperature(Double.parseDouble(st.nextToken()),TUNITS);
//									++counter;
//								}
//								if (counter != numT) {
//									System.out.println("Warning in TRange field of PressureDependence:\n" +
//											"The stated number of temperatures is: " + numT + 
//											"but the length of the temperature list is: " + counter);
//								}
//								
//								line = ChemParser.readMeaningfulLine(reader);
//								String PUNITS = "";
//								if (line.startsWith("PRange")) {
//									st = new StringTokenizer(line);
//									temp = st.nextToken(); // Should be "PRange:"
//									PUNITS = ChemParser.removeBrace(st.nextToken());
//									int numP = Integer.parseInt(st.nextToken());
//									Pressure[] listOfPs = new Pressure[numP];
//									counter = 0;
//									while (st.hasMoreTokens()) {
//										listOfPs[counter] = new Pressure(Double.parseDouble(st.nextToken()),PUNITS);
//										++counter;
//									}
//									if (counter != numP) {
//										System.out.println("Warning in PRange field of PressureDependence:\n" +
//												"The stated number of pressures is: " + numT + 
//												"but the length of the pressure list is: " + counter);
//									}
//									FastMasterEqn.setTemperatures(listOfTs);
//									PDepRateConstant.setTemperatures(listOfTs);
//									FastMasterEqn.setPressures(listOfPs);
//									PDepRateConstant.setPressures(listOfPs);
//									PDepArrheniusKinetics.setNumPressures(numP);
//									PDepArrheniusKinetics.setPressures(listOfPs);
//								}
//								else {
//									System.err.println("RMG cannot locate PRange field for PDepArrhenius.");
//									System.exit(0);
//								}
//								
//								line = ChemParser.readMeaningfulLine(reader);
//							}
//						}
//						// 6Jul2009-MRH:
//						//	RATE mode reports p-dep rxn kinetics as: A 0.0 0.0
//						//		where A is k(T,P) evaluated at the single temperature
//						//		and pressure given in the condition.txt file
//						else if (pDepKinType.toLowerCase().equals("rate"))
//							PDepRateConstant.setMode(PDepRateConstant.Mode.RATE);
//						
//					}
//					else {
//						throw new InvalidSymbolException("condition.txt: Unknown PDepKinetics = " + pDepKinType);
//					}
				}
				else throw new InvalidSymbolException("condition.txt: can't find PDepKinetics flag!");
			}
			
        	// include species (optional)
			if (!PDepRateConstant.getMode().name().equals("CHEBYSHEV") &&
					!PDepRateConstant.getMode().name().equals("PDEPARRHENIUS"))
				line = ChemParser.readMeaningfulLine(reader);
			if (line.startsWith("IncludeSpecies")) {
				StringTokenizer st = new StringTokenizer(line);
				String iS = st.nextToken();
				String fileName = st.nextToken();
				HashSet includeSpecies = readIncludeSpecies(fileName);
				((RateBasedRME)reactionModelEnlarger).addIncludeSpecies(includeSpecies);
				line = ChemParser.readMeaningfulLine(reader);
			}
			
			// read in finish controller
        	if (line.startsWith("FinishController")) {
        		line = ChemParser.readMeaningfulLine(reader);
        		StringTokenizer st = new StringTokenizer(line);
        		String index = st.nextToken();
        		String goal = st.nextToken();
        		String type = st.nextToken();
        		TerminationTester tt;
        		if (type.startsWith("Conversion")) {
        			LinkedList spc = new LinkedList();
        			while (st.hasMoreTokens()) {
        				String name = st.nextToken();
        				Species spe = (Species)speciesSet.get(name);
        				if (spe == null) throw new InvalidConversionException("Unknown reactant: " + name);
        				String conv = st.nextToken();
        				double conversion;
        				try {
        					if (conv.endsWith("%")) {
        						conversion = Double.parseDouble(conv.substring(0,conv.length()-1))/100;
        					}
        					else {
        						conversion = Double.parseDouble(conv);
        					}
        					conversionSet[49] = conversion;
        				}
        				catch (NumberFormatException e) {
        					throw new NumberFormatException("wrong number format for conversion in initial condition file!");
        				}
        				SpeciesConversion sc = new SpeciesConversion(spe, conversion);
        				spc.add(sc);
        			}
        			tt = new ConversionTT(spc);
        		}
        		else if (type.startsWith("ReactionTime")) {
        			double time = Double.parseDouble(st.nextToken());
        			String unit = ChemParser.removeBrace(st.nextToken());
        			ReactionTime rt = new ReactionTime(time, unit);
        			tt = new ReactionTimeTT(rt);
        		}
        		else {
        			throw new InvalidSymbolException("condition.txt: Unknown FinishController = " + type);
        		}
				
        		line = ChemParser.readMeaningfulLine(reader);
        		st = new StringTokenizer(line, ":");
        		String temp = st.nextToken();
        		String tol = st.nextToken();
        		double tolerance;
        		try {
        			if (tol.endsWith("%")) {
        				tolerance = Double.parseDouble(tol.substring(0,tol.length()-1))/100;
        			}
        			else {
        				tolerance = Double.parseDouble(tol);
        			}
        		}
        		catch (NumberFormatException e) {
        			throw new NumberFormatException("wrong number format for conversion in initial condition file!");
        		}
        		ValidityTester vt = null;
        		if (reactionModelEnlarger instanceof RateBasedRME) vt = new RateBasedVT(tolerance);
        		else if (reactionModelEnlarger instanceof RateBasedPDepRME) vt = new RateBasedPDepVT(tolerance);
        		else throw new InvalidReactionModelEnlargerException();
        		finishController = new FinishController(tt, vt);
        	}
        	else throw new InvalidSymbolException("condition.txt: can't find FinishController!");
			
        	// read in dynamic simulator
        	line = ChemParser.readMeaningfulLine(reader);
        	if (line.startsWith("DynamicSimulator")) {
        		StringTokenizer st = new StringTokenizer(line,":");
        		String temp = st.nextToken();
        		String simulator = st.nextToken().trim();
                        
                        //read in non-negative option if it exists: syntax would be something like this: "DynamicSimulator: DASSL: non-negative"
                        if (st.hasMoreTokens()){
                            if (st.nextToken().trim().toLowerCase().equals("non-negative")){
                                if(simulator.toLowerCase().equals("dassl")) JDAS.nonnegative = true;
                                else{
                                    System.err.println("Non-negative option is currently only supported for DASSL. Switch to DASSL solver or remove non-negative option.");
                                    System.exit(0);
                                }
                            }
                        }
        		
				numConversions = 0;//5/6/08 gmagoon: moved declaration from initializeReactionSystem() to be an attribute so it can be accessed by modelGenerator()
        		//int numConversions = 0;
				boolean autoflag = false;//5/2/08 gmagoon: updating the following if/else-if block to consider input where we want to check model validity within the ODE solver at each time step; this will be indicated by the use of a string beginning with "AUTO" after the "TimeStep" or "Conversions" line
        		// read in time step
        		line = ChemParser.readMeaningfulLine(reader);
        		if (line.startsWith("TimeStep:") && finishController.terminationTester instanceof ReactionTimeTT) {
        			st = new StringTokenizer(line);
        			temp = st.nextToken();
        			while (st.hasMoreTokens()) {
                        temp = st.nextToken();
                        if (temp.startsWith("AUTO")){//note potential opportunity for making case insensitive by using: temp.toUpperCase().startsWith("AUTO")
                            autoflag=true;
                        }
        				else if (!autoflag){//use "else if" to make sure additional numbers are not read in case numbers are erroneously used following AUTO; note that there could still be a problem if numbers come before "AUTO"
                            double tStep = Double.parseDouble(temp);
                            String unit = "sec";
                            setTimeStep(new ReactionTime(tStep, unit));
                        }
        			}      
        			((ReactionTimeTT)finishController.terminationTester).setTimeSteps(timeStep);
        		}
        		else if (line.startsWith("Conversions:") && finishController.terminationTester instanceof ConversionTT){
        			st = new StringTokenizer(line);
        			temp = st.nextToken();
        			int i=0;
        			SpeciesConversion sc = (SpeciesConversion)((ConversionTT)finishController.terminationTester).speciesGoalConversionSet.get(0);
        			Species convSpecies = sc.species;
        			Iterator iter = ((InitialStatus)(initialStatusList.get(0))).getSpeciesStatus();//10/23/07 gmagoon: changed to use first element of initialStatusList, as subsequent operations should not be affected by which one is chosen
        			double initialConc = 0;
        			while (iter.hasNext()){
        				SpeciesStatus sps = (SpeciesStatus)iter.next();
        				if (sps.species.equals(convSpecies)) initialConc = sps.concentration;
        			}
        			while (st.hasMoreTokens()){
						temp=st.nextToken();
						if (temp.startsWith("AUTO")){
							autoflag=true;
						}
						else if (!autoflag){
        					double conv = Double.parseDouble(temp);
							conversionSet[i] = (1-conv) * initialConc;
							i++;
						}
        			}
        			conversionSet[i] = (1 - conversionSet[49])* initialConc;
        			numConversions = i+1;
        		}
        		else throw new InvalidSymbolException("condition.txt: can't find time step for dynamic simulator!");
				
        		// read in atol
        		
        		line = ChemParser.readMeaningfulLine(reader);
        		if (line.startsWith("Atol:")) {
        			st = new StringTokenizer(line);
        			temp = st.nextToken();
        			atol = Double.parseDouble(st.nextToken());
        		}
        		else throw new InvalidSymbolException("condition.txt: can't find Atol for dynamic simulator!");
				
        		// read in rtol
        		
        		line = ChemParser.readMeaningfulLine(reader);
        		if (line.startsWith("Rtol:")) {
        			st = new StringTokenizer(line);
        			temp = st.nextToken();
        			String rel_tol = st.nextToken();
        			if (rel_tol.endsWith("%"))
        				rtol = Double.parseDouble(rel_tol.substring(0,rel_tol.length()-1));
        			else
        				rtol = Double.parseDouble(rel_tol);
        		}
        		else throw new InvalidSymbolException("condition.txt: can't find Rtol for dynamic simulator!");
				
				if (simulator.equals("DASPK")) {
					paraInfor = 0;//svp
					// read in SA
					line = ChemParser.readMeaningfulLine(reader);
					if (line.startsWith("Error bars")) {//svp
						st = new StringTokenizer(line,":");
						temp = st.nextToken();
						String sa = st.nextToken().trim();
						if (sa.compareToIgnoreCase("on")==0) {
							paraInfor = 1;
							error = true;
						}
						else if (sa.compareToIgnoreCase("off")==0) {
							paraInfor = 0;
							error = false;
						}
						else throw new InvalidSymbolException("condition.txt: can't find error on/off information!");
						
						
					}
					
					else throw new InvalidSymbolException("condition.txt: can't find SA information!");
					line = ChemParser.readMeaningfulLine(reader);
					if (line.startsWith("Display sensitivity coefficients")){//svp
						st = new StringTokenizer(line,":");
						temp = st.nextToken();
						String sa = st.nextToken().trim();
						if (sa.compareToIgnoreCase("on")==0){
							paraInfor = 1;
							sensitivity = true;
						}
						else if (sa.compareToIgnoreCase("off")==0){
							if (paraInfor != 1){
								paraInfor = 0;
							}
							sensitivity = false;
						}
						else throw new InvalidSymbolException("condition.txt: can't find SA on/off information!");
						//10/23/07 gmagoon: changed below from dynamicSimulator to dynamicSimulatorList
						//6/25/08 gmagoon: changed loop to use i index, and updated DASPK constructor to pass i (mirroring changes to DASSL
						//6/25/08 gmagoon: updated to pass autoflag and validity tester; this requires FinishController block of input file to be present before DynamicSimulator block, but this requirement may have already existed anyway, particularly in construction of conversion/time step lists; *perhaps we should formalize this requirement by checking to make sure validityTester is not null?
						for (int i = 0;i < initialStatusList.size();i++) {
							dynamicSimulatorList.add(new JDASPK(rtol, atol, 0, (InitialStatus)initialStatusList.get(i), i,finishController.getValidityTester(), autoflag));
						}
					}
					species = new LinkedList();
					line = ChemParser.readMeaningfulLine(reader);
					
					if (line.startsWith("Display sensitivity information") ){
						line = ChemParser.readMeaningfulLine(reader);
						System.out.println(line);
						while (!line.equals("END")){
							st = new StringTokenizer(line);
							String name = st.nextToken();
							if (name.toUpperCase().equals("ALL")) ReactionSystem.printAllSens = true; //gmagoon 12/22/09: if the line contains the word "all", turn on the flag to print out sensitivity information for everything
							species.add(name);
							line = ChemParser.readMeaningfulLine(reader);
						}
					}
					
				}
				
        		else if (simulator.equals("DASSL")) {
        			//10/23/07 gmagoon: changed below from dynamicSimulator to dynamicSimulatorList
					// for (Iterator iter = initialStatusList.iterator(); iter.hasNext(); ) {
					//      dynamicSimulatorList.add(new JDASSL(rtol, atol, 0, (InitialStatus)iter.next()));
					// }
					//11/1/07 gmagoon: changed loop to use i index, and updated DASSL constructor to pass i
					//5/5/08 gmagoon: updated to pass autoflag and validity tester; this requires FinishController block of input file to be present before DynamicSimulator block, but this requirement may have already existed anyway, particularly in construction of conversion/time step lists; *perhaps we should formalize this requirement by checking to make sure validityTester is not null?
					for (int i = 0;i < initialStatusList.size();i++) {
						dynamicSimulatorList.add(new JDASSL(rtol, atol, 0, (InitialStatus)initialStatusList.get(i), i, finishController.getValidityTester(), autoflag));
					}
        		}
        		else if (simulator.equals("Chemkin")) {
        			line = ChemParser.readMeaningfulLine(reader);
        			if (line.startsWith("ReactorType")) {
        				st = new StringTokenizer(line, ":");
        				temp = st.nextToken();
        				String reactorType = st.nextToken().trim();
        				//10/23/07 gmagoon: changed below from dynamicSimulator to dynamicSimulatorList
						for (Iterator iter = initialStatusList.iterator(); iter.hasNext(); ) {
							//dynamicSimulatorList.add(new JDASPK(rtol, atol, 0, (InitialStatus)iter.next()));
							dynamicSimulatorList.add(new Chemkin(rtol, atol, reactorType));//11/4/07 gmagoon: fixing apparent cut/paste error
						}
        			}
        		}
        		else throw new InvalidSymbolException("condition.txt: Unknown DynamicSimulator = " + simulator);
                //10/23/07 gmagoon: changed below from dynamicSimulator to dynamicSimulatorList; note: although conversionSet should actually be different for each T,P condition, it will be modified in isTPCconsistent within ReactionSystem
                for (Iterator iter = dynamicSimulatorList.iterator(); iter.hasNext(); ) {
                    double [] cs = conversionSet.clone();//11/1/07 gmagoon: trying to make sure multiple instances of conversionSet are used
                    ((DynamicSimulator)(iter.next())).addConversion(cs, numConversions);
				}
			}
        	else throw new InvalidSymbolException("condition.txt: can't find DynamicSimulator!");
			
        	// read in reaction model enlarger
			
        	/* Read in the Primary Reaction Library
        	 * MRH 12-Jun-2009
        	 * 
        	 * I've made minor changes to this piece of code.  In particular, I've
        	 * 	eliminated the on/off flag.  The user can specify as many PRLs,
        	 * 	including none, as they like.
        	 */        	
        	line = ChemParser.readMeaningfulLine(reader);
			if (line.startsWith("PrimaryReactionLibrary:")) {
				/*
				 * MRH 27Feb2010:
				 * Changing the "read in Primary Reaction Library information" code
				 * 	into it's own method.
				 * 
				 * Other modules (e.g. PopulateReactions) will be utilizing the exact code.
				 * 	Rather than copying and pasting code into other modules, just have
				 * 	everything call this new method: readAndMakePRL
				 */
				readAndMakePRL(reader);
//				// GJB modified to allow multiple primary reaction libraries
//				int Ilib = 0;
//				line = ChemParser.readMeaningfulLine(reader);
//				while (!line.equals("END")) {
//					String[] tempString = line.split("Name: ");
//					String name = tempString[tempString.length-1].trim();
//					line = ChemParser.readMeaningfulLine(reader);
//					tempString = line.split("Location: ");
//					String path = tempString[tempString.length-1].trim();
//					if (Ilib==0) {
//						//primaryReactionLibrary = new PrimaryReactionLibrary(name, path);
//						setPrimaryReactionLibrary(new PrimaryReactionLibrary(name, path));//10/14/07 gmagoon: changed to use setPrimaryReactionLibrary
//						Ilib++; 	
//					}
//					else {
//						//primaryReactionLibrary.appendPrimaryReactionLibrary(name,path);
//						getPrimaryReactionLibrary().appendPrimaryReactionLibrary(name,path);//10/14/07 gmagoon: changed to use getPrimaryReactionLibrary; check to make sure this is valid
//						Ilib++;//just in case anybody wants to track how many are processed
//					}
//					line = ChemParser.readMeaningfulLine(reader);
//				}
//				//                	System.out.println("Primary Reaction Libraries in use: " +getPrimaryReactionLibrary().getName());//10/14/07 gmagoon: changed to use getPrimaryReactionLibrary
//				if (Ilib==0) {
//					//primaryReactionLibrary = null;
//					setPrimaryReactionLibrary(null);//10/14/07 gmagoon: changed to use setPrimaryReactionLibrary; check to make sure this is valid
//				}
//				else System.out.println("Primary Reaction Libraries in use: " + getPrimaryReactionLibrary().getName());
				
			} else throw new InvalidSymbolException("condition.txt: can't find PrimaryReactionLibrary!");
			
			/*
			 * Added by MRH 12-Jun-2009
			 * 
			 * The SeedMechanism acts almost exactly as the old
			 * 	PrimaryReactionLibrary did.  Whatever is in the SeedMechanism
			 * 	will be placed in the core at the beginning of the simulation.
			 * 	The user can specify as many seed mechanisms as they like, with
			 * 	the priority (in the case of duplicates) given to the first
			 * 	instance.  There is no on/off flag.
			 */
			line = ChemParser.readMeaningfulLine(reader);
			if (line.startsWith("SeedMechanism:")) {
				int numMechs = 0;
				line = ChemParser.readMeaningfulLine(reader);
				while (!line.equals("END")) {                     		
					String[] tempString = line.split("Name: ");
					String name = tempString[tempString.length-1].trim();
					line = ChemParser.readMeaningfulLine(reader);
					tempString = line.split("Location: ");
					String path = tempString[tempString.length-1].trim();
					line = ChemParser.readMeaningfulLine(reader);
					tempString = line.split("GenerateReactions: ");
					String generateStr = tempString[tempString.length-1].trim();
					boolean generate = true;
					if (generateStr.equalsIgnoreCase("yes") ||
						generateStr.equalsIgnoreCase("on") ||
						generateStr.equalsIgnoreCase("true")){
						generate = true;
						System.out.println("Will generate cross-reactions between species in seed mechanism " + name);
					} else if(generateStr.equalsIgnoreCase("no") ||
							  generateStr.equalsIgnoreCase("off") ||
							  generateStr.equalsIgnoreCase("false")) {
						generate = false;
						System.out.println("Will NOT initially generate cross-reactions between species in seed mechanism "+ name);
						System.out.println("This may have unintended consequences");			   
					}
					else {
						System.err.println("Input file invalid");
						System.err.println("Please include a 'GenerateReactions: yes/no' line for seed mechanism "+name);
						System.exit(0);
					}
										   
					if (numMechs==0) {
						setSeedMechanism(new SeedMechanism(name, path, generate));
						++numMechs; 	
					}
					else {
						getSeedMechanism().appendSeedMechanism(name, path, generate);
						++numMechs;
					}
					line = ChemParser.readMeaningfulLine(reader);
				}
				if (numMechs != 0)	System.out.println("Seed Mechanisms in use: " + getSeedMechanism().getName());
				else setSeedMechanism(null);
			} else throw new InvalidSymbolException("Error reading condition.txt file: "
													+ "Could not locate SeedMechanism field");
			
			line = ChemParser.readMeaningfulLine(reader);
			if (line.startsWith("ChemkinUnits")) {
				line = ChemParser.readMeaningfulLine(reader);
				if (line.startsWith("Verbose:")) {
					StringTokenizer st = new StringTokenizer(line);
					String dummyString = st.nextToken();
					String OnOff = st.nextToken().toLowerCase();
					if (OnOff.equals("off")) {
						ArrheniusKinetics.setVerbose(false);
					} else if (OnOff.equals("on")) {
						ArrheniusKinetics.setVerbose(true);
					}
					line = ChemParser.readMeaningfulLine(reader);
				}
				/*
				 * MRH 3MAR2010:
				 * Adding user option regarding chemkin file
				 * 
				 * New field: If user would like the empty SMILES string 
				 * printed with each species in the thermochemistry portion 
				 * of the generated chem.inp file
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
						 * MRH 9MAR2010:
						 * MRH decided not to generate an InChI for every new species
						 * 	during an RMG simulation (especially since it is not used
						 * 	for anything).  Instead, they will only be generated in the
						 * 	post-processing, if the user asked for InChIs.
						 */
						//Species.useInChI = true;
					}
					line = ChemParser.readMeaningfulLine(reader);
				}
				if (line.startsWith("A")) {
					StringTokenizer st = new StringTokenizer(line);
					String dummyString = st.nextToken(); // Should be "A:"
					String units = st.nextToken();
					if (units.equals("moles") || units.equals("molecules"))
						ArrheniusKinetics.setAUnits(units);
					else {
						System.err.println("Units for A were not recognized: " + units);
						System.exit(0);
					}
				} else throw new InvalidSymbolException("Error reading condition.txt file: "
														+ "Could not locate Chemkin units A field.");
				line = ChemParser.readMeaningfulLine(reader);
				if (line.startsWith("Ea")) {
					StringTokenizer st = new StringTokenizer(line);
					String dummyString = st.nextToken(); // Should be "Ea:"
					String units = st.nextToken();
					if (units.equals("kcal/mol") || units.equals("cal/mol") ||
						units.equals("kJ/mol") || units.equals("J/mol") || units.equals("Kelvins"))
						ArrheniusKinetics.setEaUnits(units);
					else {
						System.err.println("Units for Ea were not recognized: " + units);
						System.exit(0);
					}
				} else throw new InvalidSymbolException("Error reading condition.txt file: "
														+ "Could not locate Chemkin units Ea field.");
			} else throw new InvalidSymbolException("Error reading condition.txt file: "
													+ "Could not locate ChemkinUnits field.");
			
        	in.close();
			
			//11/6/07 gmagoon: initializing temperatureArray and pressureArray before libraryReactionGenerator is initialized (initialization calls PDepNetwork and performs initializekLeak); UPDATE: moved after initialStatusList initialization (in case primaryReactionLibrary calls the similar pdep functions
			//                LinkedList temperatureArray = new LinkedList();
			//                LinkedList pressureArray = new LinkedList();
			//                Iterator iterIS = initialStatusList.iterator();
			//                for (Iterator iter = tempList.iterator(); iter.hasNext(); ) {
			//                      TemperatureModel tm = (TemperatureModel)iter.next();
			//                      for (Iterator iter2 = presList.iterator(); iter2.hasNext(); ){
			//                        PressureModel pm = (PressureModel)iter2.next();
			//                        InitialStatus is = (InitialStatus)iterIS.next();
			//                        temperatureArray.add(tm.getTemperature(is.getTime()));
			//                        pressureArray.add(pm.getPressure(is.getTime()));
			//                      }
			//                }
			//                PDepNetwork.setTemperatureArray(temperatureArray);
			//                PDepNetwork.setPressureArray(pressureArray);
			
			
			
			//10/4/07 gmagoon: moved to modelGeneration()
			//ReactionGenerator p_reactionGenerator = new TemplateReactionGenerator();//10/4/07 gmagoon: changed to p_reactionGenerator from reactionGenerator
			// setReactionGenerator(p_reactionGenerator);//10/4/07 gmagoon: added
        	/*
        	 * MRH 12-Jun-2009
        	 * A TemplateReactionGenerator now requires a Temperature be passed to it.
        	 * 	This allows RMG to determine the "best" kinetic parameters to use
        	 * 	in the mechanism generation.  For now, I choose to pass the first
        	 * 	temperature in the list of temperatures.  RMG only outputs one mechanism,
        	 * 	even for multiple temperature/pressure systems, so we can only have one
        	 * 	set of kinetics.
        	 */
        	Temperature t = new Temperature(300,"K");
        	for (Iterator iter = tempList.iterator(); iter.hasNext();) {
        		TemperatureModel tm  = (TemperatureModel)iter.next();
        		t = tm.getTemperature(new ReactionTime(0,"sec"));
        		setTemp4BestKinetics(t);
        		break;
        	}
			setReactionGenerator(new TemplateReactionGenerator()); //11/4/07 gmagoon: moved from modelGeneration; mysteriously, moving this later moves "Father" lines up in output at runtime, immediately after condition file (as in original code); previously, these Father lines were just before "Can't read primary reaction library files!"
			lrg = new LibraryReactionGenerator();//10/10/07 gmagoon: moved from modelGeneration (sequence lrg increases species id, and the different sequence was causing problems as main species id was 6 instead of 1); //10/31/07 gmagoon: restored this line from 10/10/07 backup: somehow it got lost along the way; 11/5/07 gmagoon: changed to use "lrg =" instead of setLibraryReactionGenerator
			//10/24/07 gmagoon: updated to use multiple reactionSystem variables
			reactionSystemList = new LinkedList();
			// LinkedList temperatureArray = new LinkedList();//10/30/07 gmagoon: added temperatureArray variable for passing to PDepNetwork; 11/6/07 gmagoon: moved before initialization of lrg;
			// LinkedList pressureArray = new LinkedList();//10/30/07 gmagoon: same for pressure;//UPDATE: commenting out: not needed if updateKLeak is done for one temperature/pressure at a time; 11/1-2/07 restored;11/6/07 gmagoon: moved before initialization of lrg;
			Iterator iter3 = initialStatusList.iterator();
			Iterator iter4 = dynamicSimulatorList.iterator();
			int i = 0;//10/30/07 gmagoon: added
			for (Iterator iter = tempList.iterator(); iter.hasNext(); ) {
				TemperatureModel tm = (TemperatureModel)iter.next();
				//InitialStatus is = (InitialStatus)iter3.next();//10/31/07 gmagoon: fixing apparent bug by moving these inside inner "for loop"
				//DynamicSimulator ds = (DynamicSimulator)iter4.next();
				for (Iterator iter2 = presList.iterator(); iter2.hasNext(); ){
					PressureModel pm = (PressureModel)iter2.next();
					InitialStatus is = (InitialStatus)iter3.next();//10/31/07 gmagoon: moved from outer "for loop""
					DynamicSimulator ds = (DynamicSimulator)iter4.next();
					// temperatureArray.add(tm.getTemperature(is.getTime()));//10/30/07 gmagoon: added; //10/31/07 added .getTemperature(is.getTime()); 11/6/07 gmagoon: moved before initialization of lrg;
					// pressureArray.add(pm.getPressure(is.getTime()));//10/30/07 gmagoon: added//UPDATE: commenting out: not needed if updateKLeak is done for one temperature/pressure at a time;11/1-2/07 restored with .getTemperature(is.getTime()) added;11/6/07 gmagoon: moved before initialization of lrg;
					//11/1/07 gmagoon: trying to make a deep copy of terminationTester when it is instance of ConversionTT
					//UPDATE: actually, I don't think this deep copy was necessary; original case with FinishController fc = new FinishController(finishController.getTerminationTester(), finishController.getValidityTester()) is probably OK; (in any case, this didn't do completetly deep copy (references to speciesConversion element in LinkedList were the same);
					// TerminationTester termTestCopy;
					// if (finishController.getTerminationTester() instanceof ConversionTT){
					//    ConversionTT termTest = (ConversionTT)finishController.getTerminationTester();
					//     LinkedList spcCopy = (LinkedList)(termTest.getSpeciesGoalConversionSetList().clone());
					//     termTestCopy = new ConversionTT(spcCopy);
					// }
					// else{
					//     termTestCopy = finishController.getTerminationTester();
					// }
					
					FinishController fc = new FinishController(finishController.getTerminationTester(), finishController.getValidityTester());//10/31/07 gmagoon: changed to create new finishController instance in each case (apparently, the finish controller becomes associated with reactionSystem in setFinishController within ReactionSystem); alteratively, could use clone, but might need to change FinishController to be "cloneable"
					// FinishController fc = new FinishController(termTestCopy, finishController.getValidityTester());
					reactionSystemList.add(new ReactionSystem(tm, pm, reactionModelEnlarger, fc, ds, getPrimaryReactionLibrary(), getReactionGenerator(), getSpeciesSeed(), is, getReactionModel(),lrg, i, equationOfState)); 
					i++;//10/30/07 gmagoon: added
					System.out.println("Created reaction system "+i+"\n");
				}
			}
			//    PDepNetwork.setTemperatureArray(temperatureArray);//10/30/07 gmagoon: passing temperatureArray to PDepNetwork; 11/6/07 gmagoon: moved before initialization of lrg;
			//    PDepNetwork.setPressureArray(pressureArray);//10/30/07 gmagoon: same for pressure;//UPDATE: commenting out: not needed if updateKLeak is done for one temperature/pressure at a time; 11/1-2/07 restored; 11/6/07 gmagoon: moved before initialization of lrg;
		}
        catch (IOException e) {
        	System.err.println("Error in read in reaction system initialization file!");
        	throw new IOException("Reaction System Initialization: " + e.getMessage());
        }
    }
    public void setReactionModel(ReactionModel p_ReactionModel) {
        reactionModel = p_ReactionModel;
    }
	
	
    public void modelGeneration() {
        //long begin_t = System.currentTimeMillis();
		try{
        	ChemGraph.readForbiddenStructure();
			setSpeciesSeed(new LinkedHashSet());//10/4/07 gmagoon moved from initializeCoreEdgeReactionModel
			//  setReactionGenerator(new TemplateReactionGenerator());//10/4/07 gmagoon: moved inside initializeReactionSystem; 11/3-4/07 gmagoon: probably reverted on or before 10/10/07 (although I have not investigated this change in detail); //11/4/07 gmagoon: moved inside initializeReactionSystems
			//  setLibraryReactionGenerator(new LibraryReactionGenerator());//10/10/07 gmagoon: moved after initializeReactionSystem
			//  initializeCoreEdgeReactionModel();//10/4/07 gmagoon moved from below to run initializeCoreEdgeReactionModel before initializeReactionSystem; 11/3-4/07 gmagoon: probably reverted on or before 10/10/07
         	initializeReactionSystems();
        }
        catch (IOException e) {
        	System.err.println(e.getMessage());
        	System.exit(0);
        }
        catch (InvalidSymbolException e) {
        	System.err.println(e.getMessage());
        	System.exit(0);
        }
		
        //10/31/07 gmagoon: initialize validList (to false) before initializeCoreEdgeReactionModel is called
        validList = new LinkedList();
        for (Integer i = 0; i<reactionSystemList.size();i++) {
            validList.add(false);
		}
        
        initializeCoreEdgeReactionModel();//10/4/07 gmagoon: moved before initializeReactionSystem; 11/3-4/07 gmagoon: probably reverted on or before 10/10/07
        //10/24/07 gmagoon: changed to use reactionSystemList
		//       LinkedList initList = new LinkedList();//10/25/07 gmagoon: moved these variables to apply to entire class
		//       LinkedList beginList = new LinkedList();
		//       LinkedList endList = new LinkedList();
		//       LinkedList lastTList = new LinkedList();
		//       LinkedList currentTList = new LinkedList();
		//       LinkedList lastPList = new LinkedList();
		//       LinkedList currentPList = new LinkedList();
		//       LinkedList conditionChangedList = new LinkedList();
		//       LinkedList reactionChangedList = new LinkedList();
        //5/6/08 gmagoon: determine whether there are intermediate time/conversion steps, type of termination tester is based on characteristics of 1st reaction system (it is assumed that they are all identical in terms of type of termination tester)
        boolean intermediateSteps = true;
        ReactionSystem rs0 = (ReactionSystem)reactionSystemList.get(0);
        if (rs0.finishController.terminationTester instanceof ReactionTimeTT){
            if (timeStep == null){
                intermediateSteps = false;
            }
        }
        else if (numConversions==1){ //if we get to this block, we presumably have a conversion terminationTester; this required moving numConversions to be attribute...alternative to using numConversions is to access one of the DynamicSimulators and determine conversion length
            intermediateSteps=false;
        }
        //10/24/07 gmagoon: note: each element of for loop could be done in parallel if desired; some modifications would be needed
        for (Iterator iter = reactionSystemList.iterator(); iter.hasNext(); ) {
            ReactionSystem rs = (ReactionSystem)iter.next();
            if ((reactionModelEnlarger instanceof RateBasedPDepRME)) {//1/2/09 gmagoon and rwest: only call initializePDepNetwork for P-dep cases
                rs.initializePDepNetwork();
            }
			
            ReactionTime init = rs.getInitialReactionTime();
            initList.add(init);
            ReactionTime begin = init;
            beginList.add(begin);
            ReactionTime end;
            if (rs.finishController.terminationTester instanceof ReactionTimeTT){
				//5/5/08 gmagoon: added below if statement to avoid null pointer exception in cases where there are no intermediate time steps specified
				if (!(timeStep==null)){
					end = (ReactionTime)timeStep.get(0);
				}             
				else{
					end= ((ReactionTimeTT)rs.finishController.terminationTester).finalTime;
				}
				//end = (ReactionTime)timeStep.get(0);
				endList.add(end);
            }
            else{
				end = new ReactionTime(1e6,"sec");
				endList.add(end);
            }
			//           int iterationNumber = 1;
			
            lastTList.add(rs.getTemperature(init));
            currentTList.add(rs.getTemperature(init));
            lastPList.add(rs.getPressure(init));
            currentPList.add(rs.getPressure(init));
            conditionChangedList.add(false);
            reactionChangedList.add(false);//10/31/07 gmagoon: added
			
            //Chemkin.writeChemkinInputFile(reactionSystem.getReactionModel(),reactionSystem.getPresentStatus());
        }
        int iterationNumber = 1;
        LinkedList terminatedList = new LinkedList();//10/24/07 gmagoon: this may not be necessary, as if one reactionSystem is terminated, I think all should be terminated
        //validList = new LinkedList();//10/31/07 gmagoon: moved before initializeCoreEdgeReactionModel
        //10/24/07 gmagoon: initialize allTerminated and allValid to true; these variables keep track of whether all the reactionSystem variables satisfy termination and validity, respectively
        boolean allTerminated = true;
        boolean allValid = true;
        
        // IF RESTART IS TURNED ON
        // Update the systemSnapshot for each ReactionSystem in the reactionSystemList
        if (readrestart) {
	        for (Integer i=0; i<reactionSystemList.size(); i++) {
	        	ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
	        	InitialStatus is = rs.getInitialStatus();
	        	putRestartSpeciesInInitialStatus(is,i);
				rs.appendUnreactedSpeciesStatus((InitialStatus)initialStatusList.get(i), rs.getPresentTemperature());
	        }
        }
        
        //10/24/07 gmagoon: note: each element of for loop could be done in parallel if desired; some modifications would be needed
        for (Integer i = 0; i<reactionSystemList.size();i++) {
            ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
            ReactionTime begin = (ReactionTime)beginList.get(i);
            ReactionTime end = (ReactionTime)endList.get(i);
            endList.set(i,rs.solveReactionSystem(begin, end, true, true, true, iterationNumber-1));
            Chemkin.writeChemkinInputFile(rs);//11/9/07 gmagoon:****temporarily commenting out: there is a NullPointerException in Reaction.toChemkinString when called from writeChemkinPdepReactions; occurs with pre-modified version of RMG as well; //11/12/07 gmagoon: restored; ****this appears to be source of non-Pdep bug
            boolean terminated = rs.isReactionTerminated();
            terminatedList.add(terminated);
            if(!terminated)
                allTerminated = false;
            boolean valid = rs.isModelValid();
            //validList.add(valid);
            validList.set(i, valid);//10/31/07 gmagoon: validList initialization moved before initializeCoreEdgeReactionModel
            if(!valid)
                allValid = false;
            reactionChangedList.set(i,false);
        }
        //9/1/09 gmagoon: if we are using QM, output a file with the CHEMKIN name, the RMG name, the (modified) InChI, and the (modified) InChIKey
        if (ChemGraph.useQM){
            writeInChIs(getReactionModel());      
        }
        writeDictionary(getReactionModel());
        //System.exit(0);
		
		System.out.println("The model core has " + ((CoreEdgeReactionModel)getReactionModel()).getReactedReactionSet().size() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
		System.out.println("The model edge has " + ((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSet().size() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getUnreactedSpeciesSet().size() + " species.");
		
		
		StringBuilder print_info = Global.diagnosticInfo;
		print_info.append("\nMolecule \t Flux\t\tTime\t \t\t \t Core \t \t Edge \t \t memory\n");
		
		print_info.append(" \t moleular \t characteristic \t findspecies \t moveUnreactedToReacted \t enlarger \t restart1 \t totalEnlarger \t resetSystem  \t readSolverFile\t writeSolverFile \t justSolver \t SolverIterations \t solverSpeciesStatus \t Totalsolver \t gc  \t restart+diagnosis \t chemkin thermo \t chemkin reactions \t validitytester \t Species \t Reactions\t Species\t Reactions \t memory used  \t allSpecies \t TotalTime \t findRateConstant\t identifyReactedSites \t reactChemGraph \t makespecies\t CheckReverseReaction \t makeTemplateReaction \t getReactionfromStruc \t genReverseFromReac");
		print_info.append("\t\t\t\t\t\t\t" + ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size()+ "\t" + ((CoreEdgeReactionModel)getReactionModel()).getReactedReactionSet().size() + "\t" + ((CoreEdgeReactionModel)getReactionModel()).getUnreactedSpeciesSet().size() + "\t" + ((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSetIncludingReverseSize() + "\t"+Global.makeSpecies+"\n");
		
		
		double solverMin = 0; 
		double vTester = 0;
		
		/*if (!restart){
		 writeRestartFile();
		 writeCoreReactions();
		 writeAllReactions();
		 }*/
		
		//System.exit(0);
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		System.out.println("Species dictionary size: "+dictionary.size());
		//boolean reactionChanged = false;//10/24/07 gmagoon: I don't know if this is even required, but I will change to use reactionChangedList (I put analogous line of code for list in above for loop); update: yes, it is required; I had been thinking of conditionChangedList
		
		double tAtInitialization = Global.tAtInitialization;
		
		//10/24/07: changed to use allTerminated and allValid	
        // step 2: iteratively grow reaction system
        while (!allTerminated || !allValid) {
        	while (!allValid) {
				
				//writeCoreSpecies();
				double pt = System.currentTimeMillis();
				// ENLARGE THE MODEL!!! (this is where the good stuff happens)
				enlargeReactionModel();//10/24/07 gmagoon: need to adjust this function
				double totalEnlarger = (System.currentTimeMillis() - pt)/1000/60;
				
				//PDepNetwork.completeNetwork(reactionSystem.reactionModel.getSpeciesSet());
				//10/24/07 gmagoon: changed to use reactionSystemList
				if ((reactionModelEnlarger instanceof RateBasedPDepRME)) {//1/2/09 gmagoon and rwest: only call initializePDepNetwork for P-dep cases
					for (Iterator iter = reactionSystemList.iterator(); iter.hasNext(); ) {
						ReactionSystem rs = (ReactionSystem)iter.next();
						rs.initializePDepNetwork();
					}
					//reactionSystem.initializePDepNetwork();
				}  
				
				pt = System.currentTimeMillis();
				//10/24/07 gmagoon: changed to use reactionSystemList
				for (Iterator iter = reactionSystemList.iterator(); iter.hasNext(); ) {
					ReactionSystem rs = (ReactionSystem)iter.next();
					rs.resetSystemSnapshot();
				}
				//reactionSystem.resetSystemSnapshot();
				double resetSystem = (System.currentTimeMillis() - pt)/1000/60;
				//10/24/07 gmagoon: changed to use reactionSystemList
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					//reactionChanged = true;
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					reactionChangedList.set(i,true);
					
					// begin = init;
					beginList.set(i, (ReactionTime)initList.get(i));
					if (rs.finishController.terminationTester instanceof ReactionTimeTT){
						//5/5/08 gmagoon: added below if statement to avoid null pointer exception in cases where there are no intermediate time steps specified
						if (!(timeStep==null)){
							endList.set(i,(ReactionTime)timeStep.get(0));
						}             
						else{
							endList.set(i, ((ReactionTimeTT)rs.finishController.terminationTester).finalTime);
						}
						// endList.set(i, (ReactionTime)timeStep.get(0));
						//end = (ReactionTime)timeStep.get(0);
					}
					else
						endList.set(i, new ReactionTime(1e6,"sec"));
					//end = new ReactionTime(1e6,"sec");
					//  iterationNumber = 1;//10/24/07 gmagoon: moved outside of loop
					
					currentTList.set(i,rs.getTemperature((ReactionTime)beginList.get(i)));
					currentPList.set(i,rs.getPressure((ReactionTime)beginList.get(i)));
					conditionChangedList.set(i,!(((Temperature)currentTList.get(i)).equals((Temperature)lastTList.get(i))) || !(((Pressure)currentPList.get(i)).equals((Pressure)lastPList.get(i))));
					//currentT = reactionSystem.getTemperature(begin);
					//currentP = reactionSystem.getPressure(begin);
					//conditionChanged = (!currentT.equals(lastT) || !currentP.equals(lastP));
				}
				iterationNumber = 1;
				
				double startTime = System.currentTimeMillis();
				
				//10/24/07 gmagoon: changed to use reactionSystemList
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					boolean reactionChanged = (Boolean)reactionChangedList.get(i);
					boolean conditionChanged = (Boolean)conditionChangedList.get(i);
					ReactionTime begin = (ReactionTime)beginList.get(i);
					ReactionTime end = (ReactionTime)endList.get(i);
					endList.set(i,rs.solveReactionSystem(begin, end, false, reactionChanged, conditionChanged, iterationNumber-1));
					//end = reactionSystem.solveReactionSystem(begin, end, false, reactionChanged, conditionChanged, iterationNumber-1);
				}
				solverMin = solverMin + (System.currentTimeMillis()-startTime)/1000/60;
				
				startTime = System.currentTimeMillis();
				//10/24/07 gmagoon: changed to use reactionSystemList
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					Chemkin.writeChemkinInputFile(rs);//10/25/07 gmagoon: ***I don't know if this will still work with multiple reaction systems: may want to modify to only write one chemkin input file for all reaction systems //11/9/07 gmagoon:****temporarily commenting out; cf. previous comment; //11/12/07 gmagoon: restored; ****this appears to be source of non-Pdep bug 
					//Chemkin.writeChemkinInputFile(reactionSystem);
				}
				//9/1/09 gmagoon: if we are using QM, output a file with the CHEMKIN name, the RMG name, the (modified) InChI, and the (modified) InChIKey
				if (ChemGraph.useQM){
					writeInChIs(getReactionModel());                    
				}
				writeDictionary(getReactionModel());
				double chemkint = (System.currentTimeMillis()-startTime)/1000/60;
				
				if (writerestart) {
					writeCoreSpecies();
					writeCoreReactions();
					writeEdgeSpecies();
					writeEdgeReactions();
					if (PDepNetwork.generateNetworks == true)	writePDepNetworks();
				}
				
				//10/24/07 gmagoon: changed to use reactionSystemList
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					System.out.println("For reaction system: "+(i+1)+" out of "+reactionSystemList.size());
					System.out.println("At this time: " + ((ReactionTime)endList.get(i)).toString());
					Species spe = SpeciesDictionary.getSpeciesFromID(1);
					double conv = rs.getPresentConversion(spe);
					System.out.print("Conversion of " + spe.getName()  + " is:");
					System.out.println(conv);
				}
				
			    System.out.println("Running Time is: " + String.valueOf((System.currentTimeMillis()-tAtInitialization)/1000/60) + " minutes.");
				System.out.println("The model edge has " + ((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSet().size() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getUnreactedSpeciesSet().size() + " species.");
				//10/24/07 gmagoon: note: all reaction systems should use the same core, but I will display for each reactionSystem for testing purposes:
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					System.out.println("For reaction system: "+(i+1)+" out of "+reactionSystemList.size());
					if (rs.getDynamicSimulator() instanceof JDASPK){
						JDASPK solver = (JDASPK)rs.getDynamicSimulator();
						System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
					}
					else{
						JDASSL solver = (JDASSL)rs.getDynamicSimulator();
						System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
					}
				}
				// if (reactionSystem.getDynamicSimulator() instanceof JDASPK){
				//	JDASPK solver = (JDASPK)reactionSystem.getDynamicSimulator();
				//	System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
				//}
				//else{
				//	JDASSL solver = (JDASSL)reactionSystem.getDynamicSimulator();
				//	System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
				//}
				
				startTime = System.currentTimeMillis();
				double mU = memoryUsed();
				double gc = (System.currentTimeMillis()-startTime)/1000/60;
				
				startTime = System.currentTimeMillis();
				//10/24/07 gmagoon: updating to use reactionSystemList
				allValid = true;
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					boolean valid = rs.isModelValid();
					if(!valid)
						allValid = false;
					validList.set(i,valid);
					//valid = reactionSystem.isModelValid();
				}
				vTester = vTester + (System.currentTimeMillis()-startTime)/1000/60;
				
				startTime = System.currentTimeMillis();
				
				writeDiagnosticInfo();
				writeEnlargerInfo();
				double restart2 = (System.currentTimeMillis()-startTime)/1000/60;
				
				int allSpecies, allReactions;
				allSpecies = SpeciesDictionary.getInstance().size();
				print_info.append(totalEnlarger + "\t" + resetSystem + "\t" + Global.readSolverFile + "\t" + Global.writeSolverFile + "\t" + Global.solvertime + "\t" + Global.solverIterations + "\t" + Global.speciesStatusGenerator +  "\t" + solverMin + "\t"  + gc + "\t"  + restart2 + "\t" + Global.chemkinThermo + '\t' + Global.chemkinReaction + "\t" + vTester + "\t" + ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size()+ "\t" + ((CoreEdgeReactionModel)getReactionModel()).getReactedReactionSet().size() + "\t" + ((CoreEdgeReactionModel)getReactionModel()).getUnreactedSpeciesSet().size() + "\t" + ((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSetIncludingReverseSize() + "\t" + mU + "\t" + allSpecies + "\t" + (System.currentTimeMillis()-Global.tAtInitialization)/1000/60 + "\t"+ String.valueOf(Global.RT_findRateConstant)+"\t"+Global.RT_identifyReactedSites+"\t"+Global.RT_reactChemGraph+"\t"+Global.makeSpecies+"\t"+Global.checkReactionReverse+"\t"+Global.makeTR+ "\t" + Global.getReacFromStruc + "\t" + Global.generateReverse+"\n");
				
        	}
			//5/6/08 gmagoon: in order to handle cases where no intermediate time/conversion steps are used, only evaluate the next block of code when there are intermediate time/conversion steps
			double startTime = System.currentTimeMillis();
			if(intermediateSteps){
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					reactionChangedList.set(i, false);
					//reactionChanged = false;
					Temperature currentT = (Temperature)currentTList.get(i);
					Pressure currentP = (Pressure)currentPList.get(i);
					lastTList.set(i,(Temperature)currentT.clone()) ;
					lastPList.set(i,(Pressure)currentP.clone());
					//lastT = (Temperature)currentT.clone();
					//lastP = (Pressure)currentP.clone();
					
					currentTList.set(i,rs.getTemperature((ReactionTime)beginList.get(i)));//10/24/07 gmagoon: ****I think this should actually be at end? (it shouldn't matter for isothermal/isobaric case)
					currentPList.set(i,rs.getPressure((ReactionTime)beginList.get(i)));
					conditionChangedList.set(i,!(((Temperature)currentTList.get(i)).equals((Temperature)lastTList.get(i))) || !(((Pressure)currentPList.get(i)).equals((Pressure)lastPList.get(i))));
					//currentT = reactionSystem.getTemperature(begin);//10/24/07 gmagoon: ****I think this should actually be at end? (it shouldn't matter for isothermal/isobaric case)
					//currentP = reactionSystem.getPressure(begin);
					//conditionChanged = (!currentT.equals(lastT) || !currentP.equals(lastP));
					
					beginList.set(i,((SystemSnapshot)(rs.getSystemSnapshotEnd().next())).time);
					// begin=((SystemSnapshot)(reactionSystem.getSystemSnapshotEnd().next())).time;
					if (rs.finishController.terminationTester instanceof ReactionTimeTT){
						if (iterationNumber < timeStep.size()){
							endList.set(i,(ReactionTime)timeStep.get(iterationNumber));
							//end = (ReactionTime)timeStep.get(iterationNumber);
						}
						else
							endList.set(i, ((ReactionTimeTT)rs.finishController.terminationTester).finalTime);
						//end = ((ReactionTimeTT)reactionSystem.finishController.terminationTester).finalTime;
					}
					else
						endList.set(i,new ReactionTime(1e6,"sec"));
					//end = new ReactionTime(1e6,"sec");
				}        
				iterationNumber++;
				startTime = System.currentTimeMillis();//5/6/08 gmagoon: moved declaration outside of if statement so it can be accessed in subsequent vTester line; previous steps are probably so fast that I could eliminate this line without much effect on normal operation with intermediate steps
				//double startTime = System.currentTimeMillis();
				//10/24/07 gmagoon: changed to use reactionSystemList
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					boolean reactionChanged = (Boolean)reactionChangedList.get(i);
					boolean conditionChanged = (Boolean)conditionChangedList.get(i);
					ReactionTime begin = (ReactionTime)beginList.get(i);
					ReactionTime end = (ReactionTime)endList.get(i);
					endList.set(i,rs.solveReactionSystem(begin, end, false, reactionChanged, false, iterationNumber-1));
					// end = reactionSystem.solveReactionSystem(begin, end, false, reactionChanged, false, iterationNumber-1);
				}
				solverMin = solverMin + (System.currentTimeMillis()-startTime)/1000/60;
				
				startTime = System.currentTimeMillis();
				
				//5/6/08 gmagoon: changed to separate validity and termination testing, and termination testing is done last...termination testing should be done even if there are no intermediate conversions; however, validity is guaranteed if there are no intermediate conversions based on previous conditional if statement
				allValid = true;
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					boolean valid = rs.isModelValid();     
					validList.set(i,valid);
					if(!valid)
						allValid = false;
				}
			}//5/6/08 gmagoon: end of block for intermediateSteps

			allTerminated = true;
        	for (Integer i = 0; i<reactionSystemList.size();i++) {
				ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
				boolean terminated = rs.isReactionTerminated();
				terminatedList.set(i,terminated);
				if(!terminated){
					allTerminated = false;
					System.out.println("Reaction System "+(i+1)+" has not reached its termination criterion");
					if (rs.isModelValid()&& runKillableToPreventInfiniteLoop(intermediateSteps, iterationNumber)) {
						System.out.println("although it seems to be valid (complete), so it was not interrupted for being invalid.");
						System.out.println("This probably means there was an error with the ODE solver, and we risk entering an endless loop.");
						System.out.println("Stopping.");
						throw new Error();
					}
				}
			}
			//     //10/24/07 gmagoon: changed to use reactionSystemList
			//     allTerminated = true;
			//     allValid = true;
			//  	for (Integer i = 0; i<reactionSystemList.size();i++) {
            //        ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
            //        boolean terminated = rs.isReactionTerminated();
            //        terminatedList.set(i,terminated);
            //        if(!terminated)
            //            allTerminated = false;
            //        boolean valid = rs.isModelValid();     
            //        validList.set(i,valid);
            //        if(!valid)
            //            allValid = false;
            //    }
            //    //terminated = reactionSystem.isReactionTerminated();
            //    //valid = reactionSystem.isModelValid();
			
			//10/24/07 gmagoon: changed to use reactionSystemList, allValid	
        	if (allValid) {
				//10/24/07 gmagoon: changed to use reactionSystemList
				for (Integer i = 0; i<reactionSystemList.size();i++) {
					ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
					System.out.println("For reaction system: "+(i+1)+" out of "+reactionSystemList.size());
					System.out.println("At this reaction time: " + ((ReactionTime)endList.get(i)).toString());
					Species spe = SpeciesDictionary.getSpeciesFromID(1);
					double conv = rs.getPresentConversion(spe);
					System.out.print("Conversion of " + spe.getName()  + " is:");
					System.out.println(conv);
				}
        		//System.out.println("At this time: " + end.toString());
        		//Species spe = SpeciesDictionary.getSpeciesFromID(1);
        		//double conv = reactionSystem.getPresentConversion(spe);
        		//System.out.print("current conversion = ");
        		//System.out.println(conv);
				
        		Runtime runTime = Runtime.getRuntime();
        		System.out.print("Memory used: ");
        		System.out.println(runTime.totalMemory());
        		System.out.print("Free memory: ");
        		System.out.println(runTime.freeMemory());
				
        		//runTime.gc();
				/* if we're not calling runTime.gc() then don't bother printing this:
				 System.out.println("After garbage collection:");
				 System.out.print("Memory used: ");
				 System.out.println(runTime.totalMemory());
				 System.out.print("Free memory: ");
				 System.out.println(runTime.freeMemory());
				 */
				
				//10/24/07 gmagoon: note: all reaction systems should use the same core, but I will display for each reactionSystem for testing purposes:
        		for (Integer i = 0; i<reactionSystemList.size();i++) {
        			ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
        			System.out.println("For reaction system: "+(i+1)+" out of "+reactionSystemList.size());
        			if (rs.getDynamicSimulator() instanceof JDASPK){
        				JDASPK solver = (JDASPK)rs.getDynamicSimulator();
        				System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
        			}
        			else{
        				JDASSL solver = (JDASSL)rs.getDynamicSimulator();
        				System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
                        System.out.println("(although rs.getReactionModel().getReactionNumber() returns "+rs.getReactionModel().getReactionNumber()+")");
        			}
        		}
				//        		if (reactionSystem.getDynamicSimulator() instanceof JDASPK){
				//					JDASPK solver = (JDASPK)reactionSystem.getDynamicSimulator();
				//					System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
				//        		}
				//				else{
				//					JDASSL solver = (JDASSL)reactionSystem.getDynamicSimulator();
				//					System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
				//      		}
				
        	}
			vTester = vTester + (System.currentTimeMillis()-startTime)/1000/60;//5/6/08 gmagoon: for case where intermediateSteps = false, this will use startTime declared just before intermediateSteps loop, and will only include termination testing, but no validity testing
        }
        
        //System.out.println("Performing model reduction");
        
		
        if (paraInfor != 0){
			System.out.println("Model Generation performed. Now generating sensitivity data.");
			//10/24/07 gmagoon: updated to use reactionSystemList
			LinkedList dynamicSimulator2List = new LinkedList();
			for (Integer i = 0; i<reactionSystemList.size();i++) {
                ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
                //6/25/08 gmagoon: updated to pass index i
                //6/25/08 gmagoon: updated to pass (dummy) finishController and autoflag (set to false here); 
                dynamicSimulator2List.add(new JDASPK(rtol, atol, paraInfor, (InitialStatus)initialStatusList.get(i),i));
                //DynamicSimulator dynamicSimulator2 = new JDASPK(rtol, atol, paraInfor, initialStatus);
                ((DynamicSimulator)dynamicSimulator2List.get(i)).addConversion(((JDASPK)rs.dynamicSimulator).conversionSet, ((JDASPK)rs.dynamicSimulator).conversionSet.length);
                //dynamicSimulator2.addConversion(((JDASPK)reactionSystem.dynamicSimulator).conversionSet, ((JDASPK)reactionSystem.dynamicSimulator).conversionSet.length);
                rs.setDynamicSimulator((DynamicSimulator)dynamicSimulator2List.get(i));
                //reactionSystem.setDynamicSimulator(dynamicSimulator2);
				
                int numSteps = rs.systemSnapshot.size() -1;
                rs.resetSystemSnapshot();
                beginList.set(i, (ReactionTime)initList.get(i));
                //begin = init;
                if (rs.finishController.terminationTester instanceof ReactionTimeTT){
					endList.set(i,((ReactionTimeTT)rs.finishController.terminationTester).finalTime);
					//end = ((ReactionTimeTT)reactionSystem.finishController.terminationTester).finalTime;
                }
                else{
					ReactionTime end = (ReactionTime)endList.get(i);
					endList.set(i, end.add(end));
					//end = end.add(end);
                }
                terminatedList.set(i, false);
                //terminated = false;
                ReactionTime begin = (ReactionTime)beginList.get(i);
                ReactionTime end = (ReactionTime)endList.get(i);
                rs.solveReactionSystemwithSEN(begin, end, true, false, false);
                //reactionSystem.solveReactionSystemwithSEN(begin, end, true, false, false);
			}
			
        }
        //10/24/07 gmagoon: updated to use reactionSystemList (***see previous comment regarding having one Chemkin input file)
        for (Integer i = 0; i<reactionSystemList.size();i++) {
			ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
			Chemkin.writeChemkinInputFile(getReactionModel(),rs.getPresentStatus());//11/9/07 gmagoon: temporarily commenting out; see previous comment; this line may not cause a problem because it is different instance of writeChemkinInputFile, but I am commenting out just in case//11/12/07 gmagoon: restored; ****this appears to be source of non-Pdep bug 
        }
		
        //9/1/09 gmagoon: if we are using QM, output a file with the CHEMKIN name, the RMG name, the (modified) InChI, and the (modified) InChIKey
        if (ChemGraph.useQM){
            writeInChIs(getReactionModel());    
        }
        writeDictionary(getReactionModel());
		
        System.out.println("Model Generation Completed");

        return;
    }
    
    //9/1/09 gmagoon: this function writes a "dictionary" with Chemkin name, RMG name, (modified) InChI, and InChIKey
    //this is based off of writeChemkinFile in ChemkinInputFile.java
    private void writeInChIs(ReactionModel p_reactionModel) {
        StringBuilder result=new StringBuilder();
        for (Iterator iter = ((CoreEdgeReactionModel)p_reactionModel).core.getSpecies(); iter.hasNext(); ) {
            Species species = (Species) iter.next();
            result.append(species.getChemkinName() + "\t"+species.getName() + "\t" + species.getChemGraph().getModifiedInChIAnew() + "\t" + species.getChemGraph().getModifiedInChIKeyAnew()+ "\n");
        }
		
		
		String file = "inchiDictionary.txt";
		
		try {
			FileWriter fw = new FileWriter(file);
			fw.write(result.toString());
			fw.close();
		}
		catch (Exception e) {
			System.out.println("Error in writing InChI file inchiDictionary.txt!");
			System.out.println(e.getMessage());
			System.exit(0);
		}
    }
    
    //9/14/09 gmagoon: function to write dictionary, based on code copied from RMG.java
    private void writeDictionary(ReactionModel rm){
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)rm;
        //Write core species to RMG_Dictionary.txt
		String coreSpecies ="";
		Iterator iter = cerm.getSpecies();
		
		if (Species.useInChI) {
			while (iter.hasNext()){
				int i=1;
				Species spe = (Species) iter.next();
				coreSpecies = coreSpecies + spe.getChemkinName() + " " + spe.getInChI() + "\n"+spe.getChemGraph().toString(i)+"\n\n";
			}
		} else {
			while (iter.hasNext()){
				int i=1;
				Species spe = (Species) iter.next();
				coreSpecies = coreSpecies + spe.getChemkinName() + "\n"+spe.getChemGraph().toString(i)+"\n\n";
			}
		}
		
		try{
			File rmgDictionary = new File("RMG_Dictionary.txt");
			FileWriter fw = new FileWriter(rmgDictionary);
			fw.write(coreSpecies);
			fw.close();
		}
		catch (IOException e) {
			System.out.println("Could not write RMG_Dictionary.txt");
			System.exit(0);
        }
    }
    
    private void parseRestartFiles() {
		parseAllSpecies();
		parseCoreSpecies();
		parseEdgeSpecies();
		parseAllReactions();
		parseCoreReactions();
		
	}
	
	
	
	private void parseEdgeReactions() {
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		//HasMap speciesMap = dictionary.dictionary;
		try{
			File coreReactions = new File("Restart/edgeReactions.txt");
			FileReader fr = new FileReader(coreReactions);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			boolean found = false;
			LinkedHashSet reactionSet = new LinkedHashSet();
			while (line != null){
				Reaction reaction = ChemParser.parseEdgeArrheniusReaction(dictionary,line,1,1);
				boolean added = reactionSet.add(reaction);
				if (!added){
					if (reaction.hasResonanceIsomerAsReactant()){
						//Structure reactionStructure = reaction.getStructure();
						found = getResonanceStructure(reaction,"reactants", reactionSet);
					}
					if (reaction.hasResonanceIsomerAsProduct() && !found){
						//Structure reactionStructure = reaction.getStructure();
						found = getResonanceStructure(reaction,"products", reactionSet);
					}
					
					if (!found){
						System.out.println("Cannot add reaction "+line+" to the Reaction Edge. All resonance isomers have already been added");
						System.exit(0);
					}
					else found = false;
				}
				//Reaction reverse = reaction.getReverseReaction();
				//if (reverse != null) reactionSet.add(reverse);
				line = ChemParser.readMeaningfulLine(reader);
			}
			((CoreEdgeReactionModel)getReactionModel()).addReactionSet(reactionSet);
		}
		catch (IOException e){
			System.out.println("Could not read the corespecies restart file");
        	System.exit(0);
		}
		
	}
	
	public void parseCoreReactions() {
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		int i=1;
		//HasMap speciesMap = dictionary.dictionary;
		try{
			File coreReactions = new File("Restart/coreReactions.txt");
			FileReader fr = new FileReader(coreReactions);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			boolean found = false;
			LinkedHashSet reactionSet = new LinkedHashSet();
			while (line != null){
				Reaction reaction = ChemParser.parseCoreArrheniusReaction(dictionary,line,1,1);//,((CoreEdgeReactionModel)reactionSystem.reactionModel));
				boolean added = reactionSet.add(reaction);
				if (!added){
					if (reaction.hasResonanceIsomerAsReactant()){
						//Structure reactionStructure = reaction.getStructure();
						found = getResonanceStructure(reaction,"reactants", reactionSet);
					}
					if (reaction.hasResonanceIsomerAsProduct() && !found){
						//Structure reactionStructure = reaction.getStructure();
						found = getResonanceStructure(reaction,"products", reactionSet);
					}
					
					if (!found){
						System.out.println("Cannot add reaction "+line+" to the Reaction Core. All resonance isomers have already been added");
						//System.exit(0);
					}
					else found = false;
				}
				
				Reaction reverse = reaction.getReverseReaction();
				if (reverse != null) {
					reactionSet.add(reverse);
					//System.out.println(2 + "\t " + line);
				}
				
				//else System.out.println(1 + "\t" + line);
				line = ChemParser.readMeaningfulLine(reader);
				i=i+1;
			}
			((CoreEdgeReactionModel)getReactionModel()).addReactedReactionSet(reactionSet);
		}
		catch (IOException e){
			System.out.println("Could not read the coreReactions restart file");
        	System.exit(0);
		}
		
	}
	
	private void parseAllReactions() {
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		int i=1;
		//HasMap speciesMap = dictionary.dictionary;
		try{
			File allReactions = new File("Restart/allReactions.txt");
			FileReader fr = new FileReader(allReactions);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			boolean found = false;
			LinkedHashSet reactionSet = new LinkedHashSet();
		OuterLoop:
			while (line != null){
				Reaction reaction = ChemParser.parseArrheniusReaction(dictionary,line,1,1,((CoreEdgeReactionModel)getReactionModel()));
				
				if (((CoreEdgeReactionModel)getReactionModel()).categorizeReaction(reaction)==-1){
					boolean added = reactionSet.add(reaction);
					if (!added){
						found = false;
						if (reaction.hasResonanceIsomerAsReactant()){
							//Structure reactionStructure = reaction.getStructure();
							found = getResonanceStructure(reaction,"reactants", reactionSet);
						}
						if (reaction.hasResonanceIsomerAsProduct() && !found){
							//Structure reactionStructure = reaction.getStructure();
							found = getResonanceStructure(reaction,"products", reactionSet);
						}
						
						if (!found){
							Iterator iter = reactionSet.iterator();
							while (iter.hasNext()){
								Reaction reacTemp = (Reaction)iter.next();
								if (reacTemp.equals(reaction)){
									reactionSet.remove(reacTemp);
									reactionSet.add(reaction);
									break;
								}
							}
							
							//System.out.println("Cannot add reaction "+line+" to the Reaction Core. All resonance isomers have already been added");
							//System.exit(0);
						}
						
						//else found = false;
					}
				}
				
				
				/*Reaction reverse = reaction.getReverseReaction();
				 if (reverse != null && ((CoreEdgeReactionModel)reactionSystem.reactionModel).isReactedReaction(reaction)) {
				 reactionSet.add(reverse);
				 //System.out.println(2 + "\t " + line);
				 }*/
				//else System.out.println(1 + "\t" + line);
				
				i=i+1;
				
				line = ChemParser.readMeaningfulLine(reader);
			}
			((CoreEdgeReactionModel)getReactionModel()).addReactionSet(reactionSet);
		}
		catch (IOException e){
			System.out.println("Could not read the corespecies restart file");
        	System.exit(0);
		}
		
	}
	
	
	private boolean getResonanceStructure(Reaction p_Reaction, String rOrP, LinkedHashSet reactionSet) {
		Structure reactionStructure = p_Reaction.getStructure();
		//Structure tempreactionStructure = new Structure(reactionStructure.getReactantList(),reactionStructure.getProductList());
		boolean found = false;
		if (rOrP.equals("reactants")){
			Iterator originalreactants = reactionStructure.getReactants();
			HashSet tempHashSet = new HashSet();
			while(originalreactants.hasNext()){
				tempHashSet.add(originalreactants.next());
			}
			Iterator reactants = tempHashSet.iterator();
			while(reactants.hasNext() && !found){
				ChemGraph reactant = (ChemGraph)reactants.next();
				if (reactant.getSpecies().hasResonanceIsomers()){
					Iterator chemGraphIterator = reactant.getSpecies().getResonanceIsomers();
					ChemGraph newChemGraph ;//= (ChemGraph)chemGraphIterator.next();
					while(chemGraphIterator.hasNext() && !found){
						newChemGraph = (ChemGraph)chemGraphIterator.next();
						reactionStructure.removeReactants(reactant);
						reactionStructure.addReactants(newChemGraph);
						reactant = newChemGraph;
						if (reactionSet.add(p_Reaction)){
							found = true;
							
						}
					}
				}
			}
			
		}
		else{
			Iterator originalproducts = reactionStructure.getProducts();
			HashSet tempHashSet = new HashSet();
			while(originalproducts.hasNext()){
				tempHashSet.add(originalproducts.next());
			}
			Iterator products = tempHashSet.iterator();
			
			while(products.hasNext() && !found){
				ChemGraph product = (ChemGraph)products.next();
				if (product.getSpecies().hasResonanceIsomers()){
					Iterator chemGraphIterator = product.getSpecies().getResonanceIsomers();
					ChemGraph newChemGraph ;//= (ChemGraph)chemGraphIterator.next();
					while(chemGraphIterator.hasNext() && !found){
						newChemGraph = (ChemGraph)chemGraphIterator.next();
						reactionStructure.removeProducts(product);
						reactionStructure.addProducts(newChemGraph);
						product = newChemGraph;
						if (reactionSet.add(p_Reaction)){
							found = true;
							
						}
					}
				}
			}
		}
		
		return found;
		
	}
	
	
	
	public void parseCoreSpecies() {
		//		String restartFileContent ="";
		//int speciesCount = 0;
		//boolean added;
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		try{
			File coreSpecies = new File ("Restart/coreSpecies.txt");
			FileReader fr = new FileReader(coreSpecies);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			//HashSet speciesSet = new HashSet();
			//			if (reactionSystem == null){//10/24/07 gmagoon: commenting out since contents of if was already commented out anyway
			//				//ReactionSystem reactionSystem = new ReactionSystem();
			//			}
			setReactionModel(new CoreEdgeReactionModel());//10/4/07 gmagoon:changed to setReactionModel
			while (line!=null) {
				
    			StringTokenizer st = new StringTokenizer(line);
    			String index = st.nextToken();
				int ID = Integer.parseInt(index);
				Species spe = dictionary.getSpeciesFromID(ID);
				if (spe == null)
					System.out.println("There was no species with ID "+ID +" in the species dictionary");
				
				((CoreEdgeReactionModel)getReactionModel()).addReactedSpecies(spe);
				line = ChemParser.readMeaningfulLine(reader);
			}
			
		}
		catch (IOException e){
			System.out.println("Could not read the corespecies restart file");
        	System.exit(0);
		}
		
	}
	public static void garbageCollect(){
		System.gc();
	}
	
	public static long memoryUsed(){
		garbageCollect();
		Runtime rT = Runtime.getRuntime();
		long uM, tM, fM;
		tM = rT.totalMemory();
		fM = rT.freeMemory();
		uM = tM - fM;
		System.out.println("After garbage collection:");
		System.out.print("Memory used: ");
		System.out.println(tM);
		System.out.print("Free memory: ");
		System.out.println(fM);
		
		return uM;
	}
	
	private HashSet readIncludeSpecies(String fileName) {
		HashSet speciesSet = new HashSet();
		try {
			File includeSpecies = new File (fileName);
			FileReader fr = new FileReader(includeSpecies);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			
			
			while (line!=null) {
				
				
    			StringTokenizer st = new StringTokenizer(line);
    			String index = st.nextToken();
    			String name = null;
    			if (!index.startsWith("(")) name = index;
    			else name = st.nextToken().trim();
				
    			Graph g = ChemParser.readChemGraph(reader);
				
				
    			ChemGraph cg = null;
    			try {
    				cg = ChemGraph.make(g);
    			}
    			catch (ForbiddenStructureException e) {
    				System.out.println("Forbidden Structure:\n" + e.getMessage());
    				System.exit(0);
    			}
				
    			Species species = Species.make(name,cg);
       			//speciesSet.put(name, species);
    			speciesSet.add(species);
				
				
    			line = ChemParser.readMeaningfulLine(reader);
				System.out.println(line);
				
    		}
			
			
		}
		catch (IOException e){
			System.out.println("Could not read the included species file");
        	System.exit(0);
		}
		return speciesSet;
	}
	
	public LinkedHashSet parseAllSpecies() {
		//		String restartFileContent ="";
		int speciesCount = 0;
		LinkedHashSet speciesSet = new LinkedHashSet();
		
		boolean added;
		try{
			long initialTime = System.currentTimeMillis();
			
			File coreSpecies = new File ("allSpecies.txt");
			BufferedReader reader = new BufferedReader(new FileReader(coreSpecies));
			String line = ChemParser.readMeaningfulLine(reader);
			int i=0;
			while (line!=null) {
				i++;
    			StringTokenizer st = new StringTokenizer(line);
    			String index = st.nextToken();
    			String name = null;
    			if (!index.startsWith("(")) name = index;
    			else name = st.nextToken().trim();
				int ID = getID(name);
				
				name = getName(name);
    			Graph g = ChemParser.readChemGraph(reader);
    			ChemGraph cg = null;
    			try {
    				cg = ChemGraph.make(g);
    			}
    			catch (ForbiddenStructureException e) {
    				System.out.println("Forbidden Structure:\n" + e.getMessage());
    				System.exit(0);
    			}
    			Species species;
    			if (ID == 0)
    				species = Species.make(name,cg);
    			else 
    				species = Species.make(name,cg,ID);
    			speciesSet.add(species);
    			double flux = 0;
    			int species_type = 1; 
     			line = ChemParser.readMeaningfulLine(reader);
     			System.out.println(line);
    		}
		}
		catch (IOException e){
			System.out.println("Could not read the allSpecies restart file");
        	System.exit(0);
		}
		return speciesSet;
	}
	
	private String getName(String name) {
		//int id;
		String number = "";
		int index=0;
		if (!name.endsWith(")")) return name;
		else {
			char [] nameChars = name.toCharArray();
			String temp = String.copyValueOf(nameChars);
			int i=name.length()-2;
			//char test = "(";
			while (i>0){
				if (name.charAt(i)== '(') {
					index=i;
					i=0;
				}
				else i = i-1;
			}
		}
		number = name.substring(0,index);
		return number;
	}
	
	private int getID(String name) {
		int id;
		String number = "";
		if (!name.endsWith(")")) return 0;
		else {
			char [] nameChars = name.toCharArray();
			
			int i=name.length()-2;
			//char test = "(";
			while (i>0){
				if (name.charAt(i)== '(') i=0;
				else{
					number = name.charAt(i)+number;
					i = i-1;
				}
			}
		}
		id = Integer.parseInt(number);
		return id;
	}
	
	private void parseEdgeSpecies() {
		//		String restartFileContent ="";
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		try{
			File edgeSpecies = new File ("Restart/edgeSpecies.txt");
			FileReader fr = new FileReader(edgeSpecies);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			//HashSet speciesSet = new HashSet();
			
			while (line!=null) {
				
    			StringTokenizer st = new StringTokenizer(line);
    			String index = st.nextToken();
				int ID = Integer.parseInt(index);
				Species spe = dictionary.getSpeciesFromID(ID);
				if (spe == null)
					System.out.println("There was no species with ID "+ID +" in the species dictionary");
				//reactionSystem.reactionModel = new CoreEdgeReactionModel();
				((CoreEdgeReactionModel)getReactionModel()).addUnreactedSpecies(spe);
				line = ChemParser.readMeaningfulLine(reader);
			}
			
		}
		catch (IOException e){
			System.out.println("Could not read the edgepecies restart file");
        	System.exit(0);
		}
		
	}
	
	
	/*private int calculateAllReactionsinReactionTemplate() {
	 int totalnum = 0;
	 TemplateReactionGenerator trg = (TemplateReactionGenerator)reactionSystem.reactionGenerator;
	 Iterator iter = trg.getReactionTemplateLibrary().getReactionTemplate();
	 while (iter.hasNext()){
	 ReactionTemplate rt = (ReactionTemplate)iter.next();
	 totalnum += rt.getNumberOfReactions();
	 }
	 
	 return totalnum;
	 }*/
	
	private void writeEnlargerInfo() {
		try {
	        File diagnosis = new File("enlarger.xls");
	        FileWriter fw = new FileWriter(diagnosis);
	        fw.write(Global.enlargerInfo.toString());
	        fw.close();
		}
		catch (IOException e) {
        	System.out.println("Cannot write enlarger file");
        	System.exit(0);
        }
	}
	
	private void writeDiagnosticInfo() {
		
		try {
	        File diagnosis = new File("diagnosis.xls");
	        FileWriter fw = new FileWriter(diagnosis);
	        fw.write(Global.diagnosticInfo.toString());
	        fw.close();
		}
		catch (IOException e) {
        	System.out.println("Cannot write diagnosis file");
        	System.exit(0);
        }
		
	}
	
	//10/25/07 gmagoon: I don't think this is used, but I will update to use reactionSystem and reactionTime as parameter to access temperature; commented-out usage of writeRestartFile will need to be modified
	//Is still incomplete.
    public void writeRestartFile(ReactionSystem p_rs, ReactionTime p_time ) {
		//writeCoreSpecies(p_rs);
		//writeCoreReactions(p_rs, p_time);
		//writeEdgeSpecies();
		//writeAllReactions(p_rs, p_time);
		//writeEdgeReactions(p_rs, p_time);
		
		//String restartFileName;
		//String restartFileContent="";
		
		/*File restartFile;
		 try {
		 restartFileName="Restart.txt";
		 restartFile = new File(restartFileName);
		 FileWriter fw = new FileWriter(restartFile);
		 restartFileContent = restartFileContent+ "TemperatureModel: Constant " +	reactionSystem.temperatureModel.getTemperature(reactionSystem.getInitialReactionTime()).getStandard()+ " (K)\n";
		 restartFileContent = restartFileContent+ "PressureModel: Constant " + reactionSystem.pressureModel.getPressure(reactionSystem.getInitialReactionTime()).getAtm() + " (atm)\n\n";
		 restartFileContent = restartFileContent+ "InitialStatus: \n";
		 int speciesCount = 1;
		 for(Iterator iter=reactionSystem.reactionModel.getSpecies();iter.hasNext();){
		 Species species = (Species) iter.next();
		 restartFileContent = restartFileContent + "("+ speciesCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
		 restartFileContent = restartFileContent + species.toStringWithoutH(1) + "\n";
		 speciesCount = speciesCount + 1;
		 }
		 restartFileContent = restartFileContent + "\n\n END \n\nInertGas:\n";
		 for (Iterator iter=reactionSystem.getPresentStatus().getInertGas(); iter.hasNext();){
		 String inertName= (String) iter.next();
		 restartFileContent = restartFileContent + inertName + " " + reactionSystem.getPresentStatus().getInertGas(inertName) + " mol/cm3 \n";
		 }
		 restartFileContent = restartFileContent + "END\n";
		 double tolerance;
		 if (reactionSystem.reactionModelEnlarger instanceof RateBasedPDepRME){
		 restartFileContent = restartFileContent + "ReactionModelEnlarger: RateBasedPDepModelEnlarger\n";
		 restartFileContent = restartFileContent + "FinishController: RateBasedPDepFinishController\n";
		 }
		 else {
		 restartFileContent = restartFileContent + "ReactionModelEnlarger: RateBasedModelEnlarger\n";
		 restartFileContent = restartFileContent + "FinishController: RateBasedFinishController\n";
		 }
		 if (reactionSystem.finishController.terminationTester instanceof ConversionTT){
		 restartFileContent = restartFileContent + "(1) Goal Conversion: ";
		 for (Iterator iter = ((ConversionTT)reactionSystem.finishController.terminationTester).getSpeciesGoalConversionSet();iter.hasNext();){
		 SpeciesConversion sc = (SpeciesConversion) iter.next();
		 restartFileContent = restartFileContent + sc.getSpecies().getName()+" "+sc.conversion + "\n";
		 
		 }
		 
		 }
		 else {
		 restartFileContent = restartFileContent + "(1) Goal ReactionTime: "+((ReactionTimeTT)reactionSystem.finishController.terminationTester).finalTime.toString();
		 }
		 
		 restartFileContent = restartFileContent + "Error Tolerance: " +((RateBasedVT) reactionSystem.finishController.validityTester).tolerance + "\n";
		 fw.write(restartFileContent);
		 fw.close();
		 }
		 catch (IOException e) {
		 System.out.println("Could not write the restart file");
		 System.exit(0);
		 }*/
		
		
	}
	
	/*Only write the forward reactions in the model core.
	 The reverse reactions are generated from the forward reactions.*/
	//10/25/07 gmagoon: added reaction system and reaction time as parameters and eliminated use of Global.temperature
	private void writeEdgeReactions(ReactionSystem p_rs, ReactionTime p_time) {
		StringBuilder restartFileContent =new StringBuilder();
		int reactionCount = 1;
		try{
			File coreSpecies = new File ("Restart/edgeReactions.txt");
			FileWriter fw = new FileWriter(coreSpecies);
			for(Iterator iter=((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSet().iterator();iter.hasNext();){
				
				Reaction reaction = (Reaction) iter.next();
				//if (reaction.getDirection()==1){
				//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
				restartFileContent = restartFileContent.append(reaction.toRestartString(p_rs.getTemperature(p_time)) + "\n");
				reactionCount = reactionCount + 1;
				//}
			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent.toString());
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart edgereactions file");
        	System.exit(0);
		}
		
	}
	
	//10/25/07 gmagoon: added reaction system and reaction time as parameters and eliminated use of Global.temperature
	private void writeAllReactions(ReactionSystem p_rs, ReactionTime p_time) {
		StringBuilder restartFileContent = new StringBuilder();
		int reactionCount = 1;
		try{
			File allReactions = new File ("Restart/allReactions.txt");
			FileWriter fw = new FileWriter(allReactions);
			for(Iterator iter=getReactionModel().getReaction();iter.hasNext();){
				
				Reaction reaction = (Reaction) iter.next();
				
				//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
				restartFileContent = restartFileContent.append(reaction.toRestartString(p_rs.getTemperature(p_time)) + "\n");
				
			}
			
			for(Iterator iter=((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSet().iterator();iter.hasNext();){
				
				Reaction reaction = (Reaction) iter.next();
				//if (reaction.getDirection()==1){
				//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
				restartFileContent = restartFileContent.append(reaction.toRestartString(p_rs.getTemperature(p_time)) + "\n");
				
			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent.toString());
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart edgereactions file");
        	System.exit(0);
		}
		
	}
	
	private void writeEdgeSpecies() {
		BufferedWriter bw = null;
		
        try {
            bw = new BufferedWriter(new FileWriter("Restart/edgeSpecies.txt"));
			for(Iterator iter=((CoreEdgeReactionModel)getReactionModel()).getUnreactedSpeciesSet().iterator();iter.hasNext();){
				Species species = (Species) iter.next();
				bw.write(species.getChemkinName());
				bw.newLine();
				int dummyInt = 0;
				bw.write(species.getChemGraph().toString(dummyInt));
				bw.newLine();
			}
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
	}
	
	//10/25/07 gmagoon: added reaction system and reaction time as parameters and eliminated use of Global.temperature
	private void writeCoreReactions(ReactionSystem p_rs, ReactionTime p_time) {
		StringBuilder restartFileContent = new StringBuilder();
		int reactionCount = 0;
		try{
			File coreSpecies = new File ("Restart/coreReactions.txt");
			FileWriter fw = new FileWriter(coreSpecies);
			for(Iterator iter=getReactionModel().getReaction();iter.hasNext();){
				
				Reaction reaction = (Reaction) iter.next();
				if (reaction.getDirection()==1){
					//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
					restartFileContent = restartFileContent.append(reaction.toRestartString(p_rs.getTemperature(p_time)) + "\n");
					reactionCount = reactionCount + 1;
				}
			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent.toString());
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart corereactions file");
        	System.exit(0);
		}
	}
	
	private void writeCoreSpecies() {
		BufferedWriter bw = null;
		
        try {
            bw = new BufferedWriter(new FileWriter("Restart/coreSpecies.txt"));
			for(Iterator iter=getReactionModel().getSpecies();iter.hasNext();){
				Species species = (Species) iter.next();
				bw.write(species.getChemkinName());
				bw.newLine();
				int dummyInt = 0;
				bw.write(species.getChemGraph().toString(dummyInt));
				bw.newLine();
			}
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
	}
	
	private void writeCoreReactions() {
		BufferedWriter bw_rxns = null;
		BufferedWriter bw_troe = null;
		BufferedWriter bw_lindemann = null;
		BufferedWriter bw_thirdbody = null;
		
        try {
            bw_rxns = new BufferedWriter(new FileWriter("Restart/coreReactions.txt"));
            bw_troe = new BufferedWriter(new FileWriter("Restart/troeReactions.txt"));
            bw_lindemann = new BufferedWriter(new FileWriter("Restart/lindemannReactions.txt"));
            bw_thirdbody = new BufferedWriter(new FileWriter("Restart/thirdBodyReactions.txt"));
            
    		String EaUnits = ArrheniusKinetics.getEaUnits();
    		String AUnits = ArrheniusKinetics.getAUnits();
    		bw_rxns.write("UnitsOfEa: " + EaUnits);
    		bw_rxns.newLine();
    		bw_troe.write("Unit:\nA: mol/cm3/s\nE: " + EaUnits + "\n\nReactions:");
    		bw_troe.newLine();
    		bw_lindemann.write("Unit:\nA: mol/cm3/s\nE: " + EaUnits + "\n\nReactions:");
    		bw_lindemann.newLine();
    		bw_thirdbody.write("Unit:\nA: mol/cm3/s\nE: " + EaUnits + "\n\nReactions :");
    		bw_thirdbody.newLine();
            
			CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)getReactionModel();
			LinkedHashSet allcoreRxns = cerm.core.reaction;
			for(Iterator iter=allcoreRxns.iterator(); iter.hasNext();){
				Reaction reaction = (Reaction) iter.next();
				if (reaction.isForward()) {
					if (reaction instanceof TROEReaction) {
						TROEReaction troeRxn = (TROEReaction) reaction;
						bw_troe.write(troeRxn.toRestartString(new Temperature(298,"K")));
						bw_troe.newLine();
					}
					else if (reaction instanceof LindemannReaction) {
						LindemannReaction lindeRxn = (LindemannReaction) reaction;
						bw_lindemann.write(lindeRxn.toRestartString(new Temperature(298,"K")));
						bw_lindemann.newLine();
					}
					else if (reaction instanceof ThirdBodyReaction) {
						ThirdBodyReaction tbRxn = (ThirdBodyReaction) reaction;
						bw_thirdbody.write(tbRxn.toRestartString(new Temperature(298,"K")));
						bw_thirdbody.newLine();
					}
					else {
						//bw.write(reaction.toChemkinString(new Temperature(298,"K")));
						bw_rxns.write(reaction.toRestartString(new Temperature(298,"K")));
						bw_rxns.newLine();
					}
				}
			}
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (bw_rxns != null) {
                    bw_rxns.flush();
                    bw_rxns.close();
                }
                if (bw_troe != null) {
                    bw_troe.flush();
                    bw_troe.close();
                }
                if (bw_lindemann != null) {
                    bw_lindemann.flush();
                    bw_lindemann.close();
                }
                if (bw_thirdbody != null) {
                    bw_thirdbody.flush();
                    bw_thirdbody.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
	}
	
	private void writeEdgeReactions() {
		BufferedWriter bw = null;
		
        try {
            bw = new BufferedWriter(new FileWriter("Restart/edgeReactions.txt"));
            
    		String EaUnits = ArrheniusKinetics.getEaUnits();
    		bw.write("UnitsOfEa: " + EaUnits);
    		bw.newLine();
            
			CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)getReactionModel();
			LinkedHashSet alledgeRxns = cerm.edge.reaction;
			for(Iterator iter=alledgeRxns.iterator(); iter.hasNext();){
				Reaction reaction = (Reaction) iter.next();
				if (reaction.isForward()) {
					//bw.write(reaction.toChemkinString(new Temperature(298,"K")));
					bw.write(reaction.toRestartString(new Temperature(298,"K")));
					bw.newLine();
				} else if (reaction.getReverseReaction().isForward()) {
					//bw.write(reaction.getReverseReaction().toChemkinString(new Temperature(298,"K")));
					bw.write(reaction.getReverseReaction().toRestartString(new Temperature(298,"K")));
					bw.newLine();
				} else
					System.out.println("Could not determine forward direction for following rxn: " + reaction.toString());
			}
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
	}
	
	private void writePDepNetworks() {
		BufferedWriter bw = null;
		
        try {
            bw = new BufferedWriter(new FileWriter("Restart/pdepnetworks.txt"));
    		int numFameTemps = PDepRateConstant.getTemperatures().length;
    		int numFamePress = PDepRateConstant.getPressures().length;
    		int numChebyTemps = ChebyshevPolynomials.getNT();
    		int numChebyPress = ChebyshevPolynomials.getNP();
    		int numPlog = PDepArrheniusKinetics.getNumPressures();
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
    		bw.write("NumberOfPLogs: " + numPlog);
    		bw.newLine();
    		bw.newLine();
    		
    		LinkedList allNets = PDepNetwork.getNetworks();
    		
			int netCounter = 0;
			for(Iterator iter=allNets.iterator(); iter.hasNext();){
				PDepNetwork pdepnet = (PDepNetwork) iter.next();
				
				++netCounter;
				bw.write("PDepNetwork #" + netCounter);
				bw.newLine();
				
				// Write netReactionList
				LinkedList netRxns = pdepnet.getNetReactions();
				bw.write("netReactionList:");
				bw.newLine();
				for (Iterator iter2=netRxns.iterator(); iter2.hasNext();) {
					
					PDepReaction currentPDepRxn = (PDepReaction)iter2.next();
					bw.write(currentPDepRxn.toString());
					bw.newLine();
					bw.write(writeRatesAndParameters(currentPDepRxn,numFameTemps,
													 numFamePress,numChebyTemps,numChebyPress,numPlog));
					
					PDepReaction currentPDepReverseRxn = currentPDepRxn.getReverseReaction();
					bw.write(currentPDepReverseRxn.toString());
					bw.newLine();
					bw.write(writeRatesAndParameters(currentPDepReverseRxn,numFameTemps,
													 numFamePress,numChebyTemps,numChebyPress,numPlog));
					
				}
				
				// Write nonincludedReactionList
				LinkedList nonIncludeRxns = pdepnet.getNonincludedReactions();
				bw.write("nonIncludedReactionList:");
				bw.newLine();
				for (Iterator iter2=nonIncludeRxns.iterator(); iter2.hasNext();) {
					
					PDepReaction currentPDepRxn = (PDepReaction)iter2.next();
					bw.write(currentPDepRxn.toString());
					bw.newLine();
					bw.write(writeRatesAndParameters(currentPDepRxn,numFameTemps,
													 numFamePress,numChebyTemps,numChebyPress,numPlog));
					
					PDepReaction currentPDepReverseRxn = currentPDepRxn.getReverseReaction();
					bw.write(currentPDepReverseRxn.toString());
					bw.newLine();
					bw.write(writeRatesAndParameters(currentPDepReverseRxn,numFameTemps,
													 numFamePress,numChebyTemps,numChebyPress,numPlog));
				}
				
				// Write pathReactionList
				LinkedList pathRxns = pdepnet.getPathReactions();
				bw.write("pathReactionList:");
				bw.newLine();
				for (Iterator iter2=pathRxns.iterator(); iter2.hasNext();) {
					PDepReaction currentPDepRxn = (PDepReaction)iter2.next();
					bw.write(currentPDepRxn.getDirection() + "\t" + currentPDepRxn.toChemkinString(new Temperature(298,"K")));
					bw.newLine();
				}
				
				bw.newLine();
				bw.newLine();
			}
    		
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (bw != null) {
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
	}
	
	public String writeRatesAndParameters(PDepReaction pdeprxn, int numFameTemps,
										  int numFamePress, int numChebyTemps, int numChebyPress, int numPlog) {
		StringBuilder sb = new StringBuilder();
		
		// Write the rate coefficients
		double[][] rateConstants = pdeprxn.getPDepRate().getRateConstants();
		for (int i=0; i<numFameTemps; i++) {
			for (int j=0; j<numFamePress; j++) {
				sb.append(rateConstants[i][j] + "\t");
			}
			sb.append("\n");
		}
		sb.append("\n");
		
		// If chebyshev polynomials are present, write them
		if (numChebyTemps != 0) {
			ChebyshevPolynomials chebyPolys = pdeprxn.getPDepRate().getChebyshev();
			for (int i=0; i<numChebyTemps; i++) {
				for (int j=0; j<numChebyPress; j++) {
					sb.append(chebyPolys.getAlpha(i,j) + "\t");
				}
				sb.append("\n");
			}
			sb.append("\n");
		}
		
		// If plog parameters are present, write them
		else if (numPlog != 0) {
			PDepArrheniusKinetics kinetics = pdeprxn.getPDepRate().getPDepArrheniusKinetics();
			for (int i=0; i<numPlog; i++) {
				double Hrxn = pdeprxn.calculateHrxn(new Temperature(298,"K"));
				sb.append(kinetics.pressures[i].getPa() + "\t" + kinetics.getKinetics(i).toChemkinString(Hrxn,new Temperature(298,"K"),false) + "\n");
			}
			sb.append("\n");
		}
		
		return sb.toString();
	}
	
	public LinkedList getTimeStep() {
        return timeStep;
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
	
    //svp
	public boolean getError(){
		return error;
	}
	
	//svp
	public boolean getSensitivity(){
		return sensitivity;
	}
	
	
	public LinkedList getSpeciesList() {
		return species;
	}
	
	//gmagoon 10/24/07: commented out getReactionSystem and setReactionSystem
	//    public ReactionSystem getReactionSystem() {
	//        return reactionSystem;
	//    }
	
	//11/2/07 gmagoon: adding accessor method for reactionSystemList
	public LinkedList getReactionSystemList(){
		return reactionSystemList;
	}
	
    //added by gmagoon 9/24/07
	//    public void setReactionSystem(ReactionSystem p_ReactionSystem) {
	//        reactionSystem = p_ReactionSystem;
	//    }
    
    //copied from ReactionSystem.java by gmagoon 9/24/07
    public ReactionModel getReactionModel() {
        return reactionModel;
    }
    
    public void readRestartSpecies() {    	
		// Read in core species -- NOTE code is almost duplicated in Read in edge species (second part of procedure)
		try {
			FileReader in = new FileReader("Restart/coreSpecies.txt");
			BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader);
			while (line != null) {
				// The first line of a new species is the user-defined name
				String totalSpeciesName = line;
				String[] splitString1 = totalSpeciesName.split("[(]");
				String[] splitString2 = splitString1[splitString1.length-1].split("[)]");
				// The remaining lines are the graph
				Graph g = ChemParser.readChemGraph(reader);
				// Make the ChemGraph, assuming it does not contain a forbidden structure
				ChemGraph cg = null;
				try {
					cg = ChemGraph.make(g);
				} catch (ForbiddenStructureException e) {
					System.out.println("Error in reading graph: Graph contains a forbidden structure.\n" + g.toString());
					System.exit(0);
				}
				// Make the species
				int intLocation = totalSpeciesName.indexOf("(" + splitString2[0] + ")");
				Species species = Species.make(totalSpeciesName.substring(0,intLocation),cg,Integer.parseInt(splitString2[0]));
				// Add the new species to the set of species
				restartCoreSpcs.add(species);
    			/*int species_type = 1; // reacted species
				 for (int i=0; i<numRxnSystems; i++) {
				 SpeciesStatus ss = new SpeciesStatus(species,species_type,y[i],yprime[i]);
				 speciesStatus[i].put(species, ss);
				 }*/
				line = ChemParser.readMeaningfulLine(reader);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// Read in edge species
		try {
			FileReader in = new FileReader("Restart/edgeSpecies.txt");
			BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader);
			while (line != null) {
				// The first line of a new species is the user-defined name
				String totalSpeciesName = line;
				String[] splitString1 = totalSpeciesName.split("[(]");
				String[] splitString2 = splitString1[splitString1.length-1].split("[)]"); // Change JDM to reflect MRH 2-11-2010
				// The remaining lines are the graph
				Graph g = ChemParser.readChemGraph(reader);
				// Make the ChemGraph, assuming it does not contain a forbidden structure
				ChemGraph cg = null;
				try {
					cg = ChemGraph.make(g);
				} catch (ForbiddenStructureException e) {
					System.out.println("Error in reading graph: Graph contains a forbidden structure.\n" + g.toString());
					System.exit(0);
				}
				// Make the species
				Species species = Species.make(splitString1[0],cg,Integer.parseInt(splitString2[0]));
				// Add the new species to the set of species
				restartEdgeSpcs.add(species);
				line = ChemParser.readMeaningfulLine(reader);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
    }
    
    public void readRestartReactions() {
    	// Grab the IDs from the core species
    	int[] coreSpcsIds = new int[restartCoreSpcs.size()]; 
    	int i = 0;
    	for (Iterator iter = restartCoreSpcs.iterator(); iter.hasNext();) {
    		Species spcs = (Species)iter.next();
    		coreSpcsIds[i] = spcs.getID();
    		++i;
    	}
    	
		// Read in core reactions
		try {
			FileReader in = new FileReader("Restart/coreReactions.txt");
			BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader);
            
            // Determine units of Ea
            StringTokenizer st = new StringTokenizer(line);
            String tempString = st.nextToken();
            String EaUnits = st.nextToken();
            
            line = ChemParser.readMeaningfulLine(reader);
            
			while (line != null) {
				if (!line.trim().equals("DUP")) {
					Reaction r = ChemParser.parseRestartReaction(line,coreSpcsIds,"core",EaUnits);
					
	        		Iterator rxnIter = restartCoreRxns.iterator();
	        		boolean foundRxn = false;
	        		while (rxnIter.hasNext()) {
	        			Reaction old = (Reaction)rxnIter.next();
	        			if (old.equals(r)) {
	        				old.addAdditionalKinetics(r.getKinetics(),1);
	        				foundRxn = true;
	        				break;
	        			}
	        		}
	        		if (!foundRxn) {
	        			if (r.hasReverseReaction()) r.generateReverseReaction();
	        			restartCoreRxns.add(r);
	        		}
				}
				
				line = ChemParser.readMeaningfulLine(reader);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			SeedMechanism.readThirdBodyReactions("Restart/thirdBodyReactions.txt");
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		try {
			SeedMechanism.readLindemannReactions("Restart/lindemannReactions.txt");
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		try {
			SeedMechanism.readTroeReactions("Restart/troeReactions.txt");
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		restartCoreRxns.addAll(SeedMechanism.reactionSet);
		
		// Read in edge reactions
		try {
			FileReader in = new FileReader("Restart/edgeReactions.txt");
			BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader);
            
            // Determine units of Ea
            StringTokenizer st = new StringTokenizer(line);
            String tempString = st.nextToken();
            String EaUnits = st.nextToken();
            
            line = ChemParser.readMeaningfulLine(reader);
            
			while (line != null) {
				if (!line.trim().equals("DUP")) {
					Reaction r = ChemParser.parseRestartReaction(line,coreSpcsIds,"edge",EaUnits);
					
	        		Iterator rxnIter = restartEdgeRxns.iterator();
	        		boolean foundRxn = false;
	        		while (rxnIter.hasNext()) {
	        			Reaction old = (Reaction)rxnIter.next();
	        			if (old.equals(r)) {
	        				old.addAdditionalKinetics(r.getKinetics(),1);
	        				foundRxn = true;
	        				break;
	        			}
	        		}
	        		if (!foundRxn) {
	        			r.generateReverseReaction();
	        			restartEdgeRxns.add(r);
	        		}
				}
				
				line = ChemParser.readMeaningfulLine(reader);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    public LinkedHashMap getRestartSpeciesStatus(int i) {
    	LinkedHashMap speciesStatus = new LinkedHashMap();
    	
		try {
			FileReader in = new FileReader("Restart/coreSpecies.txt");
			BufferedReader reader = new BufferedReader(in);
            Integer numRxnSystems = Integer.parseInt(ChemParser.readMeaningfulLine(reader));
            String line = ChemParser.readMeaningfulLine(reader);
			while (line != null) {
				// The first line of a new species is the user-defined name
				String totalSpeciesName = line;
				String[] splitString1 = totalSpeciesName.split("[(]");
				String[] splitString2 = splitString1[1].split("[)]");
				double y = 0.0;
				double yprime = 0.0;
				for (int j=0; j<numRxnSystems; j++) {
					StringTokenizer st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
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
					System.out.println("Error in reading graph: Graph contains a forbidden structure.\n" + g.toString());
					System.exit(0);
				}
				// Make the species
				Species species = Species.make(splitString1[0],cg);
				// Add the new species to the set of species
				//restartCoreSpcs.add(species);
    			int species_type = 1; // reacted species
    			SpeciesStatus ss = new SpeciesStatus(species,species_type,y,yprime);
    			speciesStatus.put(species, ss);
				line = ChemParser.readMeaningfulLine(reader);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return speciesStatus;
    }
    
    public void putRestartSpeciesInInitialStatus(InitialStatus is, int i) {
    	
		try {
			FileReader in = new FileReader("Restart/coreSpecies.txt");
			BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader);
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
					System.out.println("Error in reading graph: Graph contains a forbidden structure.\n" + g.toString());
					System.exit(0);
				}
				// Make the species
				Species species = Species.make(splitString1[0],cg);
				// Add the new species to the set of species
				//restartCoreSpcs.add(species);
				if (is.getSpeciesStatus(species) == null) {
	    			SpeciesStatus ss = new SpeciesStatus(species,1,0.0,0.0);
	    			is.putSpeciesStatus(ss);
				}
				line = ChemParser.readMeaningfulLine(reader);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
    }
    
    public void readPDepNetworks() {
    	SpeciesDictionary sd = SpeciesDictionary.getInstance();
    	LinkedList allNetworks = PDepNetwork.getNetworks(); 
    	
    	try {
			FileReader in = new FileReader("Restart/pdepnetworks.txt");
			BufferedReader reader = new BufferedReader(in);
            
			StringTokenizer st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
			String tempString = st.nextToken();
			String EaUnits = st.nextToken();
			st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
			tempString = st.nextToken();
			int numFameTs = Integer.parseInt(st.nextToken());
			st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
			tempString = st.nextToken();
			int numFamePs = Integer.parseInt(st.nextToken());
			st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
			tempString = st.nextToken();
			int numChebyTs = Integer.parseInt(st.nextToken());
			st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
			tempString = st.nextToken();
			int numChebyPs = Integer.parseInt(st.nextToken());
			st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
			tempString = st.nextToken();
			int numPlogs = Integer.parseInt(st.nextToken());
			
			double[][] rateCoefficients = new double[numFameTs][numFamePs];
			double[][] chebyPolys = new double[numChebyTs][numChebyPs];
			Kinetics[] plogKinetics = new Kinetics[numPlogs];
			
			String line = ChemParser.readMeaningfulLine(reader);	// line should be "PDepNetwork #"
			while (line != null) {
				line = ChemParser.readMeaningfulLine(reader);	// line should now be "netReactionList:"
				PDepNetwork newNetwork = new PDepNetwork();
				LinkedList netRxns = newNetwork.getNetReactions();
				LinkedList nonincludeRxns = newNetwork.getNonincludedReactions();
				
				line = ChemParser.readMeaningfulLine(reader);	// line is either data or "nonIncludedReactionList"
				
				// If line is "nonincludedreactionlist", we need to skip over this while loop
				if (!line.toLowerCase().startsWith("nonincludedreactionlist")) {
					while (!line.toLowerCase().startsWith("nonincludedreactionlist")) {
						
						// Read in the forward rxn
						String[] reactsANDprods = line.split("\\-->");
						
						PDepIsomer Reactants = null;
						String reacts = reactsANDprods[0].trim();
						if (reacts.contains("+")) {
							String[] indivReacts = reacts.split("[+]");
							String name = indivReacts[0].trim();
							Species spc1 = sd.getSpeciesFromChemkinName(name);
							if (spc1 == null) {
								spc1 = getSpeciesBySPCName(name,sd);
							}							
							name = indivReacts[1].trim();
							Species spc2 = sd.getSpeciesFromChemkinName(name);
							if (spc2 == null) {
								spc2 = getSpeciesBySPCName(name,sd);
							}	
							Reactants = new PDepIsomer(spc1,spc2);
						} else {
							String name = reacts.trim();
							Species spc = sd.getSpeciesFromChemkinName(name);
							if (spc == null) {
								spc = getSpeciesBySPCName(name,sd);
							}	
							Reactants = new PDepIsomer(spc);
						}
						
						PDepIsomer Products = null;
						String prods = reactsANDprods[1].trim();
						if (prods.contains("+")) {
							String[] indivProds = prods.split("[+]");
							String name = indivProds[0].trim();
							Species spc1 = sd.getSpeciesFromChemkinName(name);
							if (spc1 == null) {
								spc1 = getSpeciesBySPCName(name,sd);
							}	
							name = indivProds[1].trim();
							Species spc2 = sd.getSpeciesFromChemkinName(name);
							if (spc2 == null) {
								spc2 = getSpeciesBySPCName(name,sd);
							}	
							Products = new PDepIsomer(spc1,spc2);
						} else {
							String name = prods.trim();
							Species spc = sd.getSpeciesFromChemkinName(name);
							if (spc == null) {
								spc = getSpeciesBySPCName(name,sd);
							}	
							Products = new PDepIsomer(spc);
						}
						
						newNetwork.addIsomer(Reactants);
						newNetwork.addIsomer(Products);
						
						for (int i=0; i<numFameTs; i++) {
							st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
							for (int j=0; j<numFamePs; j++) {
								rateCoefficients[i][j] = Double.parseDouble(st.nextToken());
							}
						}
						
						PDepRateConstant pdepk = null;
						if (numChebyTs > 0) {
							for (int i=0; i<numChebyTs; i++) {
								st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
								for (int j=0; j<numChebyPs; j++) {
									chebyPolys[i][j] = Double.parseDouble(st.nextToken());
								}
							}
							ChebyshevPolynomials chebyshev = new ChebyshevPolynomials(numChebyTs,
																					  ChebyshevPolynomials.getTlow(), ChebyshevPolynomials.getTup(),
																					  numChebyPs, ChebyshevPolynomials.getPlow(), ChebyshevPolynomials.getPup(),
																					  chebyPolys);
							pdepk = new PDepRateConstant(rateCoefficients,chebyshev);
						} else if (numPlogs > 0) {
							for (int i=0; i<numPlogs; i++) {
								st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
								Pressure p = new Pressure(Double.parseDouble(st.nextToken()),"Pa");
								UncertainDouble dA = new UncertainDouble(Double.parseDouble(st.nextToken()),0.0,"A");
								UncertainDouble dn = new UncertainDouble(Double.parseDouble(st.nextToken()),0.0,"A");
								double Ea = Double.parseDouble(st.nextToken());
								if (EaUnits.equals("cal/mol"))
									Ea = Ea / 1000;
								else if (EaUnits.equals("J/mol"))
									Ea = Ea / 4.184 / 1000;
								else if (EaUnits.equals("kJ/mol"))
									Ea = Ea / 4.184;
								else if (EaUnits.equals("Kelvins"))
									Ea = Ea * 1.987;
								UncertainDouble dE = new UncertainDouble(Ea,0.0,"A");
								ArrheniusKinetics k = new ArrheniusKinetics(dA, dn, dE, "", 1, "", "");
								PDepArrheniusKinetics pdepAK = new PDepArrheniusKinetics(i);
								pdepAK.setKinetics(i, p, k);
								pdepk = new PDepRateConstant(rateCoefficients,pdepAK);
							}
						}
						
						PDepReaction forward = new PDepReaction(Reactants, Products, pdepk);
						
						// Read in the reverse reaction
						line = ChemParser.readMeaningfulLine(reader);
						
						for (int i=0; i<numFameTs; i++) {
							st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
							for (int j=0; j<numFamePs; j++) {
								rateCoefficients[i][j] = Double.parseDouble(st.nextToken());
							}
						}
						
						pdepk = null;
						if (numChebyTs > 0) {
							for (int i=0; i<numChebyTs; i++) {
								st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
								for (int j=0; j<numChebyPs; j++) {
									chebyPolys[i][j] = Double.parseDouble(st.nextToken());
								}
							}
							ChebyshevPolynomials chebyshev = new ChebyshevPolynomials(numChebyTs,
																					  ChebyshevPolynomials.getTlow(), ChebyshevPolynomials.getTup(),
																					  numChebyPs, ChebyshevPolynomials.getPlow(), ChebyshevPolynomials.getPup(),
																					  chebyPolys);
							pdepk = new PDepRateConstant(rateCoefficients,chebyshev);
						} else if (numPlogs > 0) {
							for (int i=0; i<numPlogs; i++) {
								st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
								Pressure p = new Pressure(Double.parseDouble(st.nextToken()),"Pa");
								UncertainDouble dA = new UncertainDouble(Double.parseDouble(st.nextToken()),0.0,"A");
								UncertainDouble dn = new UncertainDouble(Double.parseDouble(st.nextToken()),0.0,"A");
								double Ea = Double.parseDouble(st.nextToken());
								if (EaUnits.equals("cal/mol"))
									Ea = Ea / 1000;
								else if (EaUnits.equals("J/mol"))
									Ea = Ea / 4.184 / 1000;
								else if (EaUnits.equals("kJ/mol"))
									Ea = Ea / 4.184;
								else if (EaUnits.equals("Kelvins"))
									Ea = Ea * 1.987;
								UncertainDouble dE = new UncertainDouble(Ea,0.0,"A");
								ArrheniusKinetics k = new ArrheniusKinetics(dA, dn, dE, "", 1, "", "");
								PDepArrheniusKinetics pdepAK = new PDepArrheniusKinetics(i);
								pdepAK.setKinetics(i, p, k);
								pdepk = new PDepRateConstant(rateCoefficients,pdepAK);
							}
						}
						
						PDepReaction reverse = new PDepReaction(Products, Reactants, pdepk);
						reverse.setReverseReaction(forward);
						forward.setReverseReaction(reverse);
						
						netRxns.add(forward);
						
						line = ChemParser.readMeaningfulLine(reader);
					}
				}
				// This loop ends once line == "nonIncludedReactionList"
				
				line = ChemParser.readMeaningfulLine(reader);	// line is either data or "pathReactionList"
				
				if (!line.toLowerCase().startsWith("pathreactionList")) {
					
					while (!line.toLowerCase().startsWith("pathreactionlist")) {
						
						// Read in the forward rxn
						String[] reactsANDprods = line.split("\\-->");
						
						PDepIsomer Reactants = null;
						String reacts = reactsANDprods[0].trim();
						if (reacts.contains("+")) {
							String[] indivReacts = reacts.split("[+]");
							String name = indivReacts[0].trim();
							Species spc1 = sd.getSpeciesFromChemkinName(name);
							if (spc1 == null) {
								spc1 = getSpeciesBySPCName(name,sd);
							}	
							name = indivReacts[1].trim();
							Species spc2 = sd.getSpeciesFromChemkinName(name);
							if (spc2 == null) {
								spc2 = getSpeciesBySPCName(name,sd);
							}	
							Reactants = new PDepIsomer(spc1,spc2);
						} else {
							String name = reacts.trim();
							Species spc = sd.getSpeciesFromChemkinName(name);
							if (spc == null) {
								spc = getSpeciesBySPCName(name,sd);
							}	
							Reactants = new PDepIsomer(spc);
						}
						
						PDepIsomer Products = null;
						String prods = reactsANDprods[1].trim();
						if (prods.contains("+")) {
							String[] indivProds = prods.split("[+]");
							String name = indivProds[0].trim();
							Species spc1 = sd.getSpeciesFromChemkinName(name);
							if (spc1 == null) {
								spc1 = getSpeciesBySPCName(name,sd);
							}	
							name = indivProds[1].trim();
							Species spc2 = sd.getSpeciesFromChemkinName(name);
							if (spc2 == null) {
								spc2 = getSpeciesBySPCName(name,sd);
							}	
							Products = new PDepIsomer(spc1,spc2);
						} else {
							String name = prods.trim();
							Species spc = sd.getSpeciesFromChemkinName(name);
							if (spc == null) {
								spc = getSpeciesBySPCName(name,sd);
							}	
							Products = new PDepIsomer(spc);
						}
						
						newNetwork.addIsomer(Reactants);
						newNetwork.addIsomer(Products);
						
						for (int i=0; i<numFameTs; i++) {
							st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
							for (int j=0; j<numFamePs; j++) {
								rateCoefficients[i][j] = Double.parseDouble(st.nextToken());
							}
						}
						
						PDepRateConstant pdepk = null;
						if (numChebyTs > 0) {
							for (int i=0; i<numChebyTs; i++) {
								st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
								for (int j=0; j<numChebyPs; j++) {
									chebyPolys[i][j] = Double.parseDouble(st.nextToken());
								}
							}
							ChebyshevPolynomials chebyshev = new ChebyshevPolynomials(numChebyTs,
																					  ChebyshevPolynomials.getTlow(), ChebyshevPolynomials.getTup(),
																					  numChebyPs, ChebyshevPolynomials.getPlow(), ChebyshevPolynomials.getPup(),
																					  chebyPolys);
							pdepk = new PDepRateConstant(rateCoefficients,chebyshev);
						} else if (numPlogs > 0) {
							for (int i=0; i<numPlogs; i++) {
								st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
								Pressure p = new Pressure(Double.parseDouble(st.nextToken()),"Pa");
								UncertainDouble dA = new UncertainDouble(Double.parseDouble(st.nextToken()),0.0,"A");
								UncertainDouble dn = new UncertainDouble(Double.parseDouble(st.nextToken()),0.0,"A");
								double Ea = Double.parseDouble(st.nextToken());
								if (EaUnits.equals("cal/mol"))
									Ea = Ea / 1000;
								else if (EaUnits.equals("J/mol"))
									Ea = Ea / 4.184 / 1000;
								else if (EaUnits.equals("kJ/mol"))
									Ea = Ea / 4.184;
								else if (EaUnits.equals("Kelvins"))
									Ea = Ea * 1.987;
								UncertainDouble dE = new UncertainDouble(Ea,0.0,"A");
								ArrheniusKinetics k = new ArrheniusKinetics(dA, dn, dE, "", 1, "", "");
								PDepArrheniusKinetics pdepAK = new PDepArrheniusKinetics(i);
								pdepAK.setKinetics(i, p, k);
								pdepk = new PDepRateConstant(rateCoefficients,pdepAK);
							}
						}
						
						PDepReaction forward = new PDepReaction(Reactants, Products, pdepk);
						
						// Read in the reverse reaction
						line = ChemParser.readMeaningfulLine(reader);
						
						for (int i=0; i<numFameTs; i++) {
							st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
							for (int j=0; j<numFamePs; j++) {
								rateCoefficients[i][j] = Double.parseDouble(st.nextToken());
							}
						}
						
						pdepk = null;
						if (numChebyTs > 0) {
							for (int i=0; i<numChebyTs; i++) {
								st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
								for (int j=0; j<numChebyPs; j++) {
									chebyPolys[i][j] = Double.parseDouble(st.nextToken());
								}
							}
							ChebyshevPolynomials chebyshev = new ChebyshevPolynomials(numChebyTs,
																					  ChebyshevPolynomials.getTlow(), ChebyshevPolynomials.getTup(),
																					  numChebyPs, ChebyshevPolynomials.getPlow(), ChebyshevPolynomials.getPup(),
																					  chebyPolys);
							pdepk = new PDepRateConstant(rateCoefficients,chebyshev);
						} else if (numPlogs > 0) {
							for (int i=0; i<numPlogs; i++) {
								st = new StringTokenizer(ChemParser.readMeaningfulLine(reader));
								Pressure p = new Pressure(Double.parseDouble(st.nextToken()),"Pa");
								UncertainDouble dA = new UncertainDouble(Double.parseDouble(st.nextToken()),0.0,"A");
								UncertainDouble dn = new UncertainDouble(Double.parseDouble(st.nextToken()),0.0,"A");
								double Ea = Double.parseDouble(st.nextToken());
								if (EaUnits.equals("cal/mol"))
									Ea = Ea / 1000;
								else if (EaUnits.equals("J/mol"))
									Ea = Ea / 4.184 / 1000;
								else if (EaUnits.equals("kJ/mol"))
									Ea = Ea / 4.184;
								else if (EaUnits.equals("Kelvins"))
									Ea = Ea * 1.987;
								UncertainDouble dE = new UncertainDouble(Ea,0.0,"A");
								ArrheniusKinetics k = new ArrheniusKinetics(dA, dn, dE, "", 1, "", "");
								PDepArrheniusKinetics pdepAK = new PDepArrheniusKinetics(i);
								pdepAK.setKinetics(i, p, k);
								pdepk = new PDepRateConstant(rateCoefficients,pdepAK);
							}
						}
						
						PDepReaction reverse = new PDepReaction(Products, Reactants, pdepk);
						reverse.setReverseReaction(forward);
						forward.setReverseReaction(reverse);
						
						nonincludeRxns.add(forward);
						
						line = ChemParser.readMeaningfulLine(reader);
					}
				}
				// This loop ends once line == "pathReactionList"
				
				line = ChemParser.readMeaningfulLine(reader);	// line is either data or "PDepNetwork #_" or null (end of file)
				
				while (line != null && !line.toLowerCase().startsWith("pdepnetwork")) {
					
					st = new StringTokenizer(line);
					int direction = Integer.parseInt(st.nextToken());
					
					//	First token is the rxn structure: A+B=C+D
			    	//		Note: Up to 3 reactants/products allowed
			    	//			: Either "=" or "=>" will separate reactants and products
			    	String structure = st.nextToken();
			    	
			    	//	Separate the reactants from the products
			    	boolean generateReverse = false;
			    	String[] reactsANDprods = structure.split("\\=>");
			    	if (reactsANDprods.length == 1) {
			    		reactsANDprods = structure.split("[=]");
			    		generateReverse = true;
			    	}
			    	
			    	sd = SpeciesDictionary.getInstance();
			        LinkedList r = ChemParser.parseReactionSpecies(sd, reactsANDprods[0]);
			        LinkedList p = ChemParser.parseReactionSpecies(sd, reactsANDprods[1]);
					
			        Structure s = new Structure(r,p);
			        s.setDirection(direction);
					
			        //	Next three tokens are the modified Arrhenius parameters
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
			    	
			    	UncertainDouble uA = new UncertainDouble(rxn_A,0.0,"A");
			    	UncertainDouble un = new UncertainDouble(rxn_n,0.0,"A");
			    	UncertainDouble uE = new UncertainDouble(rxn_E,0.0,"A");    	
			    	
			    	//	The remaining tokens are comments
			    	String comments = "";
			    	while (st.hasMoreTokens()) {
			    		comments += st.nextToken();
			    	}
			    	
			        ArrheniusKinetics k = new ArrheniusKinetics(uA,un,uE,"",1,"",comments);
			        Reaction pathRxn = new Reaction();
					//				        if (direction == 1)
					//				        	pathRxn = Reaction.makeReaction(s,k,generateReverse);
					//				        else
					//				        	pathRxn = Reaction.makeReaction(s.generateReverseStructure(),k,generateReverse);
			        pathRxn = Reaction.makeReaction(s,k,generateReverse);
			        
			        PDepIsomer Reactants = new PDepIsomer(r);
			        PDepIsomer Products = new PDepIsomer(p);
			        PDepReaction pdeppathrxn = new PDepReaction(Reactants,Products,pathRxn);
			        newNetwork.addReaction(pdeppathrxn);
					
					line = ChemParser.readMeaningfulLine(reader);					
				}
				
				PDepNetwork.getNetworks().add(newNetwork);
				
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    /**
     * MRH 14Jan2010
     * 
     * getSpeciesBySPCName
     * 
     * Input:	String name - Name of species, normally chemical formula followed
     * 				by "J"s for radicals, and then (#)
     * 			SpeciesDictionary sd
     * 
     * This method was originally written as a complement to the method readPDepNetworks.
     * jdmo found a bug with the readrestart option.  The bug was that the method was
     * attempting to add a null species to the Isomer list.  The null species resulted
     * from searching the SpeciesDictionary by chemkinName (e.g. C4H8OJJ(48)), when the
     * chemkinName present in the dictionary was SPC(48).
     * 
	 */
    public Species getSpeciesBySPCName(String name, SpeciesDictionary sd) {
		String[] nameFromNumber = name.split("\\(");
		String newName = "SPC(" + nameFromNumber[1];
		return sd.getSpeciesFromChemkinName(newName);
    }
	
    /**
     * MRH 12-Jun-2009
     * 
     * Function initializes the model's core and edge.
     * 	The initial core species always consists of the species contained
     * 		in the condition.txt file.  If seed mechanisms exist, those species
     * 		(and the reactions given in the seed mechanism) are also added to
     * 		the core.
     * 	The initial edge species/reactions are determined by reacting the core
     * 		species by one full iteration.
     */
    public void initializeCoreEdgeModel() {
    	LinkedHashSet allInitialCoreSpecies = new LinkedHashSet();
    	LinkedHashSet allInitialCoreRxns = new LinkedHashSet();
    	
    	if (readrestart) {
    		readRestartReactions();
    		if (PDepNetwork.generateNetworks) readPDepNetworks();
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
    	
		CoreEdgeReactionModel cerm = new CoreEdgeReactionModel(allInitialCoreSpecies, allInitialCoreRxns);
		if (readrestart) {
			cerm.addUnreactedSpeciesSet(restartEdgeSpcs);
			cerm.addUnreactedReactionSet(restartEdgeRxns);
		}
		setReactionModel(cerm);
		
		PDepNetwork.reactionModel = getReactionModel();
		PDepNetwork.reactionSystem = (ReactionSystem) getReactionSystemList().get(0);
		
		// Determine initial set of reactions and edge species using only the
		// species enumerated in the input file and the seed mechanisms as the core
		if (!readrestart) {
			LinkedHashSet reactionSet;
			if (hasSeedMechanisms() && getSeedMechanism().shouldGenerateReactions()) {
				reactionSet = getReactionGenerator().react(allInitialCoreSpecies);
			}
			else {
				reactionSet = new LinkedHashSet();
				for (Iterator iter = speciesSeed.iterator(); iter.hasNext(); ) {
					Species spec = (Species) iter.next();
					reactionSet.addAll(getReactionGenerator().react(allInitialCoreSpecies, spec));
				}
			}
			reactionSet.addAll(getLibraryReactionGenerator().react(allInitialCoreSpecies));
			
	    	// Set initial core-edge reaction model based on above results
			if (reactionModelEnlarger instanceof RateBasedRME)	{
				Iterator iter = reactionSet.iterator();
	        	while (iter.hasNext()){
	        		Reaction r = (Reaction)iter.next();
	        		cerm.addReaction(r);
				}
			}
			else {
	    		// Only keep the reactions involving bimolecular reactants and bimolecular products
				Iterator iter = reactionSet.iterator();
	        	while (iter.hasNext()){
	        		Reaction r = (Reaction)iter.next();
	        		if (r.getReactantNumber() > 1 && r.getProductNumber() > 1){
	        			cerm.addReaction(r);
	        		}
					else {
						cerm.categorizeReaction(r.getStructure());
						PDepNetwork.addReactionToNetworks(r);
					}
				}
			}
		}
		
        for (Integer i = 0; i < reactionSystemList.size(); i++) {
			ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
			rs.setReactionModel(getReactionModel());
        }
		
		// We cannot return a system with no core reactions, so if this is a case we must add to the core
        while (getReactionModel().isEmpty() && !PDepNetwork.hasCoreReactions((CoreEdgeReactionModel) getReactionModel())) {
			for (Integer i = 0; i < reactionSystemList.size(); i++) {
				ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
				if (reactionModelEnlarger instanceof RateBasedPDepRME) 
					rs.initializePDepNetwork();
				rs.appendUnreactedSpeciesStatus((InitialStatus)initialStatusList.get(i), rs.getPresentTemperature());
			}
			enlargeReactionModel();
		}
		
        for (Integer i = 0; i<reactionSystemList.size();i++) {
            ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
            rs.setReactionModel(getReactionModel());
        }
		
        return;
    	
    }
    
    //## operation initializeCoreEdgeModelWithPRL()
    //9/24/07 gmagoon: moved from ReactionSystem.java
    public void initializeCoreEdgeModelWithPRL() {
        //#[ operation initializeCoreEdgeModelWithPRL()
        initializeCoreEdgeModelWithoutPRL();
		
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)getReactionModel();
		
        LinkedHashSet primarySpeciesSet = getPrimaryReactionLibrary().getSpeciesSet(); //10/14/07 gmagoon: changed to use getPrimaryReactionLibrary
        LinkedHashSet primaryReactionSet = getPrimaryReactionLibrary().getReactionSet();
        cerm.addReactedSpeciesSet(primarySpeciesSet);
        cerm.addPrimaryReactionSet(primaryReactionSet);
		
        LinkedHashSet newReactions = getReactionGenerator().react(cerm.getReactedSpeciesSet());
        
        if (reactionModelEnlarger instanceof RateBasedRME)	
        	cerm.addReactionSet(newReactions);
    	else {
    		
        	Iterator iter = newReactions.iterator();
        	while (iter.hasNext()){
        		Reaction r = (Reaction)iter.next();
        		if (r.getReactantNumber() == 2 && r.getProductNumber() == 2){
        			cerm.addReaction(r);
        		}
        	}
    	}
        
		
		
		
        
		
        return;
		
        //#]
    }
	
    //## operation initializeCoreEdgeModelWithoutPRL()
    //9/24/07 gmagoon: moved from ReactionSystem.java
    protected void initializeCoreEdgeModelWithoutPRL() {
        //#[ operation initializeCoreEdgeModelWithoutPRL()
		
		CoreEdgeReactionModel cerm = new CoreEdgeReactionModel(new LinkedHashSet(getSpeciesSeed()));
		setReactionModel(cerm);
		
		PDepNetwork.reactionModel = getReactionModel();
		PDepNetwork.reactionSystem = (ReactionSystem) getReactionSystemList().get(0);
		
		// Determine initial set of reactions and edge species using only the
		// species enumerated in the input file as the core
		LinkedHashSet reactionSet = getReactionGenerator().react(getSpeciesSeed());
		reactionSet.addAll(getLibraryReactionGenerator().react(getSpeciesSeed()));
		
    	// Set initial core-edge reaction model based on above results
		if (reactionModelEnlarger instanceof RateBasedRME)	{
			// Only keep the reactions involving bimolecular reactants and bimolecular products
			Iterator iter = reactionSet.iterator();
        	while (iter.hasNext()){
        		Reaction r = (Reaction)iter.next();
        		cerm.addReaction(r);
			}
		}
		else {
    		// Only keep the reactions involving bimolecular reactants and bimolecular products
			Iterator iter = reactionSet.iterator();
        	while (iter.hasNext()){
        		Reaction r = (Reaction)iter.next();
        		if (r.getReactantNumber() > 1 && r.getProductNumber() > 1){
        			cerm.addReaction(r);
        		}
				else {
					cerm.categorizeReaction(r.getStructure());
					PDepNetwork.addReactionToNetworks(r);
				}
			}
		}
        
		//10/9/07 gmagoon: copy reactionModel to reactionSystem; there may still be scope problems, particularly in above elseif statement
        //10/24/07 gmagoon: want to copy same reaction model to all reactionSystem variables; should probably also make similar modifications elsewhere; may or may not need to copy in ...WithPRL function
        for (Integer i = 0; i < reactionSystemList.size(); i++) {
			ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
			rs.setReactionModel(getReactionModel());
        }
        //reactionSystem.setReactionModel(getReactionModel());
		
		// We cannot return a system with no core reactions, so if this is a case we must add to the core
        while (getReactionModel().isEmpty()&&!PDepNetwork.hasCoreReactions((CoreEdgeReactionModel) getReactionModel())) {
			for (Integer i = 0; i < reactionSystemList.size(); i++) {
				ReactionSystem rs = (ReactionSystem) reactionSystemList.get(i);
				if (reactionModelEnlarger instanceof RateBasedPDepRME) 
					rs.initializePDepNetwork();
				rs.appendUnreactedSpeciesStatus((InitialStatus)initialStatusList.get(i), rs.getPresentTemperature());
			}
			enlargeReactionModel();
            
		}
		
        for (Integer i = 0; i<reactionSystemList.size();i++) {
            ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
            rs.setReactionModel(getReactionModel());
        }
		
        return;
		
        //#]
    }
	
    //## operation initializeCoreEdgeReactionModel()
    //9/24/07 gmagoon: moved from ReactionSystem.java
    public void initializeCoreEdgeReactionModel() {
		System.out.println("\nInitializing core-edge reaction model");
		// setSpeciesSeed(new LinkedHashSet());//10/4/07 gmagoon:moved from initializeReactionSystem; later moved to modelGeneration()
        //#[ operation initializeCoreEdgeReactionModel()
		//        if (hasPrimaryReactionLibrary()) initializeCoreEdgeModelWithPRL();
		//        else initializeCoreEdgeModelWithoutPRL();
		/*
		 * MRH 12-Jun-2009
		 * 
		 * I've lumped the initializeCoreEdgeModel w/ and w/o a seed mechanism
		 * 	(which used to be the PRL) into one function.  Before, RMG would
		 * 	complete one iteration (construct the edge species/rxns) before adding
		 * 	the seed mechanism to the rxn, thereby possibly estimating kinetic
		 * 	parameters for a rxn that exists in a seed mechanism
		 */
		initializeCoreEdgeModel();
		
        //#]
    }
    
    //9/24/07 gmagoon: copied from ReactionSystem.java
    public ReactionGenerator getReactionGenerator() {
        return reactionGenerator;
    }
    
    //10/4/07 gmagoon: moved from ReactionSystem.java
    public void setReactionGenerator(ReactionGenerator p_ReactionGenerator) {
		reactionGenerator = p_ReactionGenerator;
    }
    
    //9/25/07 gmagoon: moved from ReactionSystem.java
    //10/24/07 gmagoon: changed to use reactionSystemList
    //## operation enlargeReactionModel()
    public void enlargeReactionModel() {
        //#[ operation enlargeReactionModel()
        if (reactionModelEnlarger == null) throw new NullPointerException("ReactionModelEnlarger");
		System.out.println("\nEnlarging reaction model");
        reactionModelEnlarger.enlargeReactionModel(reactionSystemList, reactionModel, validList);
		
        return;
        //#]
    }
    
    //9/25/07 gmagoon: moved from ReactionSystem.java
	//## operation hasPrimaryReactionLibrary()
    public boolean hasPrimaryReactionLibrary() {
        //#[ operation hasPrimaryReactionLibrary()
        if (primaryReactionLibrary == null) return false;
        return (primaryReactionLibrary.size() > 0);
        //#]
    }
    
    public boolean hasSeedMechanisms() {
    	if (getSeedMechanism() == null) return false;
    	return (seedMechanism.size() > 0);
    }
    
    //9/25/07 gmagoon: moved from ReactionSystem.java
    public PrimaryReactionLibrary getPrimaryReactionLibrary() {
        return primaryReactionLibrary;
    }
	
    //9/25/07 gmagoon: moved from ReactionSystem.java
    public void setPrimaryReactionLibrary(PrimaryReactionLibrary p_PrimaryReactionLibrary) {
        primaryReactionLibrary = p_PrimaryReactionLibrary;
    }
    
    //10/4/07 gmagoon: added
    public LinkedHashSet getSpeciesSeed() {
        return speciesSeed;
    }
	
    //10/4/07 gmagoon: added
    public void setSpeciesSeed(LinkedHashSet p_speciesSeed) {
        speciesSeed = p_speciesSeed;
    }
	
    //10/4/07 gmagoon: added
    public LibraryReactionGenerator getLibraryReactionGenerator() {
        return lrg;
    }
	
    //10/4/07 gmagoon: added
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
    
    public void setPrimaryThermoLibrary(PrimaryThermoLibrary p_primaryThermoLibrary) {
    	primaryThermoLibrary = p_primaryThermoLibrary;
    }
    
    public static double getAtol(){
    	return atol;
    }

    public boolean runKillableToPreventInfiniteLoop(boolean intermediateSteps, int iterationNumber) {
	ReactionSystem rs0 = (ReactionSystem)reactionSystemList.get(0);
	if (!intermediateSteps)//if there are no intermediate steps (for example when using AUTO method), return true;
	    return true;
	//if there are intermediate steps, the run is killable if the iteration number exceeds the number of time steps / conversions
	else if (rs0.finishController.terminationTester instanceof ReactionTimeTT){
	    if (iterationNumber - 1 > timeStep.size()){ //-1 correction needed since when this is called, iteration number has been incremented
		return true;
	    }
	}
	else //the case where intermediate conversions are specified
	    if (iterationNumber - 1 > numConversions){ //see above; it is possible there is an off-by-one error here, so further testing will be needed
		return true;
	    }
	return false; //return false if none of the above criteria are met
    }
    
    public void readAndMakePRL(BufferedReader reader) throws IOException {
    	int Ilib = 0;
    	String line = ChemParser.readMeaningfulLine(reader);
    	while (!line.equals("END")) {
			String[] tempString = line.split("Name: ");
			String name = tempString[tempString.length-1].trim();
			line = ChemParser.readMeaningfulLine(reader);
			tempString = line.split("Location: ");
			String path = tempString[tempString.length-1].trim();
			if (Ilib==0) {
				setPrimaryReactionLibrary(new PrimaryReactionLibrary(name, path));
				Ilib++; 	
			}
			else {
				getPrimaryReactionLibrary().appendPrimaryReactionLibrary(name,path);
				Ilib++;
			}
			line = ChemParser.readMeaningfulLine(reader);
		}
		if (Ilib==0) {
			setPrimaryReactionLibrary(null);
		}
		else System.out.println("Primary Reaction Libraries in use: " + getPrimaryReactionLibrary().getName());
    }
    
    public void readAndMakePTL(BufferedReader reader) {
     	int numPTLs = 0;
     	String line = ChemParser.readMeaningfulLine(reader);
     	while (!line.equals("END")) {
     		String[] tempString = line.split("Name: ");
     		String name = tempString[tempString.length-1].trim();
			line = ChemParser.readMeaningfulLine(reader);
			tempString = line.split("Location: ");
			String path = tempString[tempString.length-1].trim();
			if (numPTLs==0) {
             	setPrimaryThermoLibrary(new PrimaryThermoLibrary(name,path));
             	++numPTLs; 	
			}
			else {
             	getPrimaryThermoLibrary().appendPrimaryThermoLibrary(name,path);
             	++numPTLs;
			}
			line = ChemParser.readMeaningfulLine(reader);
     	}
     	if (numPTLs == 0) setPrimaryThermoLibrary(null);
    }
    
    public void setSpectroscopicDataMode(String line) {
		StringTokenizer st = new StringTokenizer(line);
		String name = st.nextToken();
		String sdeType = st.nextToken().toLowerCase();
		if (sdeType.equals("frequencygroups") || sdeType.equals("default")) {
			SpectroscopicData.mode = SpectroscopicData.Mode.FREQUENCYGROUPS;
		}
		else if (sdeType.equals("therfit") || sdeType.equals("threefrequencymodel")) {
			SpectroscopicData.mode = SpectroscopicData.Mode.THREEFREQUENCY;
		}
		else if (sdeType.equals("off") || sdeType.equals("none")) {
			SpectroscopicData.mode = SpectroscopicData.Mode.OFF;
		}
		else throw new InvalidSymbolException("condition.txt: Unknown SpectroscopicDataEstimator = " + sdeType);
    }
    
    public void setPressureDependenceType(String line) {
		StringTokenizer st = new StringTokenizer(line);
		String name = st.nextToken();
		String pDepType = st.nextToken();
		if (pDepType.toLowerCase().equals("modifiedstrongcollision") ||
			pDepType.toLowerCase().equals("reservoirstate") ||
			pDepType.toLowerCase().equals("chemdis")) {
			
			reactionModelEnlarger = new RateBasedPDepRME();
			PDepNetwork.generateNetworks = true;
			
			if (pDepType.toLowerCase().equals("reservoirstate")) {
				((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new FastMasterEqn(FastMasterEqn.Mode.RESERVOIRSTATE));
				if (SpectroscopicData.mode == SpectroscopicData.Mode.OFF) {
					System.out.println("Warning: Spectroscopic data needed for pressure dependence; switching SpectroscopicDataEstimator to FrequencyGroups.");
					SpectroscopicData.mode = SpectroscopicData.Mode.FREQUENCYGROUPS;
				}
			}
			else if (pDepType.toLowerCase().equals("modifiedstrongcollision")) {
				((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new FastMasterEqn(FastMasterEqn.Mode.STRONGCOLLISION));
				if (SpectroscopicData.mode == SpectroscopicData.Mode.OFF) {
					System.out.println("Warning: Spectroscopic data needed for pressure dependence; switching SpectroscopicDataEstimator to FrequencyGroups.");
					SpectroscopicData.mode = SpectroscopicData.Mode.FREQUENCYGROUPS;
				}
			}
			else if (pDepType.toLowerCase().equals("chemdis")) {
				((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new Chemdis());
				if (SpectroscopicData.mode != SpectroscopicData.Mode.THREEFREQUENCY) {
					System.out.println("Warning: Switching SpectroscopicDataEstimator to three-frequency model.");
					SpectroscopicData.mode = SpectroscopicData.Mode.THREEFREQUENCY;
				}
			}
			else {
				throw new InvalidSymbolException("condition.txt: Unknown PDepKineticsEstimator = " + pDepType);
			}
			
			// Set temperatures and pressures to use in PDep kinetics estimation
			Temperature[] temperatures = new Temperature[8];
			temperatures[0] = new Temperature(300, "K");
			temperatures[1] = new Temperature(400, "K");
			temperatures[2] = new Temperature(600, "K");
			temperatures[3] = new Temperature(900, "K");
			temperatures[4] = new Temperature(1200, "K");
			temperatures[5] = new Temperature(1500, "K");
			temperatures[6] = new Temperature(1800, "K");
			temperatures[7] = new Temperature(2100, "K");
			FastMasterEqn.setTemperatures(temperatures);
			PDepRateConstant.setTemperatures(temperatures);
			ChebyshevPolynomials.setTlow(temperatures[0]);
			ChebyshevPolynomials.setTup(temperatures[7]);
			
			Pressure[] pressures = new Pressure[5];
			pressures[0] = new Pressure(0.01, "bar");
			pressures[1] = new Pressure(0.1, "bar");
			pressures[2] = new Pressure(1, "bar");
			pressures[3] = new Pressure(10, "bar");
			pressures[4] = new Pressure(100, "bar");
			FastMasterEqn.setPressures(pressures);
			PDepRateConstant.setPressures(pressures);
			ChebyshevPolynomials.setPlow(pressures[0]);
			ChebyshevPolynomials.setPup(pressures[4]);
			
		}
		else if (pDepType.toLowerCase().equals("off")) {
			// No pressure dependence
			reactionModelEnlarger = new RateBasedRME();
			PDepNetwork.generateNetworks = false;
		}
		else {
			throw new InvalidSymbolException("condition.txt: Unknown PressureDependence = " + pDepType);
		}
    }
    
    public String setPDepKineticsModel(String line, BufferedReader reader) throws InvalidSymbolException {
		StringTokenizer st = new StringTokenizer(line);
		String name = st.nextToken();
		String pDepKinType = st.nextToken();
		if (pDepKinType.toLowerCase().equals("chebyshev") ||
			pDepKinType.toLowerCase().equals("pdeparrhenius") ||
			pDepKinType.toLowerCase().equals("rate")) {
			
			if (pDepKinType.toLowerCase().equals("chebyshev")) {
				PDepRateConstant.setMode(PDepRateConstant.Mode.CHEBYSHEV);
				line = ChemParser.readMeaningfulLine(reader);
				if (line.startsWith("TRange")) {
					st = new StringTokenizer(line);
					String temp = st.nextToken(); // Should be "TRange:"
					String TUNITS = ChemParser.removeBrace(st.nextToken());
					double tLow = Double.parseDouble(st.nextToken());
					Temperature TMIN = new Temperature(tLow,TUNITS);
					ChebyshevPolynomials.setTlow(TMIN);
					double tHigh = Double.parseDouble(st.nextToken());
					if (tHigh <= tLow) {
						throw new InvalidSymbolException("Chebyshev Tmax is less than or equal to Tmin");
					}
					Temperature TMAX = new Temperature(tHigh,TUNITS);
					ChebyshevPolynomials.setTup(TMAX);
					if (TMAX.getK() < Global.highTemperature.getK()) 
						throw new InvalidSymbolException("Chebyshev Tmax is lower than the highest simulation temperature "+Global.highTemperature.toString() );
					if (TMIN.getK() > Global.lowTemperature.getK()) 
						throw new InvalidSymbolException("Chebyshev Tmin is higher than the lowest simulation temperature "+Global.lowTemperature.toString() );		
					int tResolution = Integer.parseInt(st.nextToken());
					int tbasisFuncs = Integer.parseInt(st.nextToken());
					if (tbasisFuncs > tResolution) {
						throw new InvalidSymbolException("The number of basis functions cannot exceed the number of grid points");
					}
					FastMasterEqn.setNumTBasisFuncs(tbasisFuncs);
					line = ChemParser.readMeaningfulLine(reader);
					String PUNITS = "";
					Pressure PMIN = new Pressure();
					Pressure PMAX = new Pressure();
					int pResolution = 0;
					int pbasisFuncs = 0;
					if (line.startsWith("PRange")) {
						st = new StringTokenizer(line);
						temp = st.nextToken(); // Should be "PRange:"
						PUNITS = ChemParser.removeBrace(st.nextToken());
						double pLow = Double.parseDouble(st.nextToken());
						PMIN = new Pressure(pLow,PUNITS);
						ChebyshevPolynomials.setPlow(PMIN);
						double pHigh = Double.parseDouble(st.nextToken());
						if (pHigh <= pLow) {
							System.err.println("Pmax is less than or equal to Pmin");
							System.exit(0);
						}
						PMAX = new Pressure(pHigh,PUNITS);
						if (PMAX.getPa() < Global.highPressure.getPa()) 
							throw new InvalidSymbolException("Chebyshev Pmax is lower than the highest simulation pressure "+Global.highPressure.toString() );
						if (PMIN.getPa() > Global.lowPressure.getPa()) 
							throw new InvalidSymbolException("Chebyshev Pmin is higher than the lowest simulation pressure "+Global.lowPressure.toString() );							
						ChebyshevPolynomials.setPup(PMAX);
						pResolution = Integer.parseInt(st.nextToken());
						pbasisFuncs = Integer.parseInt(st.nextToken());
						if (pbasisFuncs > pResolution) {
							throw new InvalidSymbolException("The number of basis functions cannot exceed the number of grid points");
						}
						FastMasterEqn.setNumPBasisFuncs(pbasisFuncs);
					}
					else {
						throw new InvalidSymbolException("RMG cannot locate PRange field for Chebyshev polynomials.");
					}
					
					// Set temperatures and pressures to use in PDep kinetics estimation
					Temperature[] temperatures = new Temperature[tResolution];
					for (int i=0; i<temperatures.length; i++) {
						double tempValueTilda = Math.cos(i*Math.PI/(temperatures.length-1)); // roots of a Chebyshev polynomial
						double tempValue = 2 / (tempValueTilda * ((1/TMAX.getK()) - (1/TMIN.getK())) + (1/TMIN.getK()) + (1/TMAX.getK()));
						temperatures[temperatures.length-i-1] = new Temperature(tempValue,TUNITS);
					}
					FastMasterEqn.setTemperatures(temperatures);
					PDepRateConstant.setTemperatures(temperatures);
					
					Pressure[] pressures = new Pressure[pResolution];
					for (int j=0; j<pressures.length; j++) {
						double pressValueTilda = Math.cos(j*Math.PI/(pressures.length-1)); // roots of a Chebyshev polynomial
						double pressValue = Math.pow(10,(pressValueTilda*(Math.log10(PMAX.getBar())-Math.log10(PMIN.getBar()))+Math.log10(PMIN.getBar())+Math.log10(PMAX.getBar()))/2);
						pressures[pressures.length-j-1] = new Pressure(pressValue,PUNITS);
					}
					FastMasterEqn.setPressures(pressures);
					PDepRateConstant.setPressures(pressures);
					line = ChemParser.readMeaningfulLine(reader);
				}
			}
			else if (pDepKinType.toLowerCase().equals("pdeparrhenius")) {
				//PDepRateConstant.setMode(PDepRateConstant.Mode.PDEPARRHENIUS);
				/*
				 *  Updated by MRH on 10Feb2010
				 *  	Allow user to specify # of T's/P's solved for in fame &
				 *  	# of PLOG's to report
				 */
				PDepRateConstant.setMode(PDepRateConstant.Mode.PDEPARRHENIUS);
				line = ChemParser.readMeaningfulLine(reader);
				if (line.startsWith("TRange")) {
					st = new StringTokenizer(line);
					String temp = st.nextToken(); // Should be "TRange:"
					String TUNITS = ChemParser.removeBrace(st.nextToken());
					int numT = Integer.parseInt(st.nextToken());
					Temperature[] listOfTs = new Temperature[numT];
					int counter = 0;
					while (st.hasMoreTokens()) {
						listOfTs[counter] = new Temperature(Double.parseDouble(st.nextToken()),TUNITS);
						++counter;
					}
					if (counter != numT) {
						System.out.println("Warning in TRange field of PressureDependence:\n" +
								"The stated number of temperatures is: " + numT + 
								"but the length of the temperature list is: " + counter);
					}
					
					line = ChemParser.readMeaningfulLine(reader);
					String PUNITS = "";
					if (line.startsWith("PRange")) {
						st = new StringTokenizer(line);
						temp = st.nextToken(); // Should be "PRange:"
						PUNITS = ChemParser.removeBrace(st.nextToken());
						int numP = Integer.parseInt(st.nextToken());
						Pressure[] listOfPs = new Pressure[numP];
						counter = 0;
						while (st.hasMoreTokens()) {
							listOfPs[counter] = new Pressure(Double.parseDouble(st.nextToken()),PUNITS);
							++counter;
						}
						if (counter != numP) {
							System.out.println("Warning in PRange field of PressureDependence:\n" +
									"The stated number of pressures is: " + numT + 
									"but the length of the pressure list is: " + counter);
						}
						FastMasterEqn.setTemperatures(listOfTs);
						PDepRateConstant.setTemperatures(listOfTs);
						FastMasterEqn.setPressures(listOfPs);
						PDepRateConstant.setPressures(listOfPs);
						PDepArrheniusKinetics.setNumPressures(numP);
						PDepArrheniusKinetics.setPressures(listOfPs);
					}
					else {
						throw new InvalidSymbolException("RMG cannot locate PRange field for PDepArrhenius.");
					}
					
					line = ChemParser.readMeaningfulLine(reader);
				}
			}
			// 6Jul2009-MRH:
			//	RATE mode reports p-dep rxn kinetics as: A 0.0 0.0
			//		where A is k(T,P) evaluated at the single temperature
			//		and pressure given in the condition.txt file
			else if (pDepKinType.toLowerCase().equals("rate"))
				PDepRateConstant.setMode(PDepRateConstant.Mode.RATE);
			
		}
		else {
			throw new InvalidSymbolException("condition.txt: Unknown PDepKinetics = " + pDepKinType);
		}
		return line;
    }
    
    public ReactionModelEnlarger getReactionModelEnlarger() {
    	return reactionModelEnlarger;
    }
}
/*********************************************************************
 File Path	: RMG\RMG\jing\rxnSys\ReactionModelGenerator.java
 *********************************************************************/


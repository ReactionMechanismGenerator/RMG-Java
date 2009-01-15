//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//!
//!******************************************************************************



package jing.rxnSys;


import java.io.*;
import jing.rxnSys.ReactionSystem;
import jing.rxn.*;
import jing.chem.*;

import java.util.*;

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
    protected double atol;
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

	protected boolean restart = false;
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
              if (line.startsWith("PrimaryThermoLibrary")){//svp
                line = ChemParser.readMeaningfulLine(reader);
              }
              else throw new InvalidSymbolException("Can't find primary thermo library!");


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
                                presList.add(new ConstantPM(p, unit));
                                //read remaining temperatures
        			while (st.hasMoreTokens()) {
                                    p = Double.parseDouble(st.nextToken());
                                    presList.add(new ConstantPM(p, unit));
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
                

        	// read in spectroscopic data estimator
        //	line = ChemParser.readMeaningfulLine(reader);
        	if (line.startsWith("SpectroscopicDataEstimator:")) {
        		StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		String sdeType = st.nextToken().toLowerCase();
        		if (sdeType.equals("frequencygroups") || sdeType.equals("default")) {
        			SpectroscopicData.useThreeFrequencyModel = false;
        		}
				else if (sdeType.equals("therfit") || sdeType.equals("threefrequencymodel")) {
        			SpectroscopicData.useThreeFrequencyModel = true;
        		}
				else throw new InvalidSymbolException("condition.txt: Unknown SpectroscopicDataEstimator = " + sdeType);
				
        	}
        	else throw new InvalidSymbolException("condition.txt: can't find SpectroscopicDataEstimator!");
			
			// read in reactants
        	line = ChemParser.readMeaningfulLine(reader);

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
					speciesnum ++;
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
        			else if (!unit.equals("mole/cm3") && !unit.equals("mol/cm3")) {
        				throw new InvalidUnitException("Species Concentration in condition.txt!");
        			}

        			//GJB to allow "unreactive" species that only follow user-defined library reactions.  
        			// They will not react according to RMG reaction families 
					boolean IsReactive = true;
					if (st.hasMoreTokens()) {
						String reactive = st.nextToken().trim();
						if (reactive.equalsIgnoreCase("unreactive"))
							IsReactive = false;
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
                                       // LinkedHashMap speStat = (LinkedHashMap)speciesStatus.clone();//10/31/07 gmagoon: trying creating multiple instances of speciesStatus to address issues with concentration normalization (last normalization seems to apply to all)
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
        			else if (!unit.equals("mole/cm3") && !unit.equals("mol/cm3")) {
        				throw new InvalidUnitException("Species Concentration in condition.txt!");
        			}
        			SystemSnapshot.putInertGas(name,inertConc);
        	   		line = ChemParser.readMeaningfulLine(reader);
        		}
           	}
        	else throw new InvalidSymbolException("condition.txt: can't find Inert gas concentration!");

        	// read in reaction model enlarger
        	line = ChemParser.readMeaningfulLine(reader);
        	if (line.startsWith("ReactionModelEnlarger:")) {
        		StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		String rmeType = st.nextToken();
        		if (rmeType.equals("RateBasedModelEnlarger")) {
        			reactionModelEnlarger = new RateBasedRME();
        			PDepNetwork.generateNetworks = false;
					line = ChemParser.readMeaningfulLine(reader);
					st = new StringTokenizer(line);
					String iS = st.nextToken();
					String fileName = st.nextToken();
					if (iS.startsWith("IncludeSpecies")) {
						HashSet includeSpecies = readIncludeSpecies(fileName);
						((RateBasedRME)reactionModelEnlarger).addIncludeSpecies(includeSpecies);
						line = ChemParser.readMeaningfulLine(reader);
					}
        		}
        		else if (rmeType.equals("RateBasedPDepModelEnlarger")) {
        			reactionModelEnlarger = new RateBasedPDepRME();
        			PDepNetwork.generateNetworks = true;
					line = ChemParser.readMeaningfulLine(reader);
					if (line.startsWith("PDepKineticsEstimator:")) {
						st = new StringTokenizer(line);
						name = st.nextToken();
						String pdkeType = st.nextToken().toLowerCase();
						if (!(reactionModelEnlarger instanceof RateBasedPDepRME))
							throw new InvalidSymbolException("condition.txt: PDepKineticsEstimator specified, but reaction model enlarger is not pressure-dependent!");
						if (pdkeType.equals("reservoirstate")) {
							((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new FastMasterEqn(FastMasterEqn.Mode.RESERVOIRSTATE));
						}
						else if (pdkeType.equals("modifiedstrongcollision")) {
							((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new FastMasterEqn(FastMasterEqn.Mode.STRONGCOLLISION));
						}
						else if (pdkeType.equals("chemdis")) {
							((RateBasedPDepRME) reactionModelEnlarger).setPDepKineticsEstimator(new Chemdis());
							if (!SpectroscopicData.useThreeFrequencyModel) {
								System.out.println("Warning: Switching SpectroscopicDataEstimator to three-frequency model.");
								SpectroscopicData.useThreeFrequencyModel = true;
							}
						}
						else {
							throw new InvalidSymbolException("condition.txt: Unknown PDepKineticsEstimator = " + pdkeType);
						}
					}
					else throw new InvalidSymbolException("condition.txt: can't find PDepKineticsEstimator!");
					line = ChemParser.readMeaningfulLine(reader);
        		}
        		else {
        			throw new InvalidSymbolException("condition.txt: Unknown ReactionModelEnlarger = " + rmeType);
        		}
				
        	}
        	else throw new InvalidSymbolException("condition.txt: can't find ReactionModelEnlarger!");
			
        	// read in finish controller
        	//line = ChemParser.readMeaningfulLine(reader);
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
               
        	line = ChemParser.readMeaningfulLine(reader);
                if (line.startsWith("PrimaryReactionLibrary:")) {
                        StringTokenizer st = new StringTokenizer(line);
                        String temp = st.nextToken();
                        String on = st.nextToken();
                        
                        if (on.compareToIgnoreCase("ON") == 0) {
                        	// GJB modified to allow multiple primary reaction libraries
                        	int Ilib = 0;
                        	line = ChemParser.readMeaningfulLine(reader);
                        	while (!line.equals("END")) {                     		
                                StringTokenizer nameST = new StringTokenizer(line);
                                temp = nameST.nextToken();
                                String name = nameST.nextToken();

                                line = ChemParser.readMeaningfulLine(reader);
                                StringTokenizer pathST = new StringTokenizer(line);
                                temp = pathST.nextToken();
                                String path = pathST.nextToken();
                                if (Ilib==0) {
                                	//primaryReactionLibrary = new PrimaryReactionLibrary(name, path);
                                        setPrimaryReactionLibrary(new PrimaryReactionLibrary(name, path));//10/14/07 gmagoon: changed to use setPrimaryReactionLibrary
                                	Ilib++; 	
                                }
                                else {
                                	//primaryReactionLibrary.appendPrimaryReactionLibrary(name,path);
                                        getPrimaryReactionLibrary().appendPrimaryReactionLibrary(name,path);//10/14/07 gmagoon: changed to use getPrimaryReactionLibrary; check to make sure this is valid
                                	Ilib++;//just in case anybody wants to track how many are processed
                                 }
                                	line = ChemParser.readMeaningfulLine(reader);
                        	}
                        	System.out.println("Primary Reaction Libraries in use: " +getPrimaryReactionLibrary().getName());//10/14/07 gmagoon: changed to use getPrimaryReactionLibrary
                        }	
                         else {
                                //primaryReactionLibrary = null;
                                setPrimaryReactionLibrary(null);//10/14/07 gmagoon: changed to use setPrimaryReactionLibrary; check to make sure this is valid
                        }
                }

        	else throw new InvalidSymbolException("condition.txt: can't find PrimaryReactionLibrary!");

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
        //#]
    }
    public void setReactionModel(ReactionModel p_ReactionModel) {
        reactionModel = p_ReactionModel;
    }
	 

	


	
 

	//## operation modelGeneration()
    public void modelGeneration() {
        //#[ operation modelGeneration()
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
				
				writeCoreSpecies();
				double pt = System.currentTimeMillis();
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
                                double chemkint = (System.currentTimeMillis()-startTime)/1000/60;
				
                                //10/24/07 gmagoon: changed to use reactionSystemList
                                for (Integer i = 0; i<reactionSystemList.size();i++) {
                                    ReactionSystem rs = (ReactionSystem)reactionSystemList.get(i);
                                    System.out.println("For reaction system: "+(i+1)+" out of "+reactionSystemList.size());
                                    System.out.println("At this time: " + ((ReactionTime)endList.get(i)).toString());
                                    Species spe = SpeciesDictionary.getSpeciesFromID(1);
                                    double conv = rs.getPresentConversion(spe);
                                    System.out.print("current conversion = ");
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
                    if(!terminated)
                        allTerminated = false;
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
                            System.out.println("At this time: " + ((ReactionTime)endList.get(i)).toString());
                            Species spe = SpeciesDictionary.getSpeciesFromID(1);
                            double conv = rs.getPresentConversion(spe);
                            System.out.print("current conversion = ");
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

        		System.out.println("After garbage collection:");
        		System.out.print("Memory used: ");
        		System.out.println(runTime.totalMemory());
        		System.out.print("Free memory: ");
        		System.out.println(runTime.freeMemory());
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
        
        System.out.println("Model Generation Completed");
		
		
        
        return;





        //#]
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
		writeCoreSpecies();
		//writeCoreReactions(p_rs, p_time);
		writeEdgeSpecies();
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
		StringBuilder restartFileContent = new StringBuilder();
		int speciesCount = 0;
		try{
			File edgeSpecies = new File ("Restart/edgeSpecies.txt");
			FileWriter fw = new FileWriter(edgeSpecies);
			for(Iterator iter=((CoreEdgeReactionModel)getReactionModel()).getUnreactedSpeciesSet().iterator();iter.hasNext();){
				Species species = (Species) iter.next();
				restartFileContent.append(species.getID());
				restartFileContent.append( "\n");
				/*Species species = (Species) iter.next();
				restartFileContent = restartFileContent + "("+ speciesCount + ") "+species.getChemkinName() + " \n ";// + 0 + " (mol/cm3) \n";
				restartFileContent = restartFileContent + species.toString(1) + "\n";
				speciesCount = speciesCount + 1;*/
			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent.toString());
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart edgespecies file");
        	System.exit(0);
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

		StringBuilder restartFileContent = new StringBuilder();
		int speciesCount = 0;
		try{
			File coreSpecies = new File ("Restart/coreSpecies.txt");
			FileWriter fw = new FileWriter(coreSpecies);
			for(Iterator iter=getReactionModel().getSpecies();iter.hasNext();){
				Species species = (Species) iter.next();
				restartFileContent.append(species.getID());
				restartFileContent.append( "\n");//("+ speciesCount + ") "+species.getChemkinName() + " \n ";// + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
				//restartFileContent = restartFileContent + species.toString(1) + "\n";
				//speciesCount = speciesCount + 1;
			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent.toString());
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart corespecies file");
        	System.exit(0);
		}

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
    	
		// Determine initial set of reactions and edge species using only the 
		// species enumerated in the input file as the core
		LinkedHashSet reactionSet = getReactionGenerator().react(getSpeciesSeed());
		reactionSet.addAll(getLibraryReactionGenerator().react(getSpeciesSeed()));
    	
    	// Set initial core-edge reaction model based on above results
		if (reactionModelEnlarger instanceof RateBasedRME)	
    		setReactionModel(new CoreEdgeReactionModel(new LinkedHashSet(getSpeciesSeed()),reactionSet));//10/4/07 gmagoon: changed to use setReactionModel
    	else {
    		setReactionModel(new CoreEdgeReactionModel(new LinkedHashSet(getSpeciesSeed())));//10/4/07 gmagoon: changed to use setReactionModel
        	// Only keep the reactions involving bimolecular reactants and bimolecular products
			Iterator iter = reactionSet.iterator();
        	while (iter.hasNext()){
        		Reaction r = (Reaction)iter.next();
        		if (r.getReactantNumber() == 2 && r.getProductNumber() == 2){
        			((CoreEdgeReactionModel)getReactionModel()).addReaction(r);
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
        while (getReactionModel().isEmpty()) {
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
        if (hasPrimaryReactionLibrary()) initializeCoreEdgeModelWithPRL();
        else initializeCoreEdgeModelWithoutPRL();

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
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ReactionModelGenerator.java
*********************************************************************/


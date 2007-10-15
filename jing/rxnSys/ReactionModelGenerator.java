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

    protected ReactionSystem reactionSystem;

    protected int paraInfor;//svp
    protected boolean error;//svp
    protected boolean sensitivity;//svp
    protected LinkedList species;//svp
    protected InitialStatus initialStatus;//svp
    protected double rtol;//svp
    protected double atol;
    protected PrimaryReactionLibrary primaryReactionLibrary;//9/24/07 gmagoon
    protected ReactionModelEnlarger reactionModelEnlarger;//9/24/07 gmagoon
    protected LinkedHashSet speciesSeed;//9/24/07 gmagoon;
    protected ReactionGenerator reactionGenerator;//9/24/07 gmagoon
    protected LibraryReactionGenerator lrg;// = new LibraryReactionGenerator();//9/24/07 gmagoon: moved from ReactionSystem.java;10/4/07 gmagoon: postponed initialization of lrg til later


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
    public void initializeReactionSystem() throws InvalidSymbolException, IOException {
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

        	TemperatureModel temperatureModel = null;
        	PressureModel pressureModel = null;
  //      	ReactionModelEnlarger reactionModelEnlarger = null;//10/9/07 gmagoon: commented out: unneeded now and causes scope problems
        	
        	FinishController finishController = null;
        	DynamicSimulator dynamicSimulator = null;
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
        	if (line.startsWith("TemperatureModel:")) {
        		StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		String modelType = st.nextToken();
        		String t = st.nextToken();
        		String unit = st.nextToken();
        		if (modelType.equals("Constant")) {
        			double temp = Double.parseDouble(t);
        			unit = ChemParser.removeBrace(unit);
        			temperatureModel = new ConstantTM(temp, unit);
					Global.temperature = new Temperature(temp,unit);
        		}
        		else if (modelType.equals("Curved")) {
        			// add reading curved temperature function here
        			temperatureModel = new CurvedTM(new LinkedList());
        		}
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
        		String p = st.nextToken();
        		String unit = st.nextToken();
        		if (modelType.equals("Constant")) {
        			double pressure = Double.parseDouble(p);
        			unit = ChemParser.removeBrace(unit);
        			pressureModel = new ConstantPM(pressure, unit);
        			Global.pressure = new Pressure(pressure, unit);
        		}
        		else if (modelType.equals("Curved")) {
        			// add reading curved temperature function here
        			pressureModel = new CurvedPM(new LinkedList());
        		}
        		else {
        			throw new InvalidSymbolException("condition.txt: Unknown PressureModel = " + modelType);
        		}
        	}
        	else throw new InvalidSymbolException("condition.txt: can't find PressureModel!");

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
        		initialStatus = new InitialStatus(speciesStatus,temperatureModel.getTemperature(initial),pressureModel.getPressure(initial));

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
        		
        		int numConversions = 0;
        		// read in time step
        		line = ChemParser.readMeaningfulLine(reader);
        		if (line.startsWith("TimeStep:") && finishController.terminationTester instanceof ReactionTimeTT) {
        			st = new StringTokenizer(line);
        			temp = st.nextToken();
        			while (st.hasMoreTokens()) {
        				double tStep = Double.parseDouble(st.nextToken());
            			String unit = "sec";
            			setTimeStep(new ReactionTime(tStep, unit));
        			}      
        			((ReactionTimeTT)finishController.terminationTester).setTimeSteps(timeStep);
        		}
        		else if (line.startsWith("Conversions:") && finishController.terminationTester instanceof ConversionTT){
        			st = new StringTokenizer(line);
        			temp = st.nextToken();
        			int i=0;
        			SpeciesConversion sc = (SpeciesConversion)((ConversionTT)finishController.terminationTester).speciesGoalConversionSet.get(0);
        			Species convSpecies = sc.species;
        			Iterator iter = initialStatus.getSpeciesStatus();
        			double initialConc = 0;
        			while (iter.hasNext()){
        				SpeciesStatus sps = (SpeciesStatus)iter.next();
        				if (sps.species.equals(convSpecies)) initialConc = sps.concentration;
        			}
        			while (st.hasMoreTokens()){
        				double conv = Double.parseDouble(st.nextToken());
            			conversionSet[i] = (1-conv) * initialConc;
            			i++;
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
                               
                               dynamicSimulator = new JDASPK(rtol, atol, 0, initialStatus);
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
        			dynamicSimulator = new JDASSL(rtol, atol, 0, initialStatus);
        		}
        		else if (simulator.equals("Chemkin")) {
        			line = ChemParser.readMeaningfulLine(reader);
        			if (line.startsWith("ReactorType")) {
        				st = new StringTokenizer(line, ":");
        				temp = st.nextToken();
        				String reactorType = st.nextToken().trim();
        				dynamicSimulator = new Chemkin(rtol, atol, reactorType);
        			}
        		}
        		else throw new InvalidSymbolException("condition.txt: Unknown DynamicSimulator = " + simulator);
                dynamicSimulator.addConversion(conversionSet, numConversions);
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
                //10/4/07 gmagoon: moved to modelGeneration()
               //ReactionGenerator p_reactionGenerator = new TemplateReactionGenerator();//10/4/07 gmagoon: changed to p_reactionGenerator from reactionGenerator
               // setReactionGenerator(p_reactionGenerator);//10/4/07 gmagoon: added
                setLibraryReactionGenerator(new LibraryReactionGenerator());//10/10/07 gmagoon: moved from modelGeneration (sequence lrg increases species id, and the different sequence was causing problems as main species id was 6 instead of 1)
        	reactionSystem = new ReactionSystem(temperatureModel, pressureModel, reactionModelEnlarger, finishController, dynamicSimulator, getPrimaryReactionLibrary(), getReactionGenerator(), getSpeciesSeed(), initialStatus, getReactionModel(),lrg);

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
                setReactionGenerator(new TemplateReactionGenerator());//10/4/07 gmagoon: moved inside initializeReactionSystem
              //  setLibraryReactionGenerator(new LibraryReactionGenerator());//10/10/07 gmagoon: moved after initializeReactionSystem
              //  initializeCoreEdgeReactionModel();//10/4/07 gmagoon moved from below to run initializeCoreEdgeReactionModel before initializeReactionSystem
         	initializeReactionSystem();
        }
        catch (IOException e) {
        	System.err.println(e.getMessage());
        	System.exit(0);
        }
        catch (InvalidSymbolException e) {
        	System.err.println(e.getMessage());
        	System.exit(0);
        }
		/*if (restart){
			reactionSystem.reactionModel = new CoreEdgeReactionModel();
			parseRestartFiles();
			((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactedSpeciesSet(reactionSystem.originalReactant);

			if (reactionSystem.primaryReactionLibrary != null){
				((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactedSpeciesSet(reactionSystem.primaryReactionLibrary.getSpeciesSet());								
				((CoreEdgeReactionModel)reactionSystem.reactionModel).addPrimaryReactionSet(reactionSystem.primaryReactionLibrary.getReactionSet());

			}
		}*/
		initializeCoreEdgeReactionModel();//10/4/07 gmagoon: moved before initializeReactionSystem
        reactionSystem.initializePDepNetwork();
		
		

        
        ReactionTime init = reactionSystem.getInitialReactionTime();
        ReactionTime begin = init;
        ReactionTime end;
        if (reactionSystem.finishController.terminationTester instanceof ReactionTimeTT)
        	end = (ReactionTime)timeStep.get(0);
        else
        	end = new ReactionTime(1e6,"sec");
        int iterationNumber = 1;

        Temperature lastT = reactionSystem.getTemperature(init);
        Temperature currentT = reactionSystem.getTemperature(init);
        Pressure lastP = reactionSystem.getPressure(init);
        Pressure currentP = reactionSystem.getPressure(init);
        boolean conditionChanged = false;

        //Chemkin.writeChemkinInputFile(reactionSystem.getReactionModel(),reactionSystem.getPresentStatus());
        
		
        end = reactionSystem.solveReactionSystem(begin, end, true, true, true, iterationNumber-1);
        Chemkin.writeChemkinInputFile(reactionSystem);
        //System.exit(0);
        boolean terminated = reactionSystem.isReactionTerminated();
        boolean valid = reactionSystem.isModelValid();
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
		System.out.println(dictionary.size());
        boolean reactionChanged = false;
		
		double tAtInitialization = Global.tAtInitialization;
		
		
        // step 2: iteratively grow reaction system
        while (!terminated || !valid) {
        	while (!valid) {
				
				writeCoreSpecies();
				double pt = System.currentTimeMillis();
				enlargeReactionModel();
				double totalEnlarger = (System.currentTimeMillis() - pt)/1000/60;
				
				//PDepNetwork.completeNetwork(reactionSystem.reactionModel.getSpeciesSet());
				reactionSystem.initializePDepNetwork();
				
				pt = System.currentTimeMillis();
				reactionSystem.resetSystemSnapshot();
				double resetSystem = (System.currentTimeMillis() - pt)/1000/60;
        		reactionChanged = true;
        		
        		begin = init;
        		if (reactionSystem.finishController.terminationTester instanceof ReactionTimeTT)
        			end = (ReactionTime)timeStep.get(0);
        		else
        			end = new ReactionTime(1e6,"sec");
        		iterationNumber = 1;
        		
        		currentT = reactionSystem.getTemperature(begin);
        		currentP = reactionSystem.getPressure(begin);
        		conditionChanged = (!currentT.equals(lastT) || !currentP.equals(lastP));
        		
        		
				
				double startTime = System.currentTimeMillis();
				end = reactionSystem.solveReactionSystem(begin, end, false, reactionChanged, conditionChanged, iterationNumber-1);
				solverMin = solverMin + (System.currentTimeMillis()-startTime)/1000/60;
				
				startTime = System.currentTimeMillis();
        		Chemkin.writeChemkinInputFile(reactionSystem);
        		double chemkint = (System.currentTimeMillis()-startTime)/1000/60;
				
				System.out.println("At this time: " + end.toString());
        		Species spe = SpeciesDictionary.getSpeciesFromID(1);
        		double conv = reactionSystem.getPresentConversion(spe);
        		System.out.print("current conversion = ");
        		System.out.println(conv);

			    System.out.println("Running Time is: " + String.valueOf((System.currentTimeMillis()-tAtInitialization)/1000/60) + " minutes.");
				System.out.println("The model edge has " + ((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSet().size() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getUnreactedSpeciesSet().size() + " species.");
				if (reactionSystem.getDynamicSimulator() instanceof JDASPK){
					JDASPK solver = (JDASPK)reactionSystem.getDynamicSimulator();
					System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
				}
				else{
					JDASSL solver = (JDASSL)reactionSystem.getDynamicSimulator();
					System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
				
				}
					

				
				startTime = System.currentTimeMillis();
				double mU = memoryUsed();
				double gc = (System.currentTimeMillis()-startTime)/1000/60;
				
				
				
				
				
				startTime = System.currentTimeMillis();
                valid = reactionSystem.isModelValid();
				vTester = vTester + (System.currentTimeMillis()-startTime)/1000/60;
				
				startTime = System.currentTimeMillis();

				writeDiagnosticInfo();
				writeEnlargerInfo();
				double restart2 = (System.currentTimeMillis()-startTime)/1000/60;
				
				int allSpecies, allReactions;
				allSpecies = SpeciesDictionary.getInstance().size();
				print_info.append(totalEnlarger + "\t" + resetSystem + "\t" + Global.readSolverFile + "\t" + Global.writeSolverFile + "\t" + Global.solvertime + "\t" + Global.solverIterations + "\t" + Global.speciesStatusGenerator +  "\t" + solverMin + "\t"  + gc + "\t"  + restart2 + "\t" + Global.chemkinThermo + '\t' + Global.chemkinReaction + "\t" + vTester + "\t" + ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size()+ "\t" + ((CoreEdgeReactionModel)getReactionModel()).getReactedReactionSet().size() + "\t" + ((CoreEdgeReactionModel)getReactionModel()).getUnreactedSpeciesSet().size() + "\t" + ((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSetIncludingReverseSize() + "\t" + mU + "\t" + allSpecies + "\t" + (System.currentTimeMillis()-Global.tAtInitialization)/1000/60 + "\t"+ String.valueOf(Global.RT_findRateConstant)+"\t"+Global.RT_identifyReactedSites+"\t"+Global.RT_reactChemGraph+"\t"+Global.makeSpecies+"\t"+Global.checkReactionReverse+"\t"+Global.makeTR+ "\t" + Global.getReacFromStruc + "\t" + Global.generateReverse+"\n");
				
        	}
        	reactionChanged = false;
        	
        	lastT = (Temperature)currentT.clone();
        	lastP = (Pressure)currentP.clone();
        	
        	currentT = reactionSystem.getTemperature(begin);
        	currentP = reactionSystem.getPressure(begin);
        	conditionChanged = (!currentT.equals(lastT) || !currentP.equals(lastP));
			
        	
        	begin=((SystemSnapshot)(reactionSystem.getSystemSnapshotEnd().next())).time;
        	if (reactionSystem.finishController.terminationTester instanceof ReactionTimeTT){
        		if (iterationNumber < timeStep.size()){
            		end = (ReactionTime)timeStep.get(iterationNumber);
        			
            	}
            	else
            		end = ((ReactionTimeTT)reactionSystem.finishController.terminationTester).finalTime;
        	}
        	else
        		end = new ReactionTime(1e6,"sec");
        	iterationNumber++;
			double startTime = System.currentTimeMillis();
			end = reactionSystem.solveReactionSystem(begin, end, false, reactionChanged, false, iterationNumber-1);
			solverMin = solverMin + (System.currentTimeMillis()-startTime)/1000/60;
			
			startTime = System.currentTimeMillis();
        	
			terminated = reactionSystem.isReactionTerminated();
        	valid = reactionSystem.isModelValid();
			
			
        	if (valid) {
        		System.out.println("At this time: " + end.toString());
        		Species spe = SpeciesDictionary.getSpeciesFromID(1);
        		double conv = reactionSystem.getPresentConversion(spe);
        		System.out.print("current conversion = ");
        		System.out.println(conv);

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
        		if (reactionSystem.getDynamicSimulator() instanceof JDASPK){
					JDASPK solver = (JDASPK)reactionSystem.getDynamicSimulator();
					System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
        		}
				else{
					JDASSL solver = (JDASSL)reactionSystem.getDynamicSimulator();
					System.out.println("The model core has " + solver.getReactionSize() + " reactions and "+ ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet().size() + " species.");
        		}
				
        	}
			vTester = vTester + (System.currentTimeMillis()-startTime)/1000/60;
        }
        
        //System.out.println("Performing model reduction");
        
        if (paraInfor != 0){
        	System.out.println("Model Generation performed. Now generating sensitivity data.");
          DynamicSimulator dynamicSimulator2 = new JDASPK(rtol, atol, paraInfor, initialStatus);
          dynamicSimulator2.addConversion(((JDASPK)reactionSystem.dynamicSimulator).conversionSet, ((JDASPK)reactionSystem.dynamicSimulator).conversionSet.length);
          reactionSystem.setDynamicSimulator(dynamicSimulator2);
          
          int numSteps = reactionSystem.systemSnapshot.size() -1;
          reactionSystem.resetSystemSnapshot();
          begin = init;
          if (reactionSystem.finishController.terminationTester instanceof ReactionTimeTT){
          		end = ((ReactionTimeTT)reactionSystem.finishController.terminationTester).finalTime;
          }
          else
 	     		end = end.add(end);
          terminated = false;
          reactionSystem.solveReactionSystemwithSEN(begin, end, true, false, false);
         
        }

        Chemkin.writeChemkinInputFile(getReactionModel(),reactionSystem.getPresentStatus());
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
			if (reactionSystem == null){
				//ReactionSystem reactionSystem = new ReactionSystem();
			}
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
	
	public void parseAllSpecies() {
//		String restartFileContent ="";
		int speciesCount = 0;

		boolean added;
		try{
			long initialTime = System.currentTimeMillis();
			//StringBuilder sb = new StringBuilder();
			//sb.append("\t Cumulative Time until after the event\n");
			//sb.append("# \t Read Graph \t ChemGraph \t Species \t total\t Memory Used");
			
			File coreSpecies = new File ("Restart/allSpecies.txt");
			BufferedReader reader = new BufferedReader(new FileReader(coreSpecies));
			String line = ChemParser.readMeaningfulLine(reader);
			HashSet speciesSet = new HashSet();
			int i=0;
			while (line!=null) {
				i++;
				//sb.append(i);sb.append("\t");
    			StringTokenizer st = new StringTokenizer(line);
    			String index = st.nextToken();
    			String name = null;
    			if (!index.startsWith("(")) name = index;
    			else name = st.nextToken().trim();
				int ID = getID(name);
				
				name = getName(name);
    			Graph g = ChemParser.readChemGraph(reader);
				double timeNow = (System.currentTimeMillis()-initialTime)/1e3;
				//sb.append(timeNow);sb.append("\t");
				
    			ChemGraph cg = null;
    			try {
    				cg = ChemGraph.make(g);
    			}
    			catch (ForbiddenStructureException e) {
    				System.out.println("Forbidden Structure:\n" + e.getMessage());
    				System.exit(0);
    			}
				timeNow = (System.currentTimeMillis()-initialTime)/1e3;
				//sb.append(timeNow);sb.append("\t");

    			Species species = Species.make(name,cg,ID);
       			//speciesSet.put(name, species);
    			speciesSet.add(species);
				timeNow = (System.currentTimeMillis()-initialTime)/1e3;
				//sb.append(timeNow);sb.append("\t");
				/*if (name.equals("C5H11.")){
					System.out.println(species.getChemkinName());
				}*/

    			double flux = 0;
    			int species_type = 1; // reacted species
    			//SpeciesStatus ss = new SpeciesStatus(species,species_type,concentration,flux);
    			//speciesStatus.put(species, ss);
    			line = ChemParser.readMeaningfulLine(reader);
				//sb.append((System.currentTimeMillis()-initialTime)/1e3);sb.append("\t");
				//sb.append(memoryUsed());sb.append("\t");
				//System.out.println(sb.toString());
				//sb.delete(0,sb.length());
    		}
			
			//reactionSystem.reactionModel = new CoreEdgeReactionModel();
			//((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactedSpeciesSet(speciesSet);
			//specs.addAll(speciesSet);
		}
		catch (IOException e){
			System.out.println("Could not read the allSpecies restart file");
        	System.exit(0);
		}

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

	//Is still incomplete.
    public void writeRestartFile() {
		writeCoreSpecies();
		//writeCoreReactions();
		writeEdgeSpecies();
		//writeAllReactions();
		//writeEdgeReactions();

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

	private void writeEdgeReactions() {
		StringBuilder restartFileContent =new StringBuilder();
		int reactionCount = 1;
		try{
			File coreSpecies = new File ("Restart/edgeReactions.txt");
			FileWriter fw = new FileWriter(coreSpecies);
			for(Iterator iter=((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSet().iterator();iter.hasNext();){

				Reaction reaction = (Reaction) iter.next();
				//if (reaction.getDirection()==1){
					//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
					restartFileContent = restartFileContent.append(reaction.toRestartString(Global.temperature) + "\n");
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

	private void writeAllReactions() {
		StringBuilder restartFileContent = new StringBuilder();
		int reactionCount = 1;
		try{
			File allReactions = new File ("Restart/allReactions.txt");
			FileWriter fw = new FileWriter(allReactions);
			for(Iterator iter=getReactionModel().getReaction();iter.hasNext();){

				Reaction reaction = (Reaction) iter.next();

				//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
				restartFileContent = restartFileContent.append(reaction.toRestartString(Global.temperature) + "\n");

			}

			for(Iterator iter=((CoreEdgeReactionModel)getReactionModel()).getUnreactedReactionSet().iterator();iter.hasNext();){

				Reaction reaction = (Reaction) iter.next();
				//if (reaction.getDirection()==1){
					//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
				restartFileContent = restartFileContent.append(reaction.toRestartString(Global.temperature) + "\n");

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

	private void writeCoreReactions() {
		StringBuilder restartFileContent = new StringBuilder();
		int reactionCount = 0;
		try{
			File coreSpecies = new File ("Restart/coreReactions.txt");
			FileWriter fw = new FileWriter(coreSpecies);
			for(Iterator iter=getReactionModel().getReaction();iter.hasNext();){

				Reaction reaction = (Reaction) iter.next();
				if (reaction.getDirection()==1){
					//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
					restartFileContent = restartFileContent.append(reaction.toRestartString(Global.temperature) + "\n");
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


    public ReactionSystem getReactionSystem() {
        return reactionSystem;
    }

    //added by gmagoon 9/24/07
    public void setReactionSystem(ReactionSystem p_ReactionSystem) {
        reactionSystem = p_ReactionSystem;
    }
    
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
    	LinkedHashSet reactionSet = getReactionGenerator().react(getSpeciesSeed());
       // LinkedHashSet reactionSetTwo = getLibraryReactionGenerator().react(getSpeciesSeed()); //gmagoon: added for testing around 10/5/07
    	reactionSet.addAll(getLibraryReactionGenerator().react(getSpeciesSeed()));
    	
    	if (reactionModelEnlarger instanceof RateBasedRME)	
    		setReactionModel(new CoreEdgeReactionModel(new LinkedHashSet(getSpeciesSeed()),reactionSet));//10/4/07 gmagoon: changed to use setReactionModel
    	else {
    		setReactionModel(new CoreEdgeReactionModel(new LinkedHashSet(getSpeciesSeed())));//10/4/07 gmagoon: changed to use setReactionModel
        	Iterator iter = reactionSet.iterator();
        	while (iter.hasNext()){
        		Reaction r = (Reaction)iter.next();
        		if (r.getReactantNumber() == 2 && r.getProductNumber() == 2){
        			((CoreEdgeReactionModel)getReactionModel()).addReaction(r);
        		}
        	}
    	}
        //10/9/07 gmagoon: copy reactionModel to reactionSystem; there may still be scope problems, particularly in above elseif statement
        reactionSystem.setReactionModel(getReactionModel());
        if (getReactionModel().isEmpty() && reactionModelEnlarger instanceof RateBasedRME) {
        	LinkedHashSet us = ((CoreEdgeReactionModel)getReactionModel()).getUnreactedSpeciesSet();
        	LinkedHashSet rs = ((CoreEdgeReactionModel)getReactionModel()).getReactedSpeciesSet();
        	LinkedHashSet newReactions = new LinkedHashSet();

        	for (Iterator iter = us.iterator(); iter.hasNext(); ) {
        		Species spe = (Species)iter.next();
        		rs.add(spe);
        		newReactions.addAll(getReactionGenerator().react(rs,spe));
        		newReactions.addAll(lrg.react(rs,spe));
        		iter.remove();
        	}

        	((CoreEdgeReactionModel)getReactionModel()).addReactionSet(newReactions);
        	((CoreEdgeReactionModel)getReactionModel()).moveFromUnreactedToReactedReaction();
        }
        else if (getReactionModel().isEmpty() && reactionModelEnlarger instanceof RateBasedPDepRME) {
        	while (getReactionModel().isEmpty()){
        		reactionSystem.initializePDepNetwork();//9/25/07 gmagoon: unresolved question: what to do with initializePDepNetwork? leave in ReactionSystem class or move elsewhere (e.g. ReactionModelGenerator)?;//10/14/07 gmagoon: I am leaning towards leaving it where it is since it is used in solveReactionSystem, but it probably warrants further investigation
        		reactionSystem.appendUnreactedSpeciesStatus(initialStatus, Global.temperature);
        		enlargeReactionModel();
        	}
        }
        //10/9/07 gmagoon: copy reactionModel to reactionSystem; there may still be scope problems, particularly in above elseif statement
        reactionSystem.setReactionModel(getReactionModel());
        return;

        //#]
    }

    //## operation initializeCoreEdgeReactionModel()
    //9/24/07 gmagoon: moved from ReactionSystem.java
    public void initializeCoreEdgeReactionModel() {
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
    //## operation enlargeReactionModel()
    public void enlargeReactionModel() {
        //#[ operation enlargeReactionModel()
        if (reactionModelEnlarger == null) throw new NullPointerException("ReactionModelEnlarger");

        reactionModelEnlarger.enlargeReactionModel(reactionSystem, getReactionModel());

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


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

    protected ReactionTime timeStep;		//## attribute timeStep

    protected String workingDirectory;		//## attribute workingDirectory

    protected ReactionSystem reactionSystem;

    protected int paraInfor;//svp
    protected boolean error;//svp
    protected boolean sensitivity;//svp
    protected LinkedList species;//svp
    protected InitialStatus initialStatus;//svp
    protected double rtol;//svp
    protected double atol;


	protected boolean restart;
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
        	ReactionModelEnlarger reactionModelEnlarger = null;
        	
        	FinishController finishController = null;
        	DynamicSimulator dynamicSimulator = null;
        	PrimaryReactionLibrary primaryReactionLibrary = null;
			String line = ChemParser.readMeaningfulLine(reader);
        	if (line.startsWith("Restart")){
				StringTokenizer st = new StringTokenizer(line);
				String token = st.nextToken();
				token = st.nextToken();
				if (token.equalsIgnoreCase("true")) {
					Runtime.getRuntime().exec("cp Restart/allSpecies.txt Restart/allSpecies1.txt");
					Runtime.getRuntime().exec("echo  >> allSpecies.txt");
					restart = true;
				}
				else if (token.equalsIgnoreCase("false")) {
					Runtime.getRuntime().exec("rm Restart/allSpecies.txt");
					restart = false;
				}
				else throw new InvalidSymbolException("UnIdentified Symbol "+token+" after Restart:");
        	}
        	else throw new InvalidSymbolException("Can't find Restart!");

			line = ChemParser.readMeaningfulLine(reader);

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
        	HashSet speciesSeed = new HashSet();
        	HashMap speciesSet = new HashMap();
        	HashMap speciesStatus = new HashMap();
			int speciesnum = 1;
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
           			speciesSet.put(name, species);
        			speciesSeed.add(species);
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
        		}
        		else if (rmeType.equals("RateBasedPDepModelEnlarger")) {
        			reactionModelEnlarger = new RateBasedPDepRME();
        		}
        		else {
        			throw new InvalidSymbolException("condition.txt: Unknown ReactionModelEnlarger = " + rmeType);
        		}
        	}
        	else throw new InvalidSymbolException("condition.txt: can't find ReactionModelEnlarger!");

        	// read in finish controller
        	line = ChemParser.readMeaningfulLine(reader);
        	if (line.startsWith("FinishController")) {
        		line = ChemParser.readMeaningfulLine(reader);
        		StringTokenizer st = new StringTokenizer(line);
        		String index = st.nextToken();
        		String goal = st.nextToken();
        		String type = st.nextToken();
        		TerminationTester tt;
        		if (type.startsWith("Conversion")) {
        			HashSet spc = new HashSet();
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
        		// read in time step
        		line = ChemParser.readMeaningfulLine(reader);
        		if (line.startsWith("TimeStep:")) {
        			st = new StringTokenizer(line);
        			temp = st.nextToken();
        			double tStep = Double.parseDouble(st.nextToken());
        			String unit = st.nextToken();
        			unit = ChemParser.removeBrace(unit);
        			setTimeStep(new ReactionTime(tStep, unit));
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
        			dynamicSimulator = new JDASSL();
        		}
        		else throw new InvalidSymbolException("condition.txt: Unknown DynamicSimulator = " + simulator);
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
                                	primaryReactionLibrary = new PrimaryReactionLibrary(name, path);
                                	Ilib++; 	
                                }
                                else {
                                	primaryReactionLibrary.appendPrimaryReactionLibrary(name,path);
                                	Ilib++;//just in case anybody wants to track how many are processed
                                 }
                                	line = ChemParser.readMeaningfulLine(reader);
                        	}
                        	System.out.println("Primary Reaction Libraries in use: " +primaryReactionLibrary.getName());
                        }	
                         else {
                                primaryReactionLibrary = null;
                        }
                }

        	else throw new InvalidSymbolException("condition.txt: can't find PrimaryReactionLibrary!");

        	in.close();

        	ReactionGenerator reactionGenerator = new TemplateReactionGenerator();

        	reactionSystem = new ReactionSystem(temperatureModel, pressureModel, reactionModelEnlarger, finishController, dynamicSimulator, primaryReactionLibrary, reactionGenerator, speciesSeed, initialStatus);

		}
        catch (IOException e) {
        	System.err.println("Error in read in reaction system initialization file!");
        	throw new IOException("Reaction System Initialization: " + e.getMessage());
        }
        //#]
    }

    private void parseRestartFiles() {
		parseAllSpecies();
		parseCoreSpecies();
		parseEdgeSpecies();
		parseAllReactions();
		parseCoreReactions();
		//parseCoreReactions();
		//parseEdgeReactions();
		/*Iterator iter = specs.iterator();
		for (int i=0;i<specs.size();i++){
			Species temp = (Species)iter.next();
			System.out.println(temp.getID()+"\t"+temp.getChemkinName()+"\t"+ temp.getResonanceIsomersHashSet().size());;
		}*/
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
			HashSet reactionSet = new HashSet();
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
			((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactionSet(reactionSet);
		}
		catch (IOException e){
			System.out.println("Could not read the corespecies restart file");
        	System.exit(0);
		}

	}

	private void parseCoreReactions() {
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		int i=1;
		//HasMap speciesMap = dictionary.dictionary;
		try{
			File coreReactions = new File("Restart/coreReactions.txt");
			FileReader fr = new FileReader(coreReactions);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			boolean found = false;
			HashSet reactionSet = new HashSet();
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
			((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactedReactionSet(reactionSet);
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
			HashSet reactionSet = new HashSet();
			OuterLoop:
			while (line != null){
				Reaction reaction = ChemParser.parseArrheniusReaction(dictionary,line,1,1,((CoreEdgeReactionModel)reactionSystem.reactionModel));
				/*if (reaction.hasReverseReaction() && ((CoreEdgeReactionModel)reactionSystem.reactionModel).isReactedReaction(reaction)){
					if (reactionSet.contains(reaction.getReverseReaction())){
						line = ChemParser.readMeaningfulLine(reader);
						continue;
					}
				}*/
				if (line.startsWith("CH2O(14) + CH3O2.(39)"))
					System.out.println("Teri maa");
				/*if (reaction.getStructure().getDirection()== -1 && ((CoreEdgeReactionModel)reactionSystem.reactionModel).isReactedReaction(reaction)){
					//line = ChemParser.readMeaningfulLine(reader);
					if (reactionSet.contains(reaction.getReverseReaction()) && !reactionSet.contains(reaction)){
						reactionSet.add(reaction);
						Iterator iter = reactionSet.iterator();
						while (iter.hasNext()){
							Reaction reacTemp = (Reaction)iter.next();
							if (reacTemp.equals(reaction.getReverseReaction())){
								reactionSet.remove(reacTemp);
								reactionSet.add(reaction.getReverseReaction());
								//line = ChemParser.readMeaningfulLine(reader);
								break;
							}
						}

					}
					else {
						line = ChemParser.readMeaningfulLine(reader);
						continue;//reaction = null;
					}
				}*/


				if (((CoreEdgeReactionModel)reactionSystem.reactionModel).categorizeReaction(reaction)==-1){
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
			((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactionSet(reactionSet);
		}
		catch (IOException e){
			System.out.println("Could not read the corespecies restart file");
        	System.exit(0);
		}

	}

	/*private void parseAllReactions() {
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		int i=1;
		//HasMap speciesMap = dictionary.dictionary;
		try{
			File allReactions = new File("Restart/allReactions.txt");
			FileReader fr = new FileReader(allReactions);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			boolean found = false;
			HashSet reactionSet = new HashSet();
			OuterLoop:
			while (line != null){
				Reaction reaction = ChemParser.parseArrheniusReaction(dictionary,line,1,1,((CoreEdgeReactionModel)reactionSystem.reactionModel));
				/*if (reaction.hasReverseReaction() && ((CoreEdgeReactionModel)reactionSystem.reactionModel).isReactedReaction(reaction)){
					if (reactionSet.contains(reaction.getReverseReaction())){
						line = ChemParser.readMeaningfulLine(reader);
						continue;
					}
				}
				if (line.startsWith("H(16) + CH3O2.(39)"))
					System.out.println("Teri maa");
				if (reaction.getStructure().getDirection()== -1 && ((CoreEdgeReactionModel)reactionSystem.reactionModel).isReactedReaction(reaction)){
					//line = ChemParser.readMeaningfulLine(reader);
					if (reactionSet.contains(reaction.getReverseReaction()) && !reactionSet.contains(reaction)){
						reactionSet.add(reaction);
						Iterator iter = reactionSet.iterator();
						while (iter.hasNext()){
							Reaction reacTemp = (Reaction)iter.next();
							if (reacTemp.equals(reaction.getReverseReaction())){
								reactionSet.remove(reacTemp);
								reactionSet.add(reaction.getReverseReaction());
								//line = ChemParser.readMeaningfulLine(reader);
								break;
							}
						}

					}
					else {
						line = ChemParser.readMeaningfulLine(reader);
						continue;//reaction = null;
					}
				}



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

				Reaction reverse = reaction.getReverseReaction();
				if (reverse != null && ((CoreEdgeReactionModel)reactionSystem.reactionModel).isReactedReaction(reaction)) {
					reactionSet.add(reverse);
					//System.out.println(2 + "\t " + line);
				}
									//else System.out.println(1 + "\t" + line);

				i=i+1;

				line = ChemParser.readMeaningfulLine(reader);
			}
			((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactionSet(reactionSet);
		}
		catch (IOException e){
			System.out.println("Could not read the corespecies restart file");
        	System.exit(0);
		}

	}*/

	private boolean getResonanceStructure(Reaction p_Reaction, String rOrP, HashSet reactionSet) {
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

	/*private void parseCoreSpecies() {
//		String restartFileContent ="";
		int speciesCount = 0;
		boolean added;
		try{
			File coreSpecies = new File ("Restart/coreSpecies.txt");
			FileReader fr = new FileReader(coreSpecies);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			HashSet speciesSet = new HashSet();

			while (line!=null) {

    			StringTokenizer st = new StringTokenizer(line);
    			String index = st.nextToken();
    			String name = null;
    			if (!index.startsWith("(")) name = index;
    			else name = st.nextToken().trim();
				int ID = getID(name);
				//String [] temp = name.split("(");
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


    			Species species = Species.make(name,cg,ID);
       			//speciesSet.put(name, species);
    			speciesSet.add(species);

				if (name.equals("C5H11.")){
					System.out.println(species.getChemkinName());
				}

    			double flux = 0;
    			int species_type = 1; // reacted species
    			//SpeciesStatus ss = new SpeciesStatus(species,species_type,concentration,flux);
    			//speciesStatus.put(species, ss);
    			line = ChemParser.readMeaningfulLine(reader);
    		}
			//reactionSystem.reactionModel = new CoreEdgeReactionModel();
			((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactedSpeciesSet(speciesSet);
			specs.addAll(speciesSet);
		}
		catch (IOException e){
			System.out.println("Could not read the corespecies restart file");
        	System.exit(0);
		}

	}*/


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
			reactionSystem.reactionModel = new CoreEdgeReactionModel();
			while (line!=null) {

    			StringTokenizer st = new StringTokenizer(line);
    			String index = st.nextToken();
				int ID = Integer.parseInt(index);
				Species spe = dictionary.getSpeciesFromID(ID);
				if (spe == null)
					System.out.println("There was no species with ID "+ID +" in the species dictionary");

				((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactedSpecies(spe);
				line = ChemParser.readMeaningfulLine(reader);
			}

		}
		catch (IOException e){
			System.out.println("Could not read the corespecies restart file");
        	System.exit(0);
		}

	}

	public void parseAllSpecies() {
//		String restartFileContent ="";
		int speciesCount = 0;

		boolean added;
		try{
			File coreSpecies = new File ("Restart/allSpecies1.txt");
			FileReader fr = new FileReader(coreSpecies);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			HashSet speciesSet = new HashSet();

			while (line!=null) {

    			StringTokenizer st = new StringTokenizer(line);
    			String index = st.nextToken();
    			String name = null;
    			if (!index.startsWith("(")) name = index;
    			else name = st.nextToken().trim();
				int ID = getID(name);
				//if (line.contains("C5H8O2"))
					//System.out.println("Teri maa");
				//System.out.println(ID + "\t");
				//String [] temp = name.split("(");
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


    			Species species = Species.make(name,cg,ID);
       			//speciesSet.put(name, species);
    			speciesSet.add(species);

				/*if (name.equals("C5H11.")){
					System.out.println(species.getChemkinName());
				}*/

    			double flux = 0;
    			int species_type = 1; // reacted species
    			//SpeciesStatus ss = new SpeciesStatus(species,species_type,concentration,flux);
    			//speciesStatus.put(species, ss);
    			line = ChemParser.readMeaningfulLine(reader);
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
				((CoreEdgeReactionModel)reactionSystem.reactionModel).addUnreactedSpecies(spe);
				line = ChemParser.readMeaningfulLine(reader);
			}

		}
		catch (IOException e){
			System.out.println("Could not read the edgepecies restart file");
        	System.exit(0);
		}

	}

	/*private void parseEdgeSpecies() {
//		String restartFileContent ="";
		int speciesCount = 0;
		int highestID=0;
		try{
			File coreSpecies = new File ("Restart/edgeSpecies.txt");
			FileReader fr = new FileReader(coreSpecies);
			BufferedReader reader = new BufferedReader(fr);
			String line = ChemParser.readMeaningfulLine(reader);
			HashSet speciesSet = new HashSet();

			while (line != null) {
				StringTokenizer st = new StringTokenizer(line);
    			String index = st.nextToken();
    			String name = null;
    			if (!index.startsWith("(")) name = index;
    			else name = st.nextToken();
				int ID = getID(name);
				name = getName(name);
				//specs[ID-1]=name;

				if (ID >highestID) highestID=ID;
    			Graph g = ChemParser.readChemGraph(reader);
    			ChemGraph cg = null;
    			try {
    				cg = ChemGraph.make(g);
    			}
    			catch (ForbiddenStructureException e) {
    				System.out.println("Forbidden Structure:\n" + e.getMessage());
    				System.exit(0);
    			}

    			Species species = Species.make(name,cg, ID);
       			//speciesSet.put(name, species);
				speciesSet.add(species);
    			double flux = 0;
    			int species_type = 1; // reacted species
    			//SpeciesStatus ss = new SpeciesStatus(species,species_type,concentration,flux);
    			//speciesStatus.put(species, ss);
    			line = ChemParser.readMeaningfulLine(reader);
    		}
			//reactionSystem.reactionModel = new CoreEdgeReactionModel();
			((CoreEdgeReactionModel)reactionSystem.reactionModel).addUnreactedSpeciesSet(speciesSet);
			specs.addAll(speciesSet);
		}
		catch (IOException e){
			System.out.println("Could not read the corespecies restart file");
        	System.exit(0);
		}

	}*/


	//## operation modelGeneration()
    public void modelGeneration() {
        //#[ operation modelGeneration()
        long begin_t = System.currentTimeMillis();
		try{
        	ChemGraph.readForbiddenStructure();
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
		if (restart){
			reactionSystem.reactionModel = new CoreEdgeReactionModel();
			parseRestartFiles();
			((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactedSpeciesSet(reactionSystem.originalReactant);

			if (reactionSystem.primaryReactionLibrary != null){
				((CoreEdgeReactionModel)reactionSystem.reactionModel).addReactedSpeciesSet(reactionSystem.primaryReactionLibrary.getSpeciesSet());								
				((CoreEdgeReactionModel)reactionSystem.reactionModel).addPrimaryReactionSet(reactionSystem.primaryReactionLibrary.getReactionSet());

			}
		}
		else reactionSystem.initializeCoreEdgeReactionModel();
        reactionSystem.initializePDepNetwork();
		//printRestartFile();

        ReactionTime delt = timeStep;
        ReactionTime init = reactionSystem.getInitialReactionTime();
        ReactionTime begin = init;
        ReactionTime end = begin.add(timeStep);

        Temperature lastT = reactionSystem.getTemperature(init);
        Temperature currentT = reactionSystem.getTemperature(init);
        Pressure lastP = reactionSystem.getPressure(init);
        Pressure currentP = reactionSystem.getPressure(init);
		int core_reaction_number, edge_reaction_number;
        boolean conditionChanged = false;

        reactionSystem.solveReactionSystem(begin, end, true, true, true);
        boolean terminated = reactionSystem.isReactionTerminated();
        boolean valid = reactionSystem.isModelValid();
		System.out.println("The model core has " + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getReactedReactionSet().size() + " reactions and "+ ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getReactedSpeciesSet().size() + " species.");
		System.out.println("The model edge has " + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedReactionSet().size() + " reactions and "+ ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedSpeciesSet().size() + " species.");
		//core_reaction_number=((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getReactedReactionSet().size();
		//edge_reaction_number=((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedReactionSet().size();
		String print_info = "Molecule \t Flux\t\tTime\t \t\t \t Core \t \t Edge \t \t memory\n";
		print_info += " \t moleular \t characteristic \t findspecies\t enlarger \t solver \t restart \t Species \t Reactions\t Species\t Reactions\n";
		print_info += "\t\t\t\t\t\t\t" + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getReactedSpeciesSet().size()+ "\t" + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getReactedReactionSet().size() + "\t" + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedSpeciesSet().size() + "\t" + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedReactionSet().size() + "\n";
		//printRestartFile();
		//((CoreEdgeReactionModel)reactionSystem.getReactionModel()).printPDepModel(reactionSystem.getPresentTemperature());;
		if (!restart){
			writeRestartFile();
			writeCoreReactions();
			writeAllReactions();
		}
		
		//writeCoreReactions();
		//writeEdgeReactions();

		Chemkin.writeChemkinInputFile(reactionSystem.getReactionModel(),reactionSystem.getPresentStatus());
		SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
		/*HashSet set = dictionary.getSpeciesSet();
		Iterator iter = set.iterator();
		for (int i=0;i<set.size();i++){
			Species temp = (Species)iter.next();
			System.out.println(temp.getChemkinName()+"\t"+temp.getID()+"\t"+ temp.getResonanceIsomersHashSet().size());
		}*/
		System.out.println(dictionary.size());
        boolean reactionChanged = false;

        // step 2: iteratively grow reaction system
        while (!terminated || !valid) {
        	while (!valid) {
				double end_t ;
				print_info += reactionSystem.enlargeReactionModel() + "\t";
				end_t = System.currentTimeMillis();
				double min=(end_t-begin_t)/1E3/60;
				print_info += min + "\t";
        		reactionSystem.resetSystemSnapshot();
        		reactionChanged = true;
        		delt = timeStep;
        		begin = init;
        		end = begin.add(delt);
        		currentT = reactionSystem.getTemperature(begin);
        		currentP = reactionSystem.getPressure(begin);
        		conditionChanged = (!currentT.equals(lastT) || !currentP.equals(lastP));
        		reactionSystem.solveReactionSystem(begin, end, false, reactionChanged, conditionChanged);
				System.out.println("At this time: " + end.toString());
        		Species spe = SpeciesDictionary.getSpeciesFromID(1);
        		double conv = reactionSystem.getPresentConversion(spe);
        		System.out.print("current conversion = ");
        		System.out.println(conv);

        		Runtime runTime = Runtime.getRuntime();
        		
        		runTime.gc();

        		System.out.println("After garbage collection:");
        		System.out.print("Memory used: ");
        		System.out.println(runTime.totalMemory());
        		System.out.print("Free memory: ");
        		System.out.println(runTime.freeMemory());
				end_t = System.currentTimeMillis();
			    
				min=(end_t-begin_t)/1E3/60;
			    System.out.println("Running Time is: " + String.valueOf(min) + " minutes.");
				System.out.println("The model edge has " + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedReactionSet().size() + " reactions and "+ ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedSpeciesSet().size() + " species.");
				System.out.println("The model core has " + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getReactedReactionSet().size() + " reactions and "+ ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getReactedSpeciesSet().size() + " species.");


				writeRestartFile();
				writeCoreReactions();
				
				double min_restart = (System.currentTimeMillis()-begin_t)/1e3/60;//getCpuTime()/1e9/60;
				print_info += min + "\t" + min_restart + "\t" + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getReactedSpeciesSet().size()+ "\t" + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getReactedReactionSet().size() + "\t" + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedSpeciesSet().size() + "\t" + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedReactionSet().size() + "\t" + runTime.totalMemory() + "\n";
				try {
			        File diagnosis = new File("diagnosis.xls");
			        FileWriter fw = new FileWriter(diagnosis);
			        fw.write(print_info);
			        fw.close();
				}
				catch (IOException e) {
		        	System.out.println("Cannot write diagnosis file");
		        	System.exit(0);
		        }
				Chemkin.writeChemkinInputFile(reactionSystem.getReactionModel(),reactionSystem.getPresentStatus());
                valid = reactionSystem.isModelValid();
        	}
        	reactionChanged = false;
        	//delt = reactionSystem.adjustTimeStep(delt);
        	lastT = (Temperature)currentT.clone();
        	lastP = (Pressure)currentP.clone();
        	//begin = end;
        	//end = begin.add(delt);
        	currentT = reactionSystem.getTemperature(begin);
        	currentP = reactionSystem.getPressure(begin);
        	conditionChanged = (!currentT.equals(lastT) || !currentP.equals(lastP));
			//Remove from this to the end of the comments when you figure out the exact problem with daspk...sandeep
			begin=init;
			//reactionSystem.resetSystemSnapshot();
			end = end.add(delt);
			
        	reactionSystem.solveReactionSystem(begin, end, false, reactionChanged, conditionChanged);
			
        	terminated = reactionSystem.isReactionTerminated();
        	valid = reactionSystem.isModelValid();
        	if (!valid) {
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

        		runTime.gc();

        		System.out.println("After garbage collection:");
        		System.out.print("Memory used: ");
        		System.out.println(runTime.totalMemory());
        		System.out.print("Free memory: ");
        		System.out.println(runTime.freeMemory());
				System.out.println("The model core has " + ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedReactionSet().size() + " reactions and "+ ((CoreEdgeReactionModel)reactionSystem.getReactionModel()).getUnreactedSpeciesSet().size() + " species.");

        	}
        }
        if (paraInfor != 0){
          DynamicSimulator dynamicSimulator2 = new JDASPK(rtol, atol, paraInfor, initialStatus);
          reactionSystem.setDynamicSimulator(dynamicSimulator2);
          reactionSystem.resetSystemSnapshot();
          begin = init;
          end = begin.add(delt);
          terminated = false;
          reactionSystem.solveReactionSystem(begin, end, true, false, false);
          while (!terminated){
            //begin = end;
        	  begin = init;
            //end = begin.add(delt);
        	  end = end.add(delt);
            reactionSystem.solveReactionSystem(begin, end, false, false,
                                               false);
            terminated = reactionSystem.isReactionTerminated();
          }
        }


        System.out.println(end);
        reactionSystem.cleanDynamicSimulator();

        return;





        //#]
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
		String restartFileContent ="";
		int reactionCount = 1;
		try{
			File coreSpecies = new File ("Restart/edgeReactions.txt");
			FileWriter fw = new FileWriter(coreSpecies);
			for(Iterator iter=((CoreEdgeReactionModel)reactionSystem.reactionModel).getUnreactedReactionSet().iterator();iter.hasNext();){

				Reaction reaction = (Reaction) iter.next();
				//if (reaction.getDirection()==1){
					//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
					restartFileContent = restartFileContent + reaction.toRestartString() + "\n";
					reactionCount = reactionCount + 1;
				//}
			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent);
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart edgereactions file");
        	System.exit(0);
		}

	}

	private void writeAllReactions() {
		String restartFileContent ="";
		int reactionCount = 1;
		try{
			File allReactions = new File ("Restart/allReactions.txt");
			FileWriter fw = new FileWriter(allReactions);
			for(Iterator iter=reactionSystem.reactionModel.getReaction();iter.hasNext();){

				Reaction reaction = (Reaction) iter.next();

				//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
				restartFileContent = restartFileContent + reaction.toRestartString() + "\n";

			}

			for(Iterator iter=((CoreEdgeReactionModel)reactionSystem.reactionModel).getUnreactedReactionSet().iterator();iter.hasNext();){

				Reaction reaction = (Reaction) iter.next();
				//if (reaction.getDirection()==1){
					//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
				restartFileContent = restartFileContent + reaction.toRestartString() + "\n";

			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent);
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart edgereactions file");
        	System.exit(0);
		}

	}

	private void writeEdgeSpecies() {
		String restartFileContent ="";
		int speciesCount = 0;
		try{
			File edgeSpecies = new File ("Restart/edgeSpecies.txt");
			FileWriter fw = new FileWriter(edgeSpecies);
			for(Iterator iter=((CoreEdgeReactionModel)reactionSystem.reactionModel).getUnreactedSpeciesSet().iterator();iter.hasNext();){
				Species species = (Species) iter.next();
				restartFileContent = restartFileContent + species.getID() + "\n";
				/*Species species = (Species) iter.next();
				restartFileContent = restartFileContent + "("+ speciesCount + ") "+species.getChemkinName() + " \n ";// + 0 + " (mol/cm3) \n";
				restartFileContent = restartFileContent + species.toString(1) + "\n";
				speciesCount = speciesCount + 1;*/
			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent);
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart edgespecies file");
        	System.exit(0);
		}
	}

	private void writeCoreReactions() {
		String restartFileContent ="";
		int reactionCount = 0;
		try{
			File coreSpecies = new File ("Restart/coreReactions.txt");
			FileWriter fw = new FileWriter(coreSpecies);
			for(Iterator iter=reactionSystem.reactionModel.getReaction();iter.hasNext();){

				Reaction reaction = (Reaction) iter.next();
				if (reaction.getDirection()==1){
					//restartFileContent = restartFileContent + "("+ reactionCount + ") "+species.getChemkinName() + "  " + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
					restartFileContent = restartFileContent + reaction.toRestartString() + "\n";
					reactionCount = reactionCount + 1;
				}
			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent);
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart corereactions file");
        	System.exit(0);
		}
	}

	private void writeCoreSpecies() {

		String restartFileContent ="";
		int speciesCount = 0;
		try{
			File coreSpecies = new File ("Restart/coreSpecies.txt");
			FileWriter fw = new FileWriter(coreSpecies);
			for(Iterator iter=reactionSystem.reactionModel.getSpecies();iter.hasNext();){
				Species species = (Species) iter.next();
				restartFileContent = restartFileContent + species.getID() + "\n";//("+ speciesCount + ") "+species.getChemkinName() + " \n ";// + reactionSystem.getPresentConcentration(species) + " (mol/cm3) \n";
				//restartFileContent = restartFileContent + species.toString(1) + "\n";
				//speciesCount = speciesCount + 1;
			}
			//restartFileContent += "\nEND";
			fw.write(restartFileContent);
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the restart corespecies file");
        	System.exit(0);
		}

	}



	public ReactionTime getTimeStep() {
        return timeStep;
    }

    public void setTimeStep(ReactionTime p_timeStep) {
        timeStep = p_timeStep;
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

    public void setReactionSystem(ReactionSystem p_ReactionSystem) {
        reactionSystem = p_ReactionSystem;
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ReactionModelGenerator.java
*********************************************************************/


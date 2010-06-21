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


import jing.rxn.*;
import jing.chem.*;
import java.util.*;
import jing.param.*;
import jing.chem.Species;
import jing.rxn.Reaction;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxn.ReactionGenerator;

//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\ReactionSystem.java
//----------------------------------------------------------------------------

/**
Under specification of the natrual reaction conditions, (Temperature, Pressure, Reactants and their initial concentrations, Conversion), ReactionSystem represents the whole reaction process.
*/
//## class ReactionSystem
public class ReactionSystem {

    protected LinkedHashSet originalReactant;		//## attribute originalReactant

    protected DynamicSimulator dynamicSimulator;
    protected FinishController finishController;
    protected InitialStatus initialStatus;
    protected PressureModel pressureModel;
    protected PrimaryKineticLibrary primaryKineticLibrary;//10/14/07 gmagoon: restored (I had commented out)
    protected ReactionGenerator reactionGenerator;//10/9/07 gmagoon: uncommented out for use in RateBasedPDepRME (I had commented out a week or two ago)
    protected ReactionModel reactionModel;
    protected ReactionModelEnlarger reactionModelEnlarger;
    protected LinkedList systemSnapshot;
    protected TemperatureModel temperatureModel;
    protected double [] reactionFlux;
    protected LibraryReactionGenerator lrg;//9/24/07 gmagoon: moved to ReactionModelGenerator.java; 9/25/07 variable is passed from ReactionModelGenerator
    protected int ind;//10/30/07 gmagoon: added index variable to number different reaction systems; starts at zero; 1/5/09 changed name from index to ind to avoid confusion with local variables used below
    protected String equationOfState; // rwest: could be "Liquid"
    protected static boolean printAllSens = false;
    // Constructors

    //## operation ReactionSystem(TemperatureModel,PressureModel,ReactionModelEnlarger,FinishController,DynamicSimulator,PrimaryKineticLibrary,ReactionGenerator,HashSet,InitialStatus)
    //9/24/07 gmagoon: reactionModel changed to parameter passed to class; setReactionModel method removed; 10/4/07: this was incorrect; setReactionModel restored
    //9/25/07 gmagoon: removed primaryKineticLibrary from parameters
    public  ReactionSystem(TemperatureModel p_temperatureModel, PressureModel p_pressureModel, ReactionModelEnlarger p_reactionModelEnlarger, FinishController p_finishController, DynamicSimulator p_dynamicSimulator, PrimaryKineticLibrary p_primaryKineticLibrary, ReactionGenerator p_reactionGenerator, LinkedHashSet p_speciesSeed, InitialStatus p_initialStatus, ReactionModel p_reactionModel, LibraryReactionGenerator p_libraryReactionGenerator, int p_index, String p_equationOfState) {
        {
            systemSnapshot=new LinkedList();
        }
        //#[ operation ReactionSystem(TemperatureModel,PressureModel,ReactionModelEnlarger,FinishController,DynamicSimulator,PrimaryKineticLibrary,ReactionGenerator,HashSet,InitialStatus)
        temperatureModel = p_temperatureModel;
        pressureModel = p_pressureModel;
        setFinishController(p_finishController);
        setReactionModel(p_reactionModel);//10/4/07 gmagoon: changed to use setReactionModel
        reactionModelEnlarger = p_reactionModelEnlarger;
        initialStatus = p_initialStatus;
        dynamicSimulator = p_dynamicSimulator;
        primaryKineticLibrary = p_primaryKineticLibrary;//10/14/07 gmagoon: restored (I had commented out)
        reactionGenerator = p_reactionGenerator; //10/4/07 gmagoon: no longer needed here;//10/9/07 gmagoon needed in RateBasedPDepRME
        originalReactant = p_speciesSeed;
        lrg = p_libraryReactionGenerator;
        systemSnapshot.add(initialStatus);
        ind = p_index;//10/30/07 gmagoon: added
        equationOfState=p_equationOfState;
		
        
        if (equationOfState=="Liquid") {
            System.out.println("Liquid phase: Not checking C=P/RT; assuming concentrations are specified correctly");
        }
		else if (!checkInitialConsistency()) {
        	System.out.println("Initial composition was not consistent: C = P/RT was not satisfied!");
        	System.out.println("The concentrations have been renormalized.");
        	//System.exit(-1);
        }

        //#]
    }
	
//	## operation checkInitialConsistency() 
    public boolean checkInitialConsistency() {
        //#[ operation checkInitialConsistency() 
        return initialStatus.isTPCConsistent(dynamicSimulator, finishController);
        //#]
    }
	
    public  ReactionSystem() {
        {
            systemSnapshot=new LinkedList();
        }
    }

    //## operation adjustTimeStep(ReactionTime)
    public ReactionTime adjustTimeStep(ReactionTime p_presentTS) {
        //#[ operation adjustTimeStep(ReactionTime)
        double time = p_presentTS.getTime();
        String unit = p_presentTS.getUnit();
        double maxCPS = -1;
        PresentStatus ps = getPresentStatus();
        for (Iterator iter = getOriginalReactant().iterator(); iter.hasNext(); ) {
          	Species reactant = (Species)iter.next();
            SpeciesStatus ss = ps.getSpeciesStatus(reactant);
            if (ss.getSpecies().getName().equals("C4H10")) {
        	    double convPerStep = -ss.getFlux()/ss.getConcentration();
            	maxCPS = convPerStep;
        	}
        }

        ReactionTime newStep = p_presentTS;
        if (maxCPS > 1e-2) {
         	newStep = new ReactionTime(10,unit);
        }
        else if (maxCPS > 1e-3) {
         	newStep = new ReactionTime(50,unit);
        }
        else {
         	newStep = new ReactionTime(100,unit);
        }

        return newStep;
          /*
          if (maxCPS == -1 || (maxCPS < upLim && maxCPS > lowLim)) return p_presentTS;

          if (maxCPS >= upLim) {
          	int order = (int)(Math.log(maxCPS/upLim)/Math.log(10));
          	ReactionTime newStep = new ReactionTime(time/Math.pow(10, order),unit);
          	System.out.println("Present CPS: " + String.valueOf(maxCPS));
          	System.out.println("order = " + String.valueOf(order));
          	System.out.println("new step = " + newStep);
          	return newStep;
          }
          else {
          	int order = (int)(Math.log(lowLim/maxCPS)/Math.log(10));
          	ReactionTime newStep = new ReactionTime(time*Math.pow(10, order),unit);
         	return newStep;
          }
        */
        /*double time = p_presentTS.getTime();
        String unit = p_presentTS.getUnit();
        double ratio = calculatePresentRadicalConcentration()/calculatePresentMoleculeConcentration();
        double largestRatio = 1E-4;

        if (ratio < largestRatio) return p_presentTS;
        else {
        	System.out.println("Radical/Molecule = "+ String.valueOf(ratio));
        	int order = (int)(Math.log(ratio/largestRatio)/Math.log(10));
        	double level = ratio/largestRatio;
        	double newT = 100;
        	if (level > 1000) newT /= 1.0E4;
        	else if (level > 100) newT /= 1.0E3;
        	else if (level > 10) newT /= 1.0E2;
        	else if (level > 1) newT /= 1.0E1;
        	System.out.println("New time step " + String.valueOf(newT));
            ReactionTime newStep = new ReactionTime(newT,unit);
            return newStep;
        }
        */







        //#]
    }

    //9/24/07: gmagoon: modified to include p_reactionModel as parameter; subsequently removed
    public void appendUnreactedSpeciesStatus(SystemSnapshot p_systemSnapshot, Temperature p_temperature) {
        
		double startTime = System.currentTimeMillis();
		
		if (!(reactionModel instanceof CoreEdgeReactionModel)) return;
        CoreEdgeReactionModel model = (CoreEdgeReactionModel)reactionModel;
        
		//double [] unreactedFlux = new double[SpeciesDictionary.getInstance().size()+1];
		double [] unreactedFlux = new double[model.getMaxSpeciesID()+1];//sizing based on MaxSpeciesID should be sufficient
                if (reactionModelEnlarger instanceof RateBasedPDepRME) {
                    for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
				PDepNetwork pdn = (PDepNetwork) iter.next();
				for (Iterator iter2 = pdn.getNetReactions().iterator(); iter2.hasNext(); ) {
					PDepReaction rxn = (PDepReaction) iter2.next();
					double forwardFlux = rxn.calculateForwardFlux(p_systemSnapshot);
					double reverseFlux = rxn.calculateReverseFlux(p_systemSnapshot);
					//System.out.println(rxn.toString() + ": " + forwardFlux + " " + reverseFlux);
					for (int j = 0; j < rxn.getReactantNumber(); j++) {
						Species species = (Species) rxn.getReactantList().get(j);
						if (model.containsAsUnreactedSpecies(species))
							unreactedFlux[species.getID()] += reverseFlux;
					}
					for (int j = 0; j < rxn.getProductNumber(); j++) {
						Species species = (Species) rxn.getProductList().get(j);
						if (model.containsAsUnreactedSpecies(species))
							unreactedFlux[species.getID()] += forwardFlux;
					}
				}
			}
// gmagoon 6/19/09: the method below was the original method for P-dep cases, but with the changes to the P-dep code, the iterator would only take into account forward reactions
//                    // first take all the unreacted reactions from PDepNetwork and calculate their rate
//                    // populate the reactionModel with all the unreacted species if they are already not there.
//                    for (Iterator iter=PDepNetwork.getNetworks().iterator(); iter.hasNext();){
//                            PDepNetwork pdn = (PDepNetwork) iter.next();
//                            Iterator reaction_iter = pdn.getNetReactions().listIterator();
//                            while (reaction_iter.hasNext()){
//                                    PDepReaction r = (PDepReaction) reaction_iter.next();
//                                    int rxnType = model.categorizeReaction(r);
//                                    //if (rxnType == 0) throw new InvalidReactionSetException();
//                                    if (rxnType == -1){
//                                                    SystemSnapshot ss = (SystemSnapshot) systemSnapshot.getLast();
//                                            double flux = r.calculateRate(ss.getTemperature(), ss.getPressure());
//                                                    for (Iterator rIter=r.getReactants(); rIter.hasNext();) {
//                                                    Species spe = (Species)rIter.next();
//                                                    double conc =0; 
//                                                    if (p_systemSnapshot.getSpeciesStatus(spe) != null)
//                                                             conc = (p_systemSnapshot.getSpeciesStatus(spe)).getConcentration();
//                                            if (conc<0)
//                                                    throw new NegativeConcentrationException(spe.getName() + ": " + String.valueOf(conc));
//                                        flux *= conc;
//
//                                    }
//                                            for (Iterator pIter = r.getProducts(); pIter.hasNext();){
//                                                    Species spe = (Species)pIter.next();
//                                                    unreactedFlux[spe.getID()] += flux;
//                                            }
//                                    }
//                            }
//                    }
                }
		
		LinkedHashSet ur = model.getUnreactedReactionSet();
            
		for (Iterator iur = ur.iterator(); iur.hasNext();) {
			Reaction r = (Reaction)iur.next();
			double flux = 0;
			if (r instanceof TemplateReaction) {
				//flux = ((TemplateReaction)r).calculateTotalPDepRate(p_temperature);
				flux = ((TemplateReaction)r).getRateConstant(p_temperature, p_systemSnapshot.getPressure());//10/26/07 gmagoon: changed to pass temperature and pressure; I assume systemSnapshot contains pressure at desired time; note: interestingly, this did not appear to turn up on my initial search for uses of templateReaction version of getRateConstant (I could have missed it by mistaking reaction version for the one I was searching for); also NetBeans did not indicate an error until I changed reaction version of getRateConstant
				//flux = ((TemplateReaction)r).getRateConstant();
			}
			else {
				// flux = r.calculateTotalRate(p_temperature);
				flux = r.getRateConstant(p_temperature);//10/26/07 gmagoon: changed to pass temperature
				// flux = r.getRateConstant();
			}
			if (flux > 0) {
				for (Iterator rIter=r.getReactants(); rIter.hasNext();) {
					Species spe = (Species)rIter.next();
					double conc = 0;
					if (p_systemSnapshot.getSpeciesStatus(spe) != null)
						conc = (p_systemSnapshot.getSpeciesStatus(spe)).getConcentration();
					if (conc<0) {
						double aTol = ReactionModelGenerator.getAtol();
						//if (Math.abs(conc) < aTol) conc = 0;
						//else throw new NegativeConcentrationException(spe.getName() + ": " + String.valueOf(conc));
						if (conc < -100.0 * aTol)
							throw new NegativeConcentrationException("Species " + spe.getName() + " has negative concentration: " + String.valueOf(conc));
					}
					flux *= conc;

				}

				for (Iterator rIter=r.getProducts(); rIter.hasNext();) {
					Species spe = (Species)rIter.next();
					if (model.containsAsUnreactedSpecies(spe)) {
						unreactedFlux[spe.getID()] += flux;
					}
				}

			}
			else {
				throw new NegativeRateException(r.toChemkinString(p_temperature) + ": " + String.valueOf(flux));
			}
		}
		p_systemSnapshot.unreactedSpeciesFlux = unreactedFlux;
		
        //#]
    }

    //## operation calculatePresentMoleculeConcentration()
    public double calculatePresentMoleculeConcentration() {
        //#[ operation calculatePresentMoleculeConcentration()
        double total = 0;

        PresentStatus ps = getPresentStatus();
        for (Iterator iter = ps.getSpeciesStatus(); iter.hasNext(); ) {
        	SpeciesStatus ss = (SpeciesStatus)iter.next();
        	Species spe = ss.getSpecies();
        	if (!spe.isRadical()) {
        		total += ss.getConcentration();
        	}
        }

        return total;

        //#]
    }

    //## operation calculatePresentRadicalConcentration()
    public double calculatePresentRadicalConcentration() {
        //#[ operation calculatePresentRadicalConcentration()
        double total = 0;

        PresentStatus ps = getPresentStatus();
        for (Iterator iter = ps.getSpeciesStatus(); iter.hasNext(); ) {
        	SpeciesStatus ss = (SpeciesStatus)iter.next();
        	Species spe = ss.getSpecies();
        	if (spe.isRadical()) {
        		total += ss.getConcentration();
        	}
        }

        return total;

        //#]
    }

    //9/25/07 gmagoon: moved to ReactionModelGenerator.java
    //## operation enlargeReactionModel()
//    public void enlargeReactionModel() {
//        //#[ operation enlargeReactionModel()
//        if (reactionModelEnlarger == null) throw new NullPointerException("ReactionModelEnlarger");
//
//       reactionModelEnlarger.enlargeReactionModel(this);
//
//        return;
//        //#]
//    }

    //## operation findMainChannel(Species,SystemSnapshot)
    public Reaction findMainChannel(Species p_species, SystemSnapshot p_systemSnapshot) {
        //#[ operation findMainChannel(Species,SystemSnapshot)
        ReactionTime rt = p_systemSnapshot.getTime();
        Temperature temp = getTemperature(rt);
        double maxFlux = 0;
        Reaction maxReaction = null;

        LinkedHashSet rs = getReactionModel().getReactionSet();
        for (Iterator iter = rs.iterator(); iter.hasNext(); ) {
        	Reaction rxn = (Reaction)iter.next();
        	double rflux = 0;
        	for (Iterator pIter = rxn.getProducts(); pIter.hasNext(); ) {
        		Species spe = ((ChemGraph)pIter.next()).getSpecies();
        		if (spe.equals(p_species)) {
        			double flux = rxn.calculateTotalRate(temp);
        			if (rxn instanceof TemplateReaction) {
        				flux = ((TemplateReaction)rxn).calculateTotalPDepRate(temp, getPressure(rt));//10/25/07 gmagoon: added pressure
        			}
        			else if (rxn instanceof ThirdBodyReaction) {
        				flux *= ((ThirdBodyReaction)rxn).calculateThirdBodyCoefficient(p_systemSnapshot);
        			}
        			for (Iterator rIter = rxn.getReactants(); rIter.hasNext(); ) {
        				Species reactant = ((ChemGraph)rIter.next()).getSpecies();
        				double concentration = p_systemSnapshot.getSpeciesStatus(reactant).getConcentration();
        				flux *= concentration;
        			}
        			rflux += flux;
        		}
        	}
        	if (rflux > maxFlux) {
        		maxFlux = rflux;
        		maxReaction = rxn;
        	}
        }

        System.out.println("The main pathway to generate " + p_species.getName() + " is ");
        System.out.println(maxReaction);
        System.out.println("The max flux is " + String.valueOf(maxFlux));

        return maxReaction;




        //#]
    }

    //## operation getColliders()
    public void getColliders() {
        //#[ operation getColliders()
        //#]
    }

    //## operation getInitialConcentration(Species)
    public double getInitialConcentration(Species p_species) {
        //#[ operation getInitialConcentration(Species)
        if (p_species == null) throw new NullPointerException();

        InitialStatus is = getInitialStatus();
        SpeciesStatus ss = is.getSpeciesStatus(p_species);

        return ss.getConcentration();
        //#]
    }

    //## operation getInitialReactionTime()
    public ReactionTime getInitialReactionTime() {
        //#[ operation getInitialReactionTime()
        return getInitialStatus().getTime();
        //#]
    }

    //## operation getPresentConcentration(Species)
    public double getPresentConcentration(Species p_species) {
        //#[ operation getPresentConcentration(Species)
        if (p_species == null) throw new NullPointerException();

        PresentStatus ps = getPresentStatus();
        SpeciesStatus ss = ps.getSpeciesStatus(p_species);

        return ss.getConcentration();
        //#]
    }

    //## operation getPresentConversion(Species)
    public double getPresentConversion(Species p_species) {
        //#[ operation getPresentConversion(Species)
        if (p_species == null) throw new NullPointerException();

        double C0 = getInitialConcentration(p_species);
        double C = getPresentConcentration(p_species);

       // if (C>C0) return -1;

        return (1-C/C0);
        //#]
    }

    //## operation getPresentPressure()
    public Pressure getPresentPressure() {
        //#[ operation getPresentPressure()
        return getPressure(getPresentStatus().getTime());
        //#]
    }

    //## operation getPresentStatus()
    public PresentStatus getPresentStatus() {
        //#[ operation getPresentStatus()
        SystemSnapshot ss = (SystemSnapshot)systemSnapshot.getLast();
        return new PresentStatus(ss);






        //#]
    }

    //## operation getPresentTemperature()
    public Temperature getPresentTemperature() {
        //#[ operation getPresentTemperature()
        return getTemperature(getPresentStatus().getTime());
        //#]
    }

    //## operation getPressure(ReactionTime)
    public Pressure getPressure(ReactionTime p_time) {
        //#[ operation getPressure(ReactionTime)
        return getPressureModel().getPressure(p_time);
        //#]
    }

    //## operation getRmin()
    public double getRmin() {
        //#[ operation getRmin()
        RateBasedVT vt = (RateBasedVT)(getFinishController().getValidityTester());
        double Rmin = vt.calculateRmin(getPresentStatus());
        return Rmin;
        //#]
    }

    //## operation getTemperature(ReactionTime)
    public Temperature getTemperature(ReactionTime p_time) {
        //#[ operation getTemperature(ReactionTime)
        return getTemperatureModel().getTemperature(p_time);
        //#]
    }

//9/25/07 gmagoon: moved to ReactionModelGenerator.java    
//    //## operation hasPrimaryKineticLibrary()
//    public boolean hasPrimaryKineticLibrary() {
//        //#[ operation hasPrimaryKineticLibrary()
//        if (primaryKineticLibrary == null) return false;
//        return (primaryKineticLibrary.size() > 0);
//        //#]
//    }

    //## operation identifyColliders()
    public HashMap identifyColliders() {
        //#[ operation identifyColliders()
        return getInitialStatus().identifyColliders();
        //#]
    }


    //## operation initializePDepNetwork()
    public void initializePDepNetwork() {
        if (!(reactionModelEnlarger instanceof RateBasedPDepRME)) {
			System.out.println("ERROR: Reaction model enlarger is not pressure-dependent!");
			System.exit(0);
		}
		
		CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) reactionModel;
		
		PDepKineticsEstimator pDepKineticsEstimator = 
				((RateBasedPDepRME) reactionModelEnlarger).getPDepKineticsEstimator();
		
		LinkedList pdnList = new LinkedList(PDepNetwork.getNetworks());
		for (Iterator iter = pdnList.iterator(); iter.hasNext(); ) {
        	PDepNetwork pdn = (PDepNetwork)iter.next();
        	if (pdn.getAltered()) {
				
        		// Update the k(T, P) estimates for the network
				pDepKineticsEstimator.runPDepCalculation(pdn, this, cerm);
				
				// Each net reaction with k(T, P) > 0 can be treated as a core or edge reaction (?)
				/*for (ListIterator<PDepReaction> iter2 = pdn.getNetReactions().listIterator(); iter2.hasNext(); ) {
					PDepReaction rxn = iter2.next();
					if (rxn.isCoreReaction())
						cerm.addReactedReaction(rxn);
					else if (rxn.isEdgeReaction())
						cerm.addUnreactedReaction(rxn);	
				}*/
			}
		}	
    }

    //## operation isFinished()
    public boolean isFinished() {
        //#[ operation isFinished()
        return (isModelValid() && isReactionTerminated());
        //#]
    }

    //## operation isModelValid()
    public boolean isModelValid() {
        //#[ operation isModelValid()
        return finishController.isModelValid();
        //#]
    }

    //## operation isReactionTerminated()
    public boolean isReactionTerminated() {
        //#[ operation isReactionTerminated()
        return finishController.isReactionTerminated();
        //#]
    }

    //## operation outputAllPathways(Species,LinkedList,SystemSnapshot,Temperature)
    /*public static void outputAllPathways(Species p_species, LinkedList p_reactionList, SystemSnapshot p_systemSnapshot, Temperature p_temperature) {
        //#[ operation outputAllPathways(Species,LinkedList,SystemSnapshot,Temperature)
        ReactionTime rt = p_systemSnapshot.getTime();
        Temperature temp = p_temperature;
        double maxFlux = 0;
        Reaction maxReaction = null;

        System.out.println("the consumption paths for " + p_species.getName());
        for (Iterator iter = p_reactionList.iterator(); iter.hasNext(); ) {
        	Reaction rxn = (Reaction)iter.next();
        	if (rxn.containsAsReactant(p_species)) {
        		double flux;
        		if (rxn instanceof TemplateReaction) {
        			flux = ((TemplateReaction)rxn).calculateTotalPDepRate(temp);
        		}
        		else if (rxn instanceof PDepNetReaction) {
        			flux = ((PDepNetReaction)rxn).calculateRate();
        		}
        		else {
        			flux = rxn.calculateTotalRate(temp);
        		}

        		if (rxn instanceof ThirdBodyReaction) {
        			flux *= ((ThirdBodyReaction)rxn).calculateThirdBodyCoefficient(p_systemSnapshot);
        		}

        		System.out.print(rxn.getStructure().toString() + '\t' + String.valueOf(flux));
        		for (Iterator rIter = rxn.getReactants(); rIter.hasNext(); ) {
        			Species reactant = ((ChemGraph)rIter.next()).getSpecies();
        			double concentration = p_systemSnapshot.getSpeciesStatus(reactant).getConcentration();
        			System.out.print('\t' + String.valueOf(concentration));
        			flux *= concentration;
        		}
        		System.out.println('\t' + String.valueOf(-flux));

        	}
        }
        System.out.println("the formationtion paths for " + p_species.getName());
        for (Iterator iter = p_reactionList.iterator(); iter.hasNext(); ) {
        	Reaction rxn = (Reaction)iter.next();
        	if (rxn.containsAsProduct(p_species)) {
        		double flux;
        		if (rxn instanceof TemplateReaction) {
        			flux = ((TemplateReaction)rxn).calculateTotalPDepRate(temp);
        		}
        		else if (rxn instanceof PDepNetReaction) {
        			flux = ((PDepNetReaction)rxn).getRate();
        		}
        		else {
        			flux = rxn.calculateTotalRate(temp);
        		}

        		if (rxn instanceof ThirdBodyReaction) {
        			flux *= ((ThirdBodyReaction)rxn).calculateThirdBodyCoefficient(p_systemSnapshot);
        		}
        		System.out.print(rxn.getStructure().toString() + '\t' + String.valueOf(flux));
        		for (Iterator rIter = rxn.getReactants(); rIter.hasNext(); ) {
        			Species reactant = ((ChemGraph)rIter.next()).getSpecies();
        			double concentration = p_systemSnapshot.getSpeciesStatus(reactant).getConcentration();
        			System.out.print('\t' + String.valueOf(concentration));
        			flux *= concentration;
        		}
        		System.out.println('\t' + String.valueOf(flux));

        	}
        }

        return;




        //#]
    }*/

    //## operation outputReactionFlux(SystemSnapshot)
    public void outputReactionFlux(SystemSnapshot p_systemSnapshot) {
        //#[ operation outputReactionFlux(SystemSnapshot)
        ReactionTime rt = p_systemSnapshot.getTime();
        Temperature t = getTemperature(rt);

        HashSet speSet = new HashSet();
        for (Iterator iter = getReactionModel().getReaction(); iter.hasNext(); ) {
        	Reaction rxn = (Reaction)iter.next();
        	double k = rxn.calculateTotalRate(t);
        	if (rxn instanceof ThirdBodyReaction) {
        		k *= ((ThirdBodyReaction)rxn).calculateThirdBodyCoefficient(p_systemSnapshot);
        	}
        	double flux = k;
        	for (Iterator rIter = rxn.getReactants(); rIter.hasNext(); ) {
        		Species spe = ((ChemGraph)rIter.next()).getSpecies();
        		double concentration = p_systemSnapshot.getSpeciesStatus(spe).getConcentration();
        		flux *= concentration;
        	}
        	if (flux>1E-7) {
        		System.out.println(rxn.toString() + "\t" + String.valueOf(flux) + '\t' + String.valueOf(rxn.calculateHrxn(t)) + '\t' + String.valueOf(rxn.calculateKeq(t)));
        		for (double temp = 400; temp<1200; temp = temp + 50) {
        			double rate = rxn.calculateTotalRate(new Temperature(temp,"K"));
        			System.out.println("Temp = " + String.valueOf(temp) + "\tRate = " + String.valueOf(rate));
        		}
        		for (Iterator rIter = rxn.getReactants(); rIter.hasNext(); ) {
        			Species spe = ((ChemGraph)rIter.next()).getSpecies();
        			speSet.add(spe);
        		}
                for (Iterator rIter = rxn.getProducts(); rIter.hasNext(); ) {
        			Species spe = ((ChemGraph)rIter.next()).getSpecies();
        			speSet.add(spe);
        		}
        	}
        }

        for (Iterator iter = speSet.iterator(); iter.hasNext(); ) {
        	Species spe = (Species)iter.next();
        	System.out.println(spe.getName()+"("+String.valueOf(spe.getID())+"): " + spe.getThermoData().toString());
        }






        //#]
    }

    //## operation printConcentrationProfile(LinkedList)
    public String printConcentrationProfile(LinkedList p_speciesList) {
        //#[ operation printConcentrationProfile(LinkedList)
        if (p_speciesList == null) throw new NullPointerException();


        if (p_speciesList.isEmpty()) return "EMPTY species list";

        // check the validity of p_speciesList and print the title line
        System.out.print("Time");
        int size = p_speciesList.size();
        for (int i=0; i<size; i++) {
        	Species spe = (Species)p_speciesList.get(i);
        	if (!spe.repOk()) throw new InvalidSpeciesException();
         	String name = spe.getChemkinName();
         	System.out.print('\t' + name);
        }
        System.out.println();

        Iterator iter = getSystemSnapshot();
        while (iter.hasNext()) {
        	SystemSnapshot ss = (SystemSnapshot)iter.next();
        	System.out.print(String.valueOf(ss.getTime().getTime()));
        	for (int i=0; i<size; i++) {
        		Species spe = (Species)p_speciesList.get(i);

         		if (spe != null) {
         			SpeciesStatus speSta = ss.getSpeciesStatus(spe);
         			double conc = 0;
         			if (speSta != null) conc = speSta.getConcentration();
         			System.out.print('\t' + String.valueOf(conc));
         		}
         	}
         	System.out.println();
        }

        return "END";
        //#]
    }

//	## operation printConcentrationProfile(LinkedList)
    public String returnConcentrationProfile(LinkedList p_speciesList) {
        //#[ operation printConcentrationProfile(LinkedList)
        if (p_speciesList == null) throw new NullPointerException();


        if (p_speciesList.isEmpty()) return "EMPTY species list";
		String output = "";
        // check the validity of p_speciesList and print the title line
        output = output + "Time";
        int size = p_speciesList.size();
        for (int i=0; i<size; i++) {
        	Species spe = (Species)p_speciesList.get(i);
        	if (!spe.repOk()) throw new InvalidSpeciesException();
         	String name = spe.getChemkinName();
			output = output + '\t' + name;
        }
		output = output + "\n";
        //System.out.println();

        Iterator iter = getSystemSnapshot();
        while (iter.hasNext()) {
        	SystemSnapshot ss = (SystemSnapshot)iter.next();
			output = output + String.valueOf(ss.getTime().getTime());
        	for (int i=0; i<size; i++) {
        		Species spe = (Species)p_speciesList.get(i);

         		if (spe != null) {
         			SpeciesStatus speSta = ss.getSpeciesStatus(spe);
         			double conc = 0;
         			if (speSta != null) conc = speSta.getConcentration();
					output = output + '\t' + String.valueOf(conc) ;
         		}
         	}
			output = output + "\n";
         	//System.out.println();
        }

        return output;
        //#]
    }
	
    //## operation printLowerBoundConcentrations(LinkedList)
//svp
public String printLowerBoundConcentrations(LinkedList p_speciesList) {
  //#[ operation printLowerBoundConcentrations(LinkedList)
  String result = "Time";
  int size = p_speciesList.size();
  for (int i = 0; i < size; i++) {
    Species spe = (Species) p_speciesList.get(i);
    String name = spe.getName();
    result += '\t' + name;
  }
  result += '\n';
  Iterator iter = getSystemSnapshot();
  while (iter.hasNext()) {
    SystemSnapshot ss = (SystemSnapshot) iter.next();
    result += String.valueOf(ss.getTime().getTime());
    for (int i = 0; i < size; i++) {
      Species spe = (Species) p_speciesList.get(i);
      if (spe != null) {
        SpeciesStatus speSta = ss.getSpeciesStatus(spe);
        double conc = 0;
        if (speSta != null) {
          conc = speSta.getConcentration();
        }
        double uncertainty = 0;
        if (ss.getTime().getTime() != 0) {
          int I = ss.getRealID(spe);
          LinkedList reactionList = ss.getReactionList();
          for (int j = 0; j < reactionList.size(); j++) {
            Reaction r = (Reaction) reactionList.get(j);
            ReactionTime rt = ss.getTime();
            Temperature t = getTemperature(rt);
            double k;
            double k_lowerbound;
            double k_upperbound;
            if (r instanceof ThirdBodyReaction){
              k = ((ThirdBodyReaction)r).calculateRate(ss);
              k_lowerbound = r.calculateLowerBoundRate(t)*((ThirdBodyReaction)r).calculateThirdBodyCoefficient(ss);
              k_upperbound = r.calculateUpperBoundRate(t)*((ThirdBodyReaction)r).calculateThirdBodyCoefficient(ss);
            }
            else if (r instanceof TROEReaction){
              k = ((TROEReaction)r).calculateRate(ss);
              k_lowerbound = r.calculateLowerBoundRate(t)*((TROEReaction)r).calculateTroeFallOff(ss);
              k_upperbound = r.calculateUpperBoundRate(t)*((TROEReaction)r).calculateTroeFallOff(ss);
            }
            else if (r instanceof LindemannReaction) {
            	k = ((LindemannReaction)r).calculateRate(ss);
            	k_lowerbound = r.calculateLowerBoundRate(t)*((LindemannReaction)r).calculateLindemannFallOff(ss);
            	k_upperbound = r.calculateUpperBoundRate(t)*((LindemannReaction)r).calculateLindemannFallOff(ss);
            }
           if (r instanceof TemplateReaction){
             k = ((TemplateReaction)r).calculateTotalPDepRate(ss.getTemperature(), ss.getPressure());//10/25/07 gmagoon: added pressure
             k_lowerbound = k/2.0;
             k_upperbound = k*2.0;
           }
           else if (r instanceof PDepReaction){
             k = ((PDepReaction) r).calculateRate(ss.getTemperature(),ss.getPressure());
             k_lowerbound = k/2.0;
             k_upperbound = k*2.0;
           }
           else {
             k = r.calculateTotalRate(ss.getTemperature());
             k_lowerbound = r.getLowerBoundRate(t);
             k_upperbound = r.getUpperBoundRate(t);
           }

            int index = I + j * getReactionModel().getSpeciesNumber() - 1;
            
              double sens = ss.getSensitivityStatus(index);
              double delta_k;
              if (sens > 0){
                delta_k = k - k_lowerbound;
              }
              else{
                delta_k = k_upperbound-k;
              }
              uncertainty -= Math.abs(sens*delta_k);
            
          }
        }
        if (conc != 0) {
          double lower_bound = conc * Math.exp(uncertainty / conc);
          result += '\t' + String.valueOf(lower_bound);
        }
        else {
          double lower_bound = 0;
          result += '\t' + String.valueOf(lower_bound);
        }
      }
    }
    result += '\n';
  }
  result += "END\n";
  return result;
  //#]
}


    //## operation printMoleFractionProfile(LinkedList)
    public String printMoleFractionProfile(LinkedList p_speciesList) {
        //#[ operation printMoleFractionProfile(LinkedList)
        if (p_speciesList == null) throw new NullPointerException();

        if (p_speciesList.isEmpty()) return "EMPTY species list";

        // check the validity of p_speciesList and print the title line
        System.out.print("Time");
        int size = p_speciesList.size();
        for (int i=0; i<size; i++) {
        	Species spe = (Species)p_speciesList.get(i);
        	if (spe != null) {
        		if (!spe.repOk()) throw new InvalidSpeciesException();
        	 	String name = spe.getChemkinName();
        	 	System.out.print('\t' + name);
        	 }
        }
        System.out.println();

        Iterator iter = getSystemSnapshot();
        while (iter.hasNext()) {
        	SystemSnapshot ss = (SystemSnapshot)iter.next();
        	double totalMole = ss.getTotalMole();
        	System.out.print(String.valueOf(ss.getTime().getTime()));
        	for (int i=0; i<size; i++) {
        		Species spe = (Species)p_speciesList.get(i);
        		if (spe != null) {
        	 		SpeciesStatus speSta = ss.getSpeciesStatus(spe);
        	 		double mf;
        	 		if (speSta==null) mf = 0;
        	 		else mf = speSta.getConcentration()/totalMole;
                                System.out.print('\t' + String.valueOf(mf));
        	 	}
         	}
         	System.out.println();
        }

        return "END";
        //#]
    }

	 //## operation printMoleFractionProfile(LinkedList)
    public String returnMoleFractionProfile(LinkedList p_speciesList) {
        //#[ operation printMoleFractionProfile(LinkedList)
        if (p_speciesList == null) throw new NullPointerException();

        if (p_speciesList.isEmpty()) return "EMPTY species list";
		String output = "";
        // check the validity of p_speciesList and print the title line
		output = output + "Time";
        int size = p_speciesList.size();
        for (int i=0; i<size; i++) {
        	Species spe = (Species)p_speciesList.get(i);
        	if (spe != null) {
        		if (!spe.repOk()) throw new InvalidSpeciesException();
        	 	String name = spe.getChemkinName();
				output = output + '\t' + name ;
        	 }
        }
		output = output +"\n";
        //System.out.println();

        Iterator iter = getSystemSnapshot();
        while (iter.hasNext()) {
        	SystemSnapshot ss = (SystemSnapshot)iter.next();
        	double totalMole = ss.getTotalMole();
			output = output + String.valueOf(ss.getTime().getTime());
        	for (int i=0; i<size; i++) {
        		Species spe = (Species)p_speciesList.get(i);
        		if (spe != null) {
        	 		SpeciesStatus speSta = ss.getSpeciesStatus(spe);
        	 		double mf;
        	 		if (speSta==null) mf = 0;
        	 		else mf = speSta.getConcentration()/totalMole;
					output = output + '\t' + String.valueOf(mf);
        	 	}
         	}
			output = output + "\n";
         	//System.out.println();
        }

        return output;
        //#]
    }
	
    //## operation printMoleFractionProfile(LinkedList)
    public String returnReactionFlux() {
        //#[ operation printMoleFractionProfile(LinkedList)
       
		StringBuilder output = new StringBuilder("");
        // check the validity of p_speciesList and print the title line
		output.append("///////ReactionFlux//////// \n \t");
      
		Iterator iter = getSystemSnapshot();
		while(iter.hasNext()){
			SystemSnapshot ss = (SystemSnapshot)iter.next();
			if (ss.getTime().time == 0.0) continue;
			output.append(ss.getTime() + "\t");
		}
		output.append("\n");	
		LinkedList reactionSet = ((SystemSnapshot)getSystemSnapshotEnd().next()).reactionList;
		LinkedList uniqueReactions = ((SystemSnapshot)getSystemSnapshotEnd().next()).getUniqueReactionList();
		
		for (int i=0; i<uniqueReactions.size(); i++) {
			Reaction rxn = (Reaction) uniqueReactions.get(i);
			int index = reactionSet.indexOf(rxn);
			if (index < 0 || index >= reactionSet.size())
				continue;
			
			iter = getSystemSnapshot();
			output.append("reaction " + (i+1) + '\t' );
			while (iter.hasNext()) {
				SystemSnapshot ss = (SystemSnapshot)iter.next();
				if (ss.getTime().time == 0.0) continue;
        	
				if (i < ss.reactionFlux.length)
					output.append( ss.reactionFlux[index] + "\t");
				else {
					System.out.println("Warning: Size of reaction set does not match number of reaction fluxes. Expected reaction flux missing.");
					output.append("\n");
					return output.toString();
				}
				
         	}
			output.append("\n");
        }

        return output.toString();
        //#]
    }
    
	 //## operation printMoleFractionProfile(LinkedList)
    /*public String returnReactionFlux() {
        //#[ operation printMoleFractionProfile(LinkedList)
       
		StringBuilder output = new StringBuilder("");
        // check the validity of p_speciesList and print the title line
		output.append("///////ReactionFlux//////// \n");
       

        Iterator iter = getSystemSnapshot();
        while (iter.hasNext()) {
        	SystemSnapshot ss = (SystemSnapshot)iter.next();
        	if (ss.getTime().time == 0.0) continue;
        	
        	
			output.append("\n Time: " + String.valueOf(ss.getTime().getTime())+"\n");
        	for (int i=0; i<ss.reactionList.size(); i++) {
        		
        		output.append("reaction " + (i+1) + '\t' + ss.reactionFlux[i] + "\n");
        	 	
         	}
			
        }

        return output.toString();
        //#]
    }*/
    
  //## operation printMostUncertainReactions(LinkedList, LinkedList)
//svp
      public String printMostUncertainReactions(LinkedList p_speciesList, LinkedList p_importantSpecies){
        //#[ operation printMostUncertainReactions(LinkedList, LinkedList)
        if (printAllSens){//12/22/09 gmagoon: make the important species equal to all the species if printAllSens is turned on
              for (int i = 0; i < p_speciesList.size(); i++) {
                Species spe = (Species) p_speciesList.get(i);; 
                p_importantSpecies.add(spe.getName());
              }
        }
        String result = "Reactions contributing most to uncertainty:\n";
        int size = p_speciesList.size();
        int n = 0;
        Iterator iter = getSystemSnapshot();
        while (iter.hasNext()) {
          SystemSnapshot ss = (SystemSnapshot) iter.next();
          if (n == 0) {
            if (ss.getTime().getTime() != 0){
            LinkedList reactionList = ss.getReactionList();
            n = 1;
          }
        }
          if (ss.getTime().getTime() != 0) {
            result += '\n';
            result += "Time:" + String.valueOf(ss.getTime().getTime())+'\n';
          }
          for (int i = 0; i < size; i++) {
            Species spe = (Species) p_speciesList.get(i);
            for (int s = 0; s <p_importantSpecies.size();s++){
              String name = (String)p_importantSpecies.get(s);

            if (spe.getName().equalsIgnoreCase(name)) {
              SpeciesStatus speSta = ss.getSpeciesStatus(spe);
              double conc = 0;
              if (speSta != null) {
                conc = speSta.getConcentration();
              }
              if (ss.getTime().getTime() != 0) {
                result += '\n';
                result += spe.getName();
                LinkedList Uncertainties = new LinkedList();
                int I = ss.getRealID(spe);
                LinkedList reactionList = ss.getReactionList();
                for (int j = 0; j < reactionList.size(); j++) {
                  Reaction r = (Reaction) reactionList.get(j);
                  ReactionTime rt = ss.getTime();
                  Temperature t = getTemperature(rt);
                  double k;
                  double k_upperbound;
                  if (r instanceof ThirdBodyReaction){
                    k = ((ThirdBodyReaction)r).calculateRate(ss);
                    k_upperbound = r.calculateUpperBoundRate(t)*((ThirdBodyReaction)r).calculateThirdBodyCoefficient(ss);
                  }
                  else if (r instanceof TROEReaction){
                    k = ((TROEReaction)r).calculateRate(ss);
                    k_upperbound = r.calculateUpperBoundRate(t)*((TROEReaction)r).calculateTroeFallOff(ss);
                  }
                  else if (r instanceof LindemannReaction) {
                	  k = ((LindemannReaction)r).calculateRate(ss);
                	  k_upperbound = r.calculateUpperBoundRate(t)*((LindemannReaction)r).calculateLindemannFallOff(ss);
                  }
                  if (r instanceof TemplateReaction){
                    k = ((TemplateReaction)r).calculateTotalPDepRate(ss.getTemperature(), ss.getPressure());//10/25/07 gmagoon: added pressure
                    k_upperbound = k*2.0;
                   }
                  else if (r instanceof PDepReaction){
                    SystemSnapshot ssnap = (SystemSnapshot) systemSnapshot.getLast();
					k = ((PDepReaction)r).calculateRate(ssnap.getTemperature(), ssnap.getPressure());//10/25/07 gmagoon: using last systemSnapshot as parameter to avoid use of Global.temperature and Global.pressure; ****is this correct handling of systemSnapshot with use of getLast?
                    k_upperbound = k*2.0;
                  }
                  else {
                    k = r.calculateTotalRate(ss.getTemperature());
                    k_upperbound = r.getUpperBoundRate(t);
                  }
                  int index = I + j * getReactionModel().getSpeciesNumber() - 1;
                  //SensitivityStatus sen_status = ss.getSensitivityStatus(index);
                  double sens=ss.getSensitivityStatus(index);
                  double delta_ln_k = 0;
                  
                    
                    int J = j+1;
                    sens /= conc;
                    sens *= k;
                    delta_ln_k = Math.abs(Math.log(k_upperbound)-Math.log(k));

                  
                  UncertainReaction u = new UncertainReaction(r,sens,delta_ln_k);
                  if (Uncertainties.size() < 5){
                    if (Uncertainties.isEmpty()){
                      Uncertainties.add(u);
                    }
                    else if (Uncertainties.size() == 1){
                      UncertainReaction r1 = (UncertainReaction)Uncertainties.getFirst();
                      double u1 = r1.getProduct();
                      if (Math.abs(sens*delta_ln_k) > u1){
                        Uncertainties.addFirst(u);
                      }
                      else {
                        Uncertainties.add(u);
                      }
                    }
                    else{
                      for (int m = 0; m < Uncertainties.size() - 1; m++){
                        UncertainReaction r1 = (UncertainReaction)Uncertainties.get(m);
                        UncertainReaction r2 = (UncertainReaction)Uncertainties.get(m+1);
                        double u1 = r1.getUncertainty();
                        double u2 = r2.getUncertainty();
                        if (Math.abs(sens*delta_ln_k) > u1 && Math.abs(sens*delta_ln_k) > u2){
                          Uncertainties.add(m,u);
                          break;
                        }
                        else if (Math.abs(sens*delta_ln_k) < u1 && Math.abs(sens*delta_ln_k) > u2){
                          Uncertainties.add(m+1,u);
                          break;
                        }

                      }
                      if (!Uncertainties.contains(u)){
                        Uncertainties.add(u);
                      }

                    }
                  }
                  else{
                    for (int p = 0; p < Uncertainties.size()-1; p++){
                      UncertainReaction r1 = (UncertainReaction)Uncertainties.get(p);
                      double u1 = r1.getProduct();
                      UncertainReaction r2 = (UncertainReaction)Uncertainties.get(p+1);
                      double u2 = r2.getProduct();
                      if (Math.abs(sens*delta_ln_k) > u1 && Math.abs(sens*delta_ln_k) > u2){
                         Uncertainties.add(p,u);
                         Uncertainties.removeLast();
                         break;
                       }
                       else if (Math.abs(sens*delta_ln_k) < u1 && Math.abs(sens*delta_ln_k) > u2){
                         Uncertainties.add(p+1,u);
                         Uncertainties.removeLast();
                         break;
                       }

                    }
                  }
                }
                for (int p = 0; p < Uncertainties.size(); p++){
                  UncertainReaction r1 = (UncertainReaction)Uncertainties.get(p);
                  double u1 = r1.getProduct();
                  double s1 = r1.getSensitivity();
                  double k1 = r1.getUncertainty();
                  Reaction rxn1 = r1.getReaction();
                  if (rxn1 instanceof PDepReaction) {
                    result += rxn1.getStructure().toString()+'\n';
                  }
                  else{
                    result += rxn1.toString() +'\n';
                  }
                  result += "dln["+name+"]/d(lnkj): "+s1+'\n';
                  result += "delta(ln kj): "+k1+'\n';
                  result += "delta(ln kj)*dln["+name+"]/d(ln kj): "+u1+'\n';
                }
              }
            }
          }
          }
        }
       return result;
       //#]
      }
      
      public String printOrderedReactions() {
    	  StringBuilder result = new StringBuilder("Reactions: \n");
    	  Iterator iter = getSystemSnapshot();
    	  SystemSnapshot ss = (SystemSnapshot) iter.next();//10/25/07 gmagoon (?) why is this done twice
    	  ss = (SystemSnapshot) iter.next();
    	  LinkedList reactionList = ss.getUniqueReactionList();
    	  for (int j = 0; j < reactionList.size(); j++) {
		    Reaction r = (Reaction) reactionList.get(j);
    		  if (r instanceof PDepReaction) {
    			  int J = j+1;
    			  result.append(J + ". " + r.getStructure().toString()+'\n');
    		  }
    		  else {
    			  int J = j+1;
                           result.append(J + ". " + r.toString(ss.getTemperature())+'\n');
    			  //result.append(J + ". " + r.toString(Global.temperature)+'\n');//10/25/07 gmagoon: changed from using global.temperature
    		  }
			}
    	return result.toString();
      }

    //## operation printSensitivityCoefficients(LinkedList, LinkedList)
      //svp
      public String printSensitivityCoefficients(LinkedList p_speciesList, LinkedList p_importantSpecies) {
        //#[ operation printSensitivityCoefficients(LinkedList, LinkedList)
          if (printAllSens){//12/22/09 gmagoon: make the important species equal to all the species if printAllSens is turned on
              for (int i = 0; i < p_speciesList.size(); i++) {
                Species spe = (Species) p_speciesList.get(i);; 
                p_importantSpecies.add(spe.getName());
              }
          }
        
        int size = p_speciesList.size();
        int n = 0;
        Iterator iter = getSystemSnapshot();
        StringBuilder result = new StringBuilder("\n");
        while (iter.hasNext()) {
        	SystemSnapshot ss = (SystemSnapshot) iter.next();
        	
        	result.append("\nSensitivities:\n");
        	
        	if (ss.getTime().getTime() != 0) {
        		result.append('\n');
        		result.append("Time:" + String.valueOf(ss.getTime().getTime())+'\n');
        	}
        	for (int i = 0; i < size; i++) {
        		Species spe = (Species) p_speciesList.get(i);
        		for (int s = 0; s < p_importantSpecies.size();s++){
        			String name = (String)p_importantSpecies.get(s);
        			if (spe.getName().equalsIgnoreCase(name)) {
        				SpeciesStatus speSta = ss.getSpeciesStatus(spe);
        				double conc = 0;
        				if (speSta != null) {
        					conc = speSta.getConcentration();
        				}
        				if (ss.getTime().getTime() != 0) {
        					int I = ss.getRealID(spe);
                 
        					LinkedList reactionList = ss.getReactionList();
        					for (int j = 0; j < reactionList.size(); j++) {
        						if (j < reactionList.size()){
        							Reaction r = (Reaction) reactionList.get(j);
        							ReactionTime rt = ss.getTime();
        							Temperature t = getTemperature(rt);
        							double k;
        							if (r instanceof TemplateReaction) {
        								k = ( (TemplateReaction) r).calculateTotalPDepRate(ss.getTemperature(), ss.getPressure());//10/25/07 gmagoon: added pressure
        							}
        							else if (r instanceof PDepReaction) {
        								SystemSnapshot ssnap = (SystemSnapshot)systemSnapshot.getLast();
        								k = ( (PDepReaction) r).calculateRate(ssnap.getTemperature(), ssnap.getPressure());//10/25/07 gmagoon: using last systemSnapshot as parameter to avoid use of Global.temperature and Global.pressure; ****is this correct handling of systemSnapshot with use of getLast?
        							}
        							else if (r instanceof ThirdBodyReaction) {
        								k = ( (ThirdBodyReaction) r).calculateRate(ss);
        							}
        							else if (r instanceof TROEReaction) {
        								k = ( (TROEReaction) r).calculateRate(ss);
        							}
        							else if (r instanceof LindemannReaction) {
        								k = ((LindemannReaction)r).calculateRate(ss);
        							}
        							else {
        								k = r.calculateTotalRate(ss.getTemperature());
        							}
        							int index = I + j * getReactionModel().getSpeciesNumber() - 1;
        							double sens = ss.getSensitivityStatus(index);
        							int J = j + 1;
        							sens /= conc;
        							sens *= k;
        							result.append("d(ln[" + spe.getName() + "])/d(lnk" + J +"): \t" +sens +'\n');
        						}
                   
        					}
        				}
        			}
        		}
        	}
        }
       return result.toString();
       //#]
      }

//## operation printSensitivityToThermo(LinkedList, LinkedList)
      //svp
      public String printSensitivityToThermo(LinkedList p_speciesList, LinkedList p_importantSpecies) {
        //#[ operation printSensitivityToThermo(LinkedList, LinkedList)
          if (printAllSens){//12/22/09 gmagoon: make the important species equal to all the species if printAllSens is turned on
              for (int i = 0; i < p_speciesList.size(); i++) {
                Species spe = (Species) p_speciesList.get(i);; 
                p_importantSpecies.add(spe.getName());
              }
          }
          
    	  int size = p_speciesList.size();
    	  String result = "Sensitivity to thermo:\n";
    	  int n = 0;
    	  Iterator iter = getSystemSnapshot();
    	  while (iter.hasNext()) {
    		  SystemSnapshot ss = (SystemSnapshot) iter.next();
    		  if (ss.getTime().getTime() != 0) {
    			  result += '\n';
    			  result += "Time:" + String.valueOf(ss.getTime().getTime())+'\n';
    		  }
    		  for (int i = 0; i < size; i++) {
    			  Species spe = (Species) p_speciesList.get(i);
    			  if (!p_importantSpecies.contains(spe.getName()))
    				  continue;
    			  SpeciesStatus speSta = ss.getSpeciesStatus(spe);
    			  for (int x = 0; x < size; x++){
    				  Species spe2 = (Species)p_speciesList.get(x);
    				  
    				  double conc = 0;
    				  if (speSta != null) {
    					  conc = speSta.getConcentration();
    				  }
    				  if (ss.getTime().getTime() != 0) {              	  
    					  int I = ss.getRealID(spe);
    					  int j = ss.getRealID(spe2);
    					  int index = I + (j-1+ss.reactionList.size()) * getReactionModel().getSpeciesNumber() -1;
                     
    					  double sens = ss.getSensitivityStatus(index);   					                       
    				  
    					  result += "d(ln["+spe.getName()+"])/d(delta_Gf("+spe2.getChemkinName()+")): "+sens/conc+'\n';
    				  }
    			  }
    		  }        
    	  }
       return result;
       //#]
      }

      //## operation printUpperBoundConcentrations(LinkedList)
      //svp
      public String printUpperBoundConcentrations(LinkedList p_speciesList) {
        //#[ operation printUpperBoundConcentrations(LinkedList)
        //System.out.print("Time");
        String result = "Time";//svp
        int size = p_speciesList.size();
        for (int i = 0; i < size; i++) {
          Species spe = (Species) p_speciesList.get(i);
          String name = spe.getName();
          result += '\t' + name;
        }
        result += '\n';
        Iterator iter = getSystemSnapshot();
        while (iter.hasNext()) {
          SystemSnapshot ss = (SystemSnapshot) iter.next();
          result += String.valueOf(ss.getTime().getTime());
          for (int i = 0; i < size; i++) {
            Species spe = (Species) p_speciesList.get(i);

            if (spe != null) {
              SpeciesStatus speSta = ss.getSpeciesStatus(spe);
              double conc = 0;
              if (speSta != null) {
                conc = speSta.getConcentration();
              }
              double uncertainty = 0;
              if (ss.getTime().getTime() != 0) {
                int I = ss.getRealID(spe);
                LinkedList reactionList = ss.getReactionList();
                for (int j = 0; j < reactionList.size(); j++) {
                  Reaction r = (Reaction) reactionList.get(j);
                  ReactionTime rt = ss.getTime();
                  Temperature t = getTemperature(rt);
                  double k;
                  double k_upperbound;
                  double k_lowerbound;
                  if (r instanceof TROEReaction){
                    k = ((TROEReaction)r).calculateRate(ss);
                    k_upperbound = r.calculateUpperBoundRate(t)*((TROEReaction)r).calculateTroeFallOff(ss);
                    k_lowerbound = r.calculateLowerBoundRate(t)*((TROEReaction)r).calculateTroeFallOff(ss);
                  }
                  else if (r instanceof LindemannReaction) {
                	  k = ((LindemannReaction)r).calculateRate(ss);
                	  k_upperbound = r.calculateUpperBoundRate(t)*((LindemannReaction)r).calculateLindemannFallOff(ss);
                	  k_lowerbound = r.calculateLowerBoundRate(t)*((LindemannReaction)r).calculateLindemannFallOff(ss);
                  }
                  else if (r instanceof ThirdBodyReaction){
                    k = ((ThirdBodyReaction)r).calculateRate(ss);
                    k_upperbound = r.calculateUpperBoundRate(t)*((ThirdBodyReaction)r).calculateThirdBodyCoefficient(ss);
                    k_lowerbound = r.calculateLowerBoundRate(t)*((ThirdBodyReaction)r).calculateThirdBodyCoefficient(ss);
                  }
                  if (r instanceof TemplateReaction){
                    k = ((TemplateReaction)r).calculateTotalPDepRate(ss.getTemperature(), ss.getPressure());//10/25/07 gmagoon: added pressure
                    k_upperbound = k*2.0;
                    k_lowerbound = k/2.0;
                  }
                  else if (r instanceof PDepReaction){
                    k = ((PDepReaction)r).calculateRate(ss.getTemperature(),ss.getPressure());
                    k_upperbound = k*2.0;
                    k_lowerbound = k/2.0;
                  }
                  else {
                    k = r.calculateTotalRate(ss.getTemperature());
                    k_upperbound = r.getUpperBoundRate(t);
                    k_lowerbound = r.getLowerBoundRate(t);
                  }
                  int index = I + j * getReactionModel().getSpeciesNumber() - 1;
                  
                    double sens = ss.getSensitivityStatus(index);
                    double delta_k;
                    if (sens >0){
                      delta_k = k_upperbound - k;
                    }
                    else {
                      delta_k = k_lowerbound-k;
                    }
                    uncertainty += Math.abs(sens * delta_k);
                  
                }
              }
              if (conc != 0) {
                double upper_bound = conc * Math.exp(uncertainty / conc);
                result += '\t' + String.valueOf(upper_bound);
              }
              else {
                double upper_bound = 0;
                result += '\t' + String.valueOf(upper_bound);
              }
            }
          }
          result += '\n';
        }
        result += "END\n";
        return result;
        //*]
      }




    //## operation resetSystemSnapshot()
    public void resetSystemSnapshot() {
        //#[ operation resetSystemSnapshot()
        systemSnapshot.clear();
        systemSnapshot.add(initialStatus);
        //#]
    }

	public void reduceModel() {
		double [] yMax;
		//yMax = getDynamicSimulator().getHighestConcentrations();
		
		
	}

    //## operation solveReactionSystem(ReactionTime,ReactionTime,boolean,boolean,boolean)
    //9/24/07 gmagoon: added p_reactionModel as parameter; subsequently removed
    public ReactionTime solveReactionSystem(ReactionTime p_beginTime, ReactionTime p_endTime, boolean p_initialization, boolean p_reactionChanged, boolean p_conditionChanged, int iterationNum) {
 
        //#[ operation solveReactionSystem(ReactionTime,ReactionTime,boolean,boolean,boolean)
        Temperature t = getTemperatureModel().getTemperature(p_beginTime);
        Pressure p = getPressureModel().getPressure(p_beginTime);

        SystemSnapshot beginStatus = (SystemSnapshot)(getSystemSnapshotEnd().next());

        if (p_reactionChanged || p_initialization || p_conditionChanged) {
             if ((reactionModelEnlarger instanceof RateBasedPDepRME)) {//1/2/09 gmagoon and rwest: only call initializePDepNetwork for P-dep cases
        	initializePDepNetwork();
             }
            p_reactionChanged = true;
            beginStatus = getInitialStatus();
        }


        if (!beginStatus.getTime().equals(p_beginTime)) throw new InvalidBeginStatusException();

        SystemSnapshot present = getDynamicSimulator().solve(p_initialization, getReactionModel(), p_reactionChanged, beginStatus, p_beginTime, p_endTime,t,p, p_conditionChanged, finishController.terminationTester, iterationNum);

        appendUnreactedSpeciesStatus(present, t);

        systemSnapshot.add(present);
        return present.time;
        //#]
    }

    public void solveReactionSystemwithSEN(ReactionTime p_beginTime, ReactionTime p_endTime, boolean p_initialization, boolean p_reactionChanged, boolean p_conditionChanged) {
    	Temperature t = getTemperatureModel().getTemperature(p_beginTime);
        Pressure p = getPressureModel().getPressure(p_beginTime);

        SystemSnapshot beginStatus = (SystemSnapshot)(getSystemSnapshotEnd().next());

        if (p_reactionChanged || p_initialization || p_conditionChanged) {
            if ((reactionModelEnlarger instanceof RateBasedPDepRME)) {//1/2/09 gmagoon and rwest: only call initializePDepNetwork for P-dep cases
        	initializePDepNetwork();
            }
            p_reactionChanged = true;
            beginStatus = getInitialStatus();
        }


        if (!beginStatus.getTime().equals(p_beginTime)) throw new InvalidBeginStatusException();

        LinkedList sS = ((JDASPK)getDynamicSimulator()).solveSEN(p_initialization, getReactionModel(), p_reactionChanged, beginStatus, p_beginTime, p_endTime,t,p, p_conditionChanged, finishController.terminationTester);

        for (int i=0; i< sS.size(); i++){
        	systemSnapshot.add(sS.get(i));
        }
        return;
		
	}
    //## operation solveReactionSystem(ReactionTime,ReactionTime,boolean,boolean,boolean)
  
    
    //## operation toString()
    public String toString() {
        //#[ operation toString()
        return reactionModel.toString();
        //#]
    }

    public LinkedHashSet getOriginalReactant() {
        return originalReactant;
    }
    
    public void setOriginalReactant(LinkedHashSet p_originalReactant){
    	originalReactant = p_originalReactant;
    }

    public DynamicSimulator getDynamicSimulator() {
        return dynamicSimulator;
    }

    public void setDynamicSimulator(DynamicSimulator p_DynamicSimulator) {
        dynamicSimulator = p_DynamicSimulator;
    }

    public FinishController getFinishController() {
        return finishController;
    }

    public void __setFinishController(FinishController p_FinishController) {
        finishController = p_FinishController;
    }

    public void _setFinishController(FinishController p_FinishController) {
        if(finishController != null)
            finishController.__setReactionSystem(null);
        __setFinishController(p_FinishController);
    }

    public void setFinishController(FinishController p_FinishController) {
        if(p_FinishController != null)
            p_FinishController._setReactionSystem(this);
        _setFinishController(p_FinishController);
    }

    public void _clearFinishController() {
        finishController = null;
    }

    public InitialStatus getInitialStatus() {
        return initialStatus;
    }

    public void setInitialStatus(InitialStatus p_InitialStatus) {
        initialStatus = p_InitialStatus;
    }

    public PressureModel getPressureModel() {
        return pressureModel;
    }

    public void setPressureModel(PressureModel p_PressureModel) {
        pressureModel = p_PressureModel;
    }


     public LibraryReactionGenerator getLibraryReactionGenerator() {
        return lrg;
     }
////10/4/07 gmagoon: used in ReactionModelGenerator.java;10/9/07: restored for use by RateBasedPDepRME
    public ReactionGenerator getReactionGenerator() {
        return reactionGenerator;
    }


    public ReactionModel getReactionModel() {
        return reactionModel;
    }
    
    //gmagoon 10/4/07: restored setReactionModel
    public void setReactionModel(ReactionModel p_ReactionModel) {
        reactionModel = p_ReactionModel;
    }

    public ReactionModelEnlarger getReactionModelEnlarger() {
        return reactionModelEnlarger;
    }

    public void setReactionModelEnlarger(ReactionModelEnlarger p_ReactionModelEnlarger) {
        reactionModelEnlarger = p_ReactionModelEnlarger;
    }

    public ListIterator getSystemSnapshot() {
        ListIterator iter=systemSnapshot.listIterator(0);
        return iter;
    }

    public ListIterator getSystemSnapshotEnd() {
        return systemSnapshot.listIterator(systemSnapshot.lastIndexOf(systemSnapshot.getLast()));
    }

    public SystemSnapshot newSystemSnapshot() {
        SystemSnapshot newSystemSnapshot = new SystemSnapshot();
        systemSnapshot.add(newSystemSnapshot);
        return newSystemSnapshot;
    }

    public void deleteSystemSnapshot(SystemSnapshot p_SystemSnapshot) {
        systemSnapshot.remove(p_SystemSnapshot);
        p_SystemSnapshot=null;
    }

    public TemperatureModel getTemperatureModel() {
        return temperatureModel;
    }

    public void setTemperatureModel(TemperatureModel p_TemperatureModel) {
        temperatureModel = p_TemperatureModel;
    }

    //10/30/07 gmagoon: added accessor method for index
    public int getIndex() {
        return ind;
    }
	

}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ReactionSystem.java
*********************************************************************/



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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.ListIterator;
import jing.chem.Species;
import jing.chem.SpeciesDictionary;
import jing.param.Global;
import jing.param.ParameterInfor;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxn.NegativeRateException;
import jing.rxn.PDepNetwork;
import jing.rxn.PDepReaction;
import jing.rxn.Reaction;
import jing.rxn.TROEReaction;
import jing.rxn.TemplateReaction;
import jing.rxn.ThirdBodyReaction;

/**
 * Common base class for the DASSL and DASPK solvers, which share a lot of
 * code.
 */
public abstract class JDAS implements DAESolver {

	protected LinkedHashMap IDTranslator = new LinkedHashMap();		//## attribute IDTranslator

    protected double atol;		//## attribute atol

    protected int parameterInfor;//svp

    protected ParameterInfor [] parameterInforArray = null;		//## attribute parameterInfor

    protected double rtol;		//## attribute rtol

    protected InitialStatus initialStatus;//svp
	protected int nState =  3 ;
    protected int neq = 3;
	protected int nParameter =0;
	protected double [] y;
	protected double [] yprime;
	protected int [] info = new int[30];
	protected LinkedList rList ;
	protected LinkedList duplicates ;
    protected LinkedList thirdBodyList ;
    protected LinkedList troeList ;
	protected StringBuilder outputString ;
	protected StringBuilder rString ;
	protected StringBuilder tbrString;
	protected StringBuilder troeString;
	protected int index; //11/1/07 gmagoon: adding index to allow appropriate naming of RWORK, IWORK****may need to make similar modification for DASPK?
	protected ValidityTester validityTester; //5/5/08 gmagoon: adding validityTester and autoflag as attributes needed for "automatic" time stepping
	protected boolean autoflag;
	protected double [] reactionFlux;
	protected double [] conversionSet;
	protected double endTime;
        protected StringBuilder thermoString = new StringBuilder();
	
	protected JDAS() {
        
    }
	
    public JDAS(double p_rtol, double p_atol, int p_parameterInfor, 
			InitialStatus p_initialStatus, int p_index, ValidityTester p_vt, 
			boolean p_autoflag) {
        
		rtol = p_rtol;
        atol = p_atol;
        index = p_index;
        validityTester = p_vt;
        autoflag = p_autoflag;

        parameterInfor = p_parameterInfor;
        initialStatus = p_initialStatus;

    }
	
	public void addSA(){
      if (parameterInfor == 0) {
        parameterInfor = 1;
      }
    }
	
	public void addConversion(double [] p_conversions, int numConversion) {
    	conversionSet = new double[numConversion];
    	for (int i=0; i< numConversion; i++)
    		conversionSet[i] = p_conversions[i];
    }

    public double[] getConversion(){
		return conversionSet;
	}
    
    public StringBuilder generatePDepODEReactionList(ReactionModel p_reactionModel, 
			SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
    
		StringBuilder rString = new StringBuilder();
        StringBuilder arrayString = new StringBuilder();
        StringBuilder rateString = new StringBuilder();
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;
        
        rList = new LinkedList();
        duplicates = new LinkedList();
    	LinkedList nonPDepList = new LinkedList();
        LinkedList pDepList = new LinkedList();

		generatePDepReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure, nonPDepList, pDepList);

		int size = nonPDepList.size() + pDepList.size() + duplicates.size();

        for (Iterator iter = nonPDepList.iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();

			if (!(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction)){
				rList.add(r);
				ODEReaction or = transferReaction(r, p_beginStatus, p_temperature, p_pressure);
				arrayString.append(or.rNum+" "+or.pNum+" ");
				for (int i=0;i<3;i++){
					if (i<or.rNum)
						arrayString.append(or.rID[i]+" ");
					else
						arrayString.append(0+" ");
				}
				for (int i=0;i<3;i++){
					if (i<or.pNum)
						arrayString.append(or.pID[i]+" ");
					else
						arrayString.append(0+" ");
				}
				// Original DASSL has these lines uncommented, while DASPK is as given (should they be different?)
				if (r.hasReverseReaction())
					arrayString.append(1 + " ");
				else
					arrayString.append(0 + " ");
				rateString.append(or.rate + " " + or.A + " " + or.n + " " + or.E + " "+r.calculateKeq(p_temperature)+ " ");
					
			}

        }

        for (Iterator iter = pDepList.iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	
			if (r instanceof PDepReaction) {
			
				rList.add(r);

				ODEReaction or = transferReaction(r, p_beginStatus, p_temperature, p_pressure);
				arrayString.append(or.rNum+" "+or.pNum+" ");
				for (int i=0;i<3;i++){
					if (i<or.rNum)
						arrayString.append(or.rID[i]+" ");
					else
						arrayString.append(0+" ");
				}
				for (int i=0;i<3;i++){
					if (i<or.pNum)
						arrayString.append(or.pID[i]+" ");
					else
						arrayString.append(0+" ");
				}

				arrayString.append(1 + " ");

				rateString.append(or.rate + " " + or.A + " " + or.n + " " + or.E + " "+r.calculateKeq(p_temperature)+" ");
			}
        }
        
		for (Iterator iter = duplicates.iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();

			//if (!(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction)){
			if (r instanceof PDepReaction) {
				rList.add(r);
				ODEReaction or = transferReaction(r, p_beginStatus, p_temperature, p_pressure);
				arrayString.append(or.rNum+" "+or.pNum+" ");
				for (int i=0;i<3;i++){
					if (i<or.rNum)
						arrayString.append(or.rID[i]+" ");
					else
						arrayString.append(0+" ");
				}
				for (int i=0;i<3;i++){
					if (i<or.pNum)
						arrayString.append(or.pID[i]+" ");
					else
						arrayString.append(0+" ");
				}
				//if (r.hasReverseReaction())
					arrayString.append(1 + " ");
				//else
					//arrayString.append(0 + " ");
				rateString.append(or.rate + " " + or.A + " " + or.n + " " + or.E + " "+r.calculateKeq(p_temperature)+ " ");
					
			}

        }
        
        rString.append(arrayString.toString()+"\n"+rateString.toString());
        return rString;
    }
	
	public void generatePDepReactionList(ReactionModel p_reactionModel, 
			SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure,
			LinkedList nonPDepList, LinkedList pDepList) {
    
		CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;
        
		for (Iterator iter = PDepNetwork.getCoreReactions(cerm).iterator(); iter.hasNext(); ) {
			PDepReaction rxn = (PDepReaction) iter.next();
			if (cerm.categorizeReaction(rxn) != 1) continue;
        		//check if this reaction is not already in the list and also check if this reaction has a reverse reaction
			// which is already present in the list.
			if (rxn.getReverseReaction() == null)
				rxn.generateReverseReaction();

			if (!rxn.reactantEqualsProduct()  && !troeList.contains(rxn) && !troeList.contains(rxn.getReverseReaction()) && !thirdBodyList.contains(rxn) && !thirdBodyList.contains(rxn.getReverseReaction()) ) {
				if (!pDepList.contains(rxn) && !pDepList.contains(rxn.getReverseReaction())){
					pDepList.add(rxn);
				}
				else if (pDepList.contains(rxn) && !pDepList.contains(rxn.getReverseReaction()))
					continue;
				else if (!pDepList.contains(rxn) && pDepList.contains(rxn.getReverseReaction())){
					Temperature T = new Temperature(298, "K");
					if (rxn.calculateKeq(T)>0.999) {
						pDepList.remove(rxn.getReverseReaction());
						pDepList.add(rxn);
					}
				}

			}
		}
		
		for (Iterator iter = p_reactionModel.getReactionSet().iterator(); iter.hasNext(); ) {
			Reaction r = (Reaction)iter.next();
			if (r.isForward() && !(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction))
				nonPDepList.add(r);
		}
			
		duplicates.clear();

	}
	
	protected LinkedHashMap generateSpeciesStatus(ReactionModel p_reactionModel, 
			double [] p_y, double [] p_yprime, int p_paraNum) {
        
		int neq = p_reactionModel.getSpeciesNumber()*(p_paraNum+1);
        if (p_y.length != neq) throw new DynamicSimulatorException();
        if (p_yprime.length != neq) throw new DynamicSimulatorException();

        LinkedHashMap speStatus = new LinkedHashMap();

        for (Iterator iter = p_reactionModel.getSpecies(); iter.hasNext(); ) {
        	Species spe = (Species)iter.next();
        	int id = getRealID(spe);
        	if (id>p_y.length) throw new UnknownReactedSpeciesException(spe.getName());
        	double conc = p_y[id-1];
        	double flux = p_yprime[id-1];

        	System.out.println(String.valueOf(spe.getID()) + '\t' + spe.getName() + '\t' + String.valueOf(conc) + '\t' + String.valueOf(flux));

        	if (conc < 0) {
        		if (Math.abs(conc) < 1.0E-19) conc = 0;
        		else throw new NegativeConcentrationException("species " + spe.getName() + " has negative conc: " + String.valueOf(conc));
        	}
        	SpeciesStatus ss = new SpeciesStatus(spe, 1, conc, flux);
        	speStatus.put(spe,ss);
        }

        return speStatus;
    }
	
	public StringBuilder generateThirdBodyReactionList(ReactionModel p_reactionModel, 
			SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
    
		int size = p_reactionModel.getReactionSet().size();
        StringBuilder arrayString = new StringBuilder();
        StringBuilder rateString = new StringBuilder();
        StringBuilder tbrString = new StringBuilder();
        Iterator iter = p_reactionModel.getReactionSet().iterator();
        thirdBodyList = new LinkedList();
        
        while (iter.hasNext()) {
        	Reaction r = (Reaction)iter.next();
            if ((r.isForward()) && (r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction)) {
            	
				ThirdBodyODEReaction or = (ThirdBodyODEReaction)transferReaction(r, p_beginStatus, p_temperature, p_pressure);
                thirdBodyList.add((ThirdBodyReaction)r);
                
            	arrayString.append(or.rNum+" "+or.pNum+" ");
    			for (int i=0;i<3;i++){
    				if (i<or.rNum)
    					arrayString.append(or.rID[i]+" ");
    				else
    					arrayString.append(0+" ");
    			}
    			for (int i=0;i<3;i++){
    				if (i<or.pNum)
    					arrayString.append(or.pID[i]+" ");
    				else
    					arrayString.append(0+" ");
    			}
    			if (r.hasReverseReaction())
    				arrayString.append(1 + " ");
    			else
    				arrayString.append(0 + " ");
    			arrayString.append(or.numCollider+" ");
    			for (int i=0; i<10; i++){
    				if (i <  or.numCollider)
    					arrayString.append(or.colliders[i] + " ");
    				else
    					arrayString.append(0 + " ");
    			}
    			rateString.append(or.rate + " " + or.A + " " + or.n + " " + or.E + " "+r.calculateKeq(p_temperature)+ " " +or.inertColliderEfficiency+" ");
    			for (int i=0; i<10; i++){
    				if (i < or.numCollider)
    					rateString.append(or.efficiency[i] + " ");
    				else
    					rateString.append(0 + " ");
    			}
                
                
             }
        }

        tbrString.append(arrayString.toString()+"\n"+rateString.toString());
        return tbrString;
    }
	
	protected StringBuilder generateTROEReactionList(ReactionModel p_reactionModel, 
			SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
		 
		int size = p_reactionModel.getReactionSet().size();
		 StringBuilder arrayString = new StringBuilder();
		 StringBuilder rateString = new StringBuilder();
		 StringBuilder troeString = new StringBuilder(); 
		 Iterator iter = p_reactionModel.getReactionSet().iterator();
		 troeList = new LinkedList();

		 while (iter.hasNext()) {
			 Reaction r = (Reaction)iter.next();

			 if (r.isForward() && r instanceof TROEReaction) {
				 TROEODEReaction or = (TROEODEReaction)transferReaction(r, p_beginStatus, p_temperature, p_pressure);
				 troeList.add((TROEReaction)r);
				 arrayString.append(or.rNum+" "+or.pNum+" ");
					for (int i=0;i<3;i++){
						if (i<or.rNum)
							arrayString.append(or.rID[i]+" ");
						else
							arrayString.append(0+" ");
					}
					for (int i=0;i<3;i++){
						if (i<or.pNum)
							arrayString.append(or.pID[i]+" ");
						else
							arrayString.append(0+" ");
					}

					if (r.hasReverseReaction())
						arrayString.append(1 + " ");
					else
						arrayString.append(0 + " ");

					arrayString.append(or.numCollider+" ");
					for (int i=0; i<10; i++){
						if (i <  or.numCollider)
							arrayString.append(or.colliders[i] + " ");
						else
							arrayString.append(0 + " ");
					}
					if (or.troe7)
						arrayString.append(0 + " ");
					else
						arrayString.append(1 + " ");
					rateString.append(or.highRate + " " + or.A + " " + or.n + " " + or.E + " "+r.calculateKeq(p_temperature)+ " "+or.inertColliderEfficiency+" ");
					for (int i=0; i<10; i++){
						if (i < or.numCollider)
							rateString.append(or.efficiency[i] + " ");
						else
							rateString.append(0 + " ");
					}
					rateString.append(or.a + " " + or.Tstar + " " + or.T2star + " " + or.T3star + " " + or.lowRate+ " ");
			 }
		 }
		 troeString.append(arrayString.toString()+"\n"+rateString.toString());
		 return troeString;

	 }
	
	public int getRealID(Species p_species) {
        Integer id = (Integer)IDTranslator.get(p_species);
        if (id == null) {
        	id = new Integer(IDTranslator.size()+1);
	        IDTranslator.put(p_species, id);
                thermoString.append(p_species.calculateG(initialStatus.getTemperature()) + " ");//10/26/07 gmagoon: changed to avoid use of Global.temperature;****ideally, current temperature would be used, but initial temperature is simplest to pass in current implementation
        }
        return id.intValue();
    }
	
	protected void initializeWorkSpace() {
		for (int i=0; i<30; i++)
			info[i] = 0;
		info[2] = 1; //print out the time steps
	}
	
	protected void initializeConcentrations(SystemSnapshot p_beginStatus, 
			ReactionModel p_reactionModel, ReactionTime p_beginTime, 
			ReactionTime p_endTime, LinkedList initialSpecies) {
		y = new double[neq];
		yprime = new double[neq];
		for (Iterator iter = p_beginStatus.getSpeciesStatus(); iter.hasNext(); ) {
			SpeciesStatus ss = (SpeciesStatus)iter.next();
			double conc = ss.getConcentration();
			double flux = ss.getFlux();
			if (ss.isReactedSpecies()) {

				Species spe = ss.getSpecies();
				int id = getRealID(spe);
				//System.out.println(String.valueOf(spe.getID()) + '\t' + spe.getName() + '\t' + String.valueOf(conc) + '\t' + String.valueOf(flux));

				y[id-1] = conc;
				yprime[id-1] = flux;
			}
		}
	}
	
	public ODEReaction transferReaction(Reaction p_reaction, SystemSnapshot p_beginStatus, 
			Temperature p_temperature, Pressure p_pressure) {
    	
    	//System.out.println(p_reaction.getStructure().toString()+"\t"+p_reaction.calculateTotalRate(Global.temperature));
        double startTime = System.currentTimeMillis();
		double dT = 1;
        Temperature Tup = new Temperature(p_temperature.getStandard()+dT, Temperature.getStandardUnit());
        Temperature Tlow = new Temperature(p_temperature.getStandard()-dT, Temperature.getStandardUnit());

        int rnum = p_reaction.getReactantNumber();
        int pnum = p_reaction.getProductNumber();

        int [] rid = new int[rnum];
        int index = 0;
        for (Iterator r_iter = p_reaction.getReactants(); r_iter.hasNext(); ) {
        	Species s = (Species)r_iter.next();
        	rid[index] = getRealID(s);
        	index++;
        }

        int [] pid = new int[pnum];
        index = 0;
        for (Iterator p_iter = p_reaction.getProducts(); p_iter.hasNext(); ) {
			Species s = (Species)p_iter.next();
        	pid[index] = getRealID(s);
        	index++;
        }

		//Global.transferReaction = Global.transferReaction + (System.currentTimeMillis() - startTime)/1000/60;
		//ODEReaction or;
        if (p_reaction instanceof PDepReaction) {
			double rate = ((PDepReaction)p_reaction).calculateRate(p_temperature, p_pressure);
				if (String.valueOf(rate).equals("NaN")){
				System.err.println(p_reaction.toChemkinString(p_temperature) + "Has bad rate probably due to Ea<DH");
				rate = 0;
			}
            ODEReaction or = new ODEReaction(rnum, pnum, rid, pid, rate);
			//Global.transferReaction = Global.transferReaction + (System.currentTimeMillis() - startTime)/1000/60;
        	return or;
        }
        else {
        	double rate = 0;
        	if (p_reaction instanceof TemplateReaction) {
				//startTime = System.currentTimeMillis();
        		//rate = ((TemplateReaction)p_reaction).getRateConstant();
				
				rate = ((TemplateReaction)p_reaction).calculateTotalRate(p_beginStatus.temperature);
				ODEReaction or = new ODEReaction(rnum, pnum, rid, pid, rate);
				//Global.transferReaction = Global.transferReaction + (System.currentTimeMillis() - startTime)/1000/60;
				
				return or;
			}
			else if (p_reaction instanceof TROEReaction){//svp
				startTime = System.currentTimeMillis();
				HashMap weightMap = ((ThirdBodyReaction)p_reaction).getWeightMap();
				int weightMapSize = weightMap.size();
				int [] colliders = new int[weightMapSize];
				double [] efficiency = new double[weightMapSize];
				Iterator colliderIter = weightMap.keySet().iterator();
				int numCollider =0;
				for (int i=0; i<weightMapSize; i++){
					String name = (String)colliderIter.next();
					Species spe = SpeciesDictionary.getInstance().getSpeciesFromName(name);
					if (spe != null){
						colliders[numCollider] = getRealID(spe);
						efficiency[numCollider] = ((Double)weightMap.get(name)).doubleValue();
						numCollider++;
					}

				}
				Global.transferReaction = Global.transferReaction + (System.currentTimeMillis() - startTime)/1000/60;
				double T2star, T3star, Tstar, a;
				T2star = ((TROEReaction)p_reaction).getT2star();
				T3star = ((TROEReaction)p_reaction).getT3star();
				Tstar = ((TROEReaction)p_reaction).getTstar();
				a = ((TROEReaction)p_reaction).getA();
				int direction = p_reaction.getDirection();
				double Keq = p_reaction.calculateKeq(p_temperature);
				double lowRate = ((TROEReaction)p_reaction).getLow().calculateRate(p_temperature, -1);
				double highRate = ((TROEReaction)p_reaction).getKinetics().calculateRate(p_temperature, -1);
				double inertColliderEfficiency = ((ThirdBodyReaction)p_reaction).calculateThirdBodyCoefficientForInerts(p_beginStatus);
				boolean troe7 = ((TROEReaction)p_reaction).getTroe7();
				TROEODEReaction or = new TROEODEReaction(rnum, pnum, rid, pid, direction, Keq, colliders, efficiency, numCollider, inertColliderEfficiency, T2star, T3star, Tstar, a, highRate, lowRate, troe7);
				return or;
			}
			else if (p_reaction instanceof ThirdBodyReaction){//svp
				startTime = System.currentTimeMillis();
				HashMap weightMap = ((ThirdBodyReaction)p_reaction).getWeightMap();
				int weightMapSize = weightMap.size();
				int [] colliders = new int[weightMapSize];
				double [] efficiency = new double[weightMapSize];
				Iterator colliderIter = weightMap.keySet().iterator();
				int numCollider =0;
				for (int i=0; i<weightMapSize; i++){
					String name = (String)colliderIter.next();
					Species spe = SpeciesDictionary.getInstance().getSpeciesFromName(name);
					if (spe != null){
						colliders[numCollider] = getRealID(spe);
						efficiency[numCollider] = ((Double)weightMap.get(name)).doubleValue();
						numCollider++;
					}

				}
				Global.transferReaction = Global.transferReaction + (System.currentTimeMillis() - startTime)/1000/60;
				
				rate = p_reaction.calculateTotalRate(p_beginStatus.temperature);
				
				double inertColliderEfficiency = ((ThirdBodyReaction)p_reaction).calculateThirdBodyCoefficientForInerts(p_beginStatus);
				//rate = p_reaction.getRateConstant();
				ThirdBodyODEReaction or = new ThirdBodyODEReaction(rnum, pnum, rid, pid, rate, colliders, efficiency,numCollider, inertColliderEfficiency);
				return or;
			}

			else{
				rate = p_reaction.calculateTotalRate(p_beginStatus.temperature);
				//startTime = System.currentTimeMillis();
				//rate = p_reaction.getRateConstant();
				ODEReaction or = new ODEReaction(rnum, pnum, rid, pid, rate);
				//Global.transferReaction = Global.transferReaction + (System.currentTimeMillis() - startTime)/1000/60;
				
				return or;
			}


        }
    
    }
	
	public double getAtol() {
        return atol;
    }

    public int getReactionSize(){
    	return rList.size()+troeList.size()+thirdBodyList.size();
    }
    public int getMaxSpeciesNumber() {
        return IDTranslator.size() - 1;
    }

    public double getRtol() {
        return rtol;
    }

	public String getEdgeReactionString(CoreEdgeReactionModel model, HashMap edgeID,
			Reaction r, Temperature temperature, Pressure pressure) {
		
		int edgeSpeciesCounter = edgeID.size();
		
		// Find the rate coefficient
		double k;
		if (r instanceof TemplateReaction)
			k = ((TemplateReaction) r).getRateConstant(temperature, pressure);
		else if (r instanceof PDepReaction)
			k = ((PDepReaction) r).calculateRate(temperature, pressure);
		else
			k = r.getRateConstant(temperature);
		
		if (k > 0) {
			int reacCount = 0;
			int prodCount = 0;
			int[] tempReacArray = {0, 0, 0};
			int[] tempProdArray = {0, 0, 0};
			//iterate over the reactants, counting and storing IDs in tempReacArray, up to a maximum of 3 reactants
			for (Iterator rIter = r.getReactants(); rIter.hasNext();) {
				reacCount++;
				Species spe = (Species) rIter.next();
				tempReacArray[reacCount - 1] = getRealID(spe);
			}
			//iterate over the products, selecting products which are not already in the core, counting and storing ID's (created sequentially in a HashMap, similar to getRealID) in tempProdArray, up to a maximum of 3
			for (Iterator pIter = r.getProducts(); pIter.hasNext();) {
				Species spe = (Species) pIter.next();
				if (model.containsAsUnreactedSpecies(spe)) {
					prodCount++;
					Integer id = (Integer) edgeID.get(spe);
					if (id == null) {
						edgeSpeciesCounter++;
						id = new Integer(edgeSpeciesCounter);
						edgeID.put(spe, id);
					}
					tempProdArray[prodCount - 1] = id;
				}
			}
			//update the output string with info for one reaction
			String str = reacCount + " " + prodCount + " " + tempReacArray[0] + " " + tempReacArray[1] + " " + tempReacArray[2] + " " + tempProdArray[0] + " " + tempProdArray[1] + " " + tempProdArray[2] + " " + k;
			return str;
		} else {
			throw new NegativeRateException(r.toChemkinString(temperature) + ": " + String.valueOf(k));
		}
	}
	
	public void getAutoEdgeReactionInfo(CoreEdgeReactionModel model, Temperature p_temperature,
			Pressure p_pressure) {
		
                //IMPORTANT: this code should pass the information needed to perform the same checks as done by the validity testing in the Java code
		//much of code below is taken or based off of code from appendUnreactedSpeciesStatus in ReactionSystem.java
		StringBuilder edgeReacInfoString = new StringBuilder();
		int edgeReactionCounter = 0;
		int edgeSpeciesCounter = 0;

		// First use reactions in unreacted reaction set, which is valid for both RateBasedRME and RateBasedPDepRME
		HashMap edgeID = new HashMap();
		LinkedHashSet ur = model.getUnreactedReactionSet();
		for (Iterator iur = ur.iterator(); iur.hasNext();) {
			edgeReactionCounter++;
			Reaction r = (Reaction) iur.next();
			String str = getEdgeReactionString(model, edgeID, r, p_temperature, p_pressure);
			edgeReacInfoString.append("\n" + str);
		}

		// For the case where validityTester is RateBasedPDepVT (assumed to also be directly associated with use of RateBasedPDepRME), consider two additional types of reactions
		if (validityTester instanceof RateBasedPDepVT) {
                        //first consider NetReactions (formerly known as PDepNetReactionList)
			for (Iterator iter0 = PDepNetwork.getNetworks().iterator(); iter0.hasNext();) {
				PDepNetwork pdn = (PDepNetwork) iter0.next();
				for (ListIterator iter = pdn.getNetReactions().listIterator(); iter.hasNext(); ) {
					PDepReaction rxn = (PDepReaction) iter.next();
					if (rxn.isEdgeReaction(model)) {
						edgeReactionCounter++;
						String str = getEdgeReactionString(model, edgeID, rxn, p_temperature, p_pressure);
						edgeReacInfoString.append("\n" + str);
					}
				}
			}
                        //second, consider kLeak of each reaction network so that the validity of each reaction network may be tested
                        //in the original CHEMDIS approach, we included a reaction and pseudospecies for each kleak/P-dep network
                        //with the FAME approach we still consider each P-dep network as a pseudospecies, but we have multiple reactions contributing to this pseudo-species, with each reaction having different reactants
                        for (Iterator iter1 = PDepNetwork.getNetworks().iterator(); iter1.hasNext();) {
                                PDepNetwork pdn = (PDepNetwork)iter1.next();
                                //account for pseudo-edge species product by incrementing the edgeSpeciesCounter and storing the ID in the tempProdArray; each of these ID's will occur once and only once; thus, note that the corresponding PDepNetwork is NOT stored to the HashMap
                                int prodCount=1;//prodCount will not be modified as each PDepNetwork will be treated as a pseudo-edge species product
                                int[] tempProdArray = {0, 0, 0};
                                edgeSpeciesCounter++;
                                tempProdArray[0]=edgeSpeciesCounter;//note that if there are no non-included reactions that have all core reactants for a particular P-dep network, then the ID will be allocated, but not used...hopefully this would not cause problems with the Fortran code
                                double k = 0.0;
                             //   double k = pdn.getKLeak(index);//index in DASSL should correspond to the same (reactionSystem/TPcondition) index as used by kLeak
                             //   if (!pdn.isActive() && pdn.getIsChemAct()) {
                             //           k = pdn.getEntryReaction().calculateTotalRate(p_temperature);
                             //   }
                                for (ListIterator<PDepReaction> iter = pdn.getNonincludedReactions().listIterator(); iter.hasNext(); ) {//cf. getLeakFlux in PDepNetwork
                                    PDepReaction rxn = iter.next();
                                    int reacCount=0;
                                    int[] tempReacArray = {0, 0, 0};
                                    boolean allCoreReac=false;//allCoreReac will be used to track whether all reactant species are in the core
                                    if (rxn.getReactant().getIncluded() && !rxn.getProduct().getIncluded()){
                                        k = rxn.calculateRate(p_temperature, p_pressure);
         
                                        //iterate over the reactants, counting and storing IDs in tempReacArray, up to a maximum of 3 reactants
                                        for (ListIterator<Species> rIter = rxn.getReactant().getSpeciesListIterator(); rIter.hasNext(); ) {
                                            allCoreReac=true;
                                            reacCount++;
                                            Species spe = (Species)rIter.next();
                                            if(model.containsAsReactedSpecies(spe)){
                                                tempReacArray[reacCount-1]=getRealID(spe);
                                            }
                                            else{
                                                allCoreReac=false;
                                            }
                                        }
                                    }
                                    else if (!rxn.getReactant().getIncluded() && rxn.getProduct().getIncluded()){
                                        PDepReaction rxnReverse = (PDepReaction)rxn.getReverseReaction();
                                        k = rxnReverse.calculateRate(p_temperature, p_pressure);
                                        //iterate over the products, counting and storing IDs in tempReacArray, up to a maximum of 3 reactants
                                        for (ListIterator<Species> rIter = rxn.getProduct().getSpeciesListIterator(); rIter.hasNext(); ) {
                                            allCoreReac=true;
                                            reacCount++;
                                            Species spe = (Species)rIter.next();
                                            if(model.containsAsReactedSpecies(spe)){
                                                tempReacArray[reacCount-1]=getRealID(spe);
                                            }
                                            else{
                                                allCoreReac=false;
                                            }
                                        }
                                    }
                                    if(allCoreReac){//only consider cases where all reactants are in the core
                                        edgeReactionCounter++;
                                        //update the output string with info for kLeak for one PDepNetwork
                                        edgeReacInfoString.append("\n" + reacCount + " " + prodCount + " " + tempReacArray[0] + " " + tempReacArray[1] + " " + tempReacArray[2] + " " + tempProdArray[0] + " " + tempProdArray[1] + " " + tempProdArray[2] + " " + k);
                                    }
                                    
                                    
                                }

                        }
		}

		//edgeSpeciesCounter = edgeID.size();
		outputString.append("\n" + ((RateBasedVT)validityTester).getTolerance() + "\n" + edgeSpeciesCounter + " " + edgeReactionCounter);//***thresh needs to be defined
		outputString.append(edgeReacInfoString);
	}
	
	public void getConcentractionFlags(ReactionModel p_reactionModel) {
		
		// Add list of flags for constantConcentration
        // one for each species, and a final one for the volume
        // if 1: DASSL will not change the number of moles of that species (or the volume)
        // if 0: DASSL will integrate the ODE as normal
        // eg. liquid phase calculations with a constant concentration of O2 (the solubility limit - replenished from the gas phase)
        // for normal use, this will be a sequence of '0 's
        outputString.append("\n");
		for (Iterator iter = p_reactionModel.getSpecies(); iter.hasNext(); ) {
        	Species spe = (Species)iter.next();
            if (spe.isConstantConcentration())
                outputString.append("1 ");
            else 
                outputString.append("0 ");
        }
        outputString.append("0 \n"); // for liquid EOS or constant volume this should be 1 
		
	}
	
	protected void solveDAE(String execName) {
		
		String workingDirectory = System.getProperty("RMG.workingDirectory");
		
		// write the input file
		File SolverInput = new File("ODESolver/SolverInput.dat");
		try {
			FileWriter fw = new FileWriter(SolverInput);
			fw.write(outputString.toString());
			fw.close();
		} catch (IOException e) {
			System.err.println("Problem writing Solver Input File!");
			e.printStackTrace();
		}
		
		// Rename RWORK and IWORK files if they exist
		File f = new File("ODESolver/RWORK_"+index+".dat");
		File newFile = new File("ODESolver/RWORK.dat");
		if(f.exists()){
			if(newFile.exists())
				newFile.delete();
			f.renameTo(newFile);
			f = new File("ODESolver/IWORK_"+index+".dat");
			newFile = new File("ODESolver/IWORK.dat");
			if(newFile.exists())
				newFile.delete();
			f.renameTo(newFile);
		}
		
		//run the solver on the input file
		boolean error = false;
        try {
        	
        	String[] command = {workingDirectory +  "/software/ODESolver/" + execName};//5/5/08 gmagoon: changed to call dasslAUTO.exe
			File runningDir = new File("ODESolver");
			
			Process solver = Runtime.getRuntime().exec(command, null, runningDir);
			InputStream is = solver.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				line = line.trim();
				if (!(line.contains("ODESOLVER SUCCESSFUL"))) {
					System.err.println("Error running the ODESolver: "+line);
					error = true;
				}          
			}
        	int exitValue = solver.waitFor();
        }
        catch (Exception e) {
        	String err = "Error in running ODESolver \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        
        //11/1/07 gmagoon: renaming RWORK and IWORK files
        renameIntermediateFiles();
	}
	
	private void renameIntermediateFiles() {
		File f = new File("ODESolver/RWORK.dat");
        File newFile = new File("ODESolver/RWORK_"+index+".dat");
        if(newFile.exists())
            newFile.delete();
        f.renameTo(newFile);
        f = new File("ODESolver/IWORK.dat");
        newFile = new File("ODESolver/IWORK_"+index+".dat");
        if(newFile.exists())
            newFile.delete();
        f.renameTo(newFile);
	}
	
}

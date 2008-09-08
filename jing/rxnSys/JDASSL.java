
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\JDASSL.java
*********************************************************************/

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


import java.awt.image.IndexColorModel;
import jing.rxn.*;
import jing.chem.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

import jing.chem.Species;
import jing.rxn.PDepNetReaction;
import jing.rxn.PDepNetwork;
import jing.rxn.Reaction;
import jing.rxn.Structure;
import jing.rxn.TROEReaction;
import jing.rxn.ThirdBodyReaction;
import jing.param.Global;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.param.ParameterInfor;

//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\JDASPK.java
//----------------------------------------------------------------------------

//## class JDASPK
public class JDASSL implements ODESolver{

    protected LinkedHashMap IDTranslator = new LinkedHashMap();		//## attribute IDTranslator

    protected double atol;		//## attribute atol

    protected int maxSpeciesNumber = 0;		//## attribute maxSpeciesNumber

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
	
    //## operation JDASPK()
    private  JDASSL() {
        //#[ operation JDASPK()
        //#]
    }
    //11/1/07 gmagoon: added p_index parameter
    //5/5/08 gmagoon: adding validityTester and autoflag as attributes needed for "automatic" time stepping
    //## operation JDASPK(double,double,int, InitialStatus)
    public  JDASSL(double p_rtol, double p_atol, int p_parameterInfor, InitialStatus p_initialStatus, int p_index, ValidityTester p_vt, boolean p_autoflag) {
        //#[ operation JDASPK(double,double,int, InitialStatus)
        rtol = p_rtol;
        atol = p_atol;
        index = p_index;
        validityTester = p_vt;
        autoflag=p_autoflag;

        parameterInfor = p_parameterInfor;
        initialStatus = p_initialStatus;


        //#]
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
    
    //## operation generatePDepODEReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
    public StringBuilder generatePDepODEReactionList(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
    	  //#[ operation generatePDepODEReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
        StringBuilder rString = new StringBuilder();
        StringBuilder arrayString = new StringBuilder();
        StringBuilder rateString = new StringBuilder();
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;
        
        rList = new LinkedList();
        duplicates = new LinkedList();
    	LinkedList nonPDepList = new LinkedList();
        LinkedList pDepList = new LinkedList();

        LinkedHashSet pDepStructureSet = new LinkedHashSet();
        for (Iterator iter = PDepNetwork.getDictionary().values().iterator(); iter.hasNext(); ) {
        	PDepNetwork pdn = (PDepNetwork)iter.next();
        	for (Iterator pdniter = pdn.getPDepNetReactionList(); pdniter.hasNext();) {
        		PDepNetReaction pdnr = (PDepNetReaction)pdniter.next();
        		if (cerm.categorizeReaction(pdnr) != 1) continue;
        		
        		//check if this reaction is not already in the list and also check if this reaction has a reverse reaction
        		// which is already present in the list.
        		if (pdnr.getReverseReaction() == null)
        			pdnr.generateReverseReaction();
        		
        		if (!pdnr.reactantEqualsProduct()  && !troeList.contains(pdnr) && !troeList.contains(pdnr.getReverseReaction()) && !thirdBodyList.contains(pdnr) && !thirdBodyList.contains(pdnr.getReverseReaction()) ) {
        			if (!pDepList.contains(pdnr) && !pDepList.contains(pdnr.getReverseReaction())){
        				pDepList.add(pdnr);
        			}
        			else if (pDepList.contains(pdnr) && !pDepList.contains(pdnr.getReverseReaction()))
        				continue;
        			else if (!pDepList.contains(pdnr) && pDepList.contains(pdnr.getReverseReaction())){
        				Temperature T = new Temperature(298, "K");
        				if (pdnr.calculateKeq(T)>0.999) {
        	      			pDepList.remove(pdnr.getReverseReaction());
        	      			pDepList.add(pdnr);
        	      		}
        			}
        				
        		}
        	}
        }
        
        /*for (Iterator iter = p_reactionModel.getReactionSet().iterator(); iter.hasNext(); ) {
    	Reaction r = (Reaction)iter.next();
    	if (r.isForward() && !(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction)) {
    		
    		nonPDepList.add(r);
    	}
    }*/
    
    LinkedList removeReactions = new LinkedList();
    for (Iterator iter = p_reactionModel.getReactionSet().iterator(); iter.hasNext(); ) {
    	Reaction r = (Reaction)iter.next();
    	
    	boolean presentInPDep = false;
    	if (r.isForward() && !(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction)) {
    		Iterator r_iter = pDepList.iterator();
    		while (r_iter.hasNext()){
    			Reaction pDepr = (Reaction)r_iter.next();
    			if (pDepr.equals(r)){
    				removeReactions.add(pDepr);
    				duplicates.add(pDepr);
    				if (!r.hasAdditionalKinetics()){
    					duplicates.add(r);
    					presentInPDep = true;
    				}
    			}
    		}
    		if (!presentInPDep)
    			nonPDepList.add(r);
    	}
    }
    for (Iterator iter = removeReactions.iterator(); iter.hasNext();){
  	  Reaction r = (Reaction)iter.next();
  	  pDepList.remove(r);
    }
    
    
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
				if (r.hasReverseReaction())
					arrayString.append(1 + " ");
				else
					arrayString.append(0 + " ");
				rateString.append(or.rate + " " + or.A + " " + or.n + " " + or.E + " "+r.calculateKeq(p_temperature)+ " ");
					
			}

        }

       
        
        for (Iterator iter = pDepList.iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
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
        
        for (Iterator iter = duplicates.iterator(); iter.hasNext(); ) {
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
				//if (r.hasReverseReaction())
					arrayString.append(1 + " ");
				//else
					//arrayString.append(0 + " ");
				rateString.append(or.rate + " " + or.A + " " + or.n + " " + or.E + " "+r.calculateKeq(p_temperature)+ " ");
					
			}

        }
        
        rString.append(arrayString.toString()+"\n"+rateString.toString());
        return rString;
        //#]
    }

    //## operation generateSpeciesStatus(ReactionModel,double [],double [],int)
    private LinkedHashMap generateSpeciesStatus(ReactionModel p_reactionModel, double [] p_y, double [] p_yprime, int p_paraNum) {
        //#[ operation generateSpeciesStatus(ReactionModel,double [],double [],int)
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
        //#]
    }

    //## operation generateSensitivityStatus(ReactionModel,double [],double [],int)
    private LinkedList generateSensitivityStatus(ReactionModel p_reactionModel, double [] p_y, double [] p_yprime, int p_paraNum) {
    	//#[ operation generateSensitivityStatus(ReactionModel,double [],double [],int)
    	int neq = p_reactionModel.getSpeciesNumber()*(p_paraNum+1);
    	if (p_y.length != neq) throw new DynamicSimulatorException();
    	if (p_yprime.length != neq) throw new DynamicSimulatorException();

    	LinkedList senStatus = new LinkedList();

    	if (p_paraNum > 0){
    		for (int i = p_reactionModel.getSpeciesNumber();i<neq;i++){
    			double sens = p_y[i];
    			double sflux = p_yprime[i];
    			int reaction_num = i/p_reactionModel.getSpeciesNumber();
    			int species_num = (i+1)%p_reactionModel.getSpeciesNumber();
    			if (species_num == 0){
    				species_num = p_reactionModel.getSpeciesNumber();
    			}
    			//String name = "dC"+String.valueOf(species_num)+"/dk"+String.valueOf(reaction_num);
    			//System.out.println(name + '\t' + String.valueOf(sens) + '\t' + String.valueOf(sflux));
    				SensitivityStatus ss = new SensitivityStatus(sens, sflux, species_num, reaction_num);
    				int index = i-p_reactionModel.getSpeciesNumber();
    				senStatus.add(ss);
    		}
    	}
    	return senStatus;
    	//#]
    }


    //## operation generateThirdBodyReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
    public StringBuilder generateThirdBodyReactionList(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
        //#[ operation generateThirdBodyReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
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
        //#]
    }

	 private StringBuilder generateTROEReactionList(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
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
	 
	 
    //## operation getRealID(Species)
    public int getRealID(Species p_species) {
        //#[ operation getRealID(Species)
        Integer id = (Integer)IDTranslator.get(p_species);
        if (id == null) {
        	maxSpeciesNumber++;
        	id = new Integer(maxSpeciesNumber);
        	IDTranslator.put(p_species, id);
        }

        return id.intValue();
        //#]
    }

    //## operation solve(boolean,ReactionModel,boolean,SystemSnapshot,ReactionTime,ReactionTime,Temperature,Pressure,boolean)
    public SystemSnapshot solve(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, ReactionTime p_beginTime, ReactionTime p_endTime, Temperature p_temperature, Pressure p_pressure, boolean p_conditionChanged, TerminationTester tt, int p_iterationNum) {
        //5/13/08 gmagoon: adding timing for entire solve step
     //   double JDASSL_timer = System.currentTimeMillis();
        
    	//first generate an id for all the species
    	Iterator spe_iter = p_reactionModel.getSpecies();
    	while (spe_iter.hasNext()){
    		Species spe = (Species)spe_iter.next();
    		int id = getRealID(spe);
    	}
    	double startTime = System.currentTimeMillis();
		
		ReactionTime rt = p_beginStatus.getTime();
        if (!rt.equals(p_beginTime)) throw new InvalidBeginStatusException();

        double tBegin = p_beginTime.getStandardTime();
        double tEnd = p_endTime.getStandardTime();

		double T = p_temperature.getK();
		double P = p_pressure.getAtm();
		
		LinkedList initialSpecies = new LinkedList();
		outputString = new StringBuilder();
        // set reaction set
        if (p_initialization || p_reactionChanged || p_conditionChanged) {
			
        	nState = p_reactionModel.getSpeciesNumber();
			nParameter = 0;
			if (parameterInfor != 0) {
				nParameter = p_reactionModel.getReactionNumber(); //svp
				if (initialStatus == null) System.out.println("initialStatus = null");
				spe_iter = initialStatus.getSpeciesStatus();
				while (spe_iter.hasNext()) {
					SpeciesStatus ss = (SpeciesStatus) spe_iter.next();
					String name = ss.getSpecies().getName();
					initialSpecies.add(name);
					nParameter++;
				}
			}
			neq = nState*(nParameter + 1);
			
//			tbrString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10  (16 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, rev(=1 or -1), ncollider, c1, c2,..c10 (20 elements)
			tbrString = generateThirdBodyReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			//troeString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10, alpha, Tstar, T2star, T3star, lowRate  (21 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, rev(=1 or -1), ncollider, c1, c2,..c10, troe(0=T or 1=F) (21 elements)
        	troeString = generateTROEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			//rString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, rev(=1 or -1)
			rString = generatePDepODEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			
			initializeWorkSpace();
			initializeConcentrations(p_beginStatus, p_reactionModel, p_beginTime, p_endTime, initialSpecies);
			
			
        }
        else
        	info[0] = 1;
        //5/5/08 gmagoon: (next two lines) for autoflag, get binary 0 or 1 corresponding to boolean true/false
        int af = 0;
        if (autoflag) af = 1;
        if (tt instanceof ConversionTT){
                SpeciesConversion sc = (SpeciesConversion)((ConversionTT)tt).speciesGoalConversionSet.get(0);    
                outputString.append(nState + "\t" + neq + "\t" +  getRealID(sc.species) + "\t" +conversionSet[p_iterationNum]+"\t"+ af+"\n");//5/5/08 gmagoon: added autoflag, needed when using dasslAUTO.exe
    		
        }
        else{
        	outputString.append(nState + "\t" + neq + "\t" +  -1 + "\t" +0+ "\t"+ af+"\n");//5/5/08 gmagoon: added autoflag, needed when using dasslAUTO.exe
    		
        }
        for (int i=0; i<neq; i++)
			outputString.append(y[i]+" ");
		outputString.append("\n");
		for (int i=0; i<neq; i++)
			outputString.append(yprime[i]+" ");
		outputString.append("\n" + tBegin+" "+tEnd+"\n");
		for (int i=0; i<30; i++)
			outputString.append(info[i]+" ");
		outputString.append("\n"+ rtol + " "+atol);
		outputString.append("\n" + p_temperature.getK() + " " + p_pressure.getPa() + "\n" + rList.size() + "\n" + rString.toString() + "\n" + thirdBodyList.size() + "\n"+tbrString.toString() + "\n" + troeList.size() + "\n" + troeString.toString()+"\n");
	
        //4/30/08 gmagoon: code for providing edge reaction info to DASSL in cases if the automatic time stepping flag is set to true
        if(autoflag){
            StringBuilder edgeReacInfoString = new StringBuilder();
            int edgeReactionCounter=0;
            int edgeSpeciesCounter=0;
            
            //much of code below is taken or based off of code from appendUnreactedSpeciesStatus in ReactionSystem.java
            CoreEdgeReactionModel model = (CoreEdgeReactionModel)p_reactionModel;
            //we assume the reaction model enlarger is either RateBasedRME or RateBasedPdepRME
            //first use reactions in unreacted reaction set, which is valid for both RateBasedRME and RateBasedPdepRME
            HashMap edgeID = new HashMap();
            LinkedHashSet ur = model.getUnreactedReactionSet();
            for (Iterator iur = ur.iterator(); iur.hasNext();) {
                edgeReactionCounter++;
            	Reaction r = (Reaction)iur.next();
                //find k, the rate coefficient
            	double k;
            	if (r instanceof TemplateReaction){
                        k = ((TemplateReaction)r).getRateConstant(p_temperature, p_pressure);
            	}
            	else {
                        k = r.getRateConstant(p_temperature);
            	}
            	if (k > 0) {
                        int reacCount=0;
                        int prodCount=0;
                        int[] tempReacArray = {0, 0, 0};
                        int[] tempProdArray = {0, 0, 0};
                        //iterate over the reactants, counting and storing IDs in tempReacArray, up to a maximum of 3 reactants
            		for (Iterator rIter=r.getReactants(); rIter.hasNext();) {
                                reacCount++;
                                Species spe = (Species)rIter.next();
                                tempReacArray[reacCount-1]=getRealID(spe);
            		}
                        //iterate over the products, selecting products which are not already in the core, counting and storing ID's (created sequentially in a HashMap, similar to getRealID) in tempProdArray, up to a maximum of 3
            		for (Iterator pIter=r.getProducts(); pIter.hasNext();) {
                                Species spe = (Species)pIter.next();
            			if (model.containsAsUnreactedSpecies(spe)) {
                                    prodCount++;
                                    Integer id = (Integer)edgeID.get(spe);
                                    if (id == null) {
                                        edgeSpeciesCounter++;
                                        id = new Integer(edgeSpeciesCounter);
                                        edgeID.put(spe, id);
                                    }
                                    tempProdArray[prodCount-1]=id;
            			}
            		}
                        //update the output string with info for one reaction
                        edgeReacInfoString.append("\n"+ reacCount + " " + prodCount + " " + tempReacArray[0] + " " + tempReacArray[1] + " " + tempReacArray[2] + " " + tempProdArray[0] + " " + tempProdArray[1] + " " + tempProdArray[2] + " " + k);
            	}
            	else {
            		throw new NegativeRateException(r.toChemkinString(p_temperature) + ": " + String.valueOf(k));
            	}
            }
            //for the case where validityTester is RateBasedPDepVT (assumed to also be directly associated with use of RateBasedPDepRME), consider two additional types of reactions
            if (validityTester instanceof RateBasedPDepVT){
                //first consider PDepNetReactionList
        	// populate the reactionModel with all the unreacted species if they are already not there.***? does this refer to categorize reaction?
        	for (Iterator iter=PDepNetwork.getDictionary().values().iterator(); iter.hasNext();){
        		PDepNetwork pdn = (PDepNetwork)iter.next();
        		Iterator reaction_iter = pdn.getPDepNetReactionList();
        		while (reaction_iter.hasNext()){
        			PDepNetReaction r = (PDepNetReaction)reaction_iter.next();
        			int rxnType = model.categorizeReaction(r);
        			if (rxnType == -1){
                                        edgeReactionCounter++;
                                        //4/30/08 gmagoon: temporary SystemSnapshot is created to pass temperature and pressure to calculateRate
                                        //ideally, calculateRate in PDepNetReaction.java would be changed to take temperature and pressure, since these are the only variables it uses (I think)
                                        SystemSnapshot dummySS=new SystemSnapshot();
                                        dummySS.setTemperature(p_temperature);
                                        dummySS.setPressure(p_pressure);
        				double k = r.calculateRate(dummySS);
                                        int reacCount=0;
                                        int prodCount=0;
                                        int[] tempReacArray = {0, 0, 0};
                                        int[] tempProdArray = {0, 0, 0};
                                        //iterate over the reactants, counting and storing IDs in tempReacArray, up to a maximum of 3 reactants
        				for (Iterator rIter=r.getReactants(); rIter.hasNext();) {
        					reacCount++;
                                                Species spe = (Species)rIter.next();
                                                tempReacArray[reacCount-1]=getRealID(spe);
                                        }
                                        //iterate over the products counting and storing ID's (created sequentially in a HashMap, similar to getRealID(except using HashMap rather than LinkedHashMap)) in tempProdArray, up to a maximum of 3
        				for (Iterator pIter = r.getProducts(); pIter.hasNext();){
        					Species spe = (Species)pIter.next();
                                                if (model.containsAsUnreactedSpecies(spe)) {
                                                    prodCount++;
                                                    Integer id = (Integer)edgeID.get(spe);
                                                    if (id == null) {
                                                         edgeSpeciesCounter++;
                                                         id = new Integer(edgeSpeciesCounter);
                                                         edgeID.put(spe, id);
                                                    }
                                                    tempProdArray[prodCount-1]=id;
                                                }
                                        }
                                        //update the output string with info for one reaction
                                        edgeReacInfoString.append("\n"+ reacCount + " " + prodCount + " " + tempReacArray[0] + " " + tempReacArray[1] + " " + tempReacArray[2] + " " + tempProdArray[0] + " " + tempProdArray[1] + " " + tempProdArray[2] + " " + k);
        			}
        		}
        	}
                //second, consider kLeak of each reaction network so that the validity of each reaction network may be tested
                for (Iterator iter = PDepNetwork.getDictionary().values().iterator(); iter.hasNext();) {
                        PDepNetwork pdn = (PDepNetwork)iter.next();
                        double k = pdn.getKLeak(index);//index in DASSL should correspond to the same (reactionSystem/TPcondition) index as used by kLeak
                        if (!pdn.isActive() && pdn.getIsChemAct()) {
                                k = pdn.getEntryReaction().calculateTotalRate(p_temperature);
                        }
                        int reacCount=0;
                        int prodCount=1;//prodCount will not be modified as each PDepNetwork will be treated as a pseudo-edge species product
                        boolean allCoreReac=true; //allCoreReac will be used to track whether all reactant species are in the core
                        int[] tempReacArray = {0, 0, 0};
                        int[] tempProdArray = {0, 0, 0};
                        //iterate over the reactants, counting and storing IDs in tempReacArray, up to a maximum of 3 reactants
                        for (Iterator rIter = pdn.getReactant().iterator(); rIter.hasNext(); ) {
                            reacCount++;
                            Species spe = (Species)rIter.next();
                            if(model.containsAsReactedSpecies(spe)){
                                tempReacArray[reacCount-1]=getRealID(spe);
                            }
                            else{
                                allCoreReac=false;
                            }
                        }
                        if(allCoreReac){//only consider cases where all reactants are in the core
                            edgeReactionCounter++;
                            //account for pseudo-edge species product by incrementing the edgeSpeciesCounter and storing the ID in the tempProdArray; each of these ID's will occur once and only once; thus, note that the corresponding PDepNetwork is NOT stored to the HashMap
                            edgeSpeciesCounter++;
                            tempProdArray[0]=edgeSpeciesCounter;
                            //update the output string with info for kLeak for one PDepNetwork
                            edgeReacInfoString.append("\n" + reacCount + " " + prodCount + " " + tempReacArray[0] + " " + tempReacArray[1] + " " + tempReacArray[2] + " " + tempProdArray[0] + " " + tempProdArray[1] + " " + tempProdArray[2] + " " + k);
                        }
                }
            }
            outputString.append("\n" + ((RateBasedVT)validityTester).getTolerance() + "\n" + edgeSpeciesCounter + " " + edgeReactionCounter);//***thresh needs to be defined
            outputString.append(edgeReacInfoString);
        }
        //end of code for autoflag=true
        int idid=0;
        LinkedHashMap speStatus = new LinkedHashMap();
        LinkedList senStatus = new LinkedList();
        
		int temp = 1;
        Global.solverPrepossesor = Global.solverPrepossesor + (System.currentTimeMillis() - startTime)/1000/60;
		if (nParameter==0) {
			startTime = System.currentTimeMillis();
        	
			idid = solveDAE();
			
			if (idid !=1 && idid != 2 && idid != 3)	{
				System.out.println("The idid from DASPK was "+idid );
				throw new DynamicSimulatorException("DASPK: SA off.");
        	}
            System.out.println("After ODE: from " + String.valueOf(tBegin) + " SEC to " + String.valueOf(endTime) + "SEC");
			Global.solvertime = Global.solvertime + (System.currentTimeMillis() - startTime)/1000/60;
			startTime = System.currentTimeMillis();
        	speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, 0);
			Global.speciesStatusGenerator = Global.speciesStatusGenerator + (System.currentTimeMillis() - startTime)/1000/60;
        }
        

		SystemSnapshot sss = new SystemSnapshot(new ReactionTime(endTime,"sec"), speStatus, p_beginStatus.getTemperature(), p_beginStatus.getPressure());
		LinkedList reactionList = new LinkedList();
		reactionList.addAll(rList);
        reactionList.addAll(duplicates);
        reactionList.addAll(thirdBodyList);
        reactionList.addAll(troeList);
        sss.setReactionList(reactionList);
		sss.setReactionFlux(reactionFlux);
                
        //5/13/08 gmagoon: timing for entire JDASSL solve routine
     //   double JDASSL_timerS = (System.currentTimeMillis() - JDASSL_timer)/1000;
     //   Global.JDASSLtime=Global.JDASSLtime+JDASSL_timerS;
     //   System.out.println("JDASSL solve timing: " + JDASSL_timerS);
     //   System.out.println("Cumulative JDASSL solve timing: " + Global.JDASSLtime);
        return sss;
        //#]
    }


    
	private int solveDAE() {
		
		String workingDirectory = System.getProperty("RMG.workingDirectory");
		
		//write the input file
		File SolverInput = new File("ODESolver/SolverInput.dat");
		try {
			FileWriter fw = new FileWriter(SolverInput);
			fw.write(outputString.toString());
			fw.close();
                        //11/1/07 gmagoon: renaming RWORK and IWORK files if they exist
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
		} catch (IOException e) {
			System.err.println("Problem writing Solver Input File!");
			e.printStackTrace();
		}
		
		//run the solver on the input file
		boolean error = false;
        try {
        	
        	String[] command = {workingDirectory +  "/software/ODESolver/dasslAUTO.exe"};//5/5/08 gmagoon: changed to call dasslAUTO.exe
			File runningDir = new File("ODESolver");
			
			Process ODESolver = Runtime.getRuntime().exec(command, null, runningDir);
			InputStream is = ODESolver.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				
				line = line.trim();
				if (!(line.contains("ODESOLVER SUCCESSFUL"))) {
					System.err.println("Error running the ODESolver: "+line);
                                       //6/25/08 gmagoon: commenting out timing code
                                       //5/13/08 gmagoon: modified to read in timing; before my modifications, only remaining content in if block was error=true;
                                       // if(error){
                                       //     //read in Fortran timing
                                       //    Global.fortrantime=Global.fortrantime + Double.parseDouble(line);
                                       //    System.out.println("Cumulative Fortran (ODE) timing: " + Global.fortrantime);
                                       // }
                                       // else{
                                       //    //read in time-stepping timing
                                       //     Global.timesteppingtime=Global.timesteppingtime + Double.parseDouble(line);
                                       //     System.out.println("Cumulative time-stepping timing: " + Global.timesteppingtime);
                                       //     error = true;
                                       // }
				}
                                
			}
        	int exitValue = ODESolver.waitFor();
        	
        }
        catch (Exception e) {
        	String err = "Error in running ODESolver \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        
        //11/1/07 gmagoon: renaming RWORK and IWORK files
        File f = new File("ODESolver/RWORK.dat");
        File newFile = new File("ODESolver/RWORK_"+index+".dat");
        if(newFile.exists())//09/08/08 gmagoon: delete newFile if it already exists to avoid "false" (failed) result when renaming; if file already exists, renaming appears to fail, returning false, on Windows Vista; I didn't notice this before on XP, but I may have missed it; three other instances of this modification were also made 
            newFile.delete();
        f.renameTo(newFile);
        f = new File("ODESolver/IWORK.dat");
        newFile = new File("ODESolver/IWORK_"+index+".dat");
        if(newFile.exists())
            newFile.delete();
        f.renameTo(newFile);
        //read the result
        File SolverOutput = new File("ODESolver/SolverOutput.dat");
        try {
        	FileReader fr = new FileReader(SolverOutput);
        	BufferedReader br = new BufferedReader(fr);
        	String line = br.readLine();
        	Global.solverIterations = Integer.parseInt(line.trim());
        	line = br.readLine();
        	
        	if (Double.parseDouble(line.trim()) != neq) {
        		System.out.println("ODESolver didnt generate all species result");
        		System.exit(0);
        	}
        	endTime = Double.parseDouble(br.readLine().trim());
        	for (int i=0; i<neq; i++){
        		line = br.readLine();
        		y[i] = Double.parseDouble(line.trim());
        	}
        	
        	for (int i=0; i<neq; i++){
        		line = br.readLine();
        		yprime[i] = Double.parseDouble(line.trim());
        	}
        	reactionFlux = new double[rList.size()+thirdBodyList.size()+troeList.size()];
        	for (int i=0; i<rList.size()+thirdBodyList.size()+troeList.size(); i++){
        		line = br.readLine();
        		reactionFlux[i] = Double.parseDouble(line.trim());
        	}
        	
        }
        catch (IOException e) {
        	String err = "Error in reading Solver Output File! \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        SolverOutput.delete();

		return 1;
	}
	private void initializeWorkSpace() {
		for (int i=0; i<30; i++)
			info[i] = 0;
		info[2] = 1; //give solution at all intermediate times
		//info[4] = 1; //use analytical jacobian
		
		
	}
	private void initializeConcentrations(SystemSnapshot p_beginStatus, ReactionModel p_reactionModel, ReactionTime p_beginTime, ReactionTime p_endTime, LinkedList initialSpecies) {
    	//System.out.println("After ODE:  from " + String.valueOf(p_beginTime.time) + "SEC to "+String.valueOf(p_endTime.time)+"SEC");
		//System.out.println("End at : " + String.valueOf(p_endTime.time) + "SEC");
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
	//## operation transferReaction(Reaction,SystemSnapshot,Temperature,Pressure)
    public ODEReaction transferReaction(Reaction p_reaction, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
        //#[ operation transferReaction(Reaction,SystemSnapshot,Temperature,Pressure)
    	
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
        if(p_reaction instanceof PDepNetReaction) {
        	//10/25/07 gmagoon: updated to use calculateRate with system snapshot (to avoid use of Global.temperature and Global.pressure)
            SystemSnapshot currentTPSnapshot = new SystemSnapshot();//10/25/07 gmagoon: make currentTPsnapshot variable, which will be used to pass temperature and pressure to calculateRate
            currentTPSnapshot.setTemperature(p_temperature);
            currentTPSnapshot.setPressure(p_pressure);
            double rate = ((PDepNetReaction)p_reaction).calculateRate(currentTPSnapshot);
            if (String.valueOf(rate).equals("NaN")){
            	System.err.println(p_reaction.toChemkinString(currentTPSnapshot.getTemperature()) + "Has bad rate probably due to Ea<DH");
            	rate = 0;
            }
            currentTPSnapshot = null;
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
        //#]
    }

    public double getAtol() {
        return atol;
    }

    public int getReactionSize(){
    	return rList.size()+troeList.size()+thirdBodyList.size();
    }
    public int getMaxSpeciesNumber() {
        return maxSpeciesNumber;
    }

    public double getRtol() {
        return rtol;
    }

}

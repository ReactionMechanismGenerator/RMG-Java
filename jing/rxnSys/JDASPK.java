
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\JDASPK.java
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
import jing.rxn.Reaction;
import jing.rxn.Structure;
import jing.param.Global;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.param.ParameterInfor;

//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\JDASPK.java
//----------------------------------------------------------------------------

//## class JDASPK
public class JDASPK implements ODESolver{

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
    protected LinkedList thirdBodyList ;
    protected LinkedList troeList ;
	protected StringBuilder outputString ;
	protected StringBuilder thermoString = new StringBuilder();
	protected StringBuilder rString;
	protected StringBuilder troeString;
	protected StringBuilder tbrString;
	protected double [] reactionFlux;
	
    //## operation JDASPK()
    private  JDASPK() {
        //#[ operation JDASPK()
        //#]
    }
    //## operation JDASPK(double,double,int, InitialStatus)
    public  JDASPK(double p_rtol, double p_atol, int p_parameterInfor, InitialStatus p_initialStatus) {
        //#[ operation JDASPK(double,double,int, InitialStatus)
        rtol = p_rtol;
        atol = p_atol;

        parameterInfor = p_parameterInfor;
        initialStatus = p_initialStatus;


        //#]
    }

    public void addSA(){
      if (parameterInfor == 0) {
        parameterInfor = 1;
      }
    }



    //## operation generatePDepODEReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
    public StringBuilder generatePDepODEReactionList(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
        //#[ operation generatePDepODEReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
        StringBuilder rString = new StringBuilder();
        StringBuilder arrayString = new StringBuilder();
        StringBuilder rateString = new StringBuilder();
        
        rList = new LinkedList();
    	LinkedList nonPDepList = new LinkedList();
        LinkedList pDepList = new LinkedList();

        LinkedHashSet pDepStructureSet = new LinkedHashSet();
        for (Iterator iter = PDepNetwork.getDictionary().values().iterator(); iter.hasNext(); ) {
        	PDepNetwork pdn = (PDepNetwork)iter.next();
        	for (Iterator pdniter = pdn.getPDepNetReactionList(); pdniter.hasNext();) {
        		PDepNetReaction pdnr = (PDepNetReaction)pdniter.next();
        		if (!pdnr.reactantEqualsProduct()) {
        			if (p_reactionModel instanceof CoreEdgeReactionModel) {
        				if (((CoreEdgeReactionModel)p_reactionModel).isReactedReaction(pdnr)) {
        					pDepList.add(pdnr);
        					pDepStructureSet.add(pdnr.getStructure());
        				}
        			}
        			else {
        				pDepList.add(pdnr);
        				pDepStructureSet.add(pdnr.getStructure());
        			}
        		}
        	}
        }
        
        for (Iterator iter = p_reactionModel.getReactionSet().iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	if (r.isForward() && !(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction)) {
        		Structure s = r.getStructure();
        		if (!pDepStructureSet.contains(s)) {
        			nonPDepList.add(r);
        		}
        	}
        }
        
        
        
		int size = nonPDepList.size() + pDepList.size();
       
       
       
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
			if (r.hasReverseReaction())
				arrayString.append(1 + " ");
			else
				arrayString.append(0 + " ");
			rateString.append(or.rate + " " + or.A + " " + or.n + " " + or.E + " "+or.Keq+" ");
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
    private double [] generateSensitivityStatus(ReactionModel p_reactionModel, double [] p_y, double [] p_yprime, int p_paraNum) {
    	//#[ operation generateSensitivityStatus(ReactionModel,double [],double [],int)
    	int neq = p_reactionModel.getSpeciesNumber()*(p_paraNum+1);
    	if (p_y.length != neq) throw new DynamicSimulatorException();
    	if (p_yprime.length != neq) throw new DynamicSimulatorException();

    	double [] senStatus = new double[nParameter*nState];

    	for (int i = p_reactionModel.getSpeciesNumber();i<neq;i++){
    		double sens = p_y[i];
    		int index = i-p_reactionModel.getSpeciesNumber();
    		senStatus[index] = p_y[i];
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
        	thermoString.append(p_species.calculateG(Global.temperature) + " ");
        }

        return id.intValue();
        //#]
    }

//  ## operation solve(boolean,ReactionModel,boolean,SystemSnapshot,ReactionTime,ReactionTime,Temperature,Pressure,boolean)
    public LinkedList solve(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, ReactionTime p_beginTime, ReactionTime p_endTime, Temperature p_temperature, Pressure p_pressure, boolean p_conditionChanged, int p_numSteps) {
    	
    	outputString = new StringBuilder();
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
		
        // set reaction set
        if (p_initialization || p_reactionChanged || p_conditionChanged) {
			
        	nState = p_reactionModel.getSpeciesNumber();
			
			
			//rString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0)
			rString = generatePDepODEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			//tbrString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10  (16 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0), ncollider, c1, c2,..c10 (20 elements)
			tbrString = generateThirdBodyReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			//troeString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10, alpha, Tstar, T2star, T3star, lowRate  (21 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0), ncollider, c1, c2,..c10, troe(0=T or 1=F) (21 elements)
        	troeString = generateTROEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
        	nParameter = 0;
			if (parameterInfor != 0) {
				nParameter = rList.size() + thirdBodyList.size() + troeList.size() + p_reactionModel.getSpeciesNumber();				
			}
			neq = nState*(nParameter + 1);
			initializeWorkSpace();
			initializeConcentrations(p_beginStatus, p_reactionModel, p_beginTime, p_endTime, initialSpecies);
			
			
        }
        outputString.append(nState +"\t" + neq+"\n");
        outputString.append( tBegin+" "+tEnd+" " + p_numSteps + "\n");
		for (int i=0; i<nState; i++)
			outputString.append(y[i]+" ");
		outputString.append("\n");
		for (int i=0; i<nState; i++)
			outputString.append(yprime[i]+" ");
		outputString.append("\n");
		for (int i=0; i<30; i++)
			outputString.append(info[i]+" ");
		outputString.append("\n"+ rtol + " "+atol);
		
		outputString.append("\n" + thermoString.toString() + "\n" + p_temperature.getK() + " " + p_pressure.getPa() + "\n" + rList.size() + "\n" + rString.toString() + "\n" + thirdBodyList.size() + "\n"+tbrString.toString() + "\n" + troeList.size() + "\n" + troeString.toString()+"\n");
		
        int idid=0;
        
        
		int temp = 1;
        Global.solverPrepossesor = Global.solverPrepossesor + (System.currentTimeMillis() - startTime)/1000/60;
		

        LinkedList systemSnapshotList = solveSEN(p_numSteps, p_reactionModel, p_beginTime, p_endTime, p_beginStatus);
        
        return systemSnapshotList;
    }
    
    //## operation solve(boolean,ReactionModel,boolean,SystemSnapshot,ReactionTime,ReactionTime,Temperature,Pressure,boolean)
    public SystemSnapshot solve(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, ReactionTime p_beginTime, ReactionTime p_endTime, Temperature p_temperature, Pressure p_pressure, boolean p_conditionChanged) {
    	
    	outputString = new StringBuilder();
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
		
        // set reaction set
        if (p_initialization || p_reactionChanged || p_conditionChanged) {
			
        	nState = p_reactionModel.getSpeciesNumber();
			
			
			//rString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0)
			rString = generatePDepODEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			//tbrString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10  (16 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0), ncollider, c1, c2,..c10 (20 elements)
			tbrString = generateThirdBodyReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			//troeString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10, alpha, Tstar, T2star, T3star, lowRate  (21 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0), ncollider, c1, c2,..c10, troe(0=T or 1=F) (21 elements)
        	troeString = generateTROEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
        	nParameter = 0;
			if (parameterInfor != 0) {
				nParameter = rList.size() + thirdBodyList.size() + troeList.size() + p_reactionModel.getSpeciesNumber();				
			}
			neq = nState*(nParameter + 1);
			initializeWorkSpace();
			initializeConcentrations(p_beginStatus, p_reactionModel, p_beginTime, p_endTime, initialSpecies);
			
			
        }
        outputString.append(nState +"\t" + neq+"\n");
        outputString.append( tBegin+" "+tEnd+"\n" );
		for (int i=0; i<nState; i++)
			outputString.append(y[i]+" ");
		outputString.append("\n");
		for (int i=0; i<nState; i++)
			outputString.append(yprime[i]+" ");
		outputString.append("\n");
		for (int i=0; i<30; i++)
			outputString.append(info[i]+" ");
		outputString.append("\n"+ rtol + " "+atol);
		
		outputString.append("\n" + thermoString.toString() + "\n" + p_temperature.getK() + " " + p_pressure.getPa() + "\n" + rList.size() + "\n" + rString.toString() + "\n" + thirdBodyList.size() + "\n"+tbrString.toString() + "\n" + troeList.size() + "\n" + troeString.toString()+"\n");
		
        int idid=0;
        LinkedHashMap speStatus = new LinkedHashMap();
        double [] senStatus = new double[nParameter*nState];
        
		int temp = 1;
        Global.solverPrepossesor = Global.solverPrepossesor + (System.currentTimeMillis() - startTime)/1000/60;
		if (nParameter==0) {
			startTime = System.currentTimeMillis();
        	//idid = solveDAE(p_initialization, reactionList, p_reactionChanged, thirdBodyReactionList, troeReactionList, nState, y, yprime, tBegin, tEnd, this.rtol, this.atol, T, P);
        	idid = solveDAE();
			if (idid !=1 && idid != 2 && idid != 3)	{
				System.out.println("The idid from DASPK was "+idid );
				throw new DynamicSimulatorException("DASPK: SA off.");
        	}
            System.out.println("After ODE: from " + String.valueOf(tBegin) + " SEC to " + String.valueOf(tEnd) + "SEC");
			Global.solvertime = Global.solvertime + (System.currentTimeMillis() - startTime)/1000/60;
			startTime = System.currentTimeMillis();
        	speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, 0);
			Global.speciesStatusGenerator = Global.speciesStatusGenerator + (System.currentTimeMillis() - startTime)/1000/60;
        }
        else {

        	//idid = solveSEN(p_initialization, reactionList, p_reactionChanged, thirdBodyReactionList, troeReactionList, nState, nParameter, this.parameterInforArray, y, yprime, tBegin, tEnd, this.rtol, this.atol, T, P);
        	idid = solveDAE();
        	
        	//if (idid != 2 && idid != 3) throw new DynamicSimulatorException("DASPK: SA on.");
        	speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, nParameter);
                System.out.println("After ODE: from " + String.valueOf(tBegin) + " SEC to " + String.valueOf(tEnd) + "SEC");
                 speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, nParameter);
                 senStatus = generateSensitivityStatus(p_reactionModel,y,yprime,nParameter);
                 SystemSnapshot sss = new SystemSnapshot(p_endTime, speStatus, senStatus, p_beginStatus.getTemperature(), p_beginStatus.getPressure());
                 sss.setIDTranslator(IDTranslator);
                 LinkedList reactionList = new LinkedList();
                 reactionList.addAll(rList);
                 reactionList.addAll(thirdBodyList);
                 reactionList.addAll(troeList);
                 sss.setReactionList(reactionList);
                 sss.setReactionFlux(reactionFlux);
                 return sss;

        }

		SystemSnapshot sss = new SystemSnapshot(p_endTime, speStatus, p_beginStatus.getTemperature(), p_beginStatus.getPressure());
		LinkedList reactionList = new LinkedList();
        reactionList.addAll(rList);
        reactionList.addAll(thirdBodyList);
        reactionList.addAll(troeList);
        sss.setReactionList(reactionList);
        sss.setReactionFlux(reactionFlux);
		
        return sss;
        //#]
    }


    
	private int solveDAE() {
		double startTime = System.currentTimeMillis();
		String workingDirectory = System.getProperty("RMG.workingDirectory");
		
		//write the input file
		File SolverInput = new File("ODESolver/SolverInput.dat");
		try {
			FileWriter fw = new FileWriter(SolverInput);
			fw.write(outputString.toString());
			fw.close();
		} catch (IOException e) {
			System.err.println("Problem writing Solver Input File!");
			e.printStackTrace();
		}
		Global.writeSolverFile +=(System.currentTimeMillis()-startTime)/1000/60;
		//run the solver on the input file
		boolean error = false;
        try {
        	 // system call for therfit
        	String[] command = {workingDirectory +  "/software/ODESolver/daspk.exe"};
			File runningDir = new File("ODESolver");
			
			Process ODESolver = Runtime.getRuntime().exec(command, null, runningDir);
			InputStream is = ODESolver.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				//System.out.println(line);
				line = line.trim();
				if (!(line.contains("ODESOLVER SUCCESSFUL"))) {
					System.err.println("Error running the ODESolver: "+line);
					error = true;
				}
			}
        	int exitValue = 4;
        	exitValue = ODESolver.waitFor();
        	System.out.println(br.readLine() + exitValue);
        	
        }
        catch (Exception e) {
        	String err = "Error in running ODESolver \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        
        startTime = System.currentTimeMillis();
        //read the result
        File SolverOutput = new File("ODESolver/SolverOutput.dat");
        try {
        	FileReader fr = new FileReader(SolverOutput);
        	BufferedReader br = new BufferedReader(fr);
        	String line = br.readLine();
        	//StringTokenizer st = new StringTokenizer(line);
        	Global.solverIterations = Integer.parseInt(line.trim());
        	line = br.readLine();
        	if (Double.parseDouble(line.trim()) != neq) {
        		System.out.println("ODESolver didnt generate all species result");
        		System.exit(0);
        	}
        	for (int i=0; i<nParameter+1; i++){
        		for (int j=0; j<nState; j++) {
        			line = br.readLine();
            		y[i*nState + j] = Double.parseDouble(line.trim());
        		}
        		line = br.readLine();
        	}
        	
        	for (int i=0; i<nParameter+1; i++){
        		for (int j=0; j<nState; j++) {
        			line = br.readLine();
            		yprime[i*nState + j] = Double.parseDouble(line.trim());
        		}
        		line = br.readLine();
        	}
        	reactionFlux = new double[rList.size()+thirdBodyList.size()+troeList.size()];
        	for (int i=0; i<rList.size()+thirdBodyList.size()+troeList.size(); i++){
        		line = br.readLine();
        		reactionFlux[i] = Double.parseDouble(line.trim());
        	}
        	for (int i=0; i<30; i++){
        		line = br.readLine();
        		info[i] = Integer.parseInt(line.trim());
        	}
        	
        }
        catch (IOException e) {
        	String err = "Error in reading Solver Output File! \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        Global.readSolverFile += (System.currentTimeMillis() - startTime)/1000/60;
		return 1;
	}
	
	
	private LinkedList solveSEN(int p_numSteps, ReactionModel p_reactionModel, ReactionTime tBegin, ReactionTime tEnd, SystemSnapshot p_beginStatus) {
		double startTime = System.currentTimeMillis();
		String workingDirectory = System.getProperty("RMG.workingDirectory");
		ReactionTime tStep = tEnd;
		LinkedList systemSnapshotList = new LinkedList();
		//write the input file
		File SolverInput = new File("ODESolver/SolverInput.dat");
		try {
			FileWriter fw = new FileWriter(SolverInput);
			fw.write(outputString.toString());
			fw.close();
		} catch (IOException e) {
			System.err.println("Problem writing Solver Input File!");
			e.printStackTrace();
		}
		Global.writeSolverFile +=(System.currentTimeMillis()-startTime)/1000/60;
		//run the solver on the input file
		boolean error = false;
        try {
        	 // system call for therfit
        	String[] command = {workingDirectory +  "/software/ODESolver/daspk.exe"};
			File runningDir = new File("ODESolver");
			
			Process ODESolver = Runtime.getRuntime().exec(command, null, runningDir);
			InputStream is = ODESolver.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				//System.out.println(line);
				line = line.trim();
				if (!(line.contains("ODESOLVER SUCCESSFUL"))) {
					System.err.println("Error running the ODESolver: "+line);
					error = true;
				}
			}
        	int exitValue = 4;
        	exitValue = ODESolver.waitFor();
        	System.out.println(br.readLine() + exitValue);
        	
        }
        catch (Exception e) {
        	String err = "Error in running ODESolver \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        
        startTime = System.currentTimeMillis();
        //read the result
        File SolverOutput = new File("ODESolver/SolverOutput.dat");
        try {
        	FileReader fr = new FileReader(SolverOutput);
        	BufferedReader br = new BufferedReader(fr);
        	String line ;

        	for (int k=0; k<p_numSteps; k++){
        		line = br.readLine();
            	if (Double.parseDouble(line.trim()) != neq) {
            		System.out.println("ODESolver didnt generate all species result");
            		System.exit(0);
            	}
            	for (int i=0; i<nParameter+1; i++){
            		for (int j=0; j<nState; j++) {
            			line = br.readLine();
                		y[i*nState + j] = Double.parseDouble(line.trim());
            		}
            		line = br.readLine();
            	}
            	
            	for (int i=0; i<nParameter+1; i++){
            		for (int j=0; j<nState; j++) {
            			line = br.readLine();
                		yprime[i*nState + j] = Double.parseDouble(line.trim());
            		}
            		line = br.readLine();
            	}
            	reactionFlux = new double[rList.size()+thirdBodyList.size()+troeList.size()];
            	for (int i=0; i<rList.size()+thirdBodyList.size()+troeList.size(); i++){
            		line = br.readLine();
            		reactionFlux[i] = Double.parseDouble(line.trim());
            	}
            	LinkedHashMap speStatus = new LinkedHashMap();
                double [] senStatus = new double[nParameter*nState];
            	
                System.out.println("After ODE: from " + String.valueOf(tBegin) + " SEC to " + String.valueOf(tEnd) + "SEC");
                speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, nParameter);
                senStatus = generateSensitivityStatus(p_reactionModel,y,yprime,nParameter);
                SystemSnapshot sss = new SystemSnapshot(tEnd, speStatus, senStatus, p_beginStatus.getTemperature(), p_beginStatus.getPressure());
                sss.setIDTranslator(IDTranslator);
                LinkedList reactionList = new LinkedList();
                reactionList.addAll(rList);
                reactionList.addAll(thirdBodyList);
                reactionList.addAll(troeList);
                sss.setReactionList(reactionList);
                systemSnapshotList.add(sss);
                sss.setReactionFlux(reactionFlux);
            	tBegin = tEnd;
            	tEnd = tEnd.add(tStep);
        	}
        	
        	
        }
        catch (IOException e) {
        	String err = "Error in reading Solver Output File! \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        Global.readSolverFile += (System.currentTimeMillis() - startTime)/1000/60;
		return systemSnapshotList;
	}
	
	
	private void initializeWorkSpace() {
		for (int i=0; i<30; i++)
			info[i] = 0;
		info[2] = 1; //print out the time steps
		info[4] = 1; //use analytical jacobian
		if (nParameter != 0) {
			info[18] = nParameter; //the number of parameters
			info[19] = 2; //perform senstivity analysis
			info[24] = 1;//staggered corrector method is used
		}
		
		
		
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

		if (nParameter != 0){//svp
						
			double [] sensitivityStatus = new double[nState*nParameter];
			int speciesNumber = p_reactionModel.getSpeciesNumber();
			
			for (int i=0; i<nParameter*speciesNumber;i++){
				sensitivityStatus[i] = 0;
			}
			p_beginStatus.addSensitivity(sensitivityStatus);							
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
        	double rate = ((PDepNetReaction)p_reaction).calculateTotalRate(p_beginStatus.temperature);
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

    public int getMaxSpeciesNumber() {
        return maxSpeciesNumber;
    }

    public double getRtol() {
        return rtol;
    }

}


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
import java.util.*;
import jing.chem.Species;
import jing.rxn.Reaction;
import jing.rxn.Structure;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.param.ParameterInfor;

//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\JDASPK.java
//----------------------------------------------------------------------------

//## class JDASPK
public class JDASPK implements SASolver, DAESolver {

    protected HashMap IDTranslator = new HashMap();		//## attribute IDTranslator

    protected double atol;		//## attribute atol

    protected int maxSpeciesNumber = 0;		//## attribute maxSpeciesNumber

     protected int parameterInfor;//svp

    protected ParameterInfor [] parameterInforArray = null;		//## attribute parameterInfor

    protected ODEReaction [] reactionList = null;		//## attribute reactionList

    protected double rtol;		//## attribute rtol

	protected TROEODEReaction [] troeReactionList;
	
    protected ThirdBodyODEReaction [] thirdBodyReactionList;		//## attribute thirdBodyReactionList

    protected InitialStatus initialStatus;//svp

    //protected int temp =1;
    // Constructors

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


    //## operation clean()
    public void clean() {
        //#[ operation clean()
        cleanDaspk();
        //#]
    }

    //## operation cleanDaspk()
    private native void cleanDaspk();

    //## operation fixingPDepRate(Structure)
    public double fixingPDepRate(Structure p_structure) {
        //#[ operation fixingPDepRate(Structure)
        //check c4h9o. beta scission
        	if (p_structure.getReactantNumber() == 1 && p_structure.getProductNumber() == 2) {
        	    ChemGraph r = (ChemGraph)p_structure.getReactants().next();
        	    Iterator it = p_structure.getProducts();
        	    ChemGraph p1 = (ChemGraph)it.next();
        	    ChemGraph p2 = (ChemGraph)it.next();

        	    if (r.getChemicalFormula().equals("C4H9O.")) {
        	    	if (p1.getChemicalFormula().equals("C2H5.") || p2.getChemicalFormula().equals("C2H5.") ) {
               			return 1E3;
        	     	}
        	    	else if (p1.getChemicalFormula().equals("C3H7.") || p2.getChemicalFormula().equals("C3H7.") ) {
               			return 1E3;
        	     	}
        	    	else if (p1.getChemicalFormula().equals("CH3.") || p2.getChemicalFormula().equals("CH3.") ) {
               			return 1E3;
        	     	}
        	    	else if (p1.getChemicalFormula().equals("H.") || p2.getChemicalFormula().equals("H.") ) {
               			return 1E3;
        	     	}
        	    }
        	}

        return 1;
        //#]
    }

    //## operation generateAllODEReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
    /*public void generateAllODEReactionList(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
        //#[ operation generateAllODEReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
        int size = p_reactionModel.getReactionSet().size();
        ODEReaction [] result = new ODEReaction[size];
        int id = 0;
        int thirdID = 0;
        Iterator iter = p_reactionModel.getReactionSet().iterator();
        while (iter.hasNext()) {
        	Reaction r = (Reaction)iter.next();
            if (r instanceof ThirdBodyReaction) {
            	thirdID++;
            	ODEReaction or = transferReaction(r, p_beginStatus, p_temperature, p_pressure);
        		result[size - thirdID] = or;
            }
            else {
        		id++;
            	ODEReaction or = transferReaction(r, p_beginStatus, p_temperature, p_pressure);
        		result[id-1] = or;
            }
        }

        reactionList = new ODEReaction[id];
        thirdBodyReactionList = new ODEReaction[thirdID];

        if (id+thirdID != size) throw new InvalidReactionSetException("Generating ODE reaction list for daspk");

        for (int i = 0; i < id; i++) {
        	reactionList[i] = result[i];

        }
        for (int i = id; i < size; i++) {
        	thirdBodyReactionList[i-id] = result[i];
        }

        return;
        //#]
    }*/

    //## operation generatePDepODEReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
    public LinkedList generatePDepODEReactionList(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
        //#[ operation generatePDepODEReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
        LinkedList nonPDepList = new LinkedList();
        LinkedList pDepList = new LinkedList();

        HashSet pDepStructureSet = new HashSet();
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
        //System.out.println("Total Number of pressure dependent reactions are "+ pDepList.size());
        for (Iterator iter = p_reactionModel.getReactionSet().iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	if (!r.reactantEqualsProduct() && !(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction)) {
        		Structure s = r.getStructure();
        		if (!pDepStructureSet.contains(s)) {
        			nonPDepList.add(r);
        		}
        	}
        }
        //System.out.println("Total Number of non pressure dependent reactions are "+ nonPDepList.size());

		int size = nonPDepList.size() + pDepList.size();
        reactionList = new ODEReaction[size];
        LinkedList all = new LinkedList();
        int id = 0;
        //System.out.println("non p_dep reactions: " + nonPDepList.size() );
        for (Iterator iter = nonPDepList.iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	
			if (!(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction)){
				all.add(r);
				ODEReaction or = transferReaction(r, p_beginStatus, p_temperature, p_pressure);
	            reactionList[id] = or;
				id++;
			}
        	
            //double rate = r.calculateRate(p_temperature);
            //if (r instanceof TemplateReaction) rate = ((TemplateReaction)r).calculatePDepRate(p_temperature);
            //System.out.println(r.getStructure().toString()+"\t rate = \t"+ String.valueOf(rate));
            //System.out.println(r.toChemkinString());
			
        }

        //System.out.println("p_dep reactions: " + pDepList.size());
        for (Iterator iter = pDepList.iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	all.add(r);

        	ODEReaction or = transferReaction(r, p_beginStatus, p_temperature, p_pressure);
            reactionList[id] = or;
            //System.out.println(r.getStructure().toString() + "\t rate = \t" + Double.toString(or.getRate()));
            //System.out.println(r.toChemkinString());
			id++;
        }

        return all;
        //#]
    }

    //## operation generateSpeciesStatus(ReactionModel,double [],double [],int)
    private HashMap generateSpeciesStatus(ReactionModel p_reactionModel, double [] p_y, double [] p_yprime, int p_paraNum) {
        //#[ operation generateSpeciesStatus(ReactionModel,double [],double [],int)
        int neq = p_reactionModel.getSpeciesNumber()*(p_paraNum+1);
        if (p_y.length != neq) throw new DynamicSimulatorException();
        if (p_yprime.length != neq) throw new DynamicSimulatorException();

        HashMap speStatus = new HashMap();

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
    public void generateThirdBodyReactionList(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
        //#[ operation generateThirdBodyReactionList(ReactionModel,SystemSnapshot,Temperature,Pressure)
        int size = p_reactionModel.getReactionSet().size();
        ThirdBodyODEReaction [] result = new ThirdBodyODEReaction[size];
        int thirdID = 0;
        Iterator iter = p_reactionModel.getReactionSet().iterator();
        while (iter.hasNext()) {
        	Reaction r = (Reaction)iter.next();

            if ((r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction)) {
            	thirdID++;
            	ThirdBodyODEReaction or = (ThirdBodyODEReaction)transferReaction(r, p_beginStatus, p_temperature, p_pressure);
        		result[thirdID-1] = or;
             }
        }

        thirdBodyReactionList = new ThirdBodyODEReaction[thirdID];

        for (int i = 0; i < thirdID; i++) {
        	thirdBodyReactionList[i] = result[i];
        }

        return;
        //#]
    }

	 private void generateTROEReactionList(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
		 int size = p_reactionModel.getReactionSet().size();
	        TROEODEReaction [] result = new TROEODEReaction[size];
	        int thirdID = 0;
	        Iterator iter = p_reactionModel.getReactionSet().iterator();
	        while (iter.hasNext()) {
	        	Reaction r = (Reaction)iter.next();

	            if (r instanceof TROEReaction) {
	            	thirdID++;
	            	TROEODEReaction or = (TROEODEReaction)transferReaction(r, p_beginStatus, p_temperature, p_pressure);
	        		result[thirdID-1] = or;
	             }
	        }

	        troeReactionList = new TROEODEReaction[thirdID];

	        for (int i = 0; i < thirdID; i++) {
	        	troeReactionList[i] = result[i];
	        }

	        return;
			
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
    public SystemSnapshot solve(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, ReactionTime p_beginTime, ReactionTime p_endTime, Temperature p_temperature, Pressure p_pressure, boolean p_conditionChanged) {
        //#[ operation solve(boolean,ReactionModel,boolean,SystemSnapshot,ReactionTime,ReactionTime,Temperature,Pressure,boolean)
        ReactionTime rt = p_beginStatus.getTime();
        if (!rt.equals(p_beginTime)) throw new InvalidBeginStatusException();

        // set time
        double tBegin = p_beginTime.getStandardTime();
        double tEnd = p_endTime.getStandardTime();

        // set reaction set
        //if (p_initialization || p_reactionChanged || p_conditionChanged) {
		generateTROEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
        generateThirdBodyReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
        LinkedList rList = generatePDepODEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);

        	//generateAllODEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
        	//p_reactionChanged = true;
        //}

        // set numbers
		//System.out.println("Total number of reactions to Daspk is "+rList.size());
        int nState = p_reactionModel.getSpeciesNumber();
        int nParameter = 0;
        LinkedList initialSpecies = new LinkedList();
//        if (parameterInfor!=null) nParameter = parameterInfor.length;
		if (parameterInfor != 0) {
			nParameter = rList.size()+thirdBodyReactionList.length; //svp
			if (initialStatus == null) System.out.println("initialStatus = null");
			Iterator spe_iter = initialStatus.getSpeciesStatus();
			while (spe_iter.hasNext()) {
				SpeciesStatus ss = (SpeciesStatus) spe_iter.next();
				String name = ss.getSpecies().getName();
				initialSpecies.add(name);
				nParameter++;
			}

		}

        int neq = nState*(nParameter+1);

        // set temperature and pressure
        double T = p_temperature.getK();
        double P = p_pressure.getAtm();

        // set initial value of y and yprime;
        double [] y = new double[neq];
        double [] yprime = new double[neq];

		//double numberOfReactedSpecies=0;
        // get the present status at t_begin, and set y and y' accordingly
        System.out.println("Before ODE: " + String.valueOf(tBegin) + "SEC");
        System.out.println("End at : " + String.valueOf(tEnd) + "SEC");
        for (Iterator iter = p_beginStatus.getSpeciesStatus(); iter.hasNext(); ) {
        	SpeciesStatus ss = (SpeciesStatus)iter.next();
        	double conc = ss.getConcentration();
        	double flux = ss.getFlux();
        	if (ss.isReactedSpecies()) {

        		Species spe = ss.getSpecies();
        		int id = getRealID(spe);
        		System.out.println(String.valueOf(spe.getID()) + '\t' + spe.getName() + '\t' + String.valueOf(conc) + '\t' + String.valueOf(flux));
         		y[id-1] = conc;
        		yprime[id-1] = flux;
        	}
        }

        if (nParameter != 0){//svp
			
			parameterInforArray = new ParameterInfor[nParameter];
			for (int i = 1; i <= rList.size()+thirdBodyReactionList.length; i++){
				parameterInforArray[i-1] = new ParameterInfor("k",i,0.00);
			}
			for (int i=rList.size()+1+thirdBodyReactionList.length; i<=nParameter;i++){
				parameterInforArray[i -
				                    1] = new ParameterInfor("CO", nParameter, 0.00);
			}

			if (p_beginTime.getTime() == 0 || p_beginTime.getTime() == 0.00) {
				LinkedList sensitivityStatus = new LinkedList();
				int reactionNumber = rList.size()+thirdBodyReactionList.length;
				int speciesNumber = p_reactionModel.getSpeciesNumber();
//            for (int i=0;i<reactionNumber*speciesNumber;i++){
				for (int i=0; i<nParameter*speciesNumber;i++){
					sensitivityStatus.add(i,null);
				}
				p_beginStatus.addSensitivity(sensitivityStatus);
				double sflux = 1;
				Iterator species_iter = p_reactionModel.getSpecies();
				while (species_iter.hasNext()) {
					Species spe = (Species) species_iter.next();
					int m = getRealID(spe);
					for (int p=0;p<nParameter;p++){
						int k = m + (p+1) * speciesNumber - 1;
                 //if (p >= rList.size()){
						if (p >= reactionNumber){
							int speciesNum = p-rList.size()-thirdBodyReactionList.length;
							String name = (String)initialSpecies.get(speciesNum);
							sflux = 0;
							SensitivityStatus senStatus;
							if (name.equalsIgnoreCase(spe.getName())){
								senStatus = new SensitivityStatus(1, sflux, m,
										p + 1);
								y[k] = 1;
							}
							else{
								senStatus = new SensitivityStatus(0,sflux,m,p+1);
								y[k]=0;
							}
							sensitivityStatus.add(k - p_reactionModel.getSpeciesNumber(),
									senStatus);
							p_beginStatus.putSensitivityStatus(k - speciesNumber, senStatus);
							
							yprime[k] = sflux;

						}
						else{
							if (p < rList.size()){
								Reaction rxn = (Reaction)rList.get(p);
								if (rxn.containsAsProduct(spe)) {
									sflux = 1;
									Iterator new_iter = rxn.getReactants();
									while (new_iter.hasNext()) {
										ChemGraph cg = (ChemGraph) new_iter.next();
										Species reactant = cg.getSpecies();
										SpeciesStatus ss = p_beginStatus.getSpeciesStatus(reactant);
										if (ss != null) {
											sflux *= ss.getConcentration();
										}
										else {
											sflux = 0;
										}
									}
									
								}
								if (rxn.containsAsReactant(spe)) {
									sflux = -1;
									Iterator new_iter = rxn.getReactants();
									while (new_iter.hasNext()) {
										ChemGraph cg = (ChemGraph) new_iter.next();
										Species reactant = cg.getSpecies();
										SpeciesStatus ss = p_beginStatus.getSpeciesStatus(reactant);
										if (ss != null) {
											sflux *= ss.getConcentration();
										}
										else {
											sflux = 0;
										}
									}
								}
								else {
									sflux = 0;
								}
								SensitivityStatus senStatus = new SensitivityStatus(0, sflux, m,
										p+1);
								sensitivityStatus.add(k-p_reactionModel.getSpeciesNumber(),senStatus);
								p_beginStatus.putSensitivityStatus(k-speciesNumber,senStatus);
								y[k] = 0;
								yprime[k] = sflux;
							}
							else {
								sflux = 0;
								SensitivityStatus senStatus = new SensitivityStatus(0, sflux, m,
										p+1);
								sensitivityStatus.add(k-p_reactionModel.getSpeciesNumber(),senStatus);
								p_beginStatus.putSensitivityStatus(k-speciesNumber,senStatus);
								y[k] = 0;
								yprime[k] = sflux;

							}
						}
					}
				}
			}
            else {
				for (int i = p_reactionModel.getSpeciesNumber(); i < y.length; i++) {
					SensitivityStatus ss = (SensitivityStatus) p_beginStatus.
					getSensitivityStatus(i-p_reactionModel.getSpeciesNumber());
					double sens = ss.getSensitivity();
					double sflux = ss.getSFlux();
					int reactionNumber = rList.size();
					int speciesNumber = p_reactionModel.getSpeciesNumber();
					Iterator species_iter = p_reactionModel.getSpecies();
					y[i] = sens;
					yprime[i] = sflux;
				}
			}
		}


		//System.out.println("Number of Reacted Species is " + numberOfReactedSpecies);

        int idid;
        HashMap speStatus = new HashMap();
        LinkedList senStatus = new LinkedList();
        double[] tPresent = {tBegin};
		int temp = 1;
        if (nParameter==0) {
        	idid = solveDAE(p_initialization, reactionList, true, thirdBodyReactionList, troeReactionList, nState, y, yprime, tBegin, tEnd, this.rtol, this.atol, T, P);
        	if (idid !=1 && idid != 2 && idid != 3)	{
				System.out.println("The idid from DASPK was "+idid + " at time "+tPresent[0]);
				throw new DynamicSimulatorException("DASPK: SA off.");
        	}
            System.out.println("After ODE: from " + String.valueOf(tBegin) + " SEC to " + String.valueOf(tEnd) + "SEC");

        	speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, 0);
        }
        else {
        	idid = solveSEN(p_initialization, reactionList, p_reactionChanged, thirdBodyReactionList, nState, nParameter, this.parameterInforArray, y, yprime, tBegin, tEnd, this.rtol, this.atol, T, P);
        	if (idid != 2 && idid != 3) throw new DynamicSimulatorException("DASPK: SA on.");
        	//speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, nParameter);
                System.out.println("After ODE: from " + String.valueOf(tBegin) + " SEC to " + String.valueOf(tEnd) + "SEC");
                 speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, nParameter);
                 senStatus = generateSensitivityStatus(p_reactionModel,y,yprime,nParameter);
                 SystemSnapshot sss = new SystemSnapshot(p_endTime, speStatus, senStatus, p_beginStatus.getTemperature(), p_beginStatus.getPressure());
                 sss.setIDTranslator(IDTranslator);
                 sss.setReactionList(rList);
                 return sss;

        }

		SystemSnapshot sss = new SystemSnapshot(p_endTime, speStatus, p_beginStatus.getTemperature(), p_beginStatus.getPressure());
		//ReactionSystem.outputAllPathways(SpeciesDictionary.getSpeciesFromName("CH4"), rList, sss, p_temperature);
        //ReactionSystem.outputAllPathways(SpeciesDictionary.getSpeciesFromName("CO"), rList, sss, p_temperature);
        //ReactionSystem.outputAllPathways(SpeciesDictionary.getSpeciesFromName("CO2"), rList, sss, p_temperature);

        return sss;
        //#]
    }

   
	//## operation solveDAE(boolean,ODEReaction [],boolean,ODEReaction [],int,double [],double [],double,double,double,double,double,double)
    private native int solveDAE(boolean p_initialization, ODEReaction [] p_reactionSet, boolean p_reactionChanged, ThirdBodyODEReaction [] p_thirdBodyReactionList, TROEODEReaction [] p_troeReactionList, int p_nState, double [] p_y, double [] p_yprime, double p_tBegin, double p_tEnd, double p_rtol, double p_atol, double p_temperature, double p_pressure);

    //## operation solveSEN(boolean,ODEReaction [],boolean,ODEReaction [],int,int,ParameterInfor [],double [],double [],double,double,double,double,double,double)
    private native int solveSEN(boolean p_initialization, ODEReaction [] p_reactionSet, boolean p_reactionChanged, ThirdBodyODEReaction [] p_thirdBodyReactionList, int p_nState, int p_nParameter, ParameterInfor [] p_parameterInfor, double [] p_y, double [] p_yprime, double p_tBegin, double p_tEnd, double p_rtol, double p_atol, double p_temperature, double p_pressure);
    static {System.loadLibrary("daspk");}


    //## operation transferReaction(Reaction,SystemSnapshot,Temperature,Pressure)
    public ODEReaction transferReaction(Reaction p_reaction, SystemSnapshot p_beginStatus, Temperature p_temperature, Pressure p_pressure) {
        //#[ operation transferReaction(Reaction,SystemSnapshot,Temperature,Pressure)
        double dT = 1;
        Temperature Tup = new Temperature(p_temperature.getStandard()+dT, Temperature.getStandardUnit());
        Temperature Tlow = new Temperature(p_temperature.getStandard()-dT, Temperature.getStandardUnit());

        int rnum = p_reaction.getReactantNumber();
        int pnum = p_reaction.getProductNumber();

        int [] rid = new int[rnum];
        int index = 0;
        for (Iterator r_iter = p_reaction.getReactants(); r_iter.hasNext(); ) {
        	Species s = ((ChemGraph)r_iter.next()).getSpecies();
        	rid[index] = getRealID(s);
        	index++;
        }

        int [] pid = new int[pnum];
        index = 0;
        for (Iterator p_iter = p_reaction.getProducts(); p_iter.hasNext(); ) {
        	Species s = ((ChemGraph)p_iter.next()).getSpecies();
        	pid[index] = getRealID(s);
        	index++;
        }
		
		//ODEReaction or;
        if(p_reaction instanceof PDepNetReaction) {
        	double rate = ((PDepNetReaction)p_reaction).getRate();
        	ODEReaction or = new ODEReaction(rnum, pnum, rid, pid, rate);
        	return or;
        }
        else {
        	double rate = 0;
        	if (p_reaction instanceof TemplateReaction) {
        		rate = ((TemplateReaction)p_reaction).calculatePDepRate(p_temperature);
				ODEReaction or = new ODEReaction(rnum, pnum, rid, pid, rate);
				return or;
			}
			else if (p_reaction instanceof TROEReaction){//svp
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
				rate = p_reaction.calculateRate(p_temperature);
				double inertColliderEfficiency = ((ThirdBodyReaction)p_reaction).calculateThirdBodyCoefficientForInerts(p_beginStatus);
				ThirdBodyODEReaction or = new ThirdBodyODEReaction(rnum, pnum, rid, pid, rate, colliders, efficiency,numCollider, inertColliderEfficiency);
				return or;
			}
			
			else{
				rate = p_reaction.calculateRate(p_temperature);
				ODEReaction or = new ODEReaction(rnum, pnum, rid, pid, rate);
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

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



package jing.rxn;


import java.io.*;
import jing.chem.*;
import java.util.*;
import jing.param.*;
import jing.rxnSys.*;
import jing.mathTool.*;
import jing.chemParser.*;
import jing.chem.Species;
import jing.rxnSys.ReactionSystem;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\PDepNetwork.java                                                                  
//----------------------------------------------------------------------------

//## class PDepNetwork 
public class PDepNetwork {
    
    protected static int ID = 0;		//## attribute ID 
    public static boolean generateNetworks;
    
    protected static HashMap dictionary = new HashMap();		//## attribute dictionary 
    
    protected Reaction entryReaction = null;		//## attribute entryReaction 
    
    protected boolean isChemAct;		//## attribute isChemAct 
    
    protected double[] kLeak;		//## attribute kLeak 
    
    protected LinkedList product = null;		//## attribute product 
    
    protected LinkedList reactant = null;		//## attribute reactant 
    
    protected static int totalNumber = 0;		//## attribute totalNumber 
    //10/30/07 gmagoon: note that implementation would probably need to be changed if temperature or pressure change with time
    //(10/30/07) UPDATE: actually, since current T,P is passed to updateKLeak and initial T in temperatureArray is only used in initializeKLeak, perhaps it is not an issue (?); should not be an issue if initializeKLeak is only called at t = 0 and kLeak is updated before being used at subsequent times; even so, updateKLeak may not be called often enough for it to correspond to current T,P in current implementation (only updated when species is added, I think?)
    protected static LinkedList temperatureArray; //10/30/07 gmagoon: added LinkedList for temperatures as static variable: note that temperatures will be repeated in cases where there are multiple pressures
    protected static LinkedList pressureArray; //10/30/07 gmagoon: added LinkedList for pressures as static variable: note that pressures will be repeated in cases where there are multiple temperatures;//UPDATE: commenting out: not needed if updateKLeak is done for one temperature/pressure at a time; //11/1-2/07 restored
    
    protected LinkedList pDepNetReactionList;
    protected LinkedList pDepNonincludedReactionList;
    protected LinkedList pDepWellList;
    
    protected HashSet nonIncludedSpecies = new HashSet();
    protected boolean altered = true;
    
    // Constructors
    
    //## operation PDepNetwork() 
    private  PDepNetwork() {
        {
            pDepNetReactionList=new LinkedList();
        }
        {
            pDepNonincludedReactionList=new LinkedList();
        }
        {
            pDepWellList=new LinkedList();
        }
        //#[ operation PDepNetwork() 
        ID = totalNumber++;
        //#]
    }

	

	//## operation addPDepWell(PDepWell) 
    private void addPDepWell(PDepWell p_newWell) {
        //#[ operation addPDepWell(PDepWell) 
        if (pDepWellList.contains(p_newWell)) return;
        
        pDepWellList.add(p_newWell);
        
        decideWellPathType();
        
        return;
        
        
        
        //#]
    }
    
    //## operation decideWellPathType() 
    private void decideWellPathType() {
        //#[ operation decideWellPathType() 
        for (Iterator iter = getPDepWellListIterator(); iter.hasNext(); ) {
        	PDepWell pdw = (PDepWell)iter.next();
        	for (Iterator pathIter = pdw.getPathsIterator(); pathIter.hasNext(); ) {
        		PDepPathReaction pdpr = (PDepPathReaction)pathIter.next();
        		int pNum = pdpr.getProductNumber();
        		LinkedList r = pdpr.getReactantList();
        		LinkedList p = pdpr.getProductList();
        
        		if (pNum == 1) {
        			Species spe = (Species)p.iterator().next();
        			if (includeAsIsomer(spe)) pdpr.setTypeAsIsomer();
        			else {
        				pdpr.setTypeAsNonIncluded();
        				nonIncludedSpecies.add(spe);
        			}
        		}
        		else if (pNum == 2) {
        			if (isChemAct) {
        				if (MathTool.isListEquivalent(r, product) && MathTool.isListEquivalent(p,reactant)) {
        				    pdpr.setTypeAsReactant();
        				}
        				else pdpr.setTypeAsProduct();
        			}
        			else {
        				pdpr.setTypeAsProduct();
        			}
        		}
        	}
        }
        return;
        
        
        //#]
    }
    
    //## operation getEntryMass() 
    public double getEntryMass() {
        //#[ operation getEntryMass() 
        double mass = 0;
        
        for (Iterator iter = reactant.iterator(); iter.hasNext();) {
        	Species spe = (Species)iter.next();
        	mass += spe.getMolecularWeight();
        }
        
        return mass;
        //#]
    }
    
    //## operation includeAsIsomer(Species) 
    public boolean includeAsIsomer(Species p_species) {
        //#[ operation includeAsIsomer(Species) 
        for (Iterator iter = getPDepWellListIterator(); iter.hasNext(); ) {
        	PDepWell pdw = (PDepWell)iter.next();
        	if ((pdw.getIsomer()).equals(p_species)) return true;
        }
        return false;
        //#]
    }
    
    //## operation initializeKLeak() 
    //10/26/07 gmagoon: changed to take temperature as parameter as part of elimination of uses of Global.temperature
    //10/30/07 gmagoon: removed temperature
    private void initializeKLeak() {
    //private void initializeKLeak(Temperature t) {
        //#[ operation initializeKLeak() 
      
        //10/30/07 gmagoon: updating to use kLeak array
        kLeak = new double[temperatureArray.size()];//10/30/07 gmagoon: initialize size of kLeak//10/31/07: moved to setTemperatureArray;
        for (Integer i = 0; i<temperatureArray.size();i++) {
             kLeak[i] = 0;
        }
        //kleak = 0
       
       
        //Temperature t = new Temperature(715,"K");
        //10/26/07 gmagoon: changed to avoid use of Global.temperature; t is passed as parameter
        //Temperature t = Global.temperature;
        
        if (isChemAct) {
            //10/30/07 gmagoon: updating to use kLeak array
            for (Integer i = 0; i<temperatureArray.size();i++) {
                kLeak[i] = entryReaction.calculateTotalRate((Temperature)temperatureArray.get(i));
            }
            //kLeak = entryReaction.calculateTotalRate(t);
        }
        else {
        	if (pDepWellList.size()!=1) {
        		System.out.println("Dissoc is not initialized correctly" + toString());
        		System.exit(0);
        	}
        	Iterator iter = getPDepWellListIterator();
        	PDepWell pdw = (PDepWell)iter.next();
        	for (Iterator pathiter = pdw.getPathsIterator(); pathiter.hasNext(); ) {
        		PDepPathReaction pdpr = (PDepPathReaction)pathiter.next();
                        //10/30/07 gmagoon: updating to use kLeak array
                        for (Integer i = 0; i<temperatureArray.size();i++) {
                            kLeak[i] += pdpr.getReaction().calculateTotalRate((Temperature)temperatureArray.get(i));
                        }
        		//kLeak += pdpr.getReaction().calculateTotalRate(t);
        	}
        }
       
      
    }
   
    //10/25/07 gmagoon: updated to take temperature, pressure
    //10/30/07 gmagoon: switched back to not take temperature/pressure when it is called for a particular reaction system; this may lead to inefficiencies\
    //UPDATE: temperature, pressure, included, along with index to specify which kLeak to update
    //11/1-2/07 gmagoon: updating to use temperature and pressure array and update all at once; updating one at a time may be causing problems; **note that this would need to be fixed in cases where temperature and pressure are not constant, since initial temp and pressure are in arrays
     public void updateKLeak(){    
    //private void updateKLeak(Temperature p_temperature, Pressure p_pressure, Integer p_index){
    //private void updateKLeak(Temperature p_temperature, Pressure p_pressure){    	
  //      if(kLeak == null)//10/31/07 gmagoon: set size of kLeak****may need further invesitgation; ideally, initializeKLeak should set size of kLeak I think; 11/6/07 gmagoon: not needed after temperatureArray, pressureArray initialization was moved before lrg initialization
  //          kLeak = new double[temperatureArray.size()];
       
        for(int i = 0; i < temperatureArray.size(); i++){
            Temperature p_temperature = (Temperature)temperatureArray.get(i);
            Pressure p_pressure = (Pressure)pressureArray.get(i);
            if (!isActive() && getIsChemAct()){
                    double rate = entryReaction.calculateTotalRate(p_temperature);//10/25/07 gmagoon: changed to avoid use of Global.temperature
                    //double rate = entryReaction.calculateTotalRate(Global.temperature);
                    //kLeak = rate;
                    kLeak[i] = rate;//10/30/07 gmagoon: updated to change a particular element of kLeak array;
                    return;
            }

            //kLeak = 0;
            kLeak[i] = 0; //10/30/07 gmagoon: updated to change a particular element of kLeak array;
            for (Iterator iter = getPDepNonincludedReactionListIterator(); iter.hasNext();){
                    PDepNetReaction pdnr = (PDepNetReaction) iter.next();
                    //10/25/07 gmagoon: updated to use calculateRate with system snapshot (to avoid use of Global.temperature and Global.pressure)
                    SystemSnapshot currentTPSnapshot = new SystemSnapshot();//10/25/07 gmagoon: make currentTPsnapshot variable, which will be used to pass temperature and pressure to calculateRate
                    currentTPSnapshot.setTemperature(p_temperature);
                    currentTPSnapshot.setPressure(p_pressure);
                    double rate = pdnr.calculateRate(currentTPSnapshot);
                    //double rate = pdnr.calculateRate();
                    if (rate < 0 ) throw new InvalidPDepNetReactionException();
                    //kLeak += rate;
                    kLeak[i] += rate;//10/30/07 gmagoon: updated to change a particular element of kLeak array;
            }
        }
    	return;
    }
    
    //## operation initializePDepNetwork(TemplateReaction) 
    private void initializePDepNetwork(Reaction p_entryReaction) {
        //#[ operation initializePDepNetwork(TemplateReaction) 
        reactant = (LinkedList)p_entryReaction.getReactantList().clone();
        if (reactant.size() == 1) {
        	isChemAct = false;
        	entryReaction = null;
         	Species spe = (Species)p_entryReaction.getReactants().next();
         	PDepWell pdw = new PDepWell(spe, p_entryReaction);
        
         	// add a new PDepWell into PDepWell list
         	addPDepWell(pdw);
         	product = null; 
         	if (p_entryReaction.structure.products.size() == 1){
         		nonIncludedSpecies.add(p_entryReaction.structure.products.get(0));
         	}
         		
        }
        else if (reactant.size() == 2) {
        	isChemAct = true;
        	entryReaction = p_entryReaction;
        	product = (LinkedList)p_entryReaction.getProductList().clone();
        	nonIncludedSpecies.add(product.get(0));
        }
        else throw new InvalidPDepNetworkTypeException();
        
        //#]
    }
    
    //## operation isActive() 
    public boolean isActive() {
        //#[ operation isActive() 
        return (!pDepWellList.isEmpty());
        
        //#]
    }
    
    //## operation makePDepNetwork(TemplateReaction) 
    public static PDepNetwork makePDepNetwork(Reaction p_entryReaction) {
        //#[ operation makePDepNetwork(TemplateReaction) 
	
    	Object key = null;
    	
    	if (!generateNetworks)
    		return null;
    	if (p_entryReaction.getReactantNumber() == 2) {
    		Species pro = (Species)p_entryReaction.getProducts().next();
        	if (pro.isTriatomicOrSmaller()) return null;
        	key = p_entryReaction.getStructure();
        }
        else {
        	key = (Species)p_entryReaction.getReactants().next();
        	if (((Species)key).isTriatomicOrSmaller()) return null;
        }
        
        PDepNetwork obj = (PDepNetwork)dictionary.get(key);
        if (obj == null) {
        	// if it is a->b, but a doesn't have three frequency model, return null
        	if (p_entryReaction.getReactantNumber() == 1) {
         		Species spe = (Species)p_entryReaction.getReactants().next();
         		if (!spe.hasThreeFrequencyModel()) {
         			return null;
         		}
           	}
        	
        	PDepNetwork pnw = new PDepNetwork();
        	dictionary.put(key, pnw);
        	pnw.initializePDepNetwork(p_entryReaction);
        	pnw.initializeKLeak();
        	return pnw;
        
        }
        else if (p_entryReaction.getReactantNumber() == 1){
        	PDepWell pdw = (PDepWell)obj.pDepWellList.get(0);
        	pdw.addPath(p_entryReaction);
        	obj.decideWellPathType();
        	obj.initializeKLeak();
        	obj.altered = true;
        	return obj;
        }
        
        else return (PDepNetwork)obj;
         
        
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        if (entryReaction != null) {
        	return entryReaction.getStructure().toString();
        }
        else {
        	Species spe = (Species)reactant.getFirst();
        	return spe.getName()+"("+ String.valueOf(spe.getID()) +")";
        	
        }
        
        //#]
    }
    
    //## operation totalSize() 
    public static int totalSize() {
        //#[ operation totalSize() 
        return dictionary.values().size();
        //#]
    }
    
    //## operation update(Species) 
    public void update(Species p_nextIsomer) {
        //#[ operation update(Species) 
        if (includeAsIsomer(p_nextIsomer)) throw new DuplicatedIsomersException();
        if (!p_nextIsomer.hasThreeFrequencyModel()) {
        	System.out.println("Species doesn't have three frequency model: " + p_nextIsomer.toString());
        	System.exit(0);	
        }
        //System.out.println("begin to new a well");
        PDepWell pdw = new PDepWell(p_nextIsomer, p_nextIsomer.getPdepPaths());
        //System.out.println("begin to add new well to system");
        addPDepWell(pdw);
        if (nonIncludedSpecies.contains(p_nextIsomer)){
        	nonIncludedSpecies.remove(p_nextIsomer);
        }

        //System.out.println("finish adding new well to system");
        return;
        //#]
    }
    
    public static int getID() {
        return ID;
    }
    
    public static HashMap getDictionary() {
        return dictionary;
    }
    
    public Reaction getEntryReaction() {
        return entryReaction;
    }
    
    public boolean getIsChemAct() {
        return isChemAct;
    }
    
    //10/30/07 gmagoon: modified to access one kLeak out of array at various temperatures
    public double getKLeak(Integer p_index)
    {
        return kLeak[p_index];
    }
    
	public void setKLeak(int index, double rate) {
		kLeak[index] = rate;
	}
    
    
	
    public LinkedList getProduct() {
        return product;
    }
    
    public LinkedList getReactant() {
        return reactant;
    }
    
    public static int getTotalNumber() {
        return totalNumber;
    }
    
    public ListIterator getPDepNetReactionListIterator() {
        ListIterator iter=pDepNetReactionList.listIterator(0);
        return iter;
    }
    
    public ListIterator getPDepNonincludedReactionListIterator() {
        ListIterator iter=pDepNonincludedReactionList.listIterator(0);
        return iter;
    }
	
	public LinkedList getPDepNetReactionList() {
        return pDepNetReactionList;
    }
    
    public LinkedList getPDepNonincludedReactionList() {
        return pDepNonincludedReactionList;
    }
    
    public ListIterator getPDepWellListIterator() {
        ListIterator iter=pDepWellList.listIterator(0);
        return iter;
    }
	
	public LinkedList getPDepWellList() {
        return pDepWellList;
    }
    
    public boolean getAltered(){
		return altered;
	}
    
	public void setAltered(boolean alt) {
		altered = alt;
	}
	
    public Iterator getNonIncludedSpecies(){
    	return nonIncludedSpecies.iterator();
    }
	
     //10/30/07: gmagoon: added function to set temperatureArray (used with kLeak); I think it should be declared as static so it belongs to class rather than object and temperatureArray will be static
    //****note that there is some inefficiency in the current implementation since temperatures will be repeated if there are multiple pressures; if calculating kLeak takes a long time, changing implementation may be warranted; UPDATE: if kLeak depends on pressure as well, there will be no inefficiency
    public static void setTemperatureArray(LinkedList p_tempArray){
        temperatureArray = p_tempArray;
    }

     //10/30/07: gmagoon: added function to set temperatureArray (used with kLeak); see comments for previous function
     //UPDATE: commenting out: not needed if updateKLeak is done for one temperature/pressure at a time
     //11/1-2/07 gmagoon: restoring
   public static void setPressureArray(LinkedList p_pressureArray){
        pressureArray = p_pressureArray;
    }
    
       /*
	public static void completeNetwork(Species p_species, HashSet hs) {
		if (p_species == null) return;
		if (p_species.getPdepPaths() == null) return;
		
		Iterator iter = dictionary.values().iterator();
		while (iter.hasNext()){
			PDepNetwork pdn = (PDepNetwork)iter.next();
			if (pdn.nonIncludedSpecies.contains(p_species)){
				PDepWell pdw = new PDepWell(p_species, p_species.getPdepPaths());
				
				pdn.completeNetwork(pdw, hs);
				//pdn.addPDepWell(pdw);
				//pdn.nonIncludedSpecies.remove(p_species);
				
				pdn.altered = true;
			}
			
		}
		
	}
	
	private void completeNetwork(PDepWell pdw, HashSet hs) {
		Iterator iter = pdw.isomer.getPdepPaths().iterator();
		addPDepWell(pdw);
		nonIncludedSpecies.remove(pdw.isomer);
		while (iter.hasNext()){
			Reaction r = (Reaction)iter.next();
			Species spe = (Species)r.getStructure().products.get(0);
			if (r.getProductNumber() == 1 && this.nonIncludedSpecies.contains(spe) && hs.contains(spe)) {
				PDepWell pdw_new = new PDepWell(spe, spe.getPdepPaths());
				completeNetwork(pdw_new, hs);
			}
		}
		
	}

	
	public static void completeNetwork(HashSet hs) {
		Iterator iter = hs.iterator();
		while (iter.hasNext()){
			Species sp = (Species)iter.next();
			completeNetwork(sp, hs);
		}
	}
    */
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\PDepNetwork.java
*********************************************************************/


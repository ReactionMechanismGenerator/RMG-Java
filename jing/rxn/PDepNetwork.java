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
        for (Iterator iter = getPDepWellList(); iter.hasNext(); ) {
        	PDepWell pdw = (PDepWell)iter.next();
        	for (Iterator pathIter = pdw.getPaths(); pathIter.hasNext(); ) {
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
    private double getEntryMass() {
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
    private boolean includeAsIsomer(Species p_species) {
        //#[ operation includeAsIsomer(Species) 
        for (Iterator iter = getPDepWellList(); iter.hasNext(); ) {
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
        	Iterator iter = getPDepWellList();
        	PDepWell pdw = (PDepWell)iter.next();
        	for (Iterator pathiter = pdw.getPaths(); pathiter.hasNext(); ) {
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
     private void updateKLeak(){    
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
            for (Iterator iter = getPDepNonincludedReactionList(); iter.hasNext();){
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
    
    //## operation runPDepCalculation(ReactionSystem) 
    public void runPDepCalculation(ReactionSystem p_reactionSystem) {
        //#[ operation runPDepCalculation(ReactionSystem) 
        if (!isActive() && getIsChemAct()) {
        	Temperature t = p_reactionSystem.getPresentTemperature();
  		    kLeak[p_reactionSystem.getIndex()] = getEntryReaction().calculateTotalRate(t);//10/30/07 gmagoon: updated to change kLeak corresponding to reactionSystem being passed
        	return;
        }
        
        // construct the input file for chemdis
        String result = writePDepNetworkHeader(p_reactionSystem);
        int index = 0;
        for (Iterator iter = getPDepWellList(); iter.hasNext(); ) {
        	index++;
        	PDepWell pdw = (PDepWell)iter.next();
        	result += "Well " + String.valueOf(index) + '\n';
        	result += pdw.toChemDisString();
        }
        
//		 run chemdis
        String dir = System.getProperty("RMG.workingDirectory");
        
        File chemdis_input;
        
        try {
        	// prepare chemdis input file, "fort.10" is the input file name
        	chemdis_input = new File("chemdis/fort.10");
        	FileWriter fw = new FileWriter(chemdis_input);
        	fw.write(result);
        	fw.close();
        }
        catch (IOException e) {
        	System.out.println("Wrong input file for chemdis!");
        	System.exit(0);
        }
        
        try {
           	// system call for chemdis
           	String[] command = {dir + "/software/chemdis/chemdis.exe"};
           	File runningDir = new File("chemdis");
            Process chemdis = Runtime.getRuntime().exec(command, null, runningDir);                     
            InputStream ips = chemdis.getInputStream();
            InputStreamReader is = new InputStreamReader(ips);
            BufferedReader br = new BufferedReader(is);
            String line=null;
            while ( (line = br.readLine()) != null) {
            	//System.out.println(line);
            }
            int exitValue = chemdis.waitFor();
        }
        catch (Exception e) {
        	System.out.println("Error in run chemdis!");
        	System.exit(0);
        }
        
        parseChemdisOutputCP();
        updateKLeak();
        
        altered = false;
        
        
        
        //#]
    }
    
//  ## operation parseChemdisOutputCP() 
    private void parseChemdisOutputCP() {
        //#[ operation parseChemdisOutputCP() 
        try {
        	String dir = System.getProperty("RMG.workingDirectory");
        	String chemdis_output = "chemdis/chemdis-rmg.out";
        
        	FileReader in = new FileReader(chemdis_output);
        	BufferedReader data = new BufferedReader(in);
        
        	String line = ChemParser.readMeaningfulLine(data);
        	line = ChemParser.readMeaningfulLine(data);
        	line = line.trim();
        	int rNum = 0;
        	if (line.startsWith("CHEMACT")) {
        		rNum = 2;
        	}                        
        	else if (line.startsWith("DISSOC")) {
        	  	rNum = 1;
        	}
        	else {
        		System.out.println("Wrong output from chemdis: unknown type!");
        		System.out.println("Unknown key word for PDep Network: " + line);
        		System.exit(0);  	
        	}
        	LinkedList reactant = new LinkedList();
        	StringTokenizer st = new StringTokenizer(line);
        	String type = st.nextToken();
        	String temp = st.nextToken();
        	temp = st.nextToken();
        	temp = st.nextToken();
        	String r1 = st.nextToken().trim();
        	r1 = r1.substring(3, r1.length());
        	int idr1 = Integer.parseInt(r1);
        	Species sr1 = SpeciesDictionary.getInstance().getSpeciesFromID(idr1);
        	String newName = sr1.getName()+"("+String.valueOf(sr1.getID())+")";
        	reactant.add(sr1);
        	if (rNum == 2) {
        		temp = st.nextToken();
        		String r2 = st.nextToken().trim();
        		r2 = r2.substring(3, r2.length());
        		int idr2 = Integer.parseInt(r2);
        		Species sr2 = SpeciesDictionary.getInstance().getSpeciesFromID(idr2); 
        		newName += "+" + sr2.getName()+"("+String.valueOf(sr2.getID())+")";
        		reactant.add(sr2);
        	}
        	
        	double Tmax=0;
        	double Tmin=0;
        	double Pmax=0;
        	double Pmin=0;
        	
        	int nT = 7; 
        	int nP = 4;
        	
        	line = ChemParser.readMeaningfulLine(data);
        	if (line.startsWith("Temperature range")) {
        		line = ChemParser.readMeaningfulLine(data);
        		st = new StringTokenizer(line);
        		String tL = st.nextToken().trim();
        		String tH = st.nextToken().trim();
        		Tmin = Double.parseDouble(tL);
        		Tmax = Double.parseDouble(tH);
        	}
        	else {
        		System.out.println("Can't read T range from chemdis output file!");
        		System.exit(0);  	
        	}
                                                 
        	line = ChemParser.readMeaningfulLine(data);
        	if (line.startsWith("Pressure range")) {
        		line = ChemParser.readMeaningfulLine(data);
        		st = new StringTokenizer(line);
        		String pL = st.nextToken().trim();
        		String pH = st.nextToken().trim();
        		Pmin = Double.parseDouble(pL);
        		Pmax = Double.parseDouble(pH);
        	}
        	else {
        		System.out.println("Can't read P range from chemdis output file!");
        		System.exit(0);  	
        	}
        	
        	pDepNetReactionList.clear();
        	pDepNonincludedReactionList.clear();
        	line = ChemParser.readMeaningfulLine(data);
        	while (!line.startsWith("END")) {
        		line = line.trim();
        		LinkedList product = new LinkedList();
        		st = new StringTokenizer(line);
        		String rxntype = st.nextToken();
        		int pNum = 1;
        		String p1 = st.nextToken().trim();
        		p1 = p1.substring(3, p1.length());
        		int idp1 = Integer.parseInt(p1);
        		Species sp1 = SpeciesDictionary.getInstance().getSpeciesFromID(idp1);
        	    product.add(sp1);
        		if (st.hasMoreTokens()) {
        			String next = st.nextToken();
         			if ((next.trim()).equals("+")) {	
        				String p2 = st.nextToken().trim();
        				p2 = p2.substring(3, p2.length());
        				int idp2 = Integer.parseInt(p2);
        				Species sp2 = SpeciesDictionary.getInstance().getSpeciesFromID(idp2);
        				product.add(sp2);
        				pNum++;
        			}
        		}
        		
        		// read chebyshev polynomial
        		double [][] alpha = new double[nT][nP];
        		for (int i = 0; i < nT; i++) {
        			line = ChemParser.readMeaningfulLine(data);
        			st = new StringTokenizer(line);
        			for (int j = 0; j < nP; j++) {
        				String a = st.nextToken().trim();
        				alpha[i][j] = Double.parseDouble(a);
        			}
        		}
        		
        		Temperature tLow = new Temperature(Tmin, "K");
        		Temperature tHigh = new Temperature(Tmax, "K");
        		Pressure pLow = new Pressure(Pmin, "Atm");
        		Pressure pHigh = new Pressure(Pmax, "Atm");
        		ChebyshevPolynomials cp = new ChebyshevPolynomials(nT, tLow, tHigh, nP, pLow, pHigh, alpha);
        		
        		PDepNetReaction pdnr = new PDepNetReaction(reactant, product, cp);
        		if (rxntype.equals("ISOMER")) {
        			//if pDepNetReactionList contains this reaction simply add the chebyshev
        			// rate of this reaction to the one already present in the system.
        			if (pDepNetReactionList.contains(pdnr)){
        				int i = pDepNetReactionList.indexOf(pdnr);
        				PDepNetReaction pdnrAlreadyPresent = (PDepNetReaction)pDepNetReactionList.get(i);
        				pdnrAlreadyPresent.itsChebyshevPolynomials.addChebyshevPolynomial(cp);
        			}
        			else
        				pDepNetReactionList.add(pdnr);
        		}
        		else if (rxntype.equals("PRODUCT")) {
        			if (pNum == 2) {
//        				if pDepNetReactionList contains this reaction simply add the chebyshev
            			// rate of this reaction to the one already present in the system.
            			if (pDepNetReactionList.contains(pdnr)){
            				int i = pDepNetReactionList.indexOf(pdnr);
            				PDepNetReaction pdnrAlreadyPresent = (PDepNetReaction)pDepNetReactionList.get(i);
            				pdnrAlreadyPresent.itsChebyshevPolynomials.addChebyshevPolynomial(cp);
            			}
            			else
            				pDepNetReactionList.add(pdnr);
        			}
        			else if (pNum == 1) {
        				if (includeAsIsomer(sp1)) pDepNetReactionList.add(pdnr);
        				else pDepNonincludedReactionList.add(pdnr);
        			}
        		} 
        		line = ChemParser.readMeaningfulLine(data);
        		
        	}
        	in.close();
        	
        	File f = new File("chemdis/fort.10");
        	File newFile = new File("chemdis/"+newName+"_input");
        	f.renameTo(newFile);
        	f = new File(chemdis_output);
        	newFile = new File("chemdis/"+newName+"_output");
        	f.renameTo(newFile);
        }
        catch (Exception e) {
        	System.out.println("Wrong output from chemdis!");
        	System.out.println(e.getMessage());
        	System.exit(0);
        }
        
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
    
    
    
    //## operation writePDepNetworkHeader(ReactionSystem) 
    private String writePDepNetworkHeader(ReactionSystem p_reactionSystem) {
        //#[ operation writePDepNetworkHeader(ReactionSystem) 
        String s = "RMG-Generated Partial Network " + String.valueOf(getID()) + '\n';
        s += "TRANGE\n 300 \t 2000 \t 10 \n";
        //double temp = p_reactionSystem.getPresentTemperature().getK();
        //s += "1\t" + Double.toString(temp) + '\n';
        s += "PRANGE\n 0.01 \t 100 \t 10 \n";
        s += "CHEBYSHEV \n 7 \t 4 \n";
        //double pres = p_reactionSystem.getPresentPressure().getAtm();
        //s += "1\t" + Double.toString(pres) + '\n';
        
        if (isChemAct) {
        	s += "CHEMACT\n";
        	s += "INPUT\n";
        	if (entryReaction == null) throw new NullPointerException();
        	Kinetics k;
        	if (entryReaction.isForward()) {
        		k = entryReaction.getKinetics();
        	}
        	else {
        		k = entryReaction.getFittedReverseKinetics();
        	}
        
        	if (k==null) throw new NullPointerException();
        
        	s += Double.toString(k.getAValue()) + '\t';
        	s += Double.toString(k.getNValue()) + '\t';
        	s += "0.0\t";
        	s += Double.toString(k.getEValue()) + '\n';
        }
        else {
        	s += "DISSOC\n";
        	s += "INPWONLY\n";
        	s += "SPC" + String.valueOf(((PDepWell)pDepWellList.getFirst()).getIsomer().getID()) + '\n';
        }
        
        s += "MASS\n";
        s += MathTool.formatDouble(getEntryMass(), 10, 2) + '\n';
        s += "PARAMETERS\n";
        LennardJones lj = ((PDepWell)pDepWellList.getFirst()).getIsomer().getLJ();
        s += MathTool.formatDouble(lj.getSigma(), 10, 2) + '\t' + MathTool.formatDouble(lj.getEpsilon(), 10, 2) + '\n';
        s += "DEDOWN\n";
        s += "INT\n";
        s += " 0.20\n";
        s += "BSGS\n";
        s += " 5.0\n";
        s += "XMG\n";
        
        double total = 0;
        HashMap colliders = p_reactionSystem.identifyColliders();
        for (Iterator iter = colliders.values().iterator(); iter.hasNext();) {
        	double conc = ((Double)iter.next()).doubleValue();
        	total += conc;
        }
        double dEdown = 0;
        for (Iterator iter = colliders.keySet().iterator(); iter.hasNext();) {
        	Object key = iter.next();
                s += "COLLIDER\n";
        	if (key instanceof Species) {
        		Species spe = (Species)key;
        		s += "!" + spe.getName() + '\n';
        		double conc = ((Double)colliders.get(spe)).doubleValue();
        		double mf = conc/total;
        		s += Double.toString(mf) + '\t' + Double.toString(spe.getMolecularWeight()) + '\t';
        		lj = spe.getLJ();
        		s += Double.toString(lj.getSigma()) + '\t' + Double.toString(lj.getEpsilon()) + '\t';
        		dEdown = spe.getDeltaEDown();
        		if (dEdown == 0) {
        			System.out.println("unknown colliders's dEdown: " + spe.getName());
        			System.exit(0);
        		}
        		s += Double.toString(dEdown) + '\n';
        	}
        	else if (key instanceof String) {
        		String name = (String)key;
        		double MW = 0.0;
        		if (name.equals("Ar") || name.equals("AR")) {
        			lj = new LennardJones();
        			dEdown = 374.0;
        			MW = 39.95;
        		}
        		else if (name.equals("N2")) {
        			lj = new LennardJones();
        			dEdown = 461.0;
        			MW = 28.01;
        		}
        		else if (name.equals("He") || name.equals("HE")) {
        			lj = new LennardJones();
        			dEdown = 291.0;
        			MW = 4.00;
        		}
        		else {
        			System.out.println("unknown colliders: " + name);
        			System.exit(0);
        		}
        		s += "!" + name + '\n';
        		double conc = ((Double)colliders.get(name)).doubleValue();
        		double mf = conc/total;
        		s += Double.toString(mf) + '\t' + Double.toString(MW) + '\t';
        		s += Double.toString(lj.getSigma()) + '\t' + Double.toString(lj.getEpsilon()) + '\t';
        		s += Double.toString(dEdown) + '\n';
        	}
        	else {
        			System.out.println("unknown colliders: " + key.toString());
        			System.exit(0);
        	}
        }
        
        return s;
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
    
    public LinkedList getProduct() {
        return product;
    }
    
    public LinkedList getReactant() {
        return reactant;
    }
    
    public static int getTotalNumber() {
        return totalNumber;
    }
    
    public ListIterator getPDepNetReactionList() {
        ListIterator iter=pDepNetReactionList.listIterator(0);
        return iter;
    }
    
    public ListIterator getPDepNonincludedReactionList() {
        ListIterator iter=pDepNonincludedReactionList.listIterator(0);
        return iter;
    }
    
    public ListIterator getPDepWellList() {
        ListIterator iter=pDepWellList.listIterator(0);
        return iter;
    }
    
    public boolean getAltered(){
		return altered;
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


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
    
    protected static HashMap dictionary = new HashMap();		//## attribute dictionary 
    
    protected Reaction entryReaction = null;		//## attribute entryReaction 
    
    protected boolean isChemAct;		//## attribute isChemAct 
    
    protected double kLeak;		//## attribute kLeak 
    
    protected LinkedList product = null;		//## attribute product 
    
    protected LinkedList reactant = null;		//## attribute reactant 
    
    protected static int totalNumber = 0;		//## attribute totalNumber 
    
    protected LinkedList pDepNetReactionList;
    protected LinkedList pDepNonincludedReactionList;
    protected LinkedList pDepWellList;
    
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
    protected void addPDepWell(PDepWell p_newWell) {
        //#[ operation addPDepWell(PDepWell) 
        if (pDepWellList.contains(p_newWell)) return;
        
        pDepWellList.add(p_newWell);
        
        decideWellPathType();
        
        return;
        
        
        
        //#]
    }
    
    //## operation decideWellPathType() 
    public void decideWellPathType() {
        //#[ operation decideWellPathType() 
        for (Iterator iter = getPDepWellList(); iter.hasNext(); ) {
        	PDepWell pdw = (PDepWell)iter.next();
        	for (Iterator pathIter = pdw.getPaths(); pathIter.hasNext(); ) {
        		PDepPathReaction pdpr = (PDepPathReaction)pathIter.next();
        		int pNum = pdpr.getProductNumber();
        		LinkedList r = pdpr.getReactantList();
        		LinkedList p = pdpr.getProductList();
        
        		if (pNum == 1) {
        			Species spe = ((ChemGraph)p.iterator().next()).getSpecies();
        			if (includeAsIsomer(spe)) pdpr.setTypeAsIsomer();
        			else pdpr.setTypeAsNonIncluded();
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
        	ChemGraph cg = (ChemGraph)iter.next();
        	mass += cg.getMolecularWeight();
        }
        
        return mass;
        //#]
    }
    
    //## operation includeAsIsomer(Species) 
    public boolean includeAsIsomer(Species p_species) {
        //#[ operation includeAsIsomer(Species) 
        for (Iterator iter = getPDepWellList(); iter.hasNext(); ) {
        	PDepWell pdw = (PDepWell)iter.next();
        	if ((pdw.getIsomer()).equals(p_species)) return true;
        }
        return false;
        //#]
    }
    
    //## operation initializeKLeak() 
    public void initializeKLeak() {
        //#[ operation initializeKLeak() 
        kLeak = 0;
        Temperature t = new Temperature(715,"K");
        if (isChemAct) {
        	kLeak = entryReaction.calculateRate(t);
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
        		kLeak += pdpr.getTemplateReaction().calculateRate(t);
        	}
        }
        
        //System.out.println("present k leak = " + String.valueOf(kLeak));
        //#]
    }
    
    //## operation initializePDepNetwork(TemplateReaction) 
    public void initializePDepNetwork(TemplateReaction p_entryReaction) {
        //#[ operation initializePDepNetwork(TemplateReaction) 
        reactant = (LinkedList)p_entryReaction.getReactantList().clone();
        if (reactant.size() == 1) {
        	isChemAct = false;
        	entryReaction = null;
        	ChemGraph cg = (ChemGraph)p_entryReaction.getReactants().next();
         	Species spe = cg.getSpecies();
         	PDepWell pdw = new PDepWell(spe);
        
         	// add a new PDepWell into PDepWell list
         	addPDepWell(pdw);
         	product = null;                      
        }
        else if (reactant.size() == 2) {
        	isChemAct = true;
        	entryReaction = p_entryReaction;
        	product = (LinkedList)p_entryReaction.getProductList().clone();
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
    public static PDepNetwork makePDepNetwork(TemplateReaction p_entryReaction) {
        //#[ operation makePDepNetwork(TemplateReaction) 
	return null;
	/*Object key = null;
        if (p_entryReaction.getReactantNumber() == 2) {
        	Species pro = ((ChemGraph)p_entryReaction.getProducts().next()).getSpecies();
        	if (pro.isTriatomicOrSmaller()) return null;
        	key = p_entryReaction.getStructure();
        }
        else {
        	key = ((ChemGraph)p_entryReaction.getReactants().next()).getSpecies();
        	if (((Species)key).isTriatomicOrSmaller()) return null;
        }
        
        Object obj = dictionary.get(key);
        if (obj == null) {
        	// if it is a->b, but a doesn't have three frequency model, return null
        	if (p_entryReaction.getReactantNumber() == 1) {
        		ChemGraph cg = (ChemGraph)p_entryReaction.getReactants().next();
         		Species spe = cg.getSpecies();
         		if (!spe.hasThreeFrequencyModel()) {
         			return null;
         		}
           	}
        	//System.out.println("construct a new PDepNetwork " + key.toString());
        	PDepNetwork pnw = new PDepNetwork();
        	dictionary.put(key, pnw);
        	pnw.initializePDepNetwork(p_entryReaction);
        	pnw.initializeKLeak();
        	return pnw;
        
        }
        else return (PDepNetwork)obj;
        
      */  
        
        //#]
    }
    
    //## operation runPDepCalculation(ReactionSystem) 
    public void runPDepCalculation(ReactionSystem p_reactionSystem) {
        //#[ operation runPDepCalculation(ReactionSystem) 
        if (!isActive() && getIsChemAct()) {
        	Temperature t = p_reactionSystem.getPresentTemperature();
        	kLeak = getEntryReaction().calculateRate(t);
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
        
        // parse the output file from chemdis
        try {
			//String dir = System.getProperty("RMG.workingDirectory");
        	String chemdis_output = "chemdis/chemdis-xmg.out";
        
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
        	reactant.add(sr1.getChemGraph());
        	if (rNum == 2) {
        		temp = st.nextToken();
        		String r2 = st.nextToken().trim();
        		r2 = r2.substring(3, r2.length());
        		int idr2 = Integer.parseInt(r2);
        		Species sr2 = SpeciesDictionary.getInstance().getSpeciesFromID(idr2); 
        		newName += "+" + sr2.getName()+"("+String.valueOf(sr2.getID())+")";
        		reactant.add(sr2.getChemGraph());
        	}
        	
        	pDepNetReactionList.clear();
        	pDepNonincludedReactionList.clear();
        	line = ChemParser.readMeaningfulLine(data);
        	while (!line.startsWith("END")) {
        		line = line.trim();
        		LinkedList product = new LinkedList();
        		st = new StringTokenizer(line);
        		String temperature = st.nextToken();
        		double doubleT = Double.parseDouble(temperature); 
        		String pressure = st.nextToken();                 
        		double doubleP = Double.parseDouble(pressure); 
        		String rxntype = st.nextToken();
        		int pNum = 1;
        		String p1 = st.nextToken().trim();
        		p1 = p1.substring(3, p1.length());
        		int idp1 = Integer.parseInt(p1);
        		Species sp1 = SpeciesDictionary.getInstance().getSpeciesFromID(idp1);
        	    product.add(sp1.getChemGraph());
        		String next = st.nextToken();
        		if ((next.trim()).equals("+")) {
        			String p2 = st.nextToken().trim();
        			p2 = p2.substring(3, p2.length());
        			int idp2 = Integer.parseInt(p2);
        			Species sp2 = SpeciesDictionary.getInstance().getSpeciesFromID(idp2);
        			product.add(sp2.getChemGraph());
        			next = st.nextToken();
        			pNum++;
        		}
        		double k = Double.parseDouble(next);
        		PDepNetReaction pdnr = new PDepNetReaction(reactant, product, k, p_reactionSystem.getPresentTemperature().getK() ,p_reactionSystem.getPresentPressure().getAtm());
        		if (rxntype.equals("ISOMER")) pDepNetReactionList.add(pdnr);
        		else if (rxntype.equals("PRODUCT")) {
        			if (pNum == 2) pDepNetReactionList.add(pdnr);
        			else if (pNum == 1) {
        				if (includeAsIsomer(sp1)) pDepNetReactionList.add(pdnr);
        				else pDepNonincludedReactionList.add(pdnr);
        			}
        		} 
        		line = ChemParser.readMeaningfulLine(data);
        		
        	}
        	in.close();
        	updateKLeak();
        	File f = new File("chemdis/fort.10");
        	File newFile = new File("chemdis/"+newName+"_input");
        	f.renameTo(newFile);
        	f = new File(chemdis_output);
        	newFile = new File("chemdis/"+newName+"_output");
        	f.renameTo(newFile);
        }
        catch (Exception e) {
        	System.out.println("Wrong output from chemdis!");
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
        	Species spe = ((ChemGraph)reactant.getFirst()).getSpecies();
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
        PDepWell pdw = new PDepWell(p_nextIsomer);
        //System.out.println("begin to add new well to system");
        addPDepWell(pdw);
        //System.out.println("finish adding new well to system");
        return;
        //#]
    }
    
    //## operation updateKLeak() 
    public void updateKLeak() {
        //#[ operation updateKLeak() 
        kLeak = 0;
        for (Iterator iter = getPDepNonincludedReactionList(); iter.hasNext(); ) {
        	PDepNetReaction pdnr = (PDepNetReaction)iter.next();
        	double rate = pdnr.getRate();
        	if (rate < 0) throw new InvalidPDepNetReactionException();
        	kLeak += rate;
        }
        return;
        //#]
    }
    
    //## operation writePDepNetworkHeader(ReactionSystem) 
    public String writePDepNetworkHeader(ReactionSystem p_reactionSystem) {
        //#[ operation writePDepNetworkHeader(ReactionSystem) 
        String s = "RMG-Generated Partial Network " + String.valueOf(getID()) + '\n';
        s += "TEMP\n";
        double temp = p_reactionSystem.getPresentTemperature().getK();
        s += "1\t" + Double.toString(temp) + '\n';
        s += "PRES\n";
        double pres = p_reactionSystem.getPresentPressure().getAtm();
        s += "1\t" + Double.toString(pres) + '\n';
        
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
    
    public double getKLeak() {
        return kLeak;
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
    
    public ListIterator getPDepNetReactionListEnd() {
        return pDepNetReactionList.listIterator(pDepNetReactionList.lastIndexOf(pDepNetReactionList.getLast()));
    }
    
    public void deletePDepNetReactionList(PDepNetReaction p_PDepNetReaction) {
        pDepNetReactionList.remove(p_PDepNetReaction);
        p_PDepNetReaction=null;
    }
    
    public ListIterator getPDepNonincludedReactionList() {
        ListIterator iter=pDepNonincludedReactionList.listIterator(0);
        return iter;
    }
    
    public ListIterator getPDepNonincludedReactionListEnd() {
        return pDepNonincludedReactionList.listIterator(pDepNonincludedReactionList.lastIndexOf(pDepNonincludedReactionList.getLast()));
    }
    
    public PDepNetReaction newPDepNonincludedReactionList() {
        PDepNetReaction newPDepNetReaction = new PDepNetReaction();
        pDepNonincludedReactionList.add(newPDepNetReaction);
        return newPDepNetReaction;
    }
    
    public void deletePDepNonincludedReactionList(PDepNetReaction p_PDepNetReaction) {
        pDepNonincludedReactionList.remove(p_PDepNetReaction);
        p_PDepNetReaction=null;
    }
    
    public ListIterator getPDepWellList() {
        ListIterator iter=pDepWellList.listIterator(0);
        return iter;
    }
    
    public ListIterator getPDepWellListEnd() {
        return pDepWellList.listIterator(pDepWellList.lastIndexOf(pDepWellList.getLast()));
    }
    
    public void deletePDepWellList(PDepWell p_PDepWell) {
        pDepWellList.remove(p_PDepWell);
        p_PDepWell=null;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\PDepNetwork.java
*********************************************************************/


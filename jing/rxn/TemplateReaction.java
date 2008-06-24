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


import jing.chem.*;
import jing.chemUtil.Graph;
import jing.chemUtil.Node;

import java.util.*;
import jing.param.*;
import jing.rxnSys.SystemSnapshot;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\TemplateReaction.java                                                                  
//----------------------------------------------------------------------------

/**
Reaction generated from templates.
Immutable object.
*/
//## class TemplateReaction 
public class TemplateReaction extends Reaction {
    
    protected PDepNetwork pDepNetwork;
    protected ReactionTemplate reactionTemplate;
    
    // Constructors
    
    //## operation TemplateReaction(Structure,RateConstant,ReactionTemplate) 
    private  TemplateReaction(Structure p_structure, Kinetics p_kinetics, ReactionTemplate p_template) {
        //#[ operation TemplateReaction(Structure,RateConstant,ReactionTemplate) 
        structure = p_structure;
        kinetics = p_kinetics;
        reactionTemplate = p_template;
		
        //#]
    }
    public  TemplateReaction() {
    }
    
    //## operation calculatePDepRate(Temperature) 
   /* public double calculatePDepRate(Temperature p_temperature) {
        //#[ operation calculatePDepRate(Temperature) 
        PDepNetwork pdn = getPDepNetwork();
        if (pdn != null) {
        	
        	Iterator iter = pdn.getPDepNetReactionList(); 
        	while (iter.hasNext()) {
        		PDepNetReaction pdnr = (PDepNetReaction)iter.next();
        		if (pdnr.getStructure().equals(getStructure())) {
        			double temp = pdnr.getTemperature();
        			if (temp != p_temperature.getK()) {
        				System.out.println("Different temperature used!");
        				System.exit(0);
        			}
        			//System.out.println("for reaction " + toString() + "\t using p dep rate:" + String.valueOf(pdnr.getRate()));
        			return pdnr.getRate();
        		}
        	}
        	iter = pdn.getPDepNonincludedReactionList(); 
        	while (iter.hasNext()) {
        		PDepNetReaction pdnr = (PDepNetReaction)iter.next();
        		if (pdnr.getStructure().equals(getStructure())) {
        			double temp = pdnr.getTemperature();
        			if (temp != p_temperature.getK()) {
        				System.out.println("Different temperature used!");
        				System.exit(0);
        			}
        			//System.out.println("for reaction " + toString() + "\t using p dep rate:" + String.valueOf(pdnr.getRate()));
        			return pdnr.getRate();
        		}
        	}
        }
        
        return calculateRate(p_temperature);
        //#]
    }*/
    
//  ## operation calculatePDepRate(Temperature) 
    public double calculateTotalPDepRate(Temperature p_temperature, Pressure p_pressure) {
        //#[ operation calculatePDepRate(Temperature) 
        PDepNetwork pdn = getPDepNetwork();
        if (pdn != null) {
        	
        	Iterator iter = pdn.getPDepNetReactionList(); 
        	while (iter.hasNext()) {
        		PDepNetReaction pdnr = (PDepNetReaction)iter.next();
        		if (pdnr.getStructure().equals(getStructure())) {
        			/*double temp = pdnr.getTemperature();
        			if (temp != p_temperature.getK()) {
        				System.out.println("Different temperature used!");
        				System.exit(0);
        			}*/
        			//System.out.println("for reaction " + toString() + "\t using p dep rate:" + String.valueOf(pdnr.getRate()));
        			//10/25/07 gmagoon: updated to use calculateRate with system snapshot (to avoid use of Global.temperature and Global.pressure)
                                SystemSnapshot currentTPSnapshot = new SystemSnapshot();//10/25/07 gmagoon: make currentTPsnapshot variable, which will be used to pass temperature and pressure to calculateRate
                                currentTPSnapshot.setTemperature(p_temperature);
                                currentTPSnapshot.setPressure(p_pressure);
                                double rate = pdnr.calculateRate(currentTPSnapshot);
                                currentTPSnapshot = null;
                                return rate;
        		}
        	}
        	iter = pdn.getPDepNonincludedReactionList(); 
        	while (iter.hasNext()) {
        		PDepNetReaction pdnr = (PDepNetReaction)iter.next();
        		if (pdnr.getStructure().equals(getStructure())) {
        			/*double temp = pdnr.getTemperature();
        			if (temp != p_temperature.getK()) {
        				System.out.println("Different temperature used!");
        				System.exit(0);
        			}*/
        			//System.out.println("for reaction " + toString() + "\t using p dep rate:" + String.valueOf(pdnr.getRate()));
        			SystemSnapshot currentTPSnapshot = new SystemSnapshot();//10/25/07 gmagoon: make currentTPsnapshot variable, which will be used to pass temperature and pressure to calculateRate
                                currentTPSnapshot.setTemperature(p_temperature);
                                currentTPSnapshot.setPressure(p_pressure);
                                double rate =  pdnr.calculateRate(currentTPSnapshot);
                                currentTPSnapshot = null;
                                return rate;
        		}
        	}
        }
        
        return calculateTotalRate(p_temperature);
        //#]
    }
    
    //## operation generateReverseForBackwardReaction() 
    private TemplateReaction generateReverseForBackwardReaction(Structure fs, Structure fsSp) {
        //#[ operation generateReverseForBackwardReaction() 
        // we need to only generate reverse reaction for backward reaction, so that we wont be stuck into a self loop.
        if (!this.isBackward()) return null;
        
        ReactionTemplate fRT = getReactionTemplate();
        ReactionTemplate rRT = null;
        
        if (fRT.isForward()) return null;
        else if (fRT.isNeutral()) rRT = fRT;
        else if (fRT.isBackward()) rRT = fRT.getReverseReactionTemplate();
        else throw new InvalidReactionTemplateDirectionException();
        
        //Structure fs = getStructure();
        LinkedList freactant = fs.getReactantList();
        LinkedList fproduct = fs.getProductList();
        
        Structure rs = new Structure(fproduct, freactant, -1*this.getDirection());
		Structure rsSp = new Structure(fsSp.products, fsSp.reactants, -1*this.getDirection());
		TemplateReaction rr = rRT.getReactionFromStructure(rsSp);
		if (rr!= null) {
			rr.setReverseReaction(this);
			return rr;
		}
		int rNum = fproduct.size();
		Kinetics k = rRT.findReverseRateConstant(rs);
		if (k == null && rRT.name.equals("R_Recombination")) {
			
			ChemGraph cg = ((ChemGraph)fproduct.get(0));
			Graph g = cg.getGraph();
			Node n = (Node)g.getCentralNodeAt(2);
			if (n == null){
				cg = ((ChemGraph)fproduct.get(1));
				g = cg.getGraph();
				n = (Node)g.getCentralNodeAt(2);
			}
			g.clearCentralNode();
			g.setCentralNode(1,n);
			k = rRT.findRateConstant(rs);	
		}
		else if (k==null && rRT.name.equals("H_Abstraction")) {
			ChemGraph cg1 = ((ChemGraph)fproduct.get(0));
			Graph g1 = cg1.getGraph();
			Node n3 = (Node)g1.getCentralNodeAt(3);
			if (n3 == null){
				cg1 = ((ChemGraph)fproduct.get(1));
				g1 = cg1.getGraph();
				n3 = (Node)g1.getCentralNodeAt(3);
				Node n2 = (Node)g1.getCentralNodeAt(2);
				g1.clearCentralNode();
				g1.setCentralNode(1,n3);
				g1.setCentralNode(2,n2);
				ChemGraph cg2 = ((ChemGraph)fproduct.get(0));
				Graph g2 = cg2.getGraph();
				Node n1 = (Node)g2.getCentralNodeAt(1);
				g2.clearCentralNode();
				g2.setCentralNode(3,n1);
			}
			else {
				Node n2 = (Node)g1.getCentralNodeAt(2);
				g1.clearCentralNode();
				g1.setCentralNode(1,n3);
				g1.setCentralNode(2,n2);
				ChemGraph cg2 = ((ChemGraph)fproduct.get(1));
				Graph g2 = cg2.getGraph();
				Node n1 = (Node)g2.getCentralNodeAt(1);
				g2.clearCentralNode();
				g2.setCentralNode(3,n1);
			}
			k = rRT.findRateConstant(rs);
		}
		else if (k==null && rRT.name.equals("intra_H_migration")) {
			ChemGraph cg = ((ChemGraph)fproduct.get(0));
			rr = rRT.calculateForwardRateConstant(cg,rs);
			if (!rr.isForward()) {
				String err = "Backward:" + structure.toString() + String.valueOf(structure.calculateKeq(new Temperature(298,"K"))) + '\n';
				err = err + "Forward:" + rr.structure.toString()+ String.valueOf(rr.structure.calculateKeq(new Temperature(298,"K"))) ;
				 
				throw new InvalidReactionDirectionException(err);
			}
			rr.setReverseReaction(this);
			rRT.addReaction(rr);
			return rr;
			
		}
		if (k==null){
			System.err.println("Couldn't find the rate constant for reaction: "+rs.toChemkinString(true)+" with "+rRT.name);
			//System.exit(0);
			return null;
		}
			rr = new TemplateReaction(rsSp,k,rRT);
			if (!(rr.getReactantNumber() ==2 && rr.getProductNumber() == 2))
				rr.pDepNetwork = PDepNetwork.makePDepNetwork(rr);

		if (!rr.isForward()) {
			String err = "Backward:" + structure.toString() + String.valueOf(structure.calculateKeq(new Temperature(298,"K"))) + '\n';
			err = err + "Forward:" + rr.structure.toString()+ String.valueOf(rr.structure.calculateKeq(new Temperature(298,"K"))) ;
			 
			throw new InvalidReactionDirectionException(err);
		}
		rr.setReverseReaction(this);
		rRT.addReaction(rr);
		return rr;
		
 
        	
        //#]
    }
    
    /**
    Requires:
    Effects: return the type of itsReactionTemplate as the type of this reaction.
    Modifies:
    */
    //## operation getType() 
    public String getType() {
        //#[ operation getType() 
        return reactionTemplate.getName();
        
        
        //#]
    }
    
    //## operation makeTemplateReaction(Structure,RateConstant,ReactionTemplate) 
    public static TemplateReaction makeTemplateReaction(Structure p_structureSp, Kinetics p_kinetics, ReactionTemplate p_template, Structure p_structure) {
        //#[ operation makeTemplateReaction(Structure,RateConstant,ReactionTemplate) 
        double PT = System.currentTimeMillis();
		TemplateReaction reaction = p_template.getReactionFromStructure(p_structureSp);
		Global.getReacFromStruc = Global.getReacFromStruc + (System.currentTimeMillis() - PT)/1000/60;
			
        if (reaction == null) {
        	reaction = new TemplateReaction(p_structureSp,p_kinetics,p_template);
        	if (reaction.isBackward()) {
				double pt = System.currentTimeMillis();
        		TemplateReaction reverse = reaction.generateReverseForBackwardReaction(p_structure, p_structureSp);
        		if (reverse == null)
        			return null;
        		reaction.setReverseReaction(reverse);
        		
				//Global.generateReverse = Global.generateReverse + (System.currentTimeMillis() - pt)/1000/60;
        	}
        	else {
        		ReactionTemplate fRT = reaction.getReactionTemplate();
                ReactionTemplate rRT = null;
                
                
                if (fRT.isNeutral()) rRT = fRT;
                else rRT = fRT.getReverseReactionTemplate();
                
                if (rRT != null){
                	TemplateReaction reverse = new TemplateReaction(p_structureSp.generateReverseStructure(),p_kinetics,rRT);
                	reaction.setReverseReaction(reverse);
            		reverse.setReverseReaction(reaction);
            		rRT.addReaction(reverse);
                	if (!(reverse.getReactantNumber() == 2 && reverse.getProductNumber() ==2))
                		reverse.pDepNetwork = PDepNetwork.makePDepNetwork(reverse);
            		
                }
        		
        	}
        	p_template.addReaction(reaction);    
        	if (!(reaction.getReactantNumber() == 2 && reaction.getProductNumber() == 2) ) {
        		reaction.pDepNetwork = PDepNetwork.makePDepNetwork(reaction);     
        	}
        	else {
        		reaction.pDepNetwork = null;
        	}
        	if (!reaction.repOk()) throw new InvalidTemplateReactionException();
        }
        else {
        	// the same identity of p_structure already exists, clear p_structure to release mem
        	Structure st = reaction.getStructure();
        	if (st!=p_structure) p_structure = null;
        }
		Global.makeTR += (System.currentTimeMillis()-PT)/1000/60;
        return reaction;
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return (super.repOk() && reactionTemplate.repOk());
        //#]
    }
    
   
    
    //## operation toString() 
    public String toString(Temperature p_temperature) {
        //#[ operation toString() 
        String s = getStructure().toString()  + '\t';
        Kinetics k = getKinetics();
        String kString = k.toChemkinString(calculateHrxn(p_temperature), p_temperature, false);
       
        return s + kString;
        
        //#]
    }
    
	   //## operation toStringWithReveseReaction() 
    public String toStringWithReveseReaction(Temperature p_temperature) {
        //#[ operation toStringWithReveseReaction() 
        TemplateReaction rr = (TemplateReaction)getReverseReaction();
        if (rr == null) return getStructure().toChemkinString(false).toString() + '\t' + getReactionTemplate().getName() + '\t' + getKinetics().toChemkinString(calculateHrxn(p_temperature), p_temperature, true);
        else {
        	TemplateReaction temp = null;
        	if (isForward()) temp = this;
        	else if (isBackward()) temp = rr;
        	else throw new InvalidReactionDirectionException();
        	
        	return temp.getStructure().toChemkinString(false).toString() + '\t' + temp.getReactionTemplate().getName() + '\t' + temp.getKinetics().toChemkinString(calculateHrxn(p_temperature), p_temperature, true);
        }
        
        //#]
    }
    
    
    public PDepNetwork getPDepNetwork() {
        return pDepNetwork;
    }
    
      public double getRateConstant(Temperature p_temperature, Pressure p_pressure){
	//public double getRateConstant(){
		if (rateConstant == 0)
                        rateConstant = calculateTotalPDepRate(p_temperature, p_pressure);
			//rateConstant = calculateTotalPDepRate(Global.temperature);
		return rateConstant;
	}
	
    public ReactionTemplate getReactionTemplate() {
        return reactionTemplate;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\TemplateReaction.java
*********************************************************************/


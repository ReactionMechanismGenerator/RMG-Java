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
    private  TemplateReaction(Structure p_structure, Kinetics[] p_kinetics, ReactionTemplate p_template) {
        //#[ operation TemplateReaction(Structure,RateConstant,ReactionTemplate) 
        structure = p_structure;
        kinetics = p_kinetics;
        reactionTemplate = p_template;
        if (kinetics != null)	
        	kineticsFromPrimaryKineticLibrary = p_kinetics[0].getFromPrimaryKineticLibrary();
		
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
        	
        	ListIterator iter = pdn.getNetReactions().listIterator(); 
        	while (iter.hasNext()) {
        		PDepReaction pdnr = (PDepReaction) iter.next();
        		if (pdnr.getStructure().equals(getStructure()))
        			return pdnr.calculateRate(p_temperature, p_pressure);
        	}
        	iter = pdn.getNonincludedReactions().listIterator(); 
        	while (iter.hasNext()) {
        		PDepReaction pdnr = (PDepReaction) iter.next();
        		if (pdnr.getStructure().equals(getStructure()))
        			return pdnr.calculateRate(p_temperature, p_pressure);
        	}
        }
        
        return calculateTotalRate(p_temperature);
        //#]
    }
    
    //## operation generateReverseForBackwardReaction() 
	private TemplateReaction generateReverseForBackwardReaction(Structure fs, Structure fsSp) {
		//#[ operation generateReverseForBackwardReaction() 
		// we need to only generate reverse reaction for backward reaction, so that we wont be stuck into a self loop.
		if (!this.isBackward()) {
			return null;
		}
		ReactionTemplate fRT = getReactionTemplate();
		ReactionTemplate rRT = null;

		if (fRT.isForward()) {
			return null;
		} else if (fRT.isNeutral()) {
			rRT = fRT;
		} else if (fRT.isBackward()) {
			rRT = fRT.getReverseReactionTemplate();
		} else {
			throw new InvalidReactionTemplateDirectionException();		//Structure fs = getStructure();
		}
		LinkedList freactant = fs.getReactantList();
		LinkedList fproduct = fs.getProductList();

		Structure rs = new Structure(fproduct, freactant, -1 * this.getDirection());
		Structure rsSp = new Structure(fsSp.products, fsSp.reactants, -1 * this.getDirection());
		TemplateReaction rr = rRT.getReactionFromStructure(rsSp);
		if (rr != null) {
			rr.setReverseReaction(this);
			return rr;
		}
		int rNum = fproduct.size();
		Kinetics[] k = rRT.findReverseRateConstant(rs);
		if (k == null && rRT.name.equals("R_Recombination")) {

			ChemGraph cg = ((ChemGraph) fproduct.get(0));
			Graph g = cg.getGraph();
			Node n = (Node) g.getCentralNodeAt(2);
			if (n == null) {
				cg = ((ChemGraph) fproduct.get(1));
				g = cg.getGraph();
				n = (Node) g.getCentralNodeAt(2);
			}
			g.clearCentralNode();
			g.setCentralNode(1, n);
			k = rRT.findRateConstant(rs);
		} else if (k == null && rRT.name.equals("H_Abstraction")) {
			ChemGraph cg1 = ((ChemGraph) fproduct.get(0));
			Graph g1 = cg1.getGraph();
			Node n3 = (Node) g1.getCentralNodeAt(3);
			if (n3 == null) {
				cg1 = ((ChemGraph) fproduct.get(1));
				g1 = cg1.getGraph();
				n3 = (Node) g1.getCentralNodeAt(3);
				Node n2 = (Node) g1.getCentralNodeAt(2);
				g1.clearCentralNode();
				g1.setCentralNode(1, n3);
				g1.setCentralNode(2, n2);
				ChemGraph cg2 = ((ChemGraph) fproduct.get(0));
				Graph g2 = cg2.getGraph();
				Node n1 = (Node) g2.getCentralNodeAt(1);
				g2.clearCentralNode();
				g2.setCentralNode(3, n1);
			} else {
				Node n2 = (Node) g1.getCentralNodeAt(2);
				g1.clearCentralNode();
				g1.setCentralNode(1, n3);
				g1.setCentralNode(2, n2);
				ChemGraph cg2 = ((ChemGraph) fproduct.get(1));
				Graph g2 = cg2.getGraph();
				Node n1 = (Node) g2.getCentralNodeAt(1);
				g2.clearCentralNode();
				g2.setCentralNode(3, n1);
			}
			k = rRT.findRateConstant(rs);
		} 
		
		/*
		 * Added by MRH on 27-Aug-2009
		 * 	This hard-coding is necessary for rxn family templates
		 * 		that are labeled "thermo_consistence".
		 * 
		 * After the chemgraphs are mutated, the central nodes for 
		 * 	the products are not correct (see example below).  These
		 * 	hard-coded portions are necessary for RMG to find Kinetics
		 * 	for the structure.
		 * 
		 * Example: CH4 + H
		 * 	CH4
		 * 	1 *1 C 0 {2,S} {3,S} {4,S} {5,S}
		 * 	2 *2 H 0 {1,S}
		 * 	3    H 0 {1,S}
		 * 	4    H 0 {1,S}
		 * 	5    H 0 {1,S}
		 * 
		 * 	H
		 * 	1 *3 H 1
		 * 
		 * After RMG has "reactChemGraph" and "mutate" the chemgraphs
		 * 	of the reactants, the products would look as such:
		 * 
		 * prod1
		 *  1 *1 C 1 {2,S} {3,S} {4,S}
		 * 	2    H 0 {1,S}
		 * 	3    H 0 {1,S}
		 * 	4    H 0 {1,S}
		 * 
		 * prod2
		 *  1 *3 H 0 {2,S}
		 *  2 *2 H 0 {1,S}
		 *  
		 *  Assuming the reaction as written (CH4+H=CH3+H2) is
		 *  	endothermic at 298K, RMG will label this structure
		 *  	as direction=-1 (backward).  When attempting to find
		 *  	Kinetics for the backward reaction, RMG will try
		 *  	to match the prod1 graph against the generic graphs
		 *  	X_H and Y_rad_birad.  It cannot match Y_rad_birad
		 *  	(because there is no *3 node) and it cannot match
		 *  	X_H (because there is no *2 node).  Thus, a "null"
		 *  	Kinetics will be returned from the findReverseRateConstant
		 *  	call.  We then relabel the central nodes on prod1
		 *  	and prod2 and attempt to get Kinetics for this
		 *  	structure.
		 *  
		 *  I am adding the following bit of code to work with the
		 *  	new reaction family Aaron Vandeputte is adding to
		 *  	RMG: "".
		 * 
		 */
		
		else if (k == null && rRT.name.equals("intra_substitutionS_isomerization")) {
			ChemGraph cg1 = ((ChemGraph) fproduct.get(0));
			Graph g1 = cg1.getGraph();
			Node n1 = (Node) g1.getCentralNodeAt(1);
			Node n2 = (Node) g1.getCentralNodeAt(2);
			Node n3 = (Node) g1.getCentralNodeAt(3);
			Node n4 = (Node) g1.getCentralNodeAt(4);
			Node n5 = (Node) g1.getCentralNodeAt(5);
			Node n6 = (Node) g1.getCentralNodeAt(6);
			Node n7 = (Node) g1.getCentralNodeAt(7);
			g1.clearCentralNode();
			g1.setCentralNode(1, n1);
			g1.setCentralNode(2, n3);
			g1.setCentralNode(3, n2);
			if (n7 != null)	{
				g1.setCentralNode(7, n4);
				g1.setCentralNode(6, n5);
				g1.setCentralNode(5, n6);
				g1.setCentralNode(4, n7);
			}
			else if (n6 != null) {
				g1.setCentralNode(6, n4);
				g1.setCentralNode(5, n5);
				g1.setCentralNode(4, n6);
			}
			else if (n5 != null) {
				g1.setCentralNode(5, n4);
				g1.setCentralNode(4, n5);
			}
			else if (n4 != null)
				g1.setCentralNode(4, n4);
			k = rRT.findRateConstant(rs);
		}
		
		else if (k == null && rRT.name.equals("intra_H_migration")) {
			ChemGraph cg = ((ChemGraph) fproduct.get(0));
			rr = rRT.calculateForwardRateConstant(cg, rs);
			if (!rr.isForward()) {
				String err = "Backward:" + structure.toString() + String.valueOf(structure.calculateKeq(new Temperature(298, "K"))) + '\n';
				err = err + "Forward:" + rr.structure.toString() + String.valueOf(rr.structure.calculateKeq(new Temperature(298, "K")));

				throw new InvalidReactionDirectionException(err);
			}
			rr.setReverseReaction(this);
			rRT.addReaction(rr);
			return rr;

		}
		if (k == null) {
			System.err.println("Couldn't find the rate constant for reaction: " + rs.toChemkinString(true) + " with " + rRT.name);
			//System.exit(0);
			return null;
		}
		rr = new TemplateReaction(rsSp, k, rRT);
		if (!rr.isForward()) {
			String err = "Backward:" + structure.toString() + String.valueOf(structure.calculateKeq(new Temperature(298, "K"))) + '\n';
			err = err + "Forward:" + rr.structure.toString() + String.valueOf(rr.structure.calculateKeq(new Temperature(298, "K")));

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
    
    public static TemplateReaction makeTemplateReaction(Structure p_structureSp,
			Kinetics[] p_kinetics, ReactionTemplate p_template, Structure p_structure) {
		
		double PT = System.currentTimeMillis();
		TemplateReaction reaction = p_template.getReactionFromStructure(p_structureSp);
		Global.getReacFromStruc = Global.getReacFromStruc + (System.currentTimeMillis() - PT) / 1000 / 60;

		if (reaction == null) {
			
			reaction = new TemplateReaction(p_structureSp, p_kinetics, p_template);
			// DEBUG: Tell console I made this reaction
			System.out.println("Created new reaction: " + reaction.toString());

			if (reaction.isBackward()) {
				
				TemplateReaction reverse = reaction.generateReverseForBackwardReaction(p_structure, p_structureSp);
				if (reverse == null)
					return null;
				reaction.setReverseReaction(reverse);
			}
			else {
				
				ReactionTemplate fRT = reaction.getReactionTemplate();
				ReactionTemplate rRT = null;

				if (fRT.isNeutral())
					rRT = fRT;
				else
					rRT = fRT.getReverseReactionTemplate();
				
				if (rRT != null) {
					TemplateReaction reverse = new TemplateReaction(p_structureSp.generateReverseStructure(), p_kinetics, rRT);
					reaction.setReverseReaction(reverse);
					reverse.setReverseReaction(reaction);
					rRT.addReaction(reverse);
				}

			}
			p_template.addReaction(reaction);
			if (!reaction.repOk()) {
				throw new InvalidTemplateReactionException();
			}
		} else {
			// the same identity of p_structure already exists, clear p_structure to release mem
			Structure st = reaction.getStructure();
			if (st != p_structure) {
				p_structure = null;
			}
		}
		Global.makeTR += (System.currentTimeMillis() - PT) / 1000 / 60;
		return reaction;

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
    	String totalString = "";
        String s = getStructure().toString()  + '\t';
        Kinetics[] k = getKinetics();
        for (int i=0; i<k.length; i++) {
        	totalString += s + k[i].toChemkinString(calculateHrxn(p_temperature), p_temperature, false);
        	if (k.length > 1) totalString += "\n";
        }
        return totalString;
        
        //#]
    }
    
    /*
     * MRH 24MAR2010:
     * 	Commented out toStringWithReverseReaction method as it is never called
     */
	   //## operation toStringWithReveseReaction() 
//    public String toStringWithReveseReaction(Temperature p_temperature) {
//        //#[ operation toStringWithReveseReaction() 
//        TemplateReaction rr = (TemplateReaction)getReverseReaction();
//        if (rr == null) return getStructure().toChemkinString(false).toString() + '\t' + getReactionTemplate().getName() + '\t' + getKinetics().toChemkinString(calculateHrxn(p_temperature), p_temperature, true);
//        else {
//        	TemplateReaction temp = null;
//        	if (isForward()) temp = this;
//        	else if (isBackward()) temp = rr;
//        	else throw new InvalidReactionDirectionException();
//        	
//        	return temp.getStructure().toChemkinString(false).toString() + '\t' + temp.getReactionTemplate().getName() + '\t' + temp.getKinetics().toChemkinString(calculateHrxn(p_temperature), p_temperature, true);
//        }
//        
//        //#]
//    }
    
    
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


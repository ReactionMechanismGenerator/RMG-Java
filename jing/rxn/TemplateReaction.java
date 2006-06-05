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
import java.util.*;
import jing.param.*;
import jing.param.Temperature;

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
    public double calculatePDepRate(Temperature p_temperature) {
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
    }
    
    //## operation generateReverseForBackwardReaction() 
    private TemplateReaction generateReverseForBackwardReaction() {
        //#[ operation generateReverseForBackwardReaction() 
        // we need to only generate reverse reaction for backward reaction, so that we wont be stuck into a self loop.
        if (!this.isBackward()) return null;
        
        ReactionTemplate fRT = getReactionTemplate();
        ReactionTemplate rRT = null;
        
        if (fRT.isForward()) return null;
        else if (fRT.isNeutral()) rRT = fRT;
        else if (fRT.isBackward()) rRT = fRT.getReverseReactionTemplate();
        else throw new InvalidReactionTemplateDirectionException();
        
        Structure fs = getStructure();
        LinkedList freactant = fs.getReactantList();
        LinkedList fproduct = fs.getProductList();
        
        Structure rs = new Structure(fproduct, freactant);
        TemplateReaction rr = rRT.getReactionFromStructure(rs);
        if (rr != null) {
        	if (!rr.isForward()) {
        		double keq1 = structure.calculateKeq(new Temperature(298,"K"));
        		double keq2 = rr.structure.calculateKeq(new Temperature(298,"K")); 
        		String err = "Backward:" + structure.toString() + '\t' + "Keq = " + String.valueOf(keq1) + '\n';
        		err = err + "Forward:" + rr.structure.toString()+ '\t' + "Keq = " + String.valueOf(keq2) ;
        			 
        		throw new InvalidReactionDirectionException(err);
        	}
        	Reaction fr = rr.getReverseReaction();
        	if (fr != null && fr != this) throw new MultipleReverseReactionException();
        	rr.setReverseReaction(this);
        	return rr;
        }
        else {
        	int rNum = fproduct.size();
        	HashSet rReactionSet = new HashSet();
        	if (rNum == 1) {
        		ChemGraph cg = (ChemGraph)fproduct.getFirst();
        		rReactionSet = rRT.reactOneReactant(cg);	
        	}
        	else if (rNum == 2) {
        		ChemGraph cg1 = (ChemGraph)fproduct.getFirst();
        		ChemGraph cg2 = (ChemGraph)fproduct.getLast(); 
        		rReactionSet = rRT.reactTwoReactants(cg1,cg2);
        		rReactionSet.addAll(rRT.reactTwoReactants(cg2,cg1));
        	}
        	else throw new InvalidReactantNumberException();
        	
        	for (Iterator iter = rReactionSet.iterator(); iter.hasNext();) {
        		TemplateReaction tr = (TemplateReaction)iter.next();
        		if (tr.getStructure().equals(rs)) {
        			if (!tr.isForward()) {
        				String err = "Backward:" + structure.toString() + String.valueOf(structure.calculateKeq(new Temperature(298,"K"))) + '\n';
        				err = err + "Forward:" + tr.structure.toString()+ String.valueOf(tr.structure.calculateKeq(new Temperature(298,"K"))) ;
        				 
        				throw new InvalidReactionDirectionException(err);
        			}
          			tr.setReverseReaction(this);
        			return tr;
        		}
        	}
        	Structure s = getStructure();
        	System.out.println("can't generate reverse reaction for: " + s.toString());
        	for (Iterator iter = s.getReactants(); iter.hasNext(); ) {
        		ChemGraph cg = (ChemGraph)iter.next();
        		System.out.println(cg.getName());
        		System.out.println(cg.getGraph());
        	}
        	
        	throw new FailToGenerateReverseReactionException(getStructure().toString());
        }
        	
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
    public static TemplateReaction makeTemplateReaction(Structure p_structure, Kinetics p_kinetics, ReactionTemplate p_template) {
        //#[ operation makeTemplateReaction(Structure,RateConstant,ReactionTemplate) 
        TemplateReaction reaction = p_template.getReactionFromStructure(p_structure);
        
        if (reaction == null) {
        	reaction = new TemplateReaction(p_structure,p_kinetics,p_template);
        	if (reaction.isBackward()) {
        		TemplateReaction reverse = reaction.generateReverseForBackwardReaction();
        		reaction.setReverseReaction(reverse);
        	}
        	p_template.addReaction(reaction);    
        	if (!(reaction.getReactantNumber() == 2 && reaction.getProductNumber() == 2)) {
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
        
        return reaction;
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return (super.repOk() && reactionTemplate.repOk());
        //#]
    }
    
    //## operation toFullString() 
    /*public String toFullString() {
        //#[ operation toFullString() 
        String s = getStructure().toString() + '\n' + getReactionTemplate().toString() + '\n';
        s = s + getKinetics().toChemkinString() + '\n';
        s = s + "Detailed rate constant information: " + getRateConstant().toString() + '\n';
        s = s + "Comments for this reaction: " + getComments().toString();
        
        return s;
        //#]
    }*/
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String s = getStructure().toString() + '\t' + getReactionTemplate().getName() + '\t';
        Kinetics k = getKinetics();
        String kString = k.toChemkinString();
        if (k instanceof ArrheniusEPKinetics) {
        	double alpha = ((ArrheniusEPKinetics)k).getAlphaValue();
        	kString = kString + '\t' + String.valueOf(alpha);
        }
        return s + kString;
        
        //#]
    }
    
	   //## operation toStringWithReveseReaction() 
    public String toStringWithReveseReaction() {
        //#[ operation toStringWithReveseReaction() 
        TemplateReaction rr = (TemplateReaction)getReverseReaction();
        if (rr == null) return getStructure().toChemkinString(false) + '\t' + getReactionTemplate().getName() + '\t' + getKinetics().toChemkinString();
        else {
        	TemplateReaction temp = null;
        	if (isForward()) temp = this;
        	else if (isBackward()) temp = rr;
        	else throw new InvalidReactionDirectionException();
        	
        	return temp.getStructure().toChemkinString(false) + '\t' + temp.getReactionTemplate().getName() + '\t' + temp.getKinetics().toChemkinString();
        }
        
        //#]
    }
    
    
    public PDepNetwork getPDepNetwork() {
        return pDepNetwork;
    }
    
    public ReactionTemplate getReactionTemplate() {
        return reactionTemplate;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\TemplateReaction.java
*********************************************************************/


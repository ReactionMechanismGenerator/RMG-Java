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
import jing.param.*;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\RateBasedPDepRME.java                                                                  
//----------------------------------------------------------------------------

//## class RateBasedPDepRME 
public class RateBasedPDepRME implements ReactionModelEnlarger {
    
    
    // Constructors
    
    public  RateBasedPDepRME() {
    }
    
    //## operation enlargeReactionModel(ReactionSystem) 
    public void enlargeReactionModel(ReactionSystem p_reactionSystem) {
        //#[ operation enlargeReactionModel(ReactionSystem) 
        Object update = getNextUpdatedObject(p_reactionSystem);
        String return_string = "";
        if (update instanceof Species) {
        	Species next = (Species)update;
        	CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionSystem.getReactionModel();
        	PresentStatus ps = p_reactionSystem.getPresentStatus();
        	
        	System.out.print("\nAdd a new reacted Species:");
        	System.out.println(next.getName());
        	return_string=next.getChemkinName();
			System.out.println(next.toStringWithoutH());
        	Temperature temp = new Temperature(715, "K");
        	double H = next.calculateH(temp);
        	double S = next.calculateS(temp);
        	double G = next.calculateG(temp);
        	double Cp = next.calculateCp(temp);
        	System.out.println("Thermo\t" + String.valueOf(H) + '\t' + String.valueOf(S)+ '\t' + String.valueOf(G)+ '\t' + String.valueOf(Cp));
        	
        	// print the reaction forms this species
            System.out.println("the reactions forming this species:");
            LinkedHashSet ur = cerm.getUnreactedReactionSet();
            for (Iterator iter = ur.iterator(); iter.hasNext();) {
            	Reaction r = (Reaction)iter.next();
            	if (r.contains(next)) {
            		double rate = r.calculateTotalRate(p_reactionSystem.getPresentTemperature());
            		if (r instanceof TemplateReaction) {
            			rate = ((TemplateReaction)r).calculateTotalPDepRate(p_reactionSystem.getPresentTemperature());
            		}
            		System.out.print(r.toString() + "\trate is " + String.valueOf(rate));
            		for (Iterator rIter = r.getReactants(); rIter.hasNext();) {
            			Species reactant = ((ChemGraph)rIter.next()).getSpecies();
            			SpeciesStatus ss = ps.getSpeciesStatus(reactant);
            			rate*=ss.getConcentration(); 
            		}
            		System.out.println("\tflux is " + String.valueOf(rate));
            	}	
            }
            
            
        	if (cerm.containsAsReactedSpecies(next)) 
        		throw new InvalidNextCandidateSpeciesException();
        	else {
        		cerm.moveFromUnreactedToReactedSpecies(next);
        		cerm.moveFromUnreactedToReactedReaction();
        	}
        	// generate new reaction set
        	LinkedHashSet newReactionSet = p_reactionSystem.getReactionGenerator().react(cerm.getReactedSpeciesSet(),next);
        	
        	// partition the reaction set into reacted reaction set and unreacted reaction set
        	// update the corresponding core and edge model of CoreEdgeReactionModel
        	cerm.addReactionSet(newReactionSet);
        
        }
        else if (update instanceof PDepNetwork) {
        	PDepNetwork pnw = (PDepNetwork)update;
        	Species nextIsomer = null;
        	if (!pnw.isActive() && pnw.getIsChemAct()) {
        		nextIsomer = ((ChemGraph)pnw.getProduct().iterator().next()).getSpecies();
        	}
        	else {
        		double maxKLeak = 0;
        		PDepNetReaction path = null;
        		for (Iterator iter = pnw.getPDepNonincludedReactionList(); iter.hasNext(); ) {
        			PDepNetReaction pdnr = (PDepNetReaction)iter.next();
        			double kleak = pdnr.getRate();
        			if (maxKLeak < kleak) {
        				maxKLeak = kleak;
        				path = pdnr;
        			}
        		}
        		
        		if (path == null) throw new InvalidReactionSystemUpdateException();
        		nextIsomer = ((ChemGraph)path.getProducts().next()).getSpecies();
        	}
        	if (nextIsomer == null) throw new InvalidReactionSystemUpdateException();
        	System.out.println("adding new isomer to the PDepNetwork: " + nextIsomer.toString());
        	pnw.update(nextIsomer);
        	pnw.runPDepCalculation(p_reactionSystem);
        }
        else throw new InvalidReactionSystemUpdateException();
        
        return ;
        
        //#]
    }
    
    //## operation getNextUpdatedObject(ReactionSystem) 
    public Object getNextUpdatedObject(ReactionSystem p_reactionSystem) {
        //#[ operation getNextUpdatedObject(ReactionSystem) 
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionSystem.getReactionModel();
        PresentStatus ps = p_reactionSystem.getPresentStatus();
        double Rmin = p_reactionSystem.getRmin();
        
        Species maxSpecies = null;
        double maxFlux = 0;
        
        for (Iterator iter = cerm.getUnreactedSpeciesSet().iterator();iter.hasNext();) {
        	Species us = (Species)iter.next();
        	double thisFlux = Math.abs(ps.getSpeciesStatus(us).getFlux());
        	if (thisFlux > maxFlux) {
        		maxFlux = thisFlux;
        		maxSpecies = us;
        	}
        }
        if (maxSpecies == null) throw new NullPointerException();
        
        double maxRleak = 0;
        PDepNetwork maxPdn = null;
        for (Iterator iter = PDepNetwork.getDictionary().values().iterator(); iter.hasNext();) {
        	PDepNetwork pdn = (PDepNetwork)iter.next();
        	double rleak = pdn.getKLeak();
        	
            for (Iterator rIter = pdn.getReactant().iterator(); rIter.hasNext(); ) {
        		ChemGraph cg = (ChemGraph)rIter.next();
        		Species spe = cg.getSpecies();
        		double conc = 0;
        		if (cerm.containsAsReactedSpecies(spe)) conc = ps.getSpeciesStatus(spe).getConcentration();
        		rleak *= conc;
        	}
        	//System.out.println("PDepNetwork " + pdn.toString() + " with RLeak: " + String.valueOf(rleak));
        	if (rleak > maxRleak) {
        		maxRleak = rleak;
        		maxPdn = pdn;
        	}
        }
        
        System.out.print("Time: ");
        System.out.println(ps.getTime());
        System.out.println("Rmin: " + String.valueOf(Rmin));
        System.out.println("unreacted Spe " + maxSpecies.getName() + " with highest flux: " + String.valueOf(maxFlux));
        System.out.println("PDepNetwork " + maxPdn.toString() + " with highest rleak: " + String.valueOf(maxRleak));
        
        
        if (maxFlux > Rmin && maxRleak < Rmin) {
        	return maxSpecies;
        }
        else {
        	if (maxPdn.getIsChemAct())
        		System.out.println("entry reaction " + maxPdn.getEntryReaction().toString());
        	else {
        		Species entry = ((ChemGraph)maxPdn.getReactant().iterator().next()).getSpecies();
        		System.out.println("entry species " + entry.toString());
        	}
        		
        	return maxPdn;
        }	
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\RateBasedPDepRME.java
*********************************************************************/


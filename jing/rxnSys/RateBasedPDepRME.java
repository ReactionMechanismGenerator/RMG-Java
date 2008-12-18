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
    
    private PDepKineticsEstimator pDepKineticsEstimator;
	
    // Constructors
    
    public  RateBasedPDepRME() {
		pDepKineticsEstimator = null;
    }
    
	public PDepKineticsEstimator getPDepKineticsEstimator() {
		return pDepKineticsEstimator;
	}
	
	public void setPDepKineticsEstimator(PDepKineticsEstimator pdke) {
		pDepKineticsEstimator = pdke;
	}
	
    //9/25/07 gmagoon: added ReactionModel as parameter
    //10/30/07 gmagoon: adjusted parameters to match ReactionModelEnlarger
    //## operation enlargeReactionModel(ReactionSystem) 
    public void enlargeReactionModel(LinkedList p_reactionSystemList, ReactionModel rm, LinkedList p_validList) {
    //public void enlargeReactionModel(ReactionSystem p_reactionSystem, ReactionModel rm) {
        //#[ operation enlargeReactionModel(ReactionSystem)
        
        //10/30/07 gmagoon: updated to loop over reactionSystemList; first loop will select what to add; second loop will add objects
        LinkedList updateList = new LinkedList();
        for (Integer i = 0; i<p_reactionSystemList.size();i++) {
            if(!(Boolean)p_validList.get(i)){
                ReactionSystem p_reactionSystem = (ReactionSystem)p_reactionSystemList.get(i);
                Object update = getNextUpdatedObject(p_reactionSystem, rm);
                //11/1/07 gmagoon: compare the update against all previous updates
                boolean repeat = false;
                int repeatIndex = -1;
                for (int j=0; j<i; j++){
                    if(update == updateList.get(j)){
                        if(!repeat){
                            repeatIndex = j;
                        } 
                        repeat = true;
                    }
                }
                if (repeat)
                    updateList.add("Same as previous update (#" + repeatIndex+")");//11/1/07 gmagoon: ****do we need to run runPDepCalculation to update kLeak for this particular reaction system (or for all reactionSystems after a network is added?)?
                else
                    updateList.add(update);
            }
            else
                updateList.add("(valid)");
        }
        //10/30/07 gmagoon: second loop; ***confirm with Sandeep that no decisions regarding what to include are being made past this point
        for (Integer i = 0; i<p_reactionSystemList.size();i++) {
            if(!(Boolean)p_validList.get(i)){      
                ReactionSystem p_reactionSystem = (ReactionSystem)p_reactionSystemList.get(i);
                Object update = (Object)updateList.get(i);
                
                String return_string = "";
        
                if (update instanceof PDepNetwork) {
					PDepNetwork pnw = (PDepNetwork)update;
					Species nextIsomer = null;
					if (!pnw.isActive() && pnw.getIsChemAct()) {
							nextIsomer = (Species)pnw.getProduct().iterator().next();
					}
					else {
						double maxKLeak = 0;
						PDepNetReaction path = null;
						for (Iterator iter = pnw.getPDepNonincludedReactionListIterator(); iter.hasNext(); ) {
							PDepNetReaction pdnr = (PDepNetReaction)iter.next();

							double kleak = pdnr.calculateRate(p_reactionSystem.getPresentStatus());//10/30/07 gmagoon: changed to include systemSnapshot (presentStatus)
							//double kleak = pdnr.calculateRate();
							if (maxKLeak < kleak) {
								maxKLeak = kleak;
								path = pdnr;
							}
						}

						if (path == null) throw new InvalidReactionSystemUpdateException();
						nextIsomer = (Species)path.getProducts().next();
					}
					if (nextIsomer == null) throw new InvalidReactionSystemUpdateException();

					if (pnw.getIsChemAct())
						System.out.println("entry reaction " + pnw.getEntryReaction().reactionToString(p_reactionSystem.getPresentTemperature()));//10/26/07: gmagoon: changed to use reactionToString rather than toString with temperature passed
						//System.out.println("entry reaction " + pnw.getEntryReaction().toString());
					else {
						Species entry = (Species)pnw.getReactant().iterator().next();
						System.out.println("entry species " + entry.toString());
					}
					System.out.println("adding new isomer to the PDepNetwork: " + nextIsomer.toString());

					if (nextIsomer.getPdepPaths() == null){
						((TemplateReactionGenerator)p_reactionSystem.getReactionGenerator()).generatePdepReactions(nextIsomer);
						((LibraryReactionGenerator)p_reactionSystem.lrg).generatePdepReactions(nextIsomer);
					}
					pnw.update(nextIsomer);
					CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionSystem.getReactionModel();
					if (!cerm.contains(nextIsomer))
						cerm.addUnreactedSpecies(nextIsomer);
					/*Iterator reactionIter = pnw.getPDepNetReactionList();
					while (reactionIter.hasNext()){
							Reaction r = (Reaction)reactionIter.next();
							cerm.addReaction(r);
					}
					}*/

					pDepKineticsEstimator.runPDepCalculation(pnw, p_reactionSystem);

					// return;//10/30/07 gmagoon: I don't think this should be here anymore
                }

                else if (update instanceof Species ) {//10/30/07 gmagoon: changed from if to else if ****confirm this is OK
					Species next = (Species)update;
					//10/10/07 gmagoon: switched to use rm
					//CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionSystem.getReactionModel();
					CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)rm;
					PresentStatus ps = p_reactionSystem.getPresentStatus();

					System.out.print("\nAdd a new reacted Species: ");
					System.out.println(next.getChemkinName());
					return_string=next.getChemkinName();
					System.out.println(next.toStringWithoutH());
					Temperature temp = new Temperature(715, "K");
					double H = next.calculateH(temp);
					double S = next.calculateS(temp);
					double G = next.calculateG(temp);
					double Cp = next.calculateCp(temp);
					System.out.println("Thermo\t" + String.valueOf(H) + " \t" + String.valueOf(S)+ " \t" + String.valueOf(G)+ " \t" + String.valueOf(Cp));



					if (cerm.containsAsReactedSpecies(next)) 
						//throw new InvalidNextCandidateSpeciesException();
						System.out.println("Species is already present in reaction model");//10/31/07 gmagoon: preliminary message; should probably be refined in future
					else {
						cerm.moveFromUnreactedToReactedSpecies(next);
						cerm.moveFromUnreactedToReactedReaction();
					}
					// generate new reaction set
					LinkedHashSet newReactionSet = p_reactionSystem.getReactionGenerator().react(cerm.getReactedSpeciesSet(),next);
					p_reactionSystem.lrg.generatePdepReactions(next);
					newReactionSet.addAll(p_reactionSystem.lrg.react(cerm.getReactedSpeciesSet(),next));

					// partition the reaction set into reacted reaction set and unreacted reaction set
					// update the corresponding core and edge model of CoreEdgeReactionModel
					Iterator iter = newReactionSet.iterator();
					while (iter.hasNext()){
						Reaction r = (Reaction)iter.next();
						if (r.getReactantNumber() ==2 && r.getProductNumber() == 2)
								cerm.addReaction(r);
					}
                }
                //11/1/07 gmagoon: added case where update is instance of string (this would come into play if there is a repeated object)
                else if (update instanceof String ){
                    //do nothing
                }
                else throw new InvalidReactionSystemUpdateException();
            }
			System.out.println("");
        }
    
        return;//10/30/07 gmagoon: added return statement in case it is needed
        //#]
    }

    //## operation getNextUpdatedObject(ReactionSystem) 
    public Object getNextUpdatedObject(ReactionSystem p_reactionSystem, ReactionModel rm) {
        //#[ operation getNextUpdatedObject(ReactionSystem) 
        //10/10/07 gmagoon: changed to use rm rather than reaction system
        //CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionSystem.getReactionModel();
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)rm;
        PresentStatus ps = p_reactionSystem.getPresentStatus();
        double Rmin = p_reactionSystem.getRmin();
        Species maxSpecies = null;
        double maxFlux = 0;
       
        for (Iterator iter = cerm.getUnreactedSpeciesSet().iterator();iter.hasNext();) {
        	Species us = (Species)iter.next();
        	double thisFlux = Math.abs(Math.abs(ps.unreactedSpeciesFlux[us.getID()]));
        	if (thisFlux >= maxFlux) {
        		maxFlux = thisFlux;
        		maxSpecies = us;
        	}
        }

        if (maxSpecies == null) throw new NullPointerException();

        double maxRleak = 0;
        PDepNetwork maxPdn = null;
        for (Iterator iter = PDepNetwork.getDictionary().values().iterator(); iter.hasNext();) {
        	PDepNetwork pdn = (PDepNetwork)iter.next();
                double rleak = pdn.getKLeak(p_reactionSystem.getIndex());//10/30/07 gmagoon: changed to pass index to getKleak 
        	//double rleak = pdn.getKLeak();
                for (Iterator rIter = pdn.getReactant().iterator(); rIter.hasNext(); ) {
            	Species spe = (Species)rIter.next();
        		//Species spe = cg.getSpecies();
        		double conc = 0;
        		if (cerm.containsAsReactedSpecies(spe) && ps.getSpeciesStatus(spe) != null) conc = ps.getSpeciesStatus(spe).getConcentration();
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
        System.out.println("Unreacted species " + maxSpecies.getName() + " has highest flux: " + String.valueOf(maxFlux));
        if (maxPdn != null)
        	System.out.println("PDepNetwork " + maxPdn.toString() + " has highest rleak: " + String.valueOf(maxRleak));

        // return max pdep network only if there is no species present
        if (maxRleak >= Rmin && maxFlux <= Rmin) {
					System.out.println("Adding the PDepNetwork");
        	if (maxPdn.getIsChemAct())
        		System.out.println("entry reaction " + maxPdn.getEntryReaction().toString());
        	else {
        		Species entry = (Species)maxPdn.getReactant().iterator().next();
        		System.out.println("entry species " + entry.toString());
        	}
        	return maxPdn;
        }
        else {
					System.out.println("Adding the unreacted species");
        	return maxSpecies;
        }	
        
        // return max species if there is no incomplete pressure dependent reaction network
        /*if (maxFlux >= Rmin && maxRleak <= Rmin) {
        	return maxSpecies;
        }
        else {
        	if (maxPdn.getIsChemAct())
        		System.out.println("entry reaction " + maxPdn.getEntryReaction().toString());
        	else {
        		Species entry = (Species)maxPdn.getReactant().iterator().next();
        		System.out.println("entry species " + entry.toString());
        	}
        		
        	return maxPdn;
        }*/	
        //#]
   }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\RateBasedPDepRME.java
*********************************************************************/


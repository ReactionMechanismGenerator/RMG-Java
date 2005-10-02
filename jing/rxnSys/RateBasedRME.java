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


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import jing.param.*;
import jing.rxn.Reaction;
import jing.chem.Species;
//import RMG;
//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\RateBasedRME.java                                                                  
//----------------------------------------------------------------------------

//## class RateBasedRME 
public class RateBasedRME implements ReactionModelEnlarger {
    
    
    // Constructors
    
    public  RateBasedRME() {
    }
    
    //## operation enlargeReactionModel(ReactionSystem) 
    public String enlargeReactionModel(ReactionSystem p_reactionSystem) {
        //#[ operation enlargeReactionModel(ReactionSystem) 
        ReactionModel rm = p_reactionSystem.getReactionModel();
        if (!(rm instanceof CoreEdgeReactionModel)) throw new InvalidReactionModelTypeException();
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)rm;
        
        PresentStatus ps = p_reactionSystem.getPresentStatus();
        String maxflux="";
        Species next = getNextCandidateSpecies(cerm,ps,maxflux);
		/*File coreSpecies = new File("Restart/coreSpecies.txt");
		ChemParser.removeEnd(coreSpecies);
		try{
			FileWriter fw = new FileWriter()
		}*/
		double time = System.currentTimeMillis();//ReactionModelGenerator.getCpuTime()/1e9/60;
        System.out.print("\nAdd a new reacted Species:");
        System.out.println(next.getName());
        Temperature temp = new Temperature(298, "K");
        double H = next.calculateH(temp);
        double S = next.calculateS(temp);
        double G = next.calculateG(temp);
        double Cp = next.calculateCp(temp);
        System.out.println("Thermo\t" + String.valueOf(H) + '\t' + String.valueOf(S)+ '\t' + String.valueOf(G)+ '\t' + String.valueOf(Cp));
        
        if (cerm.containsAsReactedSpecies(next)) 
        	throw new InvalidNextCandidateSpeciesException();
        else {
        	cerm.moveFromUnreactedToReactedSpecies(next);
        	cerm.moveFromUnreactedToReactedReaction();
        	
        }
        // generate new reaction set
        HashSet newReactionSet = p_reactionSystem.getReactionGenerator().react(cerm.getReactedSpeciesSet(),next);
		
		String restartFileContent="";
		try{
			File allReactions = new File ("Restart/allReactions.txt");
			FileWriter fw = new FileWriter(allReactions, true);
			//Species species = (Species) iter.next();
			for(Iterator iter=newReactionSet.iterator();iter.hasNext();){
				
				Reaction reaction = (Reaction) iter.next();
				if (cerm.categorizeReaction(reaction)==-1)
					restartFileContent = restartFileContent + reaction.toRestartString() + "\n";
								
			}
			
			//restartFileContent += "\nEND";
			fw.write(restartFileContent);
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the added Reactions to the allReactions file");
        	System.exit(0);
		}
        
		try{
			File coreReactions = new File ("Restart/coreReactions.txt");
			FileWriter fw = new FileWriter(coreReactions, true);
			//Species species = (Species) iter.next();
			restartFileContent="";
			for(Iterator iter=newReactionSet.iterator();iter.hasNext();){
				
				Reaction reaction = (Reaction) iter.next();
				if (cerm.categorizeReaction(reaction)==1&&reaction.getDirection()==1)
					restartFileContent = restartFileContent + reaction.toRestartString() + "\n";
				else if (cerm.categorizeReaction(reaction)==1&&reaction.getDirection()==-1)
					restartFileContent = restartFileContent + reaction.getReverseReaction().toRestartString() + "\n";
			}
			
			//restartFileContent += "\nEND";
			fw.write(restartFileContent);
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the added Reactions to the allReactions file");
        	System.exit(0);
		}
		
        // partition the reaction set into reacted reaction set and unreacted reaction set
        // update the corresponding core and edge model of CoreEdgeReactionModel
        cerm.addReactionSet(newReactionSet);
        //String return_string;
        return next.getChemkinName()+"\t" + maxflux+"\t"+ ((RateBasedVT)p_reactionSystem.finishController.validityTester).Rmin+ "\t" + time;
        
        
        //#]
    }
    
    //## operation getNextCandidateSpecies(CoreEdgeReactionModel,PresentStatus) 
    public Species getNextCandidateSpecies(CoreEdgeReactionModel p_reactionModel, PresentStatus p_presentStatus, String maxflux) {
        //#[ operation getNextCandidateSpecies(CoreEdgeReactionModel,PresentStatus) 
        HashSet unreactedSpecies = p_reactionModel.getUnreactedSpeciesSet();
         
        Species maxSpecies = null;
        double maxFlux = 0;
        
        Iterator iter = unreactedSpecies.iterator();
        while (iter.hasNext()) {
        	Species us = (Species)iter.next();
        	double thisFlux = Math.abs(p_presentStatus.getSpeciesStatus(us).getFlux());
        	if (thisFlux > maxFlux) {
        		maxFlux = thisFlux;
        		maxSpecies = us;
        	}
        }
        maxflux = ""+maxFlux;
        if (maxSpecies == null) throw new NullPointerException();
        
        System.out.print("Time: ");
        System.out.println(p_presentStatus.getTime());
        System.out.println("unreacted Spe with highest flux: " + String.valueOf(maxFlux));
        
        return maxSpecies;
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\RateBasedRME.java
*********************************************************************/


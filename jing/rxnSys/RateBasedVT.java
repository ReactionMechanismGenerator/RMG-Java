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


import java.util.*;

import jing.chem.Species;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\RateBasedVT.java                                                                  
//----------------------------------------------------------------------------

/**
In rate-based model generation process, the criteria is:
Rj<e*Rchar
in which, Rj is the unreacted species flux, e is the user tolerance, and Rchar is the system characteristic flux calculated by L2 norm of reacted species fluxes.
*/
//## class RateBasedVT 
public class RateBasedVT implements ValidityTester {
    
    protected double Rmin = 0;		//## attribute Rmin 
    
    protected double tolerance;		//## attribute tolerance 
    
    
    // Constructors
    
    //## operation RateBasedVT() 
    protected  RateBasedVT() {
        //#[ operation RateBasedVT() 
        //#]
    }
    //## operation RateBasedVT(double) 
    public  RateBasedVT(double p_tolerance) {
        //#[ operation RateBasedVT(double) 
        tolerance = p_tolerance;
        //#]
    }
    
    //## operation calculateRchar(PresentStatus) 
    public double calculateRchar(PresentStatus p_presentStatus) {
        //#[ operation calculateRchar(PresentStatus) 
        double minflux = 0;
        Iterator iter = p_presentStatus.getSpeciesStatus();
        while (iter.hasNext()) {
        	SpeciesStatus ss = (SpeciesStatus)iter.next();
        	if (ss.isReactedSpecies()) {
        		double flux = ss.getFlux();
        		minflux += flux*flux;
        	}
        }
        // calculate L2 norm of reacted species
        return Math.sqrt(minflux);
        
        
        //#]
    }
    
    //## operation calculateRmin(PresentStatus) 
    public double calculateRmin(PresentStatus p_presentStatus) {
        //#[ operation calculateRmin(PresentStatus) 
        Rmin = tolerance*calculateRchar(p_presentStatus);
        return Rmin;
        //#]
    }
    
    //## operation isModelValid(ReactionSystem) 
    public boolean isModelValid(ReactionSystem p_reactionSystem) {
        //#[ operation isModelValid(ReactionSystem) 
        // check if all the unreacted species has their fluxes under the system min flux
        PresentStatus ps = p_reactionSystem.getPresentStatus();
        calculateRmin(ps);
        for (Iterator iter =((CoreEdgeReactionModel) p_reactionSystem.getReactionModel()).getUnreactedSpeciesSet().iterator(); iter.hasNext(); ) {
        	Species s = (Species)iter.next();
        	if (ps.unreactedSpeciesFlux[s.getID()] > Rmin) return false;
        	
        }
        //System.out.println("The minimum flux required is " + Rmin);
		
        return true;
        
        
        //#]
    }
    
    public double getTolerance() {
        return tolerance;
    }
    
    public void setTolerance(double p_tolerance) {
        tolerance = p_tolerance;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\RateBasedVT.java
*********************************************************************/


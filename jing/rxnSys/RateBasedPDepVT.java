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
// jing\rxnSys\RateBasedPDepVT.java                                                                  
//----------------------------------------------------------------------------

//## class RateBasedPDepVT 
public class RateBasedPDepVT extends RateBasedVT {
    
    
    // Constructors
    
    //## operation RateBasedPDepVT(double) 
    public  RateBasedPDepVT(double p_tolerance) {
        //#[ operation RateBasedPDepVT(double) 
        super(p_tolerance);
        //#]
    }
    public  RateBasedPDepVT() {
    }
    
    //## operation isModelValid(ReactionSystem) 
    public boolean isModelValid(ReactionSystem p_reactionSystem) {
        
		// Check if all the unreacted species has their fluxes under the system min flux
        if (!super.isModelValid(p_reactionSystem))
			return false;

		// check if all the networks has their leak fluxes under the system min flux
        PresentStatus ps = p_reactionSystem.getPresentStatus();
        calculateRmin(ps);

		for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
        	PDepNetwork pdn = (PDepNetwork )iter.next();
        	if (pdn.getLeakFlux(ps) > Rmin)
            {
                System.out.println("Exceeded largest permitted flux for convergence (tolerance="+tolerance+"): " + Rmin);
                return false;
            }

        }
        
        return true;


        //#]

    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\RateBasedPDepVT.java
*********************************************************************/


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
import jing.param.Pressure;
import jing.param.Temperature;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\InitialStatus.java                                                                  
//----------------------------------------------------------------------------

//## class InitialStatus 
public class InitialStatus extends SystemSnapshot {
    
    protected static double colliderLimit = 0.1;		//## attribute colliderLimit 
    
    
    // Constructors
    
    //## operation InitialStatus() 
    private  InitialStatus() {
        //#[ operation InitialStatus() 
        //#]
    }
    //## operation InitialStatus(HashMap) 
	  public  InitialStatus(LinkedHashMap p_speciesStatus, Temperature p_temperature, Pressure p_pressure) {
	        //#[ operation InitialStatus(HashMap,Temperature,Pressure) 
	        super(new ReactionTime(0,"S"),p_speciesStatus,p_temperature,p_pressure);
	        //#]
	    }
    
    //## operation identifyColliders() 
    public HashMap identifyColliders() {
        //#[ operation identifyColliders() 
        HashMap result = new HashMap();
        double totalMole = getTotalMole();
        double adjTotalMole = 0;
        for (Iterator iter = getSpeciesStatus(); iter.hasNext(); ) {
        	SpeciesStatus ss = (SpeciesStatus)iter.next();
        	double conc = ss.getConcentration();
        	if (conc > colliderLimit*totalMole) {
        		adjTotalMole += conc;
        		result.put(ss.getSpecies(), new Double(conc)); 
        	}
        }
        for (Iterator iter = inertGas.keySet().iterator(); iter.hasNext(); ) {
        	Object key = iter.next();
        	double conc = ((Double)inertGas.get(key)).doubleValue(); 
        	if (conc < 0) throw new NegativeConcentrationException("InertGas");
        	if (conc > colliderLimit*totalMole) {
        		adjTotalMole += conc;
        		result.put(key, new Double(conc)); 
        	}
        }
        
        if (result.isEmpty()) throw new InvalidSystemCompositionException();
        
        return result;
        	
        //#]
    }
    
    private static double getColliderLimit() {
        return colliderLimit;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\InitialStatus.java
*********************************************************************/


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


import java.util.*;
import jing.rxnSys.*;
import jing.rxnSys.SystemSnapshot;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ThirdBodyReaction.java                                                                  
//----------------------------------------------------------------------------

//## class ThirdBodyReaction 
public class ThirdBodyReaction extends Reaction {
    
    protected HashMap weightMap = new HashMap();		//## attribute weightMap 
    
    
    // Constructors
    
    public  ThirdBodyReaction() {
    }
    
    //## operation calculateThirdBodyCoefficient(SystemSnapshot) 
    public double calculateThirdBodyCoefficient(SystemSnapshot p_presentStatus) {
        //#[ operation calculateThirdBodyCoefficient(SystemSnapshot) 
        double coef_total = 0;
        for (Iterator iter = p_presentStatus.getSpeciesStatus(); iter.hasNext(); ) {
        	SpeciesStatus ss = (SpeciesStatus)iter.next();
        	double conc = ss.getConcentration();
        	double coef = 1;
        	
        	String name = ss.getSpecies().getName();
        	if (weightMap.containsKey(name)) {
        		coef = ((Double)weightMap.get(name)).doubleValue();
        	}
        	
        	coef_total += coef*conc;
        }
        
        for (Iterator iter = p_presentStatus.getInertGas(); iter.hasNext(); ) {
        	String name = (String)iter.next();
        	double conc = p_presentStatus.getInertGas(name);
        	double coef = 1;
        	
        	if (weightMap.containsKey(name)) {
        		coef = ((Double)weightMap.get(name)).doubleValue();
        	}
        	
        	coef_total += coef*conc;
        }
        
        return coef_total;
        //#]
    }
    
    //## operation generateReverseReaction() 
    public void generateReverseReaction() {
        //#[ operation generateReverseReaction() 
        ThirdBodyReaction r = new ThirdBodyReaction();
        r.structure = getStructure().generateReverseStructure();
        r.rateConstant = getRateConstant();
        r.comments = "Reverse reaction";
        r.weightMap = weightMap;
        	
        r.setReverseReaction(this);
        this.setReverseReaction(r);
        
        return;
        //#]
    }
    
    //## operation make(Reaction,HashMap) 
    public static ThirdBodyReaction make(Reaction p_reaction, HashMap p_thirdBodyList) {
        //#[ operation make(Reaction,HashMap) 
        ThirdBodyReaction tbr = new ThirdBodyReaction();
        tbr.structure = p_reaction.getStructure();
        tbr.rateConstant = p_reaction.getRateConstant();
        tbr.comments = p_reaction.getComments();
        tbr.generateReverseReaction(); 
        
        p_reaction = null;
        
        return tbr;
        
        
        //#]
    }
    
    //## operation putThirdBodyCoefficient(String,double) 
    public void putThirdBodyCoefficient(String p_name, double p_coefficient) {
        //#[ operation putThirdBodyCoefficient(String,double) 
        weightMap.put(p_name,new Double(p_coefficient));
        
        
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ThirdBodyReaction.java
*********************************************************************/


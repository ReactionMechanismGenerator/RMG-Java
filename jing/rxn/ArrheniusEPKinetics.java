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
import jing.param.*;
import jing.param.Temperature;
import jing.mathTool.UncertainDouble;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ArrheniusEPKinetics.java                                                                  
//----------------------------------------------------------------------------

//## class ArrheniusEPKinetics 
public class ArrheniusEPKinetics extends ArrheniusKinetics {
    
    protected UncertainDouble alpha = new UncertainDouble(0,0,"Adder");		//## attribute alpha 
    
    
    // Constructors
    
    //## operation ArrheniusEPKinetics(UncertainDouble,UncertainDouble,UncertainDouble,UncertainDouble,String,int,String,String) 
    public  ArrheniusEPKinetics(UncertainDouble p_A, UncertainDouble p_n, UncertainDouble p_alpha, UncertainDouble p_E, String p_TRange, int p_rank, String p_source, String p_comment) {
        //#[ operation ArrheniusEPKinetics(UncertainDouble,UncertainDouble,UncertainDouble,UncertainDouble,String,int,String,String) 
        super(p_A, p_n, p_E, p_TRange, p_rank, p_source, p_comment);
        alpha = p_alpha;
        
        //#]
    }
    public  ArrheniusEPKinetics() {
    }
    
    //## operation calculateRate(Temperature,double) 
    public double calculateRate(Temperature p_temperature, double p_Hrxn) {
        //#[ operation calculateRate(Temperature,double) 
        double T = p_temperature.getStandard();
        double R = GasConstant.getKcalMolK();
        
        double Ea = E.getValue() + alpha.getValue()*p_Hrxn;
        //if (Ea<0) throw new NegativeEnergyBarrierException();
        double rate = A.getValue() * Math.pow(T, n.getValue()) * Math.exp(-Ea/R/T);
        return rate;
        
        //#]
    }
    
    //## operation getAlphaUncertainty() 
    public double getAlphaUncertainty() {
        //#[ operation getAlphaUncertainty() 
        return alpha.getUncertainty();
        //#]
    }
    
    //## operation getAlphaValue() 
    public double getAlphaValue() {
        //#[ operation getAlphaValue() 
        return alpha.getValue();
        //#]
    }
    
    //## operation multiply(double) 
    public Kinetics multiply(double p_multiple) {
        //#[ operation multiply(double) 
        UncertainDouble newA = getA().multiply(p_multiple);
        Kinetics newK = new ArrheniusEPKinetics(newA,getN(),getAlpha(),getE(),getTRange(),getRank(),getSource(),getComment());
        return newK;
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String string = "A = " + A.toString() + "; n = " + n.toString() + "; E = " + E.toString() + ";" + " alpha = " + alpha.toString() + '\n';
        string = string + "Source: " + source + '\n';
        string = string + "Comments: " + comment + '\n';
        return string;
        //#]
    }
    
    public UncertainDouble getAlpha() {
        return alpha;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ArrheniusEPKinetics.java
*********************************************************************/


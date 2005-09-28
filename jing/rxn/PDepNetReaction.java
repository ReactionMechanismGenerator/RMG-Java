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

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\PDepNetReaction.java                                                                  
//----------------------------------------------------------------------------

//## class PDepNetReaction 
public class PDepNetReaction extends Reaction {
    
    protected double pressure;		//## attribute pressure 
    
    protected double rate = 0;		//## attribute rate 
    
    protected double temperature;		//## attribute temperature 
    
    
    // Constructors
    
    //## operation PDepNetReaction(LinkedList,LinkedList,double,double,double) 
    public  PDepNetReaction(LinkedList p_reactant, LinkedList p_product, double p_rate, double p_temperature, double p_pressure) {
        //#[ operation PDepNetReaction(LinkedList,LinkedList,double,double,double) 
        if (p_rate < 0 || p_temperature < 0 || p_pressure < 0) throw new InvalidPDepNetReactionException();
        Structure st = new Structure(p_reactant, p_product,1);
        structure = st;
        rate = p_rate;
        temperature = p_temperature;
        pressure = p_pressure;
        
        //#]
    }
    public  PDepNetReaction() {
    }
    
    public double getPressure() {
        return pressure;
    }
    
    public double getRate() {
        return rate;
    }
	
//	## operation formPDepSign(String) 
    public String formPDepSign(String p_string) {
        //#[ operation formPDepSign(String) 
        StringTokenizer st = new StringTokenizer(p_string, "=");
        String s1 = st.nextToken();
        s1 += "(+m)=";
        String s2 = st.nextToken();
        s2 += "(+m)";
        return (s1+s2);
        
        //#]
    }
	
//	## operation toChemkinString() 
    public String toChemkinString() {
        //#[ operation toChemkinString() 
		
        String result = getStructure().toChemkinString(true) + '\t' + rate +" 0.0 0.0" + '\n';
        //result += getItsChebyshevPolynomials().toChemkinString() + '\n';
        return result;
        
        //#]
    }
    public void setRate(double p_rate) {
        rate = p_rate;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\PDepNetReaction.java
*********************************************************************/


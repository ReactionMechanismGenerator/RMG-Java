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
import jing.rxnSys.SystemSnapshot;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\PDepNetReaction.java                                                                  
//----------------------------------------------------------------------------

//## class PDepNetReaction 
public class PDepNetReaction extends Reaction {
    
    protected ChebyshevPolynomials itsChebyshevPolynomials;
    
    // Constructors
    
    //## operation PDepNetReaction(LinkedList,LinkedList,ChebyshevPolynomials) 
    public  PDepNetReaction(LinkedList p_reactant, LinkedList p_product, final ChebyshevPolynomials p_cp) {
        //#[ operation PDepNetReaction(LinkedList,LinkedList,ChebyshevPolynomials) 
        structure = new Structure(p_reactant, p_product,1);
        itsChebyshevPolynomials = p_cp;
        if (itsChebyshevPolynomials.getPlow().getAtm() <= 0){
            Pressure P = new Pressure(0.01,"atm");
            itsChebyshevPolynomials.setPlow(P);
          }
        //#]
    }
    public  PDepNetReaction() {
    }
    
    //## operation calculateRate(SystemSnapshot) 
    public double calculateRate(SystemSnapshot p_systemSnapshot) {
        //#[ operation calculateRate(SystemSnapshot) 
        Temperature temp = p_systemSnapshot.getTemperature();
        Pressure pres = p_systemSnapshot.getPressure();
        if (itsChebyshevPolynomials == null) throw new NullPointerException();
        return itsChebyshevPolynomials.calculateRate(temp, pres);
        //#]
    }
    
    
//  ## operation calculateRate(SystemSnapshot) 
    public double calculateRate() {
        //#[ operation calculateRate(SystemSnapshot) 
        Temperature temp = Global.temperature;
        Pressure pres = Global.pressure;
        if (itsChebyshevPolynomials == null) throw new NullPointerException();
        return itsChebyshevPolynomials.calculateRate(temp, pres);
        //#]
    }
    
    //## operation formPDepSign(String) 
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
    
    //## operation toChemkinString() 
    public String toChemkinString(Temperature t) {
        //#[ operation toChemkinString() 
        String result = formPDepSign(getStructure().toChemkinString(true).toString()) + '\t' + "1.0E0 0.0 0.0" + '\n';
        result += getItsChebyshevPolynomials().toChemkinString() + '\n';
        return result;
        
        //#]
    }
    /**
     * Generates a reaction whose structure is opposite to that of the present reaction.
     * Just appends the rate constant of this reaction to the reverse reaction.
     *
     */
    //## operation generateReverseReaction()
    public void generateReverseReaction() {
        //#[ operation generateReverseReaction()
        

        
        
        PDepNetReaction r = new PDepNetReaction(structure.products, structure.reactants, getItsChebyshevPolynomials());

  	  
  	  
        
        this.setReverseReaction(r);

        return;
        //#]
    }
    
    public String toString(){
    	return toChemkinString(Global.temperature);
    }
    
    public ChebyshevPolynomials getItsChebyshevPolynomials() {
        return itsChebyshevPolynomials;
    }
    
    public void setItsChebyshevPolynomials(ChebyshevPolynomials p_ChebyshevPolynomials) {
        itsChebyshevPolynomials = p_ChebyshevPolynomials;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\PDepNetReaction.java
*********************************************************************/


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


import jing.chem.*;
import java.util.*;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\KineticsTemplate.java                                                                  
//----------------------------------------------------------------------------

//## class KineticsTemplate 
public class KineticsTemplate {
    
    protected HashSet key;		//## attribute key 
    
    protected Kinetics kinetics;
    
    // Constructors
    
    //## operation KineticsTemplate(HashSet,Kinetics) 
    protected  KineticsTemplate(HashSet p_key, Kinetics p_kinetics) {
        //#[ operation KineticsTemplate(HashSet,Kinetics) 
        key = p_key;
        kinetics = p_kinetics;
        
        
        
        //#]
    }
    public  KineticsTemplate() {
    }
    
    //## operation printKey() 
    public String printKey() throws InvalidKineticsTemplateException, InvalidFunctionalGroupException {
        //#[ operation printKey() 
        String s = "";
        for (Iterator iter = key.iterator(); iter.hasNext(); ) {
        	Object fg = iter.next();
        	if (fg instanceof String) {
        		s = s + (String)fg;
        	}
        	else if (fg instanceof Matchable) {
        		s = s + ((Matchable)fg).getName();
        	}
        	else {
        		throw new InvalidKineticsKeyException();
        	}
        	if (iter.hasNext()) s = s + " + ";
        }
        if (s.equals("")) throw new InvalidKineticsTemplateException();
        
        return s;
        	
        	
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return kinetics.repOk();
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        return "Use kinetics for: "+ printKey() + '\n' + kinetics.toString();
        //#]
    }
    
    public HashSet getKey() {
        return key;
    }
    
    public void setKey(HashSet p_key) {
        key = p_key;
    }
    
    public Kinetics getKinetics() {
        return kinetics;
    }
    
    public void setKinetics(Kinetics p_Kinetics) {
        kinetics = p_Kinetics;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\KineticsTemplate.java
*********************************************************************/


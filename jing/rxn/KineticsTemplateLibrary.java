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


import java.io.*;
import jing.chem.*;
import java.util.*;
import jing.mathTool.*;
import jing.chemUtil.*;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\KineticsTemplateLibrary.java                                                                  
//----------------------------------------------------------------------------

//## class KineticsTemplateLibrary 
public class KineticsTemplateLibrary {
    
    protected HashMap kineticsTemplate;
    
    // Constructors
    
    public  KineticsTemplateLibrary() {
        {
            kineticsTemplate=new HashMap();
        }
    }
    
    //## operation addKinetics(HashSet,Kinetics) 
    public KineticsTemplate addKinetics(HashSet p_fgc, Kinetics p_kinetics) {
        //#[ operation addKinetics(HashSet,Kinetics) 
        KineticsTemplate old = getKineticsTemplate(p_fgc);
        // if there is already a number in the library, and it is not the root node, output information
        // otherwise, add the kt number into library, or replace the root number
        if (old != null && old.getKinetics().getRank() != 0) {
        	// if rate contant already exist, output a warning message, keep the original number, don't replace
        	String s = "";
        	for (Iterator iter = p_fgc.iterator(); iter.hasNext(); ) {
        		Object fg = iter.next();
        		if (fg instanceof String) {
        			s = s + fg;
        		}
        		else if (fg instanceof Matchable) {
        			s = s + ((Matchable)fg).getName();
        		}
        		else {
        			throw new InvalidKineticsKeyException();
        		}
        		if (iter.hasNext()) s = s + " + ";
        	}
        	System.out.println("multiple value find for: " + s);
        	return old;
        }
        else {
        	KineticsTemplate kt = new KineticsTemplate(p_fgc,p_kinetics);
        	addKineticsTemplate(kt);
        	return kt;
        }
        
        
        
        //#]
    }
    
    //## operation addKineticsTemplate(KineticsTemplate) 
    private Object addKineticsTemplate(KineticsTemplate p_kineticsTemplate) {
        //#[ operation addKineticsTemplate(KineticsTemplate) 
        return kineticsTemplate.put(p_kineticsTemplate.getKey(),p_kineticsTemplate);
        //#]
    }
    
    //## operation getKinetics(HashSet) 
    public Kinetics getKinetics(HashSet p_key) {
        //#[ operation getKinetics(HashSet) 
        KineticsTemplate kt = getKineticsTemplate(p_key);
        if (kt==null) return null;
        else return kt.getKinetics();
        //#]
    }
    
    //## operation getKineticsTemplate(HashSet) 
    public KineticsTemplate getKineticsTemplate(HashSet p_key) {
        //#[ operation getKineticsTemplate(HashSet) 
        return (KineticsTemplate)(kineticsTemplate.get(p_key));
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        // check if each kinetics is fine
        Iterator iter = getKineticsTemplate();
        while (iter.hasNext()) {
        	KineticsTemplate kt = getKineticsTemplate((HashSet)iter.next());
        	if (!kt.repOk()) return false;
        }
        return true;
        
        
        
        
        //#]
    }
    
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return kineticsTemplate.size();
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String s ="";
        int index = 0;
        
        Iterator key_iter = getKineticsTemplate();
        while (key_iter.hasNext()) {
        	index++;
        	HashSet key = (HashSet)key_iter.next();
        	KineticsTemplate kt = getKineticsTemplate(key);
        	s = s + kt.toString() + '\n';
        }
        
        return s;
        
        
        
        //#]
    }
    
    public Iterator getKineticsTemplate() {
        Iterator iter=kineticsTemplate.keySet().iterator();
        return iter;
    }
    
    public void clearKineticsTemplate() {
        kineticsTemplate.clear();
    }
    
    public void removeKineticsTemplate(KineticsTemplate p_KineticsTemplate) {
        Iterator iter=kineticsTemplate.keySet().iterator();
        while(iter.hasNext()) {
          Object key = iter.next();
          if (kineticsTemplate.get(key).equals(p_KineticsTemplate)) {
          	kineticsTemplate.remove(key);
          	break;
          }
        };
    }
    
    public void removeKineticsTemplate(HashSet key) {
        KineticsTemplate p_KineticsTemplate = getKineticsTemplate(key);
        kineticsTemplate.remove(key);
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\KineticsTemplateLibrary.java
*********************************************************************/


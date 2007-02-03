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
import jing.mathTool.*;
import jing.param.Global;
import jing.param.Temperature;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\RateConstant.java                                                                  
//----------------------------------------------------------------------------

//## class RateConstant 
public class RateConstant {
    
    protected String comments = "";		//## attribute comments 
    
    protected int distance = -1;		//## attribute distance 
    
    protected KineticsTemplate kineticsTemplate;
    protected HashSet kineticsTemplateFromLibrary;
    
    // Constructors
    
    //## operation RateConstant(KineticsTemplate,int) 
    public  RateConstant(KineticsTemplate p_kineticsTemplate, int p_distance) {
        {
            kineticsTemplateFromLibrary=new HashSet();
        }
        //#[ operation RateConstant(KineticsTemplate,int) 
        kineticsTemplate = p_kineticsTemplate;
        distance = p_distance;
        
        
        //#]
    }
    //## operation RateConstant(KineticsTemplate,HashSet,int) 
    public  RateConstant(KineticsTemplate p_kineticsTemplate, HashSet p_kineticsTemplateFromLibrary, int p_distance) {
        {
            kineticsTemplateFromLibrary=new HashSet();
        }
        //#[ operation RateConstant(KineticsTemplate,HashSet,int) 
        kineticsTemplate = p_kineticsTemplate;
        kineticsTemplateFromLibrary = p_kineticsTemplateFromLibrary;
        distance = p_distance;
        
        
        //#]
    }
    public  RateConstant() {
        {
            kineticsTemplateFromLibrary=new HashSet();
        }
    }
    
    //## operation calculateRate(Temperature,double) 
    public double calculateRate(Temperature p_temperature, double p_Hrxn) {
        //#[ operation calculateRate(Temperature,double) 
        return kineticsTemplate.getKinetics().calculateRate(p_temperature,p_Hrxn);
        //#]
    }
    
    //## operation getKinetics() 
    public Kinetics getKinetics() {
        //#[ operation getKinetics() 
        return kineticsTemplate.getKinetics();
        //#]
    }
    
    //## operation getRank() 
    public int getRank() {
        //#[ operation getRank() 
        return getKinetics().getRank();
        //#]
    }
    
    //## operation printKey() 
    public String printKey() {
        //#[ operation printKey() 
        return kineticsTemplate.printKey();
        
        
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return kineticsTemplate.repOk();
        
        
        //#]
    }
    
    //## operation toChemkinString() 
    /*public String toChemkinString() {
        //#[ operation toChemkinString() 
        return getKineticsTemplate().getKinetics().toChemkinString(Global.temperature);
        //#]
    }*/
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String s = "Kinetics: " + kineticsTemplate.toString() + " Distance = " + String.valueOf(distance) + '\n';
        if (kineticsTemplateFromLibrary != null && kineticsTemplateFromLibrary.size()>0) {
        	s = s + "Use average of: " + '\n';
        	Iterator iter = kineticsTemplateFromLibrary.iterator();
        	while (iter.hasNext()) {
        		KineticsTemplate kt = (KineticsTemplate)iter.next();
        		s = s + kt.toString() + '\n';
        	}
        }
        
        return s;
        		
        //#]
    }
    
    public String getComments() {
        return comments;
    }
    
    public void setComments(String p_comments) {
        comments = p_comments;
    }
    
    public int getDistance() {
        return distance;
    }
    
    public KineticsTemplate getKineticsTemplate() {
        return kineticsTemplate;
    }
    
    public Iterator getKineticsTemplateFromLibrary() {
        Iterator iter=kineticsTemplateFromLibrary.iterator();
        return iter;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\RateConstant.java
*********************************************************************/


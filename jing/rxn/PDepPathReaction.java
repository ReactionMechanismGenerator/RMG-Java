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
// jing\rxn\PDepPathReaction.java                                                                  
//----------------------------------------------------------------------------

//## class PDepPathReaction 
public class PDepPathReaction {
    
    protected String type;		//## attribute type 
    
    protected TemplateReaction templateReaction;
    
    // Constructors
    
    //## operation PDepPathReaction(TemplateReaction) 
    protected  PDepPathReaction(TemplateReaction p_templateReaction) {
        //#[ operation PDepPathReaction(TemplateReaction) 
        templateReaction = p_templateReaction;
        type = null;
        //#]
    }
    public  PDepPathReaction() {
    }
    
    //## operation getKinetics() 
    public Kinetics getKinetics() {
        //#[ operation getKinetics() 
        Kinetics k = null;
        if (templateReaction.isForward()) {
        	k = templateReaction.getKinetics();
        }
        else {
        	k = templateReaction.getFittedReverseKinetics();
        }
        
        return k;
        
        //#]
    }
    
    //## operation getProductList() 
    public LinkedList getProductList() {
        //#[ operation getProductList() 
        return templateReaction.getProductList();
        //#]
    }
    
    //## operation getProductNumber() 
    public int getProductNumber() {
        //#[ operation getProductNumber() 
        return templateReaction.getProductNumber();
        //#]
    }
    
    //## operation getReactantList() 
    public LinkedList getReactantList() {
        //#[ operation getReactantList() 
        return templateReaction.getReactantList();
        //#]
    }
    
    //## operation getReactantNumber() 
    public int getReactantNumber() {
        //#[ operation getReactantNumber() 
        return templateReaction.getReactantNumber();
        //#]
    }
    
    //## operation isIsomer() 
    public boolean isIsomer() {
        //#[ operation isIsomer() 
        return (type.compareToIgnoreCase("Isomer")==0 || type.compareToIgnoreCase("I")==0);
        //#]
    }
    
    //## operation isNonIncluded() 
    public boolean isNonIncluded() {
        //#[ operation isNonIncluded() 
        return (type.compareToIgnoreCase("NonIncluded")==0 || type.compareToIgnoreCase("N")==0);
        //#]
    }
    
    //## operation isProduct() 
    public boolean isProduct() {
        //#[ operation isProduct() 
        return (type.compareToIgnoreCase("Product")==0 || type.compareToIgnoreCase("P")==0);
        //#]
    }
    
    //## operation isReactant() 
    public boolean isReactant() {
        //#[ operation isReactant() 
        return (type.compareToIgnoreCase("Reactant")==0 || type.compareToIgnoreCase("R")==0);
        //#]
    }
    
    //## operation setTypeAsIsomer() 
    public void setTypeAsIsomer() {
        //#[ operation setTypeAsIsomer() 
        type = "ISOMER";
        //#]
    }
    
    //## operation setTypeAsNonIncluded() 
    public void setTypeAsNonIncluded() {
        //#[ operation setTypeAsNonIncluded() 
        type = "NONINCLUDED";
        //#]
    }
    
    //## operation setTypeAsProduct() 
    public void setTypeAsProduct() {
        //#[ operation setTypeAsProduct() 
        type = "PRODUCT";
        //#]
    }
    
    //## operation setTypeAsReactant() 
    public void setTypeAsReactant() {
        //#[ operation setTypeAsReactant() 
        type = "REACTANT";
        //#]
    }
    
    //## operation toChemDisString() 
    public String toChemDisString() {
        //#[ operation toChemDisString() 
        String s = type;
        if (isNonIncluded()) s = "PRODUCT"; 
        s += '\n';
        for (Iterator iter = getProductList().iterator(); iter.hasNext(); ) {
        	ChemGraph cg = (ChemGraph)iter.next();
        	Species spe = cg.getSpecies();
        	s += spe.toChemDisString() + " + ";
        }
        s = s.substring(0, s.length()-3) + '\n'; 
        
        Kinetics k = getKinetics(); 
        if (k == null) throw new NullPointerException();
        
        s += Double.toString(k.getAValue()) + '\t';	
        s += Double.toString(k.getNValue()) + '\t';	
        s += "0.0\t";	
        s += Double.toString(k.getEValue()) + '\n';
        
        return s;	
          
        //#]
    }
    
    public String getType() {
        return type;
    }
    
    public void setType(String p_type) {
        type = p_type;
    }
    
    public TemplateReaction getTemplateReaction() {
        return templateReaction;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\PDepPathReaction.java
*********************************************************************/


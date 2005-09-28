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



package jing.mathTool;


import java.util.*;

//## package jing::mathTool 

//----------------------------------------------------------------------------
// jing\mathTool\UncertainDouble.java                                                                  
//----------------------------------------------------------------------------

//## class UncertainDouble 
public class UncertainDouble {
    
    protected String type;		//## attribute type 
    
    protected double uncertainty;		//## attribute uncertainty 
    
    protected double value;		//## attribute value 
    
    
    // Constructors
    
    //## operation UncertainDouble(double,double,String) 
    public  UncertainDouble(double p_value, double p_uncertainty, String p_type) {
        //#[ operation UncertainDouble(double,double,String) 
        value = p_value;
        uncertainty = Math.abs(p_uncertainty);
        
        if (p_type.compareToIgnoreCase("Multiplying") == 0 || p_type.compareToIgnoreCase("Multiplier") == 0 || p_type.compareToIgnoreCase("M") == 0) {
        	type = "Multiplier";
        }
        else if (p_type.compareToIgnoreCase("Adding") == 0 || p_type.compareToIgnoreCase("Adder") == 0 || p_type.compareToIgnoreCase("A") == 0) {
        	type = "Adder";
        }
        else throw new InvalidUncertaintyTypeException(p_type);
        //#]
    }
    public  UncertainDouble() {
    }
    
    //## operation getLowerBound() 
    public double getLowerBound() {
        //#[ operation getLowerBound() 
        if (isAddingUncertainty()) {
        	return getValue()-getUncertainty();
        }
        else if (isMultiplyingUncertainty()) {
        	if (getUncertainty() == 0) throw new InvalidUncertaintyException("mutiplier is zero");
        	return getValue()/getUncertainty();
        }
        else throw new InvalidUncertaintyTypeException();
        //#]
    }
    
    //## operation getUpperBound() 
    public double getUpperBound() {
        //#[ operation getUpperBound() 
        if (isAddingUncertainty()) {
        	return getValue()+getUncertainty();
        }
        else if (isMultiplyingUncertainty()) {
        	if (getUncertainty() == 0) throw new InvalidUncertaintyException("mutiplier is zero");
        	return getValue()*getUncertainty();
        }
        else throw new InvalidUncertaintyTypeException();
        //#]
    }
    
    //## operation getValue() 
    public double getValue() {
        //#[ operation getValue() 
        return value;
        //#]
    }
    
    //## operation isAddingUncertainty() 
    public boolean isAddingUncertainty() {
        //#[ operation isAddingUncertainty() 
        return type.equals("Adder");
        //#]
    }
    
    //## operation isMultiplyingUncertainty() 
    public boolean isMultiplyingUncertainty() {
        //#[ operation isMultiplyingUncertainty() 
        return type.equals("Multiplier");
        //#]
    }
    
    //## operation multiply(double) 
    public UncertainDouble multiply(double p_multiplier) {
        //#[ operation multiply(double) 
        return new UncertainDouble(value*p_multiplier,uncertainty*p_multiplier,type);
        
        
        //#]
    }
    
    //## operation plus(double) 
    public UncertainDouble plus(double p_adder) {
        //#[ operation plus(double) 
        return new UncertainDouble(value+p_adder,uncertainty,type);
        
        
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        if (isAddingUncertainty()) {
        	return String.valueOf(value) + "+|-" + String.valueOf(uncertainty);
        }
        else if (isMultiplyingUncertainty()) {
        	return String.valueOf(value) + "*|/" + String.valueOf(uncertainty);
        }
        else throw new InvalidUncertaintyTypeException(type);
        
        
        //#]
    }
    
    public String getType() {
        return type;
    }
    
    public double getUncertainty() {
        return uncertainty;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\mathTool\UncertainDouble.java
*********************************************************************/


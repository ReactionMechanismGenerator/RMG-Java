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



package jing.param;


import java.util.*;

//## package jing::param 

//----------------------------------------------------------------------------
// jing\param\Temperature.java                                                                  
//----------------------------------------------------------------------------

//## class Temperature 
public class Temperature {
    
    /**
    It is the standard temperature in kelvin.
    */
    protected double standard;		//## attribute standard 
    
    
    // Constructors
    
    //## operation Temperature(double,String) 
    public  Temperature(double p_T, String p_unit) {
        //#[ operation Temperature(double,String) 
        if (p_unit.compareToIgnoreCase("K") == 0) {	
        	standard = p_T;
        }
        else if (p_unit.compareToIgnoreCase("C") == 0) {
        	standard = p_T + 273.15;
        }               
        else if (p_unit.compareToIgnoreCase("F") == 0) {
        	standard = (p_T - 32)*5/9 + 273.15;
        }
        else {
        	throw new InvalidUnitException();
        }
        //#]
    }
    public  Temperature() {
    }
    
    //## operation clone() 
    public Object clone() {
        //#[ operation clone() 
        return new Temperature(standard, getStandardUnit());
        //#]
    }
    
    //## operation equals(Object) 
    public boolean equals(Object p_temperature) {
        //#[ operation equals(Object) 
        if (!(p_temperature instanceof Temperature)) return false;
        Temperature t = (Temperature)p_temperature;
        
        return (t.getStandard()==getStandard());
        //#]
    }
    
    //## operation getC() 
    public double getC() {
        //#[ operation getC() 
        return (standard - 273.15);
        //#]
    }
    
    //## operation getF() 
    public double getF() {
        //#[ operation getF() 
        return ((standard-273.15)*9/5+32);
        //#]
    }
    
    //## operation getK() 
    public double getK() {
        //#[ operation getK() 
        return standard;
        //#]
    }
    
    //## operation getStandardUnit() 
    public static String getStandardUnit() {
        //#[ operation getStandardUnit() 
        return "K";
        //#]
    }
    
    //## operation hashCode() 
    public int hashCode() {
        //#[ operation hashCode() 
        return (int)standard;
        //#]
    }
    
    public double getStandard() {
        return standard;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\param\Temperature.java
*********************************************************************/


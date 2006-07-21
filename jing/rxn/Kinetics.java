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
import jing.param.Temperature;
import jing.mathTool.UncertainDouble;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\Kinetics.java                                                                  
//----------------------------------------------------------------------------

/**
This interface defines the functionality of rate calculation.
k=f(T)
Note, this k is k infinity, which is the high pressure limit rate of a reaction.  therefore, it is only a function of temperature. 
*/
//## class Kinetics 
public interface Kinetics {
    
    
    // Argument Temperaturetemperature : 
    /**
    The temperature that the rate is calculated.
    */
    //## operation calculateRate(Temperature,double) 
    double calculateRate(Temperature temperature, double Hrxn);
	
	double calculateRate(Temperature temperature);
    
    //## operation getA() 
    UncertainDouble getA();
    
    //## operation getAValue() 
    double getAValue();
    
    //## operation getComment() 
    String getComment();
    
    //## operation getE() 
    UncertainDouble getE();
    
    //## operation getEValue() 
    double getEValue();
    
    //## operation getN() 
    UncertainDouble getN();
    
    //## operation getNValue() 
    double getNValue();
    
    //## operation getRank() 
    int getRank();
    
    //## operation getSource() 
    String getSource();
	
	void setSource(String p_string);
	
	void setComments(String p_string);
    
    //## operation multiply(double) 
    Kinetics multiply(double p_multiple);
    
    //## operation repOk() 
    boolean repOk();
    
    //## operation toChemkinString() 
    String toChemkinString(double Hrxn, Temperature p_temperature);
    
    //## operation toString() 
    String toString();
	
	boolean equals(Kinetics p_k);
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\Kinetics.java
*********************************************************************/


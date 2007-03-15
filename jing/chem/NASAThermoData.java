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



package jing.chem;


import java.util.*;

import jing.param.GasConstant;
import jing.param.Temperature;

//## package jing::chem 

//----------------------------------------------------------------------------
//jing\chem\NASAThermoData.java                                                                  
//----------------------------------------------------------------------------

//## class NASAThermoData 
public class NASAThermoData {
  
  protected String dataString;		//## attribute dataString 
  protected double [] highTemperatureCoefficients = new double[7];
  protected double [] lowTemperatureCoefficients = new double[7];
  protected double lowTemperature;
  protected double highTemperature;
  protected double middleTemperature;
  
  // Constructors
  
  //## operation NASAThermoData(String) 
  public  NASAThermoData(String p_s) {
      //#[ operation NASAThermoData(String) 
      dataString = p_s;
	  lowTemperature = Double.parseDouble(dataString.substring(48,54));
	  highTemperature = Double.parseDouble(dataString.substring(57,64));
	  middleTemperature = Double.parseDouble(dataString.substring(66,73));
	  
	  String secondLine, thirdLine, fourthLine;
	  StringTokenizer st = new StringTokenizer(dataString,"\n");
	  st.nextToken();
	  secondLine = st.nextToken();
	  thirdLine = st.nextToken();
	  fourthLine = st.nextToken();
	  highTemperatureCoefficients[0] = Double.parseDouble(secondLine.substring(0,15));
	  highTemperatureCoefficients[1] = Double.parseDouble(secondLine.substring(15,30));
	  highTemperatureCoefficients[2] = Double.parseDouble(secondLine.substring(30,45));
	  highTemperatureCoefficients[3] = Double.parseDouble(secondLine.substring(45,60));
	  highTemperatureCoefficients[4] = Double.parseDouble(secondLine.substring(60,75));
	  highTemperatureCoefficients[5] = Double.parseDouble(thirdLine.substring(0,15));
	  highTemperatureCoefficients[6] = Double.parseDouble(thirdLine.substring(15,30));
	  
	  lowTemperatureCoefficients[0] = Double.parseDouble(thirdLine.substring(30,45));
	  lowTemperatureCoefficients[1] = Double.parseDouble(thirdLine.substring(45,60));
	  lowTemperatureCoefficients[2] = Double.parseDouble(thirdLine.substring(60,75));
	  lowTemperatureCoefficients[3] = Double.parseDouble(fourthLine.substring(0,15));
	  lowTemperatureCoefficients[4] = Double.parseDouble(fourthLine.substring(15,30));
	  lowTemperatureCoefficients[5] = Double.parseDouble(fourthLine.substring(30,45));
	  lowTemperatureCoefficients[6] = Double.parseDouble(fourthLine.substring(45,60));
	  
      //#]
  }
  public  NASAThermoData() {
  }
  
  //## operation toString() 
  public String toString() {
      //#[ operation toString() 
      return dataString;
      //#]
  }
  
  public String getDataString() {
      return dataString;
  }
  
  public double calculateEnthalpy(Temperature temp){
	  double T = temp.getK();
	  double enthalpy;
	  if (T < 298) throw new TemperatureOutOfRangeException();
	  else if (T < middleTemperature){
		  double [] a = new double[7];
		  a=lowTemperatureCoefficients;
		  enthalpy = a[0] + a[1]*T/2 + a[2]*T*T/3 + a[3]*T*T*T/4 + a[4]*T*T*T*T/5 + a[5]/T;
		  enthalpy = enthalpy * T * GasConstant.getCalMolK()/1000;
		  return enthalpy;
	  }
	  else if (T < highTemperature){
		  double [] a = new double[7];
		  a=highTemperatureCoefficients;
		  enthalpy = a[0] + a[1]*T/2 + a[2]*T*T/3 + a[3]*T*T*T/4 + a[4]*T*T*T*T/5 + a[5]/T;
		  enthalpy = enthalpy * T * GasConstant.getCalMolK()/1000;
		  return enthalpy;
	  }
	  else throw new TemperatureOutOfRangeException();
		  
		  
  }
  
  public double calculateEntropy(Temperature temp){
	  double T = temp.getK();
	  double entropy;
	  if (T < 298) throw new TemperatureOutOfRangeException();
	  else if (T < middleTemperature){
		  double [] a = new double[7];
		  a=lowTemperatureCoefficients;
		  entropy = a[0]*Math.log(T) + a[1]*T + a[2]*T*T/2 + a[3]*T*T*T/3 + a[4]*T*T*T*T/4 + a[6];
		  entropy = entropy * GasConstant.getCalMolK();
		  return entropy;
	  }
	  else if (T < highTemperature){
		  double [] a = new double[7];
		  a=highTemperatureCoefficients;
		  entropy = a[0]*Math.log(T) + a[1]*T + a[2]*T*T/2 + a[3]*T*T*T/3 + a[4]*T*T*T*T/4 + a[6];
		  entropy = entropy * GasConstant.getCalMolK();
		  return entropy;
	  }
	  else throw new TemperatureOutOfRangeException();
		  
		  
  }
  
  public double calculateCp(Temperature temp){
	  double T = temp.getK();
	  double Cp;
	  if (T < 298) throw new TemperatureOutOfRangeException();
	  else if (T < middleTemperature){
		  double [] a = new double[7];
		  a=lowTemperatureCoefficients;
		  Cp = a[0] + a[1]*T + a[2]*T*T + a[3]*T*T*T + a[4]*T*T*T*T;
		  Cp = Cp * T;
		  return Cp;
	  }
	  else if (T < highTemperature){
		  double [] a = new double[7];
		  a=highTemperatureCoefficients;
		  Cp = a[0] + a[1]*T + a[2]*T*T + a[3]*T*T*T + a[4]*T*T*T*T;
		  Cp = Cp * T;
		  return Cp;
	  }
	  else throw new TemperatureOutOfRangeException();
		  
		  
  }
  
//## operation calculateG(Temperature)
  public double calculateFreeEnergy(Temperature p_temperature) {
      //#[ operation calculateG(Temperature)
      double T = p_temperature.getK();

      return (calculateEnthalpy(p_temperature)*1000 - T*calculateEntropy(p_temperature))/1000;
      //#]
  }
  
  public void setDataString(String p_dataString) {
      dataString = p_dataString;
  }
  
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\NASAThermoData.java
*********************************************************************/


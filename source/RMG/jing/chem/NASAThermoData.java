////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////



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

	  /*
	   * Updated by MRH on 22-Jun-2009
	   * 	Code assumed dataString would contain 4 (and only 4) lines of code.
	   * 	For species containing more than 4 unique elements, a fifth line
	   * 	will now be present.  Here are the old and new formats (roughly):
	   * 
	   * Old:
	   * Name	C x H y O z G Tmin Tmax Tint	1
	   * a1 a2 a3 a4 a5							2
	   * a6 a7 a1 a2 a3							3
	   * a4 a5 a6 a7							4
	   * 
	   * New:
	   * Name				G Tmin Tmax Tint	1&
	   * C x H y O z Si a S b
	   * a1 a2 a3 a4 a5							2
	   * a6 a7 a1 a2 a3							3
	   * a4 a5 a6 a7							4
	   * 
	   * The ampersand at the end of the first line indicates whether the
	   * 	dataString is of the old or new format.  The Tmin, Tmax, and
	   * 	Tint are in the same locations, as are the high- and 
	   * 	low-temperature coefficients.
	   * 
	   * NOTE: The new format is not recognized by Chemkin-v.2 but is
	   * 	recognized by Chemkin-v.4
	   */
	  String firstLine, secondLine, thirdLine, fourthLine, fifthLine;
	  StringTokenizer st = new StringTokenizer(dataString,"\n");
	  firstLine = st.nextToken();
	  secondLine = st.nextToken();
	  thirdLine = st.nextToken();
	  fourthLine = st.nextToken();
	  // If there are more tokens, we have a fifth line and are therefore
	  //	using the new format
	  if (st.hasMoreTokens()) {
		  fifthLine = st.nextToken();
		  
		  highTemperatureCoefficients[0] = Double.parseDouble(thirdLine.substring(0,15));
		  highTemperatureCoefficients[1] = Double.parseDouble(thirdLine.substring(15,30));
		  highTemperatureCoefficients[2] = Double.parseDouble(thirdLine.substring(30,45));
		  highTemperatureCoefficients[3] = Double.parseDouble(thirdLine.substring(45,60));
		  highTemperatureCoefficients[4] = Double.parseDouble(thirdLine.substring(60,75));
		  highTemperatureCoefficients[5] = Double.parseDouble(fourthLine.substring(0,15));
		  highTemperatureCoefficients[6] = Double.parseDouble(fourthLine.substring(15,30));
		  
		  lowTemperatureCoefficients[0] = Double.parseDouble(fourthLine.substring(30,45));
		  lowTemperatureCoefficients[1] = Double.parseDouble(fourthLine.substring(45,60));
		  lowTemperatureCoefficients[2] = Double.parseDouble(fourthLine.substring(60,75));
		  lowTemperatureCoefficients[3] = Double.parseDouble(fifthLine.substring(0,15));
		  lowTemperatureCoefficients[4] = Double.parseDouble(fifthLine.substring(15,30));
		  lowTemperatureCoefficients[5] = Double.parseDouble(fifthLine.substring(30,45));
		  lowTemperatureCoefficients[6] = Double.parseDouble(fifthLine.substring(45,60));
	  } else {
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
	  }
	  
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


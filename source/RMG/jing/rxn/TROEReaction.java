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



package jing.rxn;


import java.util.*;
import jing.param.*;
import jing.mathTool.*;
import jing.rxnSys.SystemSnapshot;

//## package jing::rxn 

//----------------------------------------------------------------------------
//jing\rxn\TROEReaction.java                                                                  
//----------------------------------------------------------------------------

//## class TROEReaction 
public class TROEReaction extends ThirdBodyReaction {
  
  protected double T2star;		//## attribute T2star 
  
  protected double T3star;		//## attribute T3star 
  
  protected double Tstar;		//## attribute Tstar 
  
  protected double a;		//## attribute a 
  
  protected ArrheniusKinetics low;		//## attribute low 
  
  protected boolean troe7 = false;		//## attribute troe7 
  
  
  // Constructors
  
  //## operation TROEReaction() 
  private  TROEReaction() {
      //#[ operation TROEReaction() 
      //#]
  }
  
  //## operation calculateRate(SystemSnapshot) 
  public double calculateRate(SystemSnapshot p_presentStatus) {
      //#[ operation calculateRate(SystemSnapshot) 
      Temperature temp = p_presentStatus.getTemperature();
      
      //Kinetics k = getKinetics();
      //double rate = k.calculateRate(temp, calculateHrxn(temp));
      double rate = super.calculateTotalRate(temp);
      rate *= calculateTroeFallOff(p_presentStatus);
      
      return rate;
      //#]
  }
  
  //## operation calculateTroeFallOff(SystemSnapshot) 
  public double calculateTroeFallOff(SystemSnapshot p_presentStatus) {
      //#[ operation calculateTroeFallOff(SystemSnapshot) 
      Temperature temp = p_presentStatus.getTemperature();
      
      double M = calculateThirdBodyCoefficient(p_presentStatus);
      double kZero = low.calculateRate(temp,-1);
      double kInf = 0.0;
      for (int i=0; i<getKinetics().length; i++) {
    	  kInf += getKinetics()[i].calculateRate(temp,-1);
      }
      
      double Pr = kZero*M/kInf;
      double T = temp.getK();
      double Fc;
      if (troe7) // 7-parmater Troe form
      	Fc = (1.0-a)*Math.exp(-T/T3star) + a*Math.exp(-T/Tstar) + Math.exp(-T2star/T); //note: this eqn is wrong in CK-III manual (1996)
      else // 6-paramter Troe form
      	Fc = (1.0-a)*Math.exp(-T/T3star) + a*Math.exp(-T/Tstar);
      
      double small = 1.0e-30;
      double logFc = MathTool.log10(Math.max(Fc,small));
      double logPr = MathTool.log10(Math.max(Pr,small));
      double c = -0.4 - 0.67*logFc;
      double n = 0.75 - 1.27*logFc; // note: this eqn is wrong in the CK-III manual (1996)
      
      double inside = (logPr + c)/(n-0.14*(logPr + c)); 
      double logF = logFc/(1.0 + Math.pow(inside,2));
      double F = Math.pow(10.0,logF);
      
      double fallOffFactor = (Pr/(1.0+Pr))*F;
      
      return fallOffFactor;
      
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
  
  //## operation generateReverseReaction() 
  public void generateReverseReaction() {
      //#[ operation generateReverseReaction() 
      TROEReaction r = new TROEReaction();
      r.structure = getStructure().generateReverseStructure();
      r.kinetics = getKinetics();
      for (int i=0; i<r.kinetics.length; i++) {
    	  r.comments = "Reverse reaction";
      }
      r.weightMap = weightMap;
	  r.low = low;
      r.a = a;
      r.Tstar = Tstar;
      r.T3star = T3star;
      r.troe7 = troe7;
      r.T2star = T2star;
	  
      r.setReverseReaction(this);
      this.setReverseReaction(r);
      
      //this.setReverseReaction(null);
      
      return;
      //#]
  }
  
  //## operation make(Reaction,HashMap,ArrheniusKinetics,double,double,double,boolean,double) 
  public static TROEReaction make(Reaction p_reaction, HashMap p_weightMap, final ArrheniusKinetics p_low, double p_a, double p_T3star, double p_Tstar, boolean p_troe7, double p_T2star) {
      //#[ operation make(Reaction,HashMap,ArrheniusKinetics,double,double,double,boolean,double) 
      TROEReaction tr = new TROEReaction();
      tr.structure = p_reaction.getStructure();
      tr.kinetics = p_reaction.getKinetics();
      tr.comments = p_reaction.getComments();
      // generate reverse rxn
      
      
      tr.weightMap = p_weightMap;
      tr.low = p_low;
      tr.a = p_a;
      tr.Tstar = p_Tstar;
      tr.T3star = p_T3star;
      tr.troe7 = p_troe7;
      tr.T2star = p_T2star;
      
	  tr.generateReverseReaction();
	  
      return tr;
      //#]
  }
  
  //## operation toChemkinString() 
  public String toChemkinString(Temperature p_temperature) {
      //#[ operation toChemkinString() 
      String s = super.toChemkinString(p_temperature)+'\n';
      
      // write pressure-dependence parameters
      s += "LOW/" + low.toChemkinString(calculateHrxn(p_temperature),p_temperature, false) + "/\n";
      s += "TROE/" + a + '\t' + T3star + '\t' + Tstar;
      if(troe7) s += "\t" + T2star;
      s = s + "/\n";
       
      return s;
      
      //#]
  }
  
  public String toRestartString(Temperature p_temperature) {
	  String s = super.toRestartString(p_temperature) + "\n";
      
      // write pressure-dependence parameters
      s += "LOW/" + low.toChemkinString(calculateHrxn(p_temperature),p_temperature, false) + "/\n";
      s += "TROE/" + a + '\t' + T3star + '\t' + Tstar;
      if(troe7) s += "\t" + T2star;
      s = s + "/\n";
       
      return s;
  }
  
  public String toChemkinString(Temperature p_temperature, Pressure p_pressure) {
	  String s = super.toChemkinString(p_temperature) + "\n";
	  s += "LOW/ " + low.toChemkinString(calculateHrxn(p_temperature), p_temperature, false) + "/\n";
	  s += "TROE/ " + a + "\t" + T3star + "\t" + Tstar;
	  if (troe7) s += "\t" + T2star;
	  s += " /\n";
	  return s;
  }
  
  //## operation toString() 
  public String toString(Temperature p_temperature) {
      //#[ operation toString() 
      String s = getStructure().toChemkinString(true).toString() + '\n';
      for (int i=0; i<getKinetics().length; i++) {
    	  s += "kInf = " + getKinetics()[i].toChemkinString(calculateHrxn(p_temperature),p_temperature, false) + '\n';
      }
      s += "kZero = " + low.toChemkinString(calculateHrxn(p_temperature),p_temperature, false) + '\n';
      s += "a = " + a + '\t' + "T*** = " + T3star + '\t' + "T* = " + Tstar + '\t';
      if(troe7) s += "T** = " + T2star;
      s += "\n";
      
      return s;
      
      //#]
  }
  
  public double getT2star() {
      return T2star;
  }
  
  public void setT2star(double p_T2star) {
      T2star = p_T2star;
  }
  
  public double getT3star() {
      return T3star;
  }
  
  public void setT3star(double p_T3star) {
      T3star = p_T3star;
  }
  
  public double getTstar() {
      return Tstar;
  }
  
  public void setTstar(double p_Tstar) {
      Tstar = p_Tstar;
  }
  
  public double getA() {
      return a;
  }
  
  public void setA(double p_a) {
      a = p_a;
  }
  
  public ArrheniusKinetics getLow() {
      return low;
  }
  
  public void setLow(ArrheniusKinetics p_low) {
      low = p_low;
  }
  
  public boolean getTroe7() {
      return troe7;
  }
  
  public void setTroe7(boolean p_troe7) {
      troe7 = p_troe7;
  }
  
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\TROEReaction.java
*********************************************************************/


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
      
      Kinetics k = getKinetics();
      double rate = k.calculateRate(temp, calculateHrxn(temp));
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
      double kInf = getKinetics().calculateRate(temp,-1);
      
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
  
  //## operation make(Reaction,HashMap,ArrheniusKinetics,double,double,double,boolean,double) 
  public static TROEReaction make(Reaction p_reaction, HashMap p_weightMap, final ArrheniusKinetics p_low, double p_a, double p_T3star, double p_Tstar, boolean p_troe7, double p_T2star) {
      //#[ operation make(Reaction,HashMap,ArrheniusKinetics,double,double,double,boolean,double) 
      TROEReaction tr = new TROEReaction();
      tr.structure = p_reaction.getStructure();
      tr.rateConstant = p_reaction.getRateConstant();
      tr.comments = p_reaction.getComments();
      // generate reverse rxn
      tr.generateReverseReaction();
      
      tr.weightMap = p_weightMap;
      tr.low = p_low;
      tr.a = p_a;
      tr.Tstar = p_Tstar;
      tr.T3star = p_T3star;
      tr.troe7 = p_troe7;
      tr.T2star = p_T2star;
      
      return tr;
      //#]
  }
  
  //## operation toChemkinString() 
  public String toChemkinString() {
      //#[ operation toChemkinString() 
      String s = super.toChemkinString()+'\n';
      
      // write pressure-dependence parameters
      s += "LOW/" + low.toChemkinString() + "/\n";
      s += "TROE/" + a + '\t' + T3star + '\t' + Tstar;
      if(troe7) s += "\t" + T2star;
      s = s + "/\n";
       
      return s;
      
      //#]
  }
  
  //## operation toString() 
  public String toString() {
      //#[ operation toString() 
      String s = getStructure().toChemkinString(true) + '\n';
      s += "kInf = " + getKinetics().toChemkinString() + '\n';
      s += "kZero = " + low.toChemkinString() + '\n';
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


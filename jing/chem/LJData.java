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
import jing.param.*;
import jing.param.Temperature;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\LJData.java
//----------------------------------------------------------------------------

/**
Mutable objects.  Becareful when use this.plus() or this.minus(), this.sum(), since this will be changed after such operations.
*/
//## class LJData
public class LJData extends LJGroupData {
    protected int na;
//this class is analogous to ThermoData, and was based off of it
    
    // Constructors

    //## operation LJData()
    public  LJData() {
        //#[ operation LJData()
        super();
        //#]
    }
    //## operation LJData(double,double,double,double,double,double,double,double,double,String)
    public  LJData(double dTc, double dPc, double dVc, double dTb,String p_comments,int p_na) {
        //#[ operation LJData(double,double,double,double,double,double,double,double,double,String)
        super(dTc, dPc, dVc, dTb, p_comments);
        na = p_na;



        //#]
    }
    //## operation LJData(ThermoGAValue)
    public  LJData(LJGroupData p_ga, int p_na) {
        //#[ operation LJData(ThermoGAValue)
        super(p_ga);
        na = p_na;
        //#]
    }

    //## operation LJData(String, LJData, String)
    //svp
    public LJData(String p_name, LJData p_td, String p_comments){
      //#[ operation LJData(String, LJData, String)
      super(p_name, p_td, p_comments);
      //#]
    }
    
    public LJData(String p_name, LJData p_td, String p_comments, String p_source, int p_na) {
    	super(p_name, p_td, p_comments, p_source);
        na = p_na;
    }



    //## operation copy()
    public LJData copy() {
        //#[ operation copy()
        if (this == null) return null;

        return new LJData(dTc,dPc,dVc,dTb,comments,na);

        //#]
    }


//    //## operation minus(ThermoGAValue)
//    public void minus(ThermoGAValue p_thermoData) {
//        //#[ operation minus(ThermoGAValue)
//        H298 -= p_thermoData.H298;
//        S298 -= p_thermoData.S298;
//        Cp300 -= p_thermoData.Cp300;
//        Cp400 -= p_thermoData.Cp400;
//        Cp500 -= p_thermoData.Cp500;
//        Cp600 -= p_thermoData.Cp600;
//        Cp800 -= p_thermoData.Cp800;
//        Cp1000 -= p_thermoData.Cp1000;
//        Cp1500 -= p_thermoData.Cp1500;
//        dH -= Math.pow(p_thermoData.dH,2);//svp
//        dS -= Math.pow(p_thermoData.dS,2);//svp
//        dCp -= Math.pow(p_thermoData.dCp,2);//svp
//
//
//
//
//        //#]
//    }
//
//    //## operation multiply(int)
//    public void multiply(int p_multiplier) {
//        //#[ operation multiply(int)
//        if (p_multiplier == 1) return;
//
//        H298 *= p_multiplier;
//        S298 *= p_multiplier;
//        Cp300 *= p_multiplier;
//        Cp400 *= p_multiplier;
//        Cp500 *= p_multiplier;
//        Cp600 *= p_multiplier;
//        Cp800 *= p_multiplier;
//        Cp1000 *= p_multiplier;
//        Cp1500 *= p_multiplier;
//        dH *= p_multiplier;//svp
//        dS *= p_multiplier;//svp
//        dCp *= p_multiplier;//svp
//
//
//        name = name + "redundancy = " + p_multiplier + '\n';
//        //#]
//    }

    //## operation plus(ThermoGAValue)
    public void plus(LJGroupData p_LJData) {
        //#[ operation plus(ThermoGAValue)
        if (p_LJData == null) return;

        dTc += p_LJData.dTc;
        dPc += p_LJData.dPc;
        dVc += p_LJData.dVc;
        dTb += p_LJData.dTb;

        if (p_LJData.getName() != null) {
                if (name == null) name = p_LJData.getName() + '\n';
                else name = name + p_LJData.getName() + '\n';
        }




        //#]
    }

    //## operation toString()
    public String toString() {
        //#[ operation toString()
        return super.toString();
        //#]
    }
    
    //calculates Tc in K
    public double calculateTc() {
        return this.calculateTb()/(0.584+0.965*dTc-dTc*dTc);
    }

    //calculates Pc in bar
    public double calculatePc() {
        return 1.0/((0.113+0.0032*na-dPc)*(0.113+0.0032*na-dPc));
    }
    
    //calculates Vc in cc/mol
    public double calculateVc() {
        return 17.5 + dVc;
    }
    
    //calculates Tb in K
    public double calculateTb() {
        return 198.0 + dTb;
    }
    
    //calculates omega (accentric factor)
    public double calculateOmega(){
        double f = calculateTb()/calculateTc();
        return (-1.0*Math.log(calculatePc()/1.01325) -5.927 + 6.096/f + 1.289*Math.log(f) - 0.169*Math.pow(f, 6) ) / ( 15.252 - 15.688/f - 13.472*Math.log(f) + 0.436*Math.pow(f,6));
    }

    //calculates sigma in angstroms
    public double calculateSigma(){
        return (2.3442*Math.exp(0.1303*calculateOmega())) * Math.pow(calculateTc()*1.01325/calculatePc(),1.0/3.0);
    }

    //calculates epsilon in K
    public double calculateEpsilon(){
        return 0.8109*Math.exp(-0.6228*calculateOmega())*calculateTc();
    }

}
/*********************************************************************
        File Path	: RMG\RMG\jing\chem\LJData.java
*********************************************************************/


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

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\ThermoGAValue.java
//----------------------------------------------------------------------------

/**
Immutable data holds all the benson's group value.
*/
//## class ThermoGAValue
public class ThermoGAValue {

    protected double Cp1000 = 0;		//## attribute Cp1000

    protected double Cp1500 = 0;		//## attribute Cp1500

    protected double Cp300 = 0;		//## attribute Cp300

    protected double Cp400 = 0;		//## attribute Cp400

    protected double Cp500 = 0;		//## attribute Cp500

    protected double Cp600 = 0;		//## attribute Cp600

    protected double Cp800 = 0;		//## attribute Cp800

    protected double H298 = 0;		//## attribute H298

    protected double S298 = 0;		//## attribute S298

    protected double dH = 0;//svp

    protected double dS = 0;//svp

    protected double dCp = 0;//svp

    protected static double T_HIGH = 1500;		//## attribute T_HIGH

    protected static double T_LOW = 300;		//## attribute T_LOW

    protected String comments = null;		//## attribute comments

    protected String name = null;		//## attribute name
    
    protected String source = null;


    // Constructors

    //## operation ThermoGAValue()
    public  ThermoGAValue() {
        //#[ operation ThermoGAValue()
        H298 = 0;
        S298 = 0;
        Cp300 = 0;
        Cp400 = 0;
        Cp500 = 0;
        Cp600 = 0;
        Cp800 = 0;
        Cp1000 = 0;
        Cp1500 = 0;
        dH = 0;
        dS = 0;
        dCp = 0;



        //#]
    }
    //## operation ThermoGAValue(double,double,double,double,double,double,double,double,double,String)
    public  ThermoGAValue(double p_H298, double p_S298, double p_Cp300, double p_Cp400, double p_Cp500, double p_Cp600, double p_Cp800, double p_Cp1000, double p_Cp1500, double p_dH,double p_dS,double p_dCp,String p_comments) {
        //#[ operation ThermoGAValue(double,double,double,double,double,double,double,double,double,String)
        H298 = p_H298;
        S298 = p_S298;
        Cp300 = p_Cp300;
        Cp400 = p_Cp400;
        Cp500 = p_Cp500;
        Cp600 = p_Cp600;
        Cp800 = p_Cp800;
        Cp1000 = p_Cp1000;
        Cp1500 = p_Cp1500;
        dH = p_dH;
        dS = p_dS;
        dCp = p_dCp;
        comments = p_comments;



        //#]
    }
    //## operation ThermoGAValue(String,ThermoGAValue,String)
    public  ThermoGAValue(String p_name, ThermoGAValue p_ga, String p_comments) {
        //#[ operation ThermoGAValue(String,ThermoGAValue,String)
        H298 = p_ga.H298;
        S298 = p_ga.S298;
        Cp300 = p_ga.Cp300;
        Cp400 = p_ga.Cp400;
        Cp500 = p_ga.Cp500;
        Cp600 = p_ga.Cp600;
        Cp800 = p_ga.Cp800;
        Cp1000 = p_ga.Cp1000;
        Cp1500 = p_ga.Cp1500;
        dH = p_ga.dH;
        dS = p_ga.dS;
        dCp = p_ga.dCp;
        comments = p_comments;
        name = p_name;



        //#]
    }
    
    public ThermoGAValue(String p_name, ThermoGAValue p_ga, String p_comments, String p_source) {
        H298 = p_ga.H298;
        S298 = p_ga.S298;
        Cp300 = p_ga.Cp300;
        Cp400 = p_ga.Cp400;
        Cp500 = p_ga.Cp500;
        Cp600 = p_ga.Cp600;
        Cp800 = p_ga.Cp800;
        Cp1000 = p_ga.Cp1000;
        Cp1500 = p_ga.Cp1500;
        dH = p_ga.dH;
        dS = p_ga.dS;
        dCp = p_ga.dCp;
        comments = p_comments;
        name = p_name;
        source = p_source;

	}    
    
    //## operation ThermoGAValue(ThermoGAValue)
    public  ThermoGAValue(ThermoGAValue p_ga) {
        //#[ operation ThermoGAValue(ThermoGAValue)
        H298 = p_ga.H298;
        S298 = p_ga.S298;
        Cp300 = p_ga.Cp300;
        Cp400 = p_ga.Cp400;
        Cp500 = p_ga.Cp500;
        Cp600 = p_ga.Cp600;
        Cp800 = p_ga.Cp800;
        Cp1000 = p_ga.Cp1000;
        Cp1500 = p_ga.Cp1500;
        dH = p_ga.dH;
        dS = p_ga.dS;
        dCp = p_ga.dCp;
        comments = p_ga.comments;
        name = p_ga.name;


        //#]
    }

    //## operation toString()
    public String toString() {
        //#[ operation toString()
        String s = "";
        s = s + String.valueOf(H298) + '\t';
        s = s + String.valueOf(S298) + '\t';
        s = s + String.valueOf(Cp300) + '\t';
        s = s + String.valueOf(Cp400) + '\t';
        s = s + String.valueOf(Cp500) + '\t';
        s = s + String.valueOf(Cp600) + '\t';
        s = s + String.valueOf(Cp800) + '\t';
        s = s + String.valueOf(Cp1000) + '\t';
        s = s + String.valueOf(Cp1500);

        return s;
        //#]
    }

    protected double getCp1000() {
        return Cp1000;
    }

    protected double getCp1500() {
        return Cp1500;
    }

    protected double getCp300() {
        return Cp300;
    }

    protected double getCp400() {
        return Cp400;
    }

    protected double getCp500() {
        return Cp500;
    }

    protected double getCp600() {
        return Cp600;
    }

    protected double getCp800() {
        return Cp800;
    }

    protected double getH298() {
        return H298;
    }

    protected double getS298() {
        return S298;
    }

    //svp
    protected double getDH() {
      return Math.pow(dH,0.5);
    }

    //svp
    protected double getDS(){
      return Math.pow(dS,0.5);
    }

    //svp
    protected double getDCp() {
      return Math.pow(dCp,0.5);
    }


    private static double getT_HIGH() {
        return T_HIGH;
    }

    private static void setT_HIGH(double p_T_HIGH) {
        T_HIGH = p_T_HIGH;
    }

    private static double getT_LOW() {
        return T_LOW;
    }

    private static void setT_LOW(double p_T_LOW) {
        T_LOW = p_T_LOW;
    }

    public String getComments() {
        return comments;
    }

    public String getName() {
        return name;
    }

}
/*********************************************************************
        File Path	: RMG\RMG\jing\chem\ThermoGAValue.java
*********************************************************************/


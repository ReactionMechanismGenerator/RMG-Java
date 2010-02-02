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
//import jing.chem.*;


import java.util.*;

//## package jing::chem

//----------------------------------------------------------------------------
// based on jing\chem\ThermoGAValue.java
//----------------------------------------------------------------------------

/**
Immutable data holds Joback group values.
*/
//## class ThermoGAValue
public class LJGroupData {

    protected double dTc = 0;		
    protected double dPc = 0;  
    protected double dVc = 0;
    protected double dTb = 0;


    protected String comments = null;		//## attribute comments

    protected String name = null;		//## attribute name
    
    protected String source = null;


    // Constructors

    //## operation ThermoGAValue()
    public  LJGroupData() {
        dTc = 0;		
        dPc = 0;  
        dVc = 0;
        dTb = 0;



        //#]
    }
    //## operation ThermoGAValue(double,double,double,double,double,double,double,double,double,String)
    public  LJGroupData(double p_dTc, double p_dPc, double p_dVc, double p_dTb,String p_comments) {
        //#[ operation ThermoGAValue(double,double,double,double,double,double,double,double,double,String)
        dTc = p_dTc;
        dPc = p_dPc;
        dVc = p_dVc;
        dTb = p_dTb;
        comments = p_comments;



        //#]
    }
    //## operation ThermoGAValue(String,ThermoGAValue,String)
    public  LJGroupData(String p_name, LJGroupData p_ga, String p_comments) {
        //#[ operation ThermoGAValue(String,ThermoGAValue,String)
        dTc = p_ga.dTc;
        dPc = p_ga.dPc;
        dVc = p_ga.dVc;
        dTb = p_ga.dTb;
        comments = p_comments;
        name = p_name;



        //#]
    }
    
    public LJGroupData(String p_name, LJGroupData p_ga, String p_comments, String p_source) {
        dTc = p_ga.dTc;
        dPc = p_ga.dPc;
        dVc = p_ga.dVc;
        dTb = p_ga.dTb;
        comments = p_comments;
        name = p_name;
        source = p_source;

	}    
    
    //## operation ThermoGAValue(ThermoGAValue)
    public  LJGroupData(LJGroupData p_ga) {
        //#[ operation ThermoGAValue(ThermoGAValue)
        dTc = p_ga.dTc;
        dPc = p_ga.dPc;
        dVc = p_ga.dVc;
        dTb = p_ga.dTb;
        comments = p_ga.comments;
        name = p_ga.name;


        //#]
    }

    //## operation toString()
    public String toString() {
        //#[ operation toString()
        String s = "";
        s = s + String.valueOf(dTc) + '\t';
        s = s + String.valueOf(dPc) + '\t';
        s = s + String.valueOf(dVc) + '\t';
        s = s + String.valueOf(dTb);

        return s;
        //#]
    }

    protected double getdTc() {
        return dTc;
    }

    protected double getdPc() {
        return dPc;
    }

    protected double getdVc() {
        return dVc;
    }

    protected double getdTb() {
        return dTb;
    }

    public String getComments() {
        return comments;
    }

    public String getName() {
        return name;
    }

}
/*********************************************************************
        File Path	: RMG\RMG\LJGroupData.java
*********************************************************************/
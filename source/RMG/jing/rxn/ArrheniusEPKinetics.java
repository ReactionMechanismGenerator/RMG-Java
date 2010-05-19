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
import jing.mathTool.UncertainDouble;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ArrheniusEPKinetics.java                                                                  
//----------------------------------------------------------------------------

//## class ArrheniusEPKinetics 
public class ArrheniusEPKinetics extends ArrheniusKinetics {
    
    protected UncertainDouble alpha = new UncertainDouble(0,0,"Adder");		//## attribute alpha 
    
    
    // Constructors
    
    //## operation ArrheniusEPKinetics(UncertainDouble,UncertainDouble,UncertainDouble,UncertainDouble,String,int,String,String) 
    public  ArrheniusEPKinetics(UncertainDouble p_A, UncertainDouble p_n, UncertainDouble p_alpha, UncertainDouble p_E, String p_TRange, int p_rank, String p_source, String p_comment) {
        //#[ operation ArrheniusEPKinetics(UncertainDouble,UncertainDouble,UncertainDouble,UncertainDouble,String,int,String,String) 
        super(p_A, p_n, p_E, p_TRange, p_rank, p_source, p_comment);
        alpha = p_alpha;
        
        //#]
    }
    public  ArrheniusEPKinetics() {
    }
	
	///////////////////////////////////////
	// Notice that we don't redefine getEValue() here but instead inherit it from ArrheniusKinetics
	// which means that if we imagine 
	//   Ea = Eo + alpha * Hrxn
	// then getEValue() returns Eo NOT Ea
	// but we DO redefine toChemkinString to return the Ea
	///////////////////////////////////////
	
	public String toChemkinString(double p_Hrxn, Temperature p_temperature, boolean includeComments){
		double Ea = getEaValue(p_Hrxn);
		// If reported Arrhenius Ea value was computed using Evans-Polanyi relationship,
		//	inform user (in chem.inp file) of what deltaHrxn(T) was used.
		if ((int)alpha.getValue() != 0) {
			comment += " (Ea computed using deltaHrxn(T=298K) of " + p_Hrxn + " kcal/mol)";
		}
		Object [] formatString = new Object[5];
		
//		formatString[0] = new Double(getAValue());
    	double tempDouble = new Double(getAValue());
    	if (AUnits.equals("moles")) {
    		formatString[0] = tempDouble;
    	} else if (AUnits.equals("molecules")) {
    		formatString[0] = tempDouble / 6.022e23;
    	}
				
		formatString[1] = new Double(getNValue());
		
		if (EaUnits.equals("kcal/mol")) {
			formatString[2] = Ea;
		}
		else if (EaUnits.equals("cal/mol")) {
			formatString[2] = Ea * 1000.0;
		}
		else if (EaUnits.equals("kJ/mol")) {
			formatString[2] = Ea * 4.184;
		}
		else if (EaUnits.equals("J/mol")) {
			formatString[2] = Ea * 4184.0;
		}
		else if (EaUnits.equals("Kelvins")) {
			formatString[2] = Ea / 1.987e-3;
		}
		
		formatString[3] = source; formatString[4] = comment;
		
		if (includeComments){
			return String.format("%1.7e \t %2.5f \t %3.5f \t !%s  %s", formatString);
		}
		else
			return String.format("%1.7e \t %2.5f \t %3.5f \t ", formatString);
		
		//return String.valueOf(getAValue()) + '\t' + String.valueOf(getNValue()) + '\t' + Ea + "\t!" + source + " "+comment;
        
	}
	
    //## operation calculateRate(Temperature,double) 
    public double calculateRate(Temperature p_temperature, double p_Hrxn) {
        //#[ operation calculateRate(Temperature,double) 
        double T = p_temperature.getStandard();
        double R = GasConstant.getKcalMolK();
        
        double Ea = E.getValue() + alpha.getValue()*p_Hrxn;
        //if (Ea<0) throw new NegativeEnergyBarrierException();
        double rate = A.getValue() * Math.pow(T, n.getValue()) * Math.exp(-Ea/R/T);
        return rate;
        
        //#]
    }
    
    //## operation getAlphaUncertainty() 
    public double getAlphaUncertainty() {
        //#[ operation getAlphaUncertainty() 
        return alpha.getUncertainty();
        //#]
    }
    
    //## operation getAlphaValue() 
    public double getAlphaValue() {
        //#[ operation getAlphaValue() 
        return alpha.getValue();
        //#]
    }
    
    //## operation multiply(double) 
    public Kinetics multiply(double p_multiple) {
        //#[ operation multiply(double) 
        UncertainDouble newA = getA().multiply(p_multiple);
        Kinetics newK = new ArrheniusEPKinetics(newA,getN(),getAlpha(),getE(),getTRange(),getRank(),getSource(),getComment());
        return newK;
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String string = "A = " + A.toString() + "; n = " + n.toString() + "; E = " + E.toString() + ";" + " alpha = " + alpha.toString() + '\n';
        string = string + "Source: " + source + '\n';
        string = string + "Comments: " + comment + '\n';
        return string;
        //#]
    }
    
    public UncertainDouble getAlpha() {
        return alpha;
    }
    //gmagoon 05/19/10: it looks like in rest of code, Hrxn at 298 K is used, so this temperature should be used elsewhere for consistency
    public double getEaValue(double p_Hrxn) {
        return E.getValue()+alpha.getValue()*p_Hrxn;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ArrheniusEPKinetics.java
*********************************************************************/


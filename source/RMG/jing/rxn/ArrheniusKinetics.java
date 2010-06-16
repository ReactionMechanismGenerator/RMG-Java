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


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import jing.param.*;
import jing.chem.GATPFitException;
import jing.mathTool.*;
import jing.param.Temperature;
import jing.mathTool.UncertainDouble;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ArrheniusKinetics.java                                                                  
//----------------------------------------------------------------------------

//## class ArrheniusKinetics 
public class ArrheniusKinetics implements Kinetics {
    
    protected UncertainDouble A;		//## attribute A 
    
    protected UncertainDouble E;		//## attribute E 
    
    protected String TRange;		//## attribute TRange 
    
    protected String comment;		//## attribute comment 
    
    protected UncertainDouble n;		//## attribute n 
    
    protected int rank;		//## attribute rank 
    
    protected String source;		//## attribute source 
    /*
     * 29Jun2009-MRH: Added fromPrimaryKineticLibrary attribute
     * 	When RMG computes total rate or ask for kinetics for a given rxn,
     * 		it first checks if the rxn/structure of interest has kinetics
     * 		from a PRL.  If so, it will use those numbers.
     */
    protected boolean fromPrimaryKineticLibrary;
    
    protected static boolean verbose = false;
    
    protected static String EaUnits;
    
    protected static String AUnits;
    
    // Constructors
    
    //## operation ArrheniusKinetics(UncertainDouble,UncertainDouble,UncertainDouble,String,int,String,String) 
    public  ArrheniusKinetics(UncertainDouble p_A, UncertainDouble p_n, UncertainDouble p_E, String p_TRange, int p_rank, String p_source, String p_comment) {
        //#[ operation ArrheniusKinetics(UncertainDouble,UncertainDouble,UncertainDouble,String,int,String,String) 
        A = p_A;
        n = p_n;
        E = p_E;
        TRange = p_TRange;
        rank = p_rank;
        source = p_source;
        comment = p_comment;
        
        
        //#]
    }
    public  ArrheniusKinetics() {
    }
    
    //## operation average(HashSet) 
    public static final Kinetics average(LinkedHashSet p_kSet) {//06/26/09 gmagoon: made p_kSet a LinkedHashSet rather than a HashSet...this seems to make the averaged values reproducible; previously, numerical-errors caused slight differences in results when averaging was done in a different order
        //#[ operation average(HashSet) 
        int size = p_kSet.size();
        if (size == 0) return null;
                
        String type = "ArrheniusKinetics";
        double sum_logA=0;
        double sum_n=0;
        double sum_alpha=0;
        double sum_E=0;
        double max_A=0;
        double min_A=0;
        double max_n=Double.MIN_VALUE;
        double min_n=Double.MAX_VALUE;
        double max_E=Double.MIN_VALUE;
        double min_E=Double.MAX_VALUE;
        double max_alpha=Double.MIN_VALUE;
        double min_alpha=Double.MAX_VALUE;
        //String source = "Average Rate Constants calculated from:\n";
        /*
         * Commented out by MRH on 11-Jun-2009
         * 	Some of the numbers reported in the chem.inp files are averages
         * 	of averages of averages of ... To condense the length of the "source"
         * 	string, I've shorted the expression to "Average of:"
         */
        String source = "Average of: (";
                
        int index = 0;
        Iterator iter = p_kSet.iterator();
        while (iter.hasNext()) {                         
        	index++;
        	ArrheniusKinetics k = (ArrheniusKinetics)iter.next();
            if (k instanceof ArrheniusEPKinetics) {
            	ArrheniusEPKinetics kep = (ArrheniusEPKinetics)k;
            	type = "ArrheniusEPKinetics";
            	UncertainDouble uAlpha = kep.getAlpha();
            	double alpha = uAlpha.getValue();
        		double alpha_upper = uAlpha.getUpperBound();
        		double alpha_lower = uAlpha.getLowerBound();
        		if (max_alpha < alpha_upper) max_alpha = alpha_upper;
        		if (min_alpha > alpha_lower) min_alpha = alpha_lower;
        		sum_alpha += alpha;
        	}
            
            UncertainDouble uA = k.getA();    	
        	double A = uA.getValue();
        	double A_upper = uA.getUpperBound();
        	double A_lower = uA.getLowerBound();
        	if (max_A < A_upper) max_A = A_upper;
        	if (min_A > A_lower) min_A = A_lower;
        	sum_logA += Math.log(A);
                	
            UncertainDouble uN = k.getN();    	
        	double n = uN.getValue();
        	double n_upper = uN.getUpperBound();
        	double n_lower = uN.getLowerBound();
        	if (max_n < n_upper) max_n = n_upper;
        	if (min_n > n_lower) min_n = n_lower;
        	sum_n += n;
                
            UncertainDouble uE = k.getE();    	
        	double E = uE.getValue();
        	double E_upper = uE.getUpperBound();
        	double E_lower = uE.getLowerBound();
        	if (max_E < E_upper) max_E = E_upper;
        	if (min_E > E_lower) min_E = E_lower;
        	sum_E += E;                                             
        	
        	/*
        	 * Update the source of this node. (MRH 11-Jun-2009)
        	 * 	The source string begins with "Average of:".  This line
        	 * 	updates that string with one of the sets of nodes used
        	 * 	in the averaging process.  The && is my way of separting
        	 * 	the sets of nodes from one another.
        	 * 
        	 * Before, RMG kept no record of what was being averaged.
        	 */
        	source += k.source + " && ";
            //source = source + "(" + String.valueOf(index) +")" + k.toChemkinString() + '\n';
        }
        /*
         * This next line removes the last " && " expression and closes
         * 	the parentheses.  The syntax is "source.length() - 4" because
         * 	the " && " expression is length 4.
         */
        source = source.substring(0,source.length()-4) + ")";
                
        double new_A = Math.exp(sum_logA/size);
        double new_dA = Math.max(max_A/new_A, new_A/min_A);
        UncertainDouble A_average = new UncertainDouble(new_A, new_dA, "Multiplier");
                
        double new_n = sum_n/size;
        double new_dn = Math.max(max_n-new_n, new_n-min_n);
        UncertainDouble n_average = new UncertainDouble(new_n, new_dn, "Adder");
                
        double new_E = sum_E/size;
        double new_dE = Math.max(max_E-new_E, new_E-min_E);
        UncertainDouble E_average = new UncertainDouble(new_E, new_dE, "Adder");
                
        if (!getVerbose()) source = "Average:";
        if (type.equals("ArrheniusKinetics")) {
        	return new ArrheniusKinetics(A_average, n_average, E_average,"Unknown",5, source, "Average");
        }	                	
        else if (type.equals("ArrheniusEPKinetics")) {
        	double new_alpha = sum_alpha/size;
        	double new_dalpha = Math.max(max_alpha-new_alpha, new_alpha-min_alpha);
        	UncertainDouble alpha_average = new UncertainDouble(new_alpha, new_dalpha, "Adder");
        	return new ArrheniusEPKinetics(A_average, n_average, alpha_average, E_average, "Unknown",5, source, "Average");
        }
        else throw new InvalidKineticsTypeException("Unknown Kinetics Type: " + type);
        
        
        //#]
    }
    
    //## operation calculateRate(Temperature,double) 
    public double calculateRate(Temperature p_temperature, double p_Hrxn) {
        //#[ operation calculateRate(Temperature,double) 
        double T = p_temperature.getStandard();
        double R = GasConstant.getKcalMolK();
        
        //if (E.getValue() < 0) throw new NegativeEnergyBarrierException();        
        double rate = A.getValue() * Math.pow(T, n.getValue()) * Math.exp(-E.getValue()/R/T);
        return rate;
        
        
        //#]
    }
	
//	## operation calculateRate(Temperature,double) 
    public double calculateRate(Temperature p_temperature) {
        //#[ operation calculateRate(Temperature,double) 
        double T = p_temperature.getStandard();
        double R = GasConstant.getKcalMolK();
        
        //if (E.getValue() < 0) throw new NegativeEnergyBarrierException();        
        double rate = A.getValue() * Math.pow(T, n.getValue()) * Math.exp(-E.getValue()/R/T);
        return rate;
        
        
        //#]
    }
	
	public boolean equals(Kinetics p_k){
		if (p_k == null) return true;
		
		if (Math.abs((p_k.getA().getValue()-A.getValue())/A.getValue()) > 0.01)
			return false;
		if (Math.abs((p_k.getE().getValue()-E.getValue())/E.getValue()) > 0.01 && E.getValue() != 0)
			return false;
		if (Math.abs((p_k.getN().getValue()-n.getValue())/n.getValue()) > 0.01 && n.getValue() != 0)
			return false;
		return true;
	}
	
    //## operation getAValue() 
    public double getAValue() {
        //#[ operation getAValue() 
        return A.getValue();
        //#]
    }
    
    //## operation getEValue() 
    public double getEValue() {
        //#[ operation getEValue() 
        return E.getValue();
        //#]
    }
    
    //## operation getNValue() 
    public double getNValue() {
        //#[ operation getNValue() 
        return n.getValue();
        //#]
    }
    
    //## operation multiply(double) 
    public Kinetics multiply(double p_multiple) {
        //#[ operation multiply(double) 
        UncertainDouble newA = getA().multiply(p_multiple);
        Kinetics newK = new ArrheniusKinetics(newA,getN(),getE(),getTRange(),getRank(),getSource(),getComment());
        return newK;
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        if (getAValue()<0 || getA().getLowerBound()<0 || getA().getUpperBound()<0) return false;
        return true;
        //#]
    }
    
    //## operation toChemkinString() 
    public String toChemkinString(double Hrxn, Temperature p_temperature, boolean includeComments) {
        //#[ operation toChemkinString() 
    	Object [] formatString = new Object[5];
		
//    	formatString[0] = new Double(getAValue());
    	double tempDouble = new Double(getAValue());
    	if (AUnits.equals("moles")) {
    		formatString[0] = tempDouble;
    	} else if (AUnits.equals("molecules")) {
    		formatString[0] = tempDouble / 6.022e23;
    	}
		
		formatString[1] = new Double(getNValue());
		
//		formatString[2] = new Double(E.getValue());
		tempDouble = E.getValue();
		if (EaUnits.equals("kcal/mol")) {
			formatString[2] = tempDouble;
		}
		else if (EaUnits.equals("cal/mol")) {
			formatString[2] = tempDouble * 1000.0;
		}
		else if (EaUnits.equals("kJ/mol")) {
			formatString[2] = tempDouble * 4.184;
		}
		else if (EaUnits.equals("J/mol")) {
			formatString[2] = tempDouble * 4184.0;
		}
		else if (EaUnits.equals("Kelvins")) {
			formatString[2] = tempDouble / 1.987e-3;
		}
		
		formatString[3] = source; formatString[4] = comment;
		
    	if (includeComments)
    		return String.format("%1.3e \t %2.2f \t %3.2f \t !%s  %s", formatString);
    	else
    		return String.format("%1.3e \t %2.2f \t %3.2f \t ", formatString);
    	
        //String.valueOf(getAValue()) + '\t' + String.valueOf(getNValue()) + '\t' + String.valueOf(getEValue() + "\t!" + source + " "+comment);
        //#]
    }
   
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        return "A = " + A.toString() + "; n = " + n.toString() + "; E = " + E.toString() + '\n' + "T Range: " + TRange + '\n' + "Source: " + source + '\n' + "Comment: " + comment + '\n' + "Rank: " + String.valueOf(rank);
        
        
        //#]
    }
    
    public UncertainDouble getA() {
        return A;
    }
    
    public UncertainDouble getE() {
        return E;
    }
    
    public String getTRange() {
        return TRange;
    }
    
    public String getComment() {
        return comment;
    }
    
    public UncertainDouble getN() {
        return n;
    }
    
    public int getRank() {
        return rank;
    }
    
    public String getSource() {
        return source;
    }
	
	public void setSource(String p_source){
		source = p_source;
		return;
	}
	
	public void setComments(String p_comments) {
		comment = p_comments;
	}
	
	public void setFromPrimaryKineticLibrary(boolean p_boolean) {
		fromPrimaryKineticLibrary = p_boolean;
	}
	
	public boolean getFromPrimaryKineticLibrary() {
		return fromPrimaryKineticLibrary;
	}
	
	public static void setVerbose(boolean p_boolean) {
		verbose = p_boolean;
	}
	
	public static boolean getVerbose() {
		return verbose;
	}
	
    public static void setEaUnits(String p_string) {
    	EaUnits = p_string;
    }
    
    public static String getEaUnits() {
    	return EaUnits;
    }
    
    public static void setAUnits(String p_string) {
    	AUnits = p_string;
    }
    
    public static String getAUnits() {
    	return AUnits;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ArrheniusKinetics.java
*********************************************************************/


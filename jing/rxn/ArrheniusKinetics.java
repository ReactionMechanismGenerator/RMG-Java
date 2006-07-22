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
    public static final Kinetics average(HashSet p_kSet) {
        //#[ operation average(HashSet) 
        int size = p_kSet.size();
        if (size == 0) return null;
                
        String type = "ArrehiusKinetics";
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
        String source = "Average Rate Constants calculated from:\n";
                
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
                    
            //source = source + "(" + String.valueOf(index) +")" + k.toChemkinString() + '\n';
        }
                
        double new_A = Math.exp(sum_logA/size);
        double new_dA = Math.max(max_A/new_A, new_A/min_A);
        UncertainDouble A_average = new UncertainDouble(new_A, new_dA, "Multiplier");
                
        double new_n = sum_n/size;
        double new_dn = Math.max(max_n-new_n, new_n-min_n);
        UncertainDouble n_average = new UncertainDouble(new_n, new_dn, "Adder");
                
        double new_E = sum_E/size;
        double new_dE = Math.max(max_E-new_E, new_E-min_E);
        UncertainDouble E_average = new UncertainDouble(new_E, new_dE, "Adder");
                
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
    public String toChemkinString(double Hrxn, Temperature p_temperature) {
        //#[ operation toChemkinString() 
        return String.valueOf(getAValue()) + '\t' + String.valueOf(getNValue()) + '\t' + String.valueOf(getEValue() + "\t!" + source + " "+comment);
        //#]
    }
    
	public String toChemkinStringNoComments(double Hrxn, Temperature p_temperature) {
		return String.valueOf(getAValue()) + '\t' + String.valueOf(getNValue()) + '\t' + String.valueOf(getEValue());
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
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ArrheniusKinetics.java
*********************************************************************/


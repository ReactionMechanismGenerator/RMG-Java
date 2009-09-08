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
import jing.mathTool.*;
import jing.param.Pressure;
import jing.param.Temperature;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ChebyshevPolynomials.java                                                                  
//----------------------------------------------------------------------------

//## class ChebyshevPolynomials 
public class ChebyshevPolynomials implements PDepKinetics {
    
    protected static int NP;		//## attribute NP 
    
    protected static int NT;		//## attribute NT 
    
    protected static Pressure Plow;		//## attribute Plow 
    
    protected static Pressure Pup;		//## attribute Pup 
    
    protected static Temperature Tlow;		//## attribute Tlow 
    
    protected static Temperature Tup;		//## attribute Tup 
    
    protected double [][] alpha;		//## attribute alpha 
    
    
    // Constructors
    
    //## operation ChebyshevPolynomials(int,Temperature,Temperature,int,Pressure,Pressure,double [][]) 
    public  ChebyshevPolynomials(int p_nT, Temperature p_Tlow, Temperature p_Tup, int p_nP, Pressure p_Plow, Pressure p_Pup, double [][] p_alpha) {
        //#[ operation ChebyshevPolynomials(int,Temperature,Temperature,int,Pressure,Pressure,double [][]) 
        if (p_alpha == null) throw new NullPointerException();
        
        if (p_nT < 0 || p_nP < 0 || p_alpha.length < p_nT || p_alpha[0].length < p_nP) throw new InvalidChebyshevPolynomialsException();
        if (p_Tup.getK() < p_Tlow.getK() || p_Pup.getAtm() < p_Plow.getAtm()) throw new InvalidChebyshevPolynomialsException();
         
        NT = p_nT;
        NP = p_nP;
        Tup = p_Tup;
        Tlow = p_Tlow;
        Pup = p_Pup;
        Plow = p_Plow;
        
        alpha = new double[NT][NP];
        
        for (int i = 0; i < NT; i++) {
        	for (int j = 0; j < NP; j++) {
        		alpha[i][j] = p_alpha[i][j];
        	}
        }
        
        
        //#]
    }
    public  ChebyshevPolynomials() {
    }
    
    //## operation calculatePAvg(Pressure) 
    private double calculatePAvg(Pressure p_pressure) {
        //#[ operation calculatePAvg(Pressure) 
        double P = p_pressure.getAtm();
        double Pmin = Plow.getAtm();
        double Pmax = Pup.getAtm();
        
        double result = (2*Math.log(P)-Math.log(Pmin)-Math.log(Pmax))/(Math.log(Pmax)-Math.log(Pmin));
        
        return result;
        
        //#]
    }
    
    public void addChebyshevPolynomial(ChebyshevPolynomials cbp){
    	if (NP != cbp.NP || NT != cbp.NT || Plow.getAtm() !=cbp.Plow.getAtm() || Pup.getAtm() != cbp.Pup.getAtm() || Tup.getK() != cbp.Tup.getK() || Tlow.getK() != cbp.Tlow.getK()){
    		System.err.println("The chebyshev polynomial: \n" + toChemkinString()+  "cannot be added directly to :\n" +cbp.toChemkinString());
    		
    		System.exit(0);
    	}
    	
    	//let us take 10 Gauss-Tchebycheff pivot points
    	int dT = 10;
    	int dP = 10;
    	double [] T_pivot = new double[10];
    	double [] P_pivot = new double[10];
    	double [][] k_total = new double[10][10];
    	
    	for (int i=1; i<=dT; i++){
    		T_pivot[i-1] = Math.cos((2*i-1)*Math.PI/2/dT);
    		T_pivot[i-1] = 2/(T_pivot[i-1]*(1/Tup.getK() - 1/Tlow.getK()) + 1/Tlow.getK() + 1/Tup.getK());
    	}
    	
    	for (int i=1; i<=dP; i++){
    		
    		P_pivot[i-1] = Math.cos((2*i - 1)*Math.PI/2/dP);
    		P_pivot[i-1] = Math.pow(10,( (P_pivot[i-1]*(Math.log10(Pup.getAtm()) -Math.log10(Plow.getAtm())) + Math.log10(Plow.getAtm()) + Math.log10(Pup.getAtm()))/2 ));
    	}
    	
    	for (int i=0; i<dT; i++){
    		for (int j=0; j<dP; j++){
    			k_total[i][j] = Math.log10(calculateRate(new Temperature(T_pivot[i],"K"),new Pressure(P_pivot[j],"atm")) + cbp.calculateRate(new Temperature(T_pivot[i],"K"),new Pressure(P_pivot[j],"atm")));
    		}
    	}
    	
    	for(int n=0; n<NT; n++){
    		for (int m=0; m<NP; m++){
    			double anm = 0;
    			for (int i=0; i<dT; i++){
    				for (int j=0; j<dP; j++){
    					anm = anm + k_total[i][j] * calculatePhi((n+1), Math.cos((i+0.5)*Math.PI/dT)) * calculatePhi((m+1), Math.cos((j+0.5)*Math.PI/dP));
    				}
    			}
    			alpha[n][m] = anm*4/dT/dP;
    		}
    	}
    	
    	for (int i=0 ;i<NT; i++){
    		alpha[i][0] = alpha[i][0]/2.0;
    	}
    	
    	for (int i=0 ;i<NP; i++){
    		alpha[0][i] = alpha[0][i]/2.0;
    	}
    	
    	return;
    	
    }
    
    //## operation calculatePhi(int,double) 
    private double calculatePhi(int p_i, double p_x) {
        //#[ operation calculatePhi(int,double) 
        if (p_i < 1) throw new InvalidChebyshevPolynomialsException();
        double result = Math.cos((p_i-1)*Math.acos(p_x));
        return result;
        //#]
    }
    
    //## operation calculateRate(Temperature,Pressure) 
    public double calculateRate(Temperature p_temperature, Pressure p_pressure) {
        //#[ operation calculateRate(Temperature,Pressure) 
        if (p_temperature.getK()>Tup.getK() || p_temperature.getK()<Tlow.getK()) throw new TOutOfRangeException();
        if (p_pressure.getAtm()>Pup.getAtm() || p_pressure.getAtm()<Plow.getAtm()) throw new POutOfRangeException();
        
        double Tavg = calculateTAvg(p_temperature);
        double Pavg = calculatePAvg(p_pressure);
        
        double result = 0;
        for (int i = 0; i < NT; i++) {
        	for (int j = 0; j < NP; j++) {
        		result += alpha[i][j]*calculatePhi(i+1, Tavg)*calculatePhi(j+1,Pavg);
        	}
        }
        
        double k = Math.pow(10,result);//changed by Sally 1/24/06
        
        return k;
        //#]
    }
    
    //## operation calculateTAvg(Temperature) 
    private double calculateTAvg(Temperature p_temperature) {
        //#[ operation calculateTAvg(Temperature) 
        double T = p_temperature.getK();
        double Tmin = Tlow.getK();
        double Tmax = Tup.getK();
        
        double result = (2.0/T-1.0/Tmin-1.0/Tmax)/(1.0/Tmax-1.0/Tmin);
        return result;
        //#]
    }
    
    //## operation toChemkinString() 
    public String toChemkinString() {
        //#[ operation toChemkinString() 
        String result = "TCHEB / " + Tlow.getK() + " " + Tup.getK() + " /";
        result += "\tPCHEB / " + Plow.getAtm() + " " + Pup.getAtm() + " /\n";
        result += "CHEB / " + NT + '\t' + NP + " /\n";
        for (int i = 0; i < NT; i++) {
        	result += "CHEB / ";
        	for (int j = 0; j < NP; j++) {
        		result += String.format("%1.7e",alpha[i][j]) + " ";
        	}
        	result += "/\n";
        }
        return result;
        
        //#]
    }
    
    public static int getNP() {
        return NP;
    }
    
    public void setNP(int p_NP) {
        NP = p_NP;
    }
    
    public static int getNT() {
        return NT;
    }
    
    public void setNT(int p_NT) {
        NT = p_NT;
    }
    
    public static Pressure getPlow() {
        return Plow;
    }
    
    public static void setPlow(Pressure p_Plow) {
        Plow = p_Plow;
    }
    
    public static Pressure getPup() {
        return Pup;
    }
    
    public static void setPup(Pressure p_Pup) {
        Pup = p_Pup;
    }
    
    public static Temperature getTlow() {
        return Tlow;
    }
    
    public static void setTlow(Temperature p_Tlow) {
        Tlow = p_Tlow;
    }
    
    public static Temperature getTup() {
        return Tup;
    }
    
    public static void setTup(Temperature p_Tup) {
        Tup = p_Tup;
    }
    
    public double getAlpha(int i2, int i1) {
        return alpha[i2][i1];
    }
    
    public void setAlpha(int i2, int i1, double p_alpha) {
        alpha[i2][i1] = p_alpha;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ChebyshevPolynomials.java
*********************************************************************/


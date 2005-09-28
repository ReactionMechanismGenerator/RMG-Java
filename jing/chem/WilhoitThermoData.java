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



package jing.chem;


import Jama.*;
import java.util.*;
import jing.param.*;
import jing.param.Temperature;

//## package jing::chem 

//----------------------------------------------------------------------------
// jing\chem\WilhoitThermoData.java                                                                  
//----------------------------------------------------------------------------

//## class WilhoitThermoData 
public class WilhoitThermoData {
    
    protected final double B = 500.0;		//## attribute B 
    
    protected double I;		//## attribute I 
    
    protected double J;		//## attribute J 
    
    protected double [] a = new double[4];		//## attribute a 
    
    protected double cpInfinity = -1;		//## attribute cpInfinity 
    
    protected double cpZero = -1;		//## attribute cpZero 
    
    
    // Constructors
    
    //## operation WilhoitThermoData() 
    private  WilhoitThermoData() {
        //#[ operation WilhoitThermoData() 
        //#]
    }
    //## operation WilhoitThermoData(double,double,double [],double,double) 
    private  WilhoitThermoData(double p_cpZero, double p_cpInfinity, double [] p_a, double p_I, double p_J) {
        //#[ operation WilhoitThermoData(double,double,double [],double,double) 
        cpZero = p_cpZero;
        cpInfinity = p_cpInfinity;
        for (int i = 0; i < 3; i++) {
        	a[i] = p_a[i];
        } 
        I = p_I;
        J = p_J;
        
        //#]
    }
    
    //## operation calculateAIJ(ChemGraph) 
    private void calculateAIJ(ChemGraph p_chemGraph) {
        //#[ operation calculateAIJ(ChemGraph) 
        // declare and initialize some variables
        double deltaCp = cpInfinity - cpZero;
        double[] T = {300., 400., 500., 600., 800., 1000., 1500.};
        double[] sum = {0., 0., 0., 0., 0., 0., 0.};
        double[] y = new double[7];
        double[] b = new double[7];
        double[] cp = new double[7];
        int nT = T.length; int nOrder = 3;
        double yi2; double Aij;
        
        // get Cp values at Benson's temperatures, H98 and S298
        ThermoData thermoData = p_chemGraph.getThermoData();
        cp[0] = thermoData.getCp300();
        cp[1] = thermoData.getCp400();
        cp[2] = thermoData.getCp500();
        cp[3] = thermoData.getCp600();
        cp[4] = thermoData.getCp800();
        cp[5] = thermoData.getCp1000();
        cp[6] = thermoData.getCp1500(); 
        
        double H298 = thermoData.getH298(); 
        double S298 = thermoData.getS298();
        
        // calculate y = T/(T+B)
        for (int i=0; i <= nT-1; i++) {
        	y[i] = T[i] / (T[i] + B);
        }
        
        // calculate b, the rhs for the linear least-squares problem
        for (int i=0; i <= nT-1; i++) {
        	yi2 = Math.pow(y[i],2.);
        	b[i] = (cp[i] - cpZero - yi2*deltaCp) / (deltaCp*yi2*(y[i]-1.0));
        }
        
        // construct JAMA matrix object, fill with rhs
        Matrix rhs = new Matrix(b, 1);  // construct r.h.s.
        rhs = rhs.transpose();          // make row into column
        
        // construct JAMA matrix object, fill w/ ones
        Matrix A = new Matrix(nT, nOrder+1, 1.);
        
        // form the matrix A for th LLS problem
        for (int i=0; i <= nT-1; i++) {
        	for (int j=0; j <= nOrder; j++) {
        		Aij = Math.pow(y[i],(double)j);
        		A.set(i, j, Aij);
        	}
        }
        
        // Use QR factorization to solve the linear least-squares problem for a0...a3
        QRDecomposition QR = A.qr();
        Matrix coeffs = QR.solve(rhs);
        
        // print out the matrix and results
        // A.print(16, 4);
        // rhs.print(25,16);
        // coeffs.print(25, 16);
        
        // convert Matrix object to a double-precision vector
        double[] coeffsArray = coeffs.getColumnPackedCopy();
        
        // save parameters in chemgraph attribute called wilhoitParameters
        for (int i=0; i<=nOrder; i++) {
        	if(p_chemGraph.isMonatomic())
        		a[i] = 0.0;
        	else
        		a[i] = coeffsArray[i];
        }
        
        // calculate I and J, the integration constants in the enthalpy and entropy eqns.
        Temperature T298 = new Temperature(298.0,"K");
        I = 0.0;
        J = 0.0;
        double Itemp = H298 - calculateH(T298);
        I = 0.0;
        J = 0.0;
        double Jtemp = S298 - calculateS(T298);
        I = Itemp;
        J = Jtemp;
        
        //#]
    }
    
    //## operation calculateCp(Temperature) 
    public double calculateCp(Temperature p_temperature) {
        //#[ operation calculateCp(Temperature) 
        return calculateCp(p_temperature.getK());
        //#]
    }
    
    //## operation calculateCp(double) 
    public double calculateCp(double p_T) {
        //#[ operation calculateCp(double) 
        double T = p_T;
        
        // calculate some intermediate variables
        double deltacp = cpInfinity - cpZero;
        double y = T/(T+B);
        double y2 = Math.pow(y,2.);
        
        // calculate sum of a(i)*y^i (3rd-order Wilhoit polynomial)
        double aySum = 0.;
        for (int i=0; i<=3; i++) {
        	aySum += a[i]*Math.pow(y,(double)i);
        }
        
        double cp = cpZero + deltacp*y2*(1 + (y-1)*aySum);
        return cp;  // cal/mol-K
        
        //#]
    }
    
    //## operation calculateCpInfinity(ChemGraph) 
    private double calculateCpInfinity(ChemGraph p_chemGraph) {
        //#[ operation calculateCpInfinity(ChemGraph) 
        double cp;
        final double R = GasConstant.getCalMolK(); // cal/mol-K
        int Natoms = p_chemGraph.getAtomNumber();
        int Nrotors = p_chemGraph.getInternalRotor();
        
        if (p_chemGraph.isMonatomic())
        	cp = 2.5*R;
        else if (p_chemGraph.isLinear())
        	cp = (3.0*Natoms - 1.5)*R;
        else
        	cp = (3.0*Natoms - (2.0 + 0.5*Nrotors))*R;
        
        cpInfinity = cp;
        return cp;
        //#]
    }
    
    //## operation calculateCpZero(ChemGraph) 
    private double calculateCpZero(ChemGraph p_chemGraph) {
        //#[ operation calculateCpZero(ChemGraph) 
        double cp;
        final double R = GasConstant.getCalMolK(); // cal/mol-K
        
        if (p_chemGraph.isMonatomic())
        	cp = 2.5*R;
        else if (p_chemGraph.isLinear())
        	cp = 3.5*R;
        else
        	cp = 4.0*R;
        
        cpZero = cp;	
        return cp;
        	
        
        //#]
    }
    
    //## operation calculateG(Temperature) 
    public double calculateG(Temperature p_temperature) {
        //#[ operation calculateG(Temperature) 
        return calculateG(p_temperature.getK());
        //#]
    }
    
    //## operation calculateG(double) 
    public double calculateG(double p_T) {
        //#[ operation calculateG(double) 
        return (calculateH(p_T)*1000-p_T*calculateS(p_T))/1000; 
        //#]
    }
    
    //## operation calculateH(Temperature) 
    public double calculateH(Temperature p_temperature) {
        //#[ operation calculateH(Temperature) 
        return calculateH(p_temperature.getK());
        //#]
    }
    
    //## operation calculateH(double) 
    public double calculateH(double p_T) {
        //#[ operation calculateH(double) 
        double T = p_T;
        
        double a0=a[0]; 
        double a1=a[1]; 
        double a2=a[2]; 
        double a3=a[3];
                
        // calculate some intermediate variables
        double deltacp = cpInfinity - cpZero;
        double y = T/(T+B);
        double y2 = Math.pow(y,2.);
        double aSum = a0 + a1 + a2 + a3;
        double[] b = new double[4];
        b[0] = 0.5*a0 + 1./6.*a1 + 1./6.*a2 + 1./6.*a3;
        b[1] = 1./3.*a1 + 1./12.*a2 + 1./12.*a3;
        b[2] = 0.25*a2 + 0.05*a3;
        b[3] = 0.2*a3;
        
        // calculate sum in the enthalpy equation
        double hSum = 0.;
        for (int i=0; i<=3; i++) {
        	hSum += b[i]*Math.pow(y,(double)i);
        }
        
        double H = I + cpZero/1000.*T - deltacp/1000.*T*((2.+aSum)*(0.5*y - 1. + (1./y - 1.)*Math.log(T/y)) + y2*hSum);
        return H;  // kcal/mol
        
        //#]
    }
    
    //## operation calculateS(Temperature) 
    public double calculateS(Temperature p_temperature) {
        //#[ operation calculateS(Temperature) 
        return calculateS(p_temperature.getK());
        //#]
    }
    
    //## operation calculateS(double) 
    public double calculateS(double p_T) {
        //#[ operation calculateS(double) 
        double T = p_T;
        
        // calculate some intermediate variables
        double deltacp = cpInfinity - cpZero;
        double y = T/(T+B);
        double y2 = Math.pow(y,2.);
        
        // calculate sum in entropy equation
        double sSum = 0.;
        for (int i=0; i<=3; i++) {
        	sSum += a[i]/(2.+(double)i)*Math.pow(y,(double)i);
        }
        
        double S = J + cpInfinity*Math.log(T) - deltacp*(Math.log(y) + y + y2*sSum);
        return S;  // cal/mol-K
        
        //#]
    }
    
    //## operation make(ChemGraph) 
    public static WilhoitThermoData make(final ChemGraph p_chemGraph) {
        //#[ operation make(ChemGraph) 
        WilhoitThermoData wtd = new WilhoitThermoData();
        
        wtd.calculateCpZero(p_chemGraph);
        wtd.calculateCpInfinity(p_chemGraph);
        wtd.calculateAIJ(p_chemGraph);
        
        return wtd;
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String s = "Cp0 = " + cpZero + '\n';
        s += "CpInf = " + cpInfinity + '\n';    
        s += "B = " + B + '\n';
        s += "I = " + I + '\n';
        s += "J = " + J + '\n';
        
        for (int i=0; i<4; i++) {
        	s += "a"+ i +" = " + a[i] + '\n';
        }
        return s; 
        
        //#]
    }
    
    public final double getB() {
        return B;
    }
    
    public double getI() {
        return I;
    }
    
    public void setI(double p_I) {
        I = p_I;
    }
    
    public double getJ() {
        return J;
    }
    
    public void setJ(double p_J) {
        J = p_J;
    }
    
    public double getA(int i1) {
        return a[i1];
    }
    
    public void setA(int i1, double p_a) {
        a[i1] = p_a;
    }
    
    public double getCpInfinity() {
        return cpInfinity;
    }
    
    public void setCpInfinity(double p_cpInfinity) {
        cpInfinity = p_cpInfinity;
    }
    
    public double getCpZero() {
        return cpZero;
    }
    
    public void setCpZero(double p_cpZero) {
        cpZero = p_cpZero;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\WilhoitThermoData.java
*********************************************************************/


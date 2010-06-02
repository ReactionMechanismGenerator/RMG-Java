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
import jing.chemUtil.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.mathTool.*;
import jing.param.Temperature;

public class GATP_Solvation implements GeneralSolvationGAPP {
	
    private static GATP_Solvation INSTANCE = new GATP_Solvation();

    // Constructors
    private  GATP_Solvation() {
    }
	
	public ThermoData generateSolvThermoData(ChemGraph p_chemGraph) {
		double r_solute=p_chemGraph.getRadius();               // Returns VdW radius in meter
		double r_solvent; r_solvent=3.498e-10;// 3.311;        // Manually assigned solvent radius [=] meter Calculated using Connolly solvent excluded volume from Chem3dPro
		double r_cavity=r_solute+r_solvent;                    // Cavity radius [=] meter
		double rho; rho=3.09e27;  //0.00381;                   // number density of solvent [=] molecules/m^3   Value here is for decane using density =0.73 g/cm3
		double parameter_y=4.1887902*rho*Math.pow(r_solvent, 3);  // Parameter y from Ashcraft Thesis Refer pg no. 60. (4/3)*pi*rho*r^3
		double parameter_ymod=parameter_y/(1-parameter_y);     // parameter_ymod= y/(1-y) Defined for convenience
		double R=8.314;                                        // Gas constant units J/mol K
		double T=298;                                          // Standard state temperature
		
		// Definitions of K0, K1 and K2 correspond to those for K0', K1' and K2' respectively from Ashcraft's Thesis (-d/dT of K0,K1,K2)
		double K0= -R*(-Math.log(1-parameter_y)+(4.5*parameter_ymod*parameter_ymod));
		double K1= (R*0.5/r_solvent)*((6*parameter_ymod)+(18*parameter_ymod*parameter_ymod));
		double K2= -(R*0.25/(r_solvent*r_solvent))*((12*parameter_ymod)+(18*parameter_ymod*parameter_ymod));
		
		// Basic definition of entropy change of solvation from Ashcfrat's Thesis
		double deltaS0;
		deltaS0=K0+(K1*r_cavity)+(K2*r_cavity*r_cavity);
		
		// Generation of Abraham Solute Parameters
		AbramData result_Abraham= new AbramData();
		result_Abraham= p_chemGraph.getAbramData();

		// Solute descriptors from the Abraham Model
		double S=result_Abraham.S;
		double B=result_Abraham.B;
		double E=result_Abraham.E;
		double L=result_Abraham.L;
		double A=result_Abraham.A;
		
        //Manually specified solvent descriptors 
		// (constants here are for dry decane, from from M.H. Abraham et al. / J. Chromatogr. A 1037 (2004) 29–47 )
        double c = 0.156;
        double s = 0;
        double b = 0; 
        double e = -0.143;
        double l = 0.989;
        double a = 0;
		
		// Dry octan-1-ol, from M.H. Abraham et al. / J. Chromatogr. A 1037 (2004) 29–47 
		// c=-0.120; e=-0.203; s=0.560; a=3.560; b=0.702; l=0.939
		
        double logK = c + s*S + b*B + e*E + l*L + a*A;    // Implementation of Abraham Model
        double deltaG0_decane = -8.314*298*logK;
		//       System.out.println("The free energy of solvation in decane at 298K w/o reference state corrections  = " + deltaG0_decane +" J/mol for " );
		
        // Calculation of enthalpy change of solvation using the data obtained above
        double deltaH0=deltaG0_decane+(T*deltaS0);
        deltaS0=deltaS0/4.18;   //unit conversion from J/mol to cal/mol
        deltaH0=deltaH0/4180;   //unit conversion from J/mol to kcal/mol
		
		// Generation of Gas Phase data to add to the solution phase quantities
		ThermoData solvationCorrection = new ThermoData(deltaH0, deltaS0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,"Solvation correction");
		
		// Now, solvationCorrection contains solution phase estimates of CORRECTION TO H298, S298 and all the gas phase heat capacities.
		// Assuming the solution phase heat capcities to be the same as that in the gas phase we wouls now want to pass on this
		// modified version of result to the kinetics codes. This might require reading in a keyword from the condition.txt file.
		// Exactly how this will be done is yet to be figured out.
		
        return solvationCorrection;
    }
	
	
	protected static GATP_Solvation getINSTANCE() {
        return INSTANCE;
    }
}

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

public class LJData extends LJGroupData {
    protected int na;

    // Constructors
    public LJData() {
        super();
    }

    public void plus(LJGroupData p_LJData) {
        if (p_LJData == null) return;

        dTc += p_LJData.dTc;
        dPc += p_LJData.dPc;
        dVc += p_LJData.dVc;
        dTb += p_LJData.dTb;
        shapeIndex += p_LJData.shapeIndex;

        if (p_LJData.getName() != null) {
        	if (name == null) name = "(" + p_LJData.getName();
        	else name = name + ", " + p_LJData.getName();
        }
    }

    public String toString() {
        //return super.toString();
    	return "LJ parameters estimated by RMG:" + 
    		" Tc=" + String.format("%8.3f", calculateTc()) + "K" +
    		" Pc=" + String.format("%8.3f", calculatePc()) + "bar" +
    		" Vc=" + String.format("%8.3f", calculateVc()) + "cm3/mol" +
    		" Tb=" + String.format("%8.3f", calculateTb()) + "K";
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
    
    //calculates omega (acentric factor)
    /*
     * Source:
     * B.I. Lee and M.G. Kesler;
     * "A Generalized Thermodynamic Correlation Based on Three-Parameter Corresponding States";
     * AIChE J., (1975), 21(3), 510-527 
     * 
     */
    public double calculateOmega(){
        double f = calculateTb()/calculateTc();
        return (-1.0*Math.log(calculatePc()/1.01325) -5.92714 + 6.09648/f + 1.28862*Math.log(f) - 0.169347*Math.pow(f, 6) ) / ( 15.2518 - 15.6875/f - 13.4721*Math.log(f) + 0.43577*Math.pow(f,6));
        // Equation only valid for f <= 0.8
    }

    //calculates sigma in angstroms
    /*
     * Source:
     * L.S. Tee, S. Gotoh, and W.E. Stewart;
     * "Molecular Parameters for Normal Fluids: The Lennard-Jones 12-6 Potential";
     * Ind. Eng. Chem., Fundamentals (1966), 5(3), 346-63
     * Table III: Correlation iii
     */
    public double calculateSigma(){
        return (2.3442*Math.exp(0.1303*calculateOmega())) * Math.pow(calculateTc()*1.01325/calculatePc(),1.0/3.0);
    }

    //calculates epsilon in K
    /*
     * Source:
     * L.S. Tee, S. Gotoh, and W.E. Stewart;
     * "Molecular Parameters for Normal Fluids: The Lennard-Jones 12-6 Potential";
     * Ind. Eng. Chem., Fundamentals (1966), 5(3), 346-63
     * Table III: Correlation iii
     */
    public double calculateEpsilon(){
        return 0.8109*Math.exp(-0.6228*calculateOmega())*calculateTc();
    }

}
/*********************************************************************
        File Path	: RMG\RMG\jing\chem\LJData.java
*********************************************************************/


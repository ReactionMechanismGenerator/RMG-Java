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

// Amrit Jalan, April 19, 2009
/**
 Replicates the functions of ThermoGAvalue for the Abraham Solvation model i.e. defines varibles like S, E, B, V and A
 Immutable data holds all the Platt's group value.
 */

public class UnifacGAValue {
    protected double R = 0;	
    protected double Q = 0;	
	
	// Constructors
	
    public UnifacGAValue() {
        R = 0;
        Q = 0;
    }
	
    public  UnifacGAValue(double p_R, double p_Q) {
        R = p_R;
        Q = p_Q;
    }
	
    public  UnifacGAValue(UnifacGAValue p_ga) {
        R = p_ga.R;
        Q = p_ga.Q;
    }
	
    public String toString() {
        String s = "";
        s = s + String.valueOf(R) + '\t';
        s = s + String.valueOf(Q);
		
        return s;
    }
	
	protected double getR() {
        return R;
    }
	
    protected double getQ() {
        return Q;
    }
	
}
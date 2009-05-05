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
 */

/**
Immutable data holds all the Platt's group value.
*/
//## class ThermoGAValue


public class AbrahamGAValue {

    protected double S = 0;		//## attribute S

    protected double B = 0;		//## attribute B

    protected double E = 0;		//## attribute E

    protected double L = 0;		//## attribute L

    protected double A= 0;      //## attribute A



// Constructors

    //## operation ThermoGAValue()
    public AbrahamGAValue() {
        //#[ operation ThermoGAValue()
        S = 0;
        B = 0;
        E = 0;
        L = 0;
        A=  0;

        //#]
    }

    //## operation ThermoGAValue(double,double,double,double,double,double,double,double,double,String)
    public  AbrahamGAValue(double p_S, double p_B, double p_E, double p_L, double p_A) {
        //#[ operation ThermoGAValue(double,double,double,double,double,double,double,double,double,String)
        S = p_S;
        B = p_B;
        E = p_E;
        L = p_L;
        A=  p_A;

        //#]
    }

    //## operation ThermoGAValue(String,ThermoGAValue,String)
    public  AbrahamGAValue(AbrahamGAValue p_ga) {
        //#[ operation ThermoGAValue(String,ThermoGAValue,String)
        S = p_ga.S;
        B = p_ga.B;
        E = p_ga.E;
        L = p_ga.L;
        A=  p_ga.A;

        //#]
    }

//## operation toString()
    public String toString() {
        //#[ operation toString()
        String s = "";
        s = s + String.valueOf(S) + '\t';
        s = s + String.valueOf(B) + '\t';
        s = s + String.valueOf(E) + '\t';
        s = s + String.valueOf(L) + '\t';
        s = s + String.valueOf(A);

        return s;
        //#]
    }

     protected double getS() {
        return S;
    }

    protected double getB() {
        return B;
    }

    protected double getE() {
        return E;
    }

    protected double getL() {
        return L;
    }

    protected double getA() {
        return A;
    }

}
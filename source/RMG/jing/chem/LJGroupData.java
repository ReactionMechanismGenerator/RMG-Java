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

public class LJGroupData {

    protected double dTc = 0;		
    protected double dPc = 0;  
    protected double dVc = 0;
    protected double dTb = 0;
    protected int shapeIndex = 0;
    protected String comments = null;
    protected String name = null;
    protected String source = null;


    // Constructors
    public  LJGroupData() {
        dTc = 0;		
        dPc = 0;  
        dVc = 0;
        dTb = 0;
        shapeIndex = 0;
    }
    
    public  LJGroupData(double p_dTc, double p_dPc, double p_dVc, double p_dTb, int p_shape, String p_comments) {
        dTc = p_dTc;
        dPc = p_dPc;
        dVc = p_dVc;
        dTb = p_dTb;
        shapeIndex = p_shape;
        comments = p_comments;
    }

    public  LJGroupData(String p_name, LJGroupData p_ga, String p_comments) {
        dTc = p_ga.dTc;
        dPc = p_ga.dPc;
        dVc = p_ga.dVc;
        dTb = p_ga.dTb;
        shapeIndex = p_ga.shapeIndex;
        comments = p_comments;
        name = p_name;
    }

    public String getName() {
        return name;
    }

}
/*********************************************************************
        File Path	: RMG\RMG\LJGroupData.java
*********************************************************************/
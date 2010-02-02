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



package jing.param;


import java.util.*;

//## package jing::param

//----------------------------------------------------------------------------
// jing\param\ParameterInfor.java
//----------------------------------------------------------------------------

//## class ParameterInfor
public class ParameterInfor {

    protected int ID = -1;		//## attribute ID

    protected String name;		//## attribute name

    protected double range;		//## attribute range

    protected int type = -1;		//## attribute type


    // Constructors

    //## operation ParameterInfor(String,int,double)
    public  ParameterInfor(String p_name, int p_ID, double p_range) {
        //#[ operation ParameterInfor(String,int,double)
        name = p_name.toUpperCase();

        if (name.equals("C0")) type = 1;
        else if (name.equals("T")) type = 2;
         else if (name.equalsIgnoreCase("k")) type = 3;//svp
         else if (name.equalsIgnoreCase("G")) type = 4; //sandeep for thermo
        else type = -1;

        ID = p_ID;
        range = p_range;


        //#]
    }
    public  ParameterInfor() {
    }

    public int getID() {
        return ID;
    }

    public String getName() {
        return name;
    }

    public double getRange() {
        return range;
    }

    public int getType() {
        return type;
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\param\ParameterInfor.java
*********************************************************************/


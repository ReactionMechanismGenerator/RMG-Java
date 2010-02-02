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



package jing.rxnSys;


import java.util.*;
import jing.param.*;

//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\DynamicSimulator.java
//----------------------------------------------------------------------------

//## class DynamicSimulator
public interface DynamicSimulator {


    //## operation clean()
    //void clean();

    //## operation solve(boolean,ReactionModel,boolean,SystemSnapshot,ReactionTime,ReactionTime,Temperature,Pressure,boolean)
    SystemSnapshot solve(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, ReactionTime p_beginTime, ReactionTime p_endTime, Temperature p_temperature, Pressure p_pressure,boolean p_conditionChanged, TerminationTester tt, int p_iterationNum);

    void addConversion(double [] p_conversion, int p_numConversiosn);
    double[] getConversion();
}
/*********************************************************************
        File Path	: RMG\RMG\jing\rxnSys\DynamicSimulator.java
*********************************************************************/




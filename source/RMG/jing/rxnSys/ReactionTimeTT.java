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

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\ReactionTimeTT.java                                                                  
//----------------------------------------------------------------------------

/**
Using reaction time to control the termination of the reaction system.  ie., if the reaction system has reacted for a given amount of time, the reaction system terminates.
*/

//## class ReactionTimeTT 
public class ReactionTimeTT implements TerminationTester {

    protected ReactionTime finalTime;		//## attribute finalTime 

    protected LinkedList timeStep;

    // Constructors
 
	//## operation ReactionTimeTT()

    private  ReactionTimeTT() {
        //#[ operation ReactionTimeTT() 
        //#]
    }

    //## operation ReactionTimeTT(ReactionTime) 
    public  ReactionTimeTT(ReactionTime p_finalTime) {
        //#[ operation ReactionTimeTT(ReactionTime) 
        finalTime = p_finalTime;
        //#]

    }

    protected void setTimeSteps(LinkedList p_timeStep) {
        //5/5/08 gmagoon: added null check to avoid null pointer exception when adding finalTime in cases where no intermediate time steps are specified
        if (p_timeStep == null){
            timeStep = new LinkedList();
        }
        else{
            timeStep = p_timeStep;
        }
    	timeStep.add(finalTime);
    }

    

    //## operation getFinalTime() 
    protected ReactionTime getFinalTime() {
        //#[ operation getFinalTime() 
        return finalTime;
        //#]
    }

    

    //## operation isReactionTerminated(InitialStatus,PresentStatus) 
	public boolean isReactionTerminated(InitialStatus p_initialStatus, PresentStatus p_presentStatus) {
        //#[ operation isReactionTerminated(InitialStatus,PresentStatus) 
        ReactionTime t = p_presentStatus.getPresentTime();
        if (t.reach(finalTime)) return true;
        else return false;
        //#]
    }
   
    //## operation setFinalTime(ReactionTime) 
    protected void setFinalTime(ReactionTime p_finalTime) {
        //#[ operation setFinalTime(ReactionTime) 
        finalTime = p_finalTime;
        //#]
    }
    
}

/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ReactionTimeTT.java
*********************************************************************/

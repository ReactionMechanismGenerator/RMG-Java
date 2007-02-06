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
    	timeStep = p_timeStep;
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


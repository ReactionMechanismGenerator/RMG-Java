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
// jing\rxnSys\FinishController.java                                                                  
//----------------------------------------------------------------------------

//## class FinishController 
public class FinishController {
    
    protected ReactionSystem reactionSystem;
    protected TerminationTester terminationTester;
    protected ValidityTester validityTester;
    
    // Constructors
    
    //## operation FinishController(TerminationTester,ValidityTester) 
    public  FinishController(TerminationTester p_terminationTester, ValidityTester p_validityTester) {
        //#[ operation FinishController(TerminationTester,ValidityTester) 
        terminationTester = p_terminationTester;
        validityTester = p_validityTester;
        
        
        //#]
    }
    public  FinishController() {
    }
    
    //## operation isFinished() 
    public boolean isFinished() {
        //#[ operation isFinished() 
        return (isReactionTerminated() && isModelValid());
        //#]
    }
    
    //## operation isModelValid() 
    public boolean isModelValid() {
        //#[ operation isModelValid() 
        return (validityTester.isModelValid(reactionSystem));
        //#]
    }
    
    //## operation isReactionTerminated() 
    public boolean isReactionTerminated() {
        //#[ operation isReactionTerminated() 
        return terminationTester.isReactionTerminated(reactionSystem.getInitialStatus(), reactionSystem.getPresentStatus());
        
        
        //#]
    }
    
    //## operation newTerminationTester(String,Object) 
    public TerminationTester newTerminationTester(String p_type, Object p_end) {
        //#[ operation newTerminationTester(String,Object) 
        if (p_type.equals("Conversion")) {
        	LinkedList scs = (LinkedList)p_end;
        	terminationTester = new ConversionTT(scs);
        }
        else if (p_type.equals("ReactionTime")) {
        	ReactionTime rt = (ReactionTime)p_end;
        	terminationTester = new ReactionTimeTT(rt);
        }
        return terminationTester;
        //#]
    }
    
    public ReactionSystem getReactionSystem() {
        return reactionSystem;
    }
    
    public void __setReactionSystem(ReactionSystem p_ReactionSystem) {
        reactionSystem = p_ReactionSystem;
    }
    
    public void _setReactionSystem(ReactionSystem p_ReactionSystem) {
        if(reactionSystem != null)
            reactionSystem.__setFinishController(null);
        __setReactionSystem(p_ReactionSystem);
    }
    
    public void setReactionSystem(ReactionSystem p_ReactionSystem) {
        if(p_ReactionSystem != null)
            p_ReactionSystem._setFinishController(this);
        _setReactionSystem(p_ReactionSystem);
    }
    
    public void _clearReactionSystem() {
        reactionSystem = null;
    }
    
    public TerminationTester getTerminationTester() {
        return terminationTester;
    }
    
    public void setTerminationTester(TerminationTester p_TerminationTester) {
        terminationTester = p_TerminationTester;
    }
    
    public ValidityTester getValidityTester() {
        return validityTester;
    }
    
    public void setValidityTester(ValidityTester p_ValidityTester) {
        validityTester = p_ValidityTester;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\FinishController.java
*********************************************************************/


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



package jing.rxn;


import java.util.*;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\LibraryReaction.java                                                                  
//----------------------------------------------------------------------------

/**
the reaction in the kinetics library.  the part inherited from super class Reaction is immutable.  but the _Key attribute part is mutable.  Say, we can change the _Key for search a reaction in a kinetics library.
*/
//## class LibraryReaction 
public class LibraryReaction extends Reaction {
    
    protected Structure key;		//## attribute key 
    
    
    // Constructors
    
    //## operation LibraryReaction(Structure,RateConstant) 
    protected  LibraryReaction(Structure p_structure, RateConstant p_rateConstant) {
        //#[ operation LibraryReaction(Structure,RateConstant) 
        structure = p_structure;
        rateConstant = p_rateConstant;
        comments = "Library Reaction";
        
        
        
        //#]
    }
    //## operation LibraryReaction(Structure,RateConstant,String) 
    private  LibraryReaction(Structure p_structure, RateConstant p_rateConstant, String p_comments) {
        //#[ operation LibraryReaction(Structure,RateConstant,String) 
        structure = p_structure;
        rateConstant = p_rateConstant;
        if (comments != null) comments = "Libarry Reaction: " + p_comments;
        else comments = "Libarry Reaction";
        
        
        
        //#]
    }
    public  LibraryReaction() {
    }
    
    //## operation getKey() 
    public Structure getKey() {
        //#[ operation getKey() 
        return key;
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\LibraryReaction.java
*********************************************************************/


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
// jing\rxn\ReactionLibrary.java                                                                  
//----------------------------------------------------------------------------

/**
This is a super user-defined reaction library.  Here, "super" means that all the library reactions whose reactant(s)/product(s) appears in the reaction system will be included in the final reaction system model.  Also, the kinetics in this library has the higher priority than the kinetics in our default reaction/kinetics template.  So, be careful when make this ReactionLibrary.  we will include things in this library without any checking.  But, if there is an alternative kinetics we find in our kinetics template, we will output it to warn user.
*/
//## class ReactionLibrary 
public class ReactionLibrary {
    
    private static ReactionLibrary INSTANCE = new ReactionLibrary();		//## attribute INSTANCE 
    
    protected HashSet libraryReaction;
    
    // Constructors
    
    /**
    Requires:
    Effects: this is the only constructor in this singleton reaction library.  it is protected, which means no user can construct it.  the construction should go through instance() method to check if there is only one instance of this class.
    Modifies: itsLibraryReaction
    */
    //## operation ReactionLibrary() 
    private  ReactionLibrary() {
        {
            libraryReaction=new HashSet();
        }
        //#[ operation ReactionLibrary() 
        libraryReaction = new HashSet();
        //#]
    }
    
    //## operation addLibraryReaction(LibraryReaction) 
    public void addLibraryReaction(LibraryReaction p_LibraryReaction) {
        //#[ operation addLibraryReaction(LibraryReaction) 
        libraryReaction.add(p_LibraryReaction);
        //#]
    }
    
    //## operation clearLibraryReaction() 
    public void clearLibraryReaction() {
        //#[ operation clearLibraryReaction() 
        libraryReaction.clear();
        //#]
    }
    
    //## operation getLibraryReaction() 
    public Iterator getLibraryReaction() {
        //#[ operation getLibraryReaction() 
        Iterator iter=libraryReaction.iterator();
        return iter;
        //#]
    }
    
    /**
    Requires:
    Effects: check if this reaction library is empty.  if it is, return true; otherwise, return false;
    Modifies:
    */
    //## operation isEmpty() 
    public boolean isEmpty() {
        //#[ operation isEmpty() 
        return (size() == 0);
        
        
        
        //#]
    }
    
    //## operation removeLibraryReaction(LibraryReaction) 
    public void removeLibraryReaction(LibraryReaction p_LibraryReaction) {
        //#[ operation removeLibraryReaction(LibraryReaction) 
        libraryReaction.remove(p_LibraryReaction);
        //#]
    }
    
    /**
    Requires:
    Effects: check if all the library reaction in this reaction library are valid
    Modifies:
    */
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        Iterator iter = getLibraryReaction();
        Reaction reaction;
        while (iter.hasNext()) {
        	reaction = (Reaction)iter.next();
        	if (!reaction.repOk()) return false;
        }
        
        return true;
        //#]
    }
    
    /**
    Requires:
    Effects: return the size of itsLibraryReaction.
    Modifies:
    */
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return libraryReaction.size();
        //#]
    }
    
    protected static ReactionLibrary getINSTANCE() {
        return INSTANCE;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ReactionLibrary.java
*********************************************************************/


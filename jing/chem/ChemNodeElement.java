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



package jing.chem;


import java.util.*;

//## package jing::chem 

//----------------------------------------------------------------------------
// jing\chem\ChemNodeElement.java                                                                  
//----------------------------------------------------------------------------

/**
Atom represents the stuff that could be put in the node site in a chem graph.  Atom is designed as an interface inherited by ChemElement, Radical, and Ion.  There are three types of Atom:  
(1) ChemElement: no free radical, no charge on the central chemical element
(2) Radical: with free radical, no chare on the central chemical element
(3) Ion: no free radical, with charge on teh central chemical element. (not implement in this package, but
users can extend to define Ion if they want.)
*/
//## class ChemNodeElement 
public interface ChemNodeElement {
    
    
    //## operation changeRadical(int,String) 
    ChemNodeElement changeRadical(int p_radical, String p_spin);
    
    //## operation getFreeElectron() 
    FreeElectron getFreeElectron();
    
    //## operation getName() 
    String getName();
    
    //## operation isAny() 
    boolean isAny();
    
    //## operation isCarbon() 
    boolean isCarbon();
    
    //## operation isHydrogen() 
    boolean isHydrogen();
    
    //## operation isNonH() 
    boolean isNonH();
    
    //## operation isOxygen() 
    boolean isOxygen();
    
    //## operation isRadical() 
    boolean isRadical();
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\ChemNodeElement.java
*********************************************************************/


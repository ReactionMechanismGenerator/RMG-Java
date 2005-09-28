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


import java.io.*;
import java.util.*;
import jing.chemUtil.*;
import jing.chemParser.*;

//## package jing::chem 

//----------------------------------------------------------------------------
// jing\chem\FGElementDictionary.java                                                                  
//----------------------------------------------------------------------------

//## class FGElementDictionary 
public class FGElementDictionary {
    
    private static FGElementDictionary INSTANCE = new FGElementDictionary();		//## attribute INSTANCE 
    
    protected HashMap dictionary;		//## attribute dictionary 
    
    
    // Constructors
    
    //## operation FGElementDictionary() 
    private  FGElementDictionary() {
        //#[ operation FGElementDictionary() 
        dictionary = new HashMap();
        
        
        //#]
    }
    
    //## operation getFGElement(String) 
    public FGElement getFGElement(String p_name) {
        //#[ operation getFGElement(String) 
        return (FGElement)(dictionary.get(p_name));
        //#]
    }
    
    //## operation getInstance() 
    public static FGElementDictionary getInstance() {
        //#[ operation getInstance() 
        return INSTANCE;
        //#]
    }
    
    //## operation putFGElement(FGElement) 
    public FGElement putFGElement(FGElement p_fGElement) {
        //#[ operation putFGElement(FGElement) 
        return (FGElement)(dictionary.put(p_fGElement.name,p_fGElement));
         
        //#]
    }
    
    public HashMap getDictionary() {
        return dictionary;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FGElementDictionary.java
*********************************************************************/


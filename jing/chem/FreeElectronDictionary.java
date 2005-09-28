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
// jing\chem\FreeElectronDictionary.java                                                                  
//----------------------------------------------------------------------------

//## class FreeElectronDictionary 
public class FreeElectronDictionary {
    
    private static FreeElectronDictionary INSTANCE = new FreeElectronDictionary();		//## attribute INSTANCE 
    
    protected HashMap dictionary;		//## attribute dictionary 
    
    
    // Constructors
    
    //## operation FreeElectronDictionary() 
    private  FreeElectronDictionary() {
        //#[ operation FreeElectronDictionary() 
        dictionary = new HashMap();
        //#]
    }
    
    //## operation getFreeElectron(String) 
    public FreeElectron getFreeElectron(String p_name) {
        //#[ operation getFreeElectron(String) 
        return (FreeElectron)(dictionary.get(p_name));
        //#]
    }
    
    //## operation getInstance() 
    protected static FreeElectronDictionary getInstance() {
        //#[ operation getInstance() 
        return INSTANCE;
        //#]
    }
    
    //## operation putFreeElectron(FreeElectron) 
    public void putFreeElectron(FreeElectron p_freeElectron) {
        //#[ operation putFreeElectron(FreeElectron) 
        dictionary.put(p_freeElectron.name, p_freeElectron);
        //#]
    }
    
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return dictionary.size();
        //#]
    }
    
    public HashMap getDictionary() {
        return dictionary;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FreeElectronDictionary.java
*********************************************************************/


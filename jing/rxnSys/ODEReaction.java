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
// jing\rxnSys\ODEReaction.java                                                                  
//----------------------------------------------------------------------------

//## class ODEReaction 
public class ODEReaction {
    
    protected double A;		//## attribute A 
    
    protected double E;		//## attribute E 
    
    protected double H;		//## attribute H 
    
    protected double Keq;		//## attribute Keq 
    
    protected double alpha;		//## attribute alpha 
    
    protected double dHdT;		//## attribute dHdT 
    
    protected double dKeqdT;		//## attribute dKeqdT 
    
    protected int direction = 0;		//## attribute direction 
    
    protected double n;		//## attribute n 
    
    protected int [] pID = null;		//## attribute pID 
    
    protected int pNum = -1;		//## attribute pNum 
    
    protected int [] rID = null;		//## attribute rID 
    
    protected int rNum = -1;		//## attribute rNum 
    
    protected double rate = 0;		//## attribute rate 
    
    protected boolean setRate = false;		//## attribute setRate 
    
    
    // Constructors
    
    //## operation ODEReaction(int,int,int [],int [],int,double,double,double,double,double,double,double,double) 
    public  ODEReaction(int p_rNum, int p_pNum, int [] p_rID, int [] p_pID, int p_direction, double p_A, double p_n, double p_E, double p_alpha, double p_H, double p_dHdT, double p_Keq, double p_dKeqdT) {
        //#[ operation ODEReaction(int,int,int [],int [],int,double,double,double,double,double,double,double,double) 
        rNum = p_rNum;
        pNum = p_pNum;
        rID = p_rID;
        pID = p_pID;
        direction = p_direction;
        A = p_A;
        n = p_n;
        E = p_E;
        alpha = p_alpha;
        H = p_H;
        dHdT = p_dHdT;
        Keq = p_Keq;
        dKeqdT = p_dKeqdT;
        setRate = false;
        rate = 0;
        
        //#]
    }
    //## operation ODEReaction(int,int,int [],int [],double) 
    public  ODEReaction(int p_rNum, int p_pNum, int [] p_rID, int [] p_pID, double p_rate) {
        //#[ operation ODEReaction(int,int,int [],int [],double) 
        rNum = p_rNum;
        pNum = p_pNum;
        rID = p_rID;
        pID = p_pID;
        direction = 1;
        A = 0;
        n = 0;
        E = 0;
        alpha = 0;
        H = 0;
        dHdT = 0;
        Keq = 0;
        dKeqdT = 0;
        setRate = true;
        rate = p_rate;
        
        //#]
    }
    public  ODEReaction() {
    }
    
    //## operation Message_21() 
    public void Message_21() {
        //#[ operation Message_21() 
        //#]
    }
    
    public double getA() {
        return A;
    }
    
    public double getE() {
        return E;
    }
    
    public double getH() {
        return H;
    }
    
    public double getKeq() {
        return Keq;
    }
    
    public double getAlpha() {
        return alpha;
    }
    
    public double getDHdT() {
        return dHdT;
    }
    
    public double getDKeqdT() {
        return dKeqdT;
    }
    
    public int getDirection() {
        return direction;
    }
    
    public double getN() {
        return n;
    }
    
    public int getPID(int i1) {
        return pID[i1];
    }
    
    public int getPNum() {
        return pNum;
    }
    
    public int getRID(int i1) {
        return rID[i1];
    }
    
    public int getRNum() {
        return rNum;
    }
    
    public double getRate() {
        return rate;
    }
    
    public boolean getSetRate() {
        return setRate;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ODEReaction.java
*********************************************************************/


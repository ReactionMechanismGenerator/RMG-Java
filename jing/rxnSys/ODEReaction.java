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
// jing\rxnSys\ODEReaction.java                                                                  
//----------------------------------------------------------------------------

//## class ODEReaction 
public class ODEReaction {
    
    protected double A;		//## attribute A 
    
    protected double E;		//## attribute E 
    
    protected double H;		//## attribute H 
    
    protected double Keq=0;		//## attribute Keq 
    
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
	
	public String toString(){
		StringBuilder s = new StringBuilder();
		for (int i=0; i<rNum; i++) {
			s.append(rID[i]+"  ");
		}
		s.append("  =  ");
		
		for (int i=0; i<pNum; i++) {
			s.append(pID[i]+"  ");
		}
		
		s.append("    "+ rate +"\n");
		
		return s.toString();
	}
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ODEReaction.java
*********************************************************************/


package jing.rxnSys;

import jing.rxn.ArrheniusKinetics;


public class TROEODEReaction extends ODEReaction {
	  
	  protected int[] colliders = null;
	  protected double[] efficiency = null;
	  protected int numCollider;
	  protected double inertColliderEfficiency;
	  protected double T2star;		//## attribute T2star 
	  
	  protected double T3star;		//## attribute T3star 
	  
	  protected double Tstar;		//## attribute Tstar 
	  
	  protected double a;		//## attribute a 
	  
	  protected double lowRate;
	  protected double highRate; 
	  protected boolean troe7 = false;		//## attribute troe7 
	
	  
	  public TROEODEReaction(int p_rNum, int p_pNum, int [] p_rID, int [] p_pID, int p_direction,  double p_Keq, int [] p_collider, double [] p_efficiency, int p_numCollider, double p_inertColliderEfficiency, double p_T2star, double p_T3star, double p_Tstar, double p_a, double p_highRate, double p_lowRate, boolean p_troe7) {
        //#[ operation ODEReaction(int,int,int [],int [],int,double,double,double,double,double,double,double,double) 
        rNum = p_rNum;
        pNum = p_pNum;
        rID = p_rID;
        pID = p_pID;
        direction = p_direction;
        Keq = p_Keq;
        setRate = false;
        highRate = p_highRate;
		lowRate = p_lowRate;
        colliders = p_collider;
		efficiency = p_efficiency;
		numCollider = p_numCollider;
		inertColliderEfficiency = p_inertColliderEfficiency;
		T2star = p_T2star;
		T3star = p_T3star;
		Tstar = p_Tstar;
		a = p_a;
		troe7 = p_troe7;
        //#]
    }
	  
	  
}

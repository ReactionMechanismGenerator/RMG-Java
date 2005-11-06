package jing.rxnSys;

public class ThirdBodyODEReaction extends ODEReaction {

		protected int[] colliders = null;
		protected double[] efficiency = null;
		protected int numCollider;
		protected double inertColliderEfficiency;
		
		
		public  ThirdBodyODEReaction(int p_rNum, int p_pNum, int [] p_rID, int [] p_pID, int p_direction, double p_A, double p_n, double p_E, double p_alpha, double p_H, double p_dHdT, double p_Keq, double p_dKeqdT, int [] p_collider, double [] p_efficiency, int p_numCollider, double p_inertColliderEfficiency) {
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
	        colliders = p_collider;
			efficiency = p_efficiency;
			numCollider = p_numCollider;
			p_inertColliderEfficiency = inertColliderEfficiency;
	        //#]
	    }
		
//		## operation ODEReaction(int,int,int [],int [],double) 
	    public  ThirdBodyODEReaction(int p_rNum, int p_pNum, int [] p_rID, int [] p_pID, double p_rate, int [] p_collider, double [] p_efficiency, int p_numCollider, double p_inertColliderEfficiency) {
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
			colliders = p_collider;
			efficiency = p_efficiency;
			numCollider = p_numCollider;
			p_inertColliderEfficiency = inertColliderEfficiency;
	        //#]
	    }
	
}

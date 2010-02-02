package jing.rxnSys;

public class LindemannODEReaction extends ODEReaction {
	protected int[] colliders = null;
	protected double[] efficiency = null;
	protected int numCollider;
	protected double inertColliderEfficiency;
	protected double lowRate;
	protected double highRate;
	
	public LindemannODEReaction(int p_rNum, int p_pNum, int [] p_rID, int [] p_pID, int p_direction,  double p_Keq, int [] p_collider, double [] p_efficiency, int p_numCollider, double p_inertColliderEfficiency, double p_highRate, double p_lowRate) {
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
	}
}

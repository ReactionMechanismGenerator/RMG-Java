////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

package jing.rxn;

import jing.rxnSys.CoreEdgeReactionModel;
import jing.rxnSys.ReactionSystem;

/**
 * Represents a generic pressure-dependent kinetics estimator. The original one
 * was Chemdis, but due to a licensing change this is being phased out in favor
 * of FastMasterEqn.
 * @author jwallen
 */
public interface PDepKineticsEstimator {

	/**
	 * Runs a pressure-dependent calculation by preparing the input file,
	 * calling the executable, parsing the output file, and updating the
	 * network/system accordingly.
	 * @param pdn The pressure-dependent reaction network of interest
	 * @param rxnSystem The reaction system of interest
	 */
	public void runPDepCalculation(PDepNetwork pdn, ReactionSystem rxnSystem,
			CoreEdgeReactionModel cerm);
	
}

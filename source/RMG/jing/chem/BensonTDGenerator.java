package jing.chem;

import jing.rxnSys.Logger;

public class BensonTDGenerator extends TDGenerator{

	public BensonTDGenerator(){
		thermoGAPP = GATP.getINSTANCE();
	}

	/**
	 * Generates a thermochemistry object based solely on Benson Group additivity.<BR><BR>
	 * In case an adequate ring correction, or polycyclic ring correction is not found. 
	 * This will return the thermochemistry object WITHOUT corrections but will warn the user for this.
	 */
	@Override
	public ThermoData generateThermo(ChemGraph chemGraph) {
		ThermoData thermo = thermoGAPP.generateThermoData(chemGraph);
		
		boolean fusedPolycyclic = chemGraph.containsFusedRingAtoms();
		
		//Check if cyclic RSCs have been found in case of a cyclic molecule that are non trivial nodes such as six-membered ring.
		if (!chemGraph.isAcyclic() && !fusedPolycyclic && ((GATP)thermoGAPP).getMonoCyclicRSCs().isImperfectMatch()) {
			Logger.warning("Could not find a non trivial ring correction!" +
					"It is advised to review the thermochemistry data of this species.");
		}
		//Check if polycyclic RSCs have been found in case of a polycyclic molecule:
		else if (! chemGraph.isAcyclic() && fusedPolycyclic && ((GATP)thermoGAPP).getPolycyclic() == null){
			Logger.info("Polycyclic ring system with fused ring atoms.");
			Logger.error("Could not find a polycyclic ring strain correction.");                		
			Logger.error("Thermochemistry of polycyclic species does not contain corrections for the" +
					"polyclic ring system. It is advised to review the thermochemistry data of this species.");
		}
	
		return thermo;
	}


}

package jing.chem;

import jing.rxnSys.Logger;

public class HybridTDGenerator extends TDGenerator {
    public HybridTDGenerator() {
        thermoGAPP = GATP.getINSTANCE();
        thermoQM = QMTP.getINSTANCE();
    }

    /**
     * Acyclic molecules are estimated through Benson group additivity. For (poly)cyclic molecules, adequate ring
     * corrections are sought for. If they are not available, fall back to QMTP. If QMTP fails, fall back to Benson GA
     * estimates without ring corrections.
     */
    @Override
    public ThermoData generateThermo(ChemGraph chemGraph) {
        if (chemGraph.isAcyclic()) {
            return thermoGAPP.generateThermoData(chemGraph);
        } else {
            ThermoData thermo = thermoGAPP.generateThermoData(chemGraph);
            boolean fusedPolycyclic = chemGraph.containsFusedRingAtoms();
            BensonRingCorrections monoCyclicRSCs = ((GATP) thermoGAPP)
                    .getMonoCyclicRSCs();
            // Check if cyclic RSCs have been found in case of a cyclic molecule that are non trivial nodes such as
// six-membered ring.
            if (chemGraph.fromprimarythermolibrary) {
                // No action is needed because the thermo comes from a library
            } else if (!fusedPolycyclic && monoCyclicRSCs.isImperfectMatch()) {
                Logger.info("Could not find a non trivial ring correction!"
                        + "Trying QMTP...");
                TDGenerator gen = new QMForCyclicsGenerator();
                return gen.generateThermo(chemGraph);
            }
            // Check if polycyclic RSCs have been found in case of a polycyclic molecule:
            else if (fusedPolycyclic
                    && ((GATP) thermoGAPP).getPolycyclic() == null) {
                Logger.error("Could not find a polycyclic ring strain correction."
                        + "Trying QMTP...");
                TDGenerator gen = new QMForCyclicsGenerator();
                return gen.generateThermo(chemGraph);
            }
            return thermo;
        }
    }
}

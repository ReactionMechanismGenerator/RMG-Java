////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

package jing.rxnSys;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import jing.chem.Species;
import jing.param.Temperature;
import jing.rxn.PDepIsomer;
import jing.rxn.PDepKineticsEstimator;
import jing.rxn.PDepNetwork;
import jing.rxn.PDepReaction;
import jing.rxn.Reaction;

/**
 *
 * @author jwallen
 */
public class RateBasedPDepRME implements ReactionModelEnlarger {

	private PDepKineticsEstimator pDepKineticsEstimator;
	
    // Constructors
    
    public  RateBasedPDepRME() {
		pDepKineticsEstimator = null;
    }
    
	public PDepKineticsEstimator getPDepKineticsEstimator() {
		return pDepKineticsEstimator;
	}
	
	public void setPDepKineticsEstimator(PDepKineticsEstimator pdke) {
		pDepKineticsEstimator = pdke;
	}

	public void enlargeReactionModel(LinkedList rxnSystemList, ReactionModel rm, LinkedList validList) {
		
		CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) rm;
		
		// Iterate over reaction systems, enlarging each individually
		LinkedList updateList = new LinkedList();
        for (int i = 0; i < rxnSystemList.size(); i++) {
			
			// Don't need to enlarge if the system is already valid
			if ((Boolean) validList.get(i))
				continue;
				
			ReactionSystem rxnSystem = (ReactionSystem) rxnSystemList.get(i);
            
			// Select what to add
			PresentStatus ps = rxnSystem.getPresentStatus();
			double Rmin = rxnSystem.getRmin();
			
			// Determine flux of all species (combining both pDep and non-pDep systems)
			int len = cerm.getMaxSpeciesID() + 1;
			double[] flux = new double[len];
			for (int n = 0; n < len; n++)
				flux[n] = 0.0;
			for (Iterator iter = cerm.getUnreactedSpeciesSet().iterator(); iter.hasNext(); ) {
				Species us = (Species) iter.next();
				flux[us.getID()] = Math.abs(ps.getUnreactedSpeciesFlux(us));
			}
			for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
				PDepNetwork pdn = (PDepNetwork) iter.next();
				for (Iterator iter2 = pdn.getNetReactions().iterator(); iter2.hasNext(); ) {
					PDepReaction rxn = (PDepReaction) iter2.next();
					for (int j = 0; j < rxn.getReactantNumber(); j++) {
						Species species = (Species) rxn.getReactantList().get(j);
						if (cerm.containsAsUnreactedSpecies(species))
							flux[species.getID()] += rxn.calculateFlux(ps);
					}
						
				}
			}
			
			// Determine species with maximum flux and its flux
			Species maxSpecies = null;
			double maxFlux = 0;
			for (Iterator iter = cerm.getUnreactedSpeciesSet().iterator(); iter.hasNext(); ) {
				Species us = (Species) iter.next();
				if (Math.abs(flux[us.getID()]) >= maxFlux) {
					maxFlux = Math.abs(flux[us.getID()]);
					maxSpecies = us;
				}
			}
			if (maxSpecies == null) throw new NullPointerException();
			
			// Determine network with maximum leak flux and its flux
			PDepNetwork maxNetwork = null;
			double maxLeak = 0;
			for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
				PDepNetwork pdn = (PDepNetwork) iter.next();
				if (Math.abs(pdn.getLeakFlux(ps)) >= maxLeak) {
					maxLeak = Math.abs(pdn.getLeakFlux(ps));
					maxNetwork = pdn;
				}
			}
			if (maxNetwork == null) throw new NullPointerException();

			// Output results of above calculations to console
			System.out.print("Time: ");
			System.out.println(ps.getTime());
			System.out.println("Rmin: " + String.valueOf(Rmin));
			System.out.println("Unreacted species " + maxSpecies.getName() + " has highest flux: " + String.valueOf(maxFlux));
			System.out.println("Network " + maxNetwork.getID() + " has highest leak flux: " + String.valueOf(maxLeak));

			if (maxFlux > Rmin && (maxFlux > maxLeak || maxLeak < Rmin)) {

				// Add a species to the core
				System.out.print("\nAdd a new reacted Species: ");
				System.out.println(maxSpecies.getChemkinName());
				System.out.println(maxSpecies.toStringWithoutH());
				Temperature temp = new Temperature(715, "K");
				double H = maxSpecies.calculateH(temp);
				double S = maxSpecies.calculateS(temp);
				double G = maxSpecies.calculateG(temp);
				double Cp = maxSpecies.calculateCp(temp);
				System.out.println("Thermo\t" + String.valueOf(H) + " \t" + String.valueOf(S)+ " \t" + String.valueOf(G)+ " \t" + String.valueOf(Cp));

				if (cerm.containsAsReactedSpecies(maxSpecies))
					System.out.println("Species " + maxSpecies.getName() + "(" +
							Integer.toString(maxSpecies.getID()) +
							") is already present in reaction model");
				else {

					// Move the species and appropriate reactions from the edge to the core
					cerm.moveFromUnreactedToReactedSpecies(maxSpecies);
					cerm.moveFromUnreactedToReactedReaction();

					// Generate new reaction set; partition into core and edge
					LinkedHashSet newReactionSet = rxnSystem.getReactionGenerator().react(cerm.getReactedSpeciesSet(),maxSpecies);
					newReactionSet.addAll(rxnSystem.getLibraryReactionGenerator().react(cerm.getReactedSpeciesSet(),maxSpecies));
					Iterator rxnIter = newReactionSet.iterator();
					while (rxnIter.hasNext()){
						Reaction r = (Reaction) rxnIter.next();
						if (r.getReactantNumber() > 1 && r.getProductNumber() > 1)
							cerm.addReaction(r);
						else {
							cerm.categorizeReaction(r.getStructure());
							PDepNetwork.addReactionToNetworks(r);
						}
					}

					// Also make species included in all PDepNetworks for which the species is a unimolecular isomer
					for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
						PDepNetwork pdn = (PDepNetwork) iter.next();
						if (pdn.contains(maxSpecies))
							pdn.makeIsomerIncluded(pdn.getIsomer(maxSpecies));
					}
				}

			}
			else if (maxLeak > Rmin) {

				PDepIsomer isomer = maxNetwork.getMaxLeakIsomer(ps);
				System.out.print("\nAdd a new included Species: " + isomer.toString() +
						" to network " + maxNetwork.getID());

				maxNetwork.makeIsomerIncluded(isomer);

			}
			
			System.out.println("");
        }
    
	}
}

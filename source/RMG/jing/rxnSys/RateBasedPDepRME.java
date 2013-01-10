// //////////////////////////////////////////////////////////////////////////////
//
// RMG - Reaction Mechanism Generator
//
// Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
// RMG Team (rmg_dev@mit.edu)
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// //////////////////////////////////////////////////////////////////////////////
package jing.rxnSys;

import java.util.AbstractSequentialList;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.StringTokenizer;
import jing.chem.Species;
import jing.param.Temperature;
import jing.rxn.FastMasterEqn;
import jing.rxn.PDepException;
import jing.rxn.PDepIsomer;
import jing.rxn.PDepKineticsEstimator;
import jing.rxn.PDepNetwork;
import jing.rxn.PDepReaction;
import jing.rxn.Reaction;
import jing.rxn.Structure;

/**
 * A rate-based reaction model enlarger for use when pressure-dependent kinetics estimation is desired. In addition to
 * checking the species fluxes, this class contains a piece of the activated species algorithm (ASA), a method for
 * ensuring that pressure-dependent networks are being explored in sufficient detail.
 * 
 * @author jwallen
 */
public class RateBasedPDepRME implements ReactionModelEnlarger {
    /**
     * The pressure-dependent kinetics estimator to use. Currently this will only hold a FastMasterEqn object, as the
     * Chemdis class has been depreciated. This is the object called when a pressure-dependent calculation is run.
     */
    private PDepKineticsEstimator pDepKineticsEstimator;

    // ==========================================================================
    //
    // Constructors
    //
    /**
     * Default constructor. Does not set the pressure-dependent kinetics estimator.
     */
    public RateBasedPDepRME() {
        pDepKineticsEstimator = null;
    }

    // ==========================================================================
    //
    // Accessors
    //
    /**
     * Returns the current pressure-dependent kinetics estimator.
     * 
     * @return The current pressure-dependent kinetics estimator
     */
    public PDepKineticsEstimator getPDepKineticsEstimator() {
        return pDepKineticsEstimator;
    }

    /**
     * Sets the pressure-dependent kinetics estimator.
     * 
     * @param pdke
     *            The new pressure-dependent kinetics estimator
     */
    public void setPDepKineticsEstimator(PDepKineticsEstimator pdke) {
        pDepKineticsEstimator = pdke;
    }

    /**
     * Enlarges the reaction model by either adding a species to the core or making a unimolecular isomer included in a
     * PDepNetwork. The action taken is based on the fluxes of species.
     * 
     * @param rxnSystemList
     *            The reaction systems in the simulation
     * @param rm
     *            The current reaction model in the simulation
     * @param validList
     *            A boolean list of the validity status of each reaction system
     */
    public void enlargeReactionModel(LinkedList rxnSystemList,
            ReactionModel rm, LinkedList validList) {
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) rm;
        // Lists of species to add to the core or to explore in PDepNetworks
        // These lists will be of the same length as rxnSystemList, and will
        // contain null where there is no corresponding species; this is done
        // so that we know which reaction system corresponds with each update
        // (although it doesn't matter very much, since they all share the same
        // model core and edge anyway)
        // They should also not contain duplicates, i.e. each candidate species
        // should only appear in one of the two lists, and only once in that
        // list; since adding to core is a superset of making included, species
        // that would be in both lists are placed in coreUpdateList
        LinkedList coreUpdateList = new LinkedList();
        LinkedList leakUpdateList = new LinkedList();
        // Iterate over reaction systems, enlarging each individually
        for (int i = 0; i < rxnSystemList.size(); i++) {
            Logger.info(String.format("Reaction system %d of %d", i + 1,
                    rxnSystemList.size()));
            // Don't need to enlarge if the system is already valid
            if ((Boolean) validList.get(i)) {
                coreUpdateList.add(null);
                leakUpdateList.add(null);
                Logger.info("Model is valid");
                continue;
            }
            ReactionSystem rxnSystem = (ReactionSystem) rxnSystemList.get(i);
            PresentStatus ps = rxnSystem.getPresentStatus();
            // Get Rmin
            double Rmin = rxnSystem.getRmin();
            // Determine flux of all species (combining both pDep and non-pDep systems)
            int len = cerm.getMaxSpeciesID() + 1;
            double[] flux = new double[len];
            for (int n = 0; n < len; n++)
                flux[n] = 0.0;
            // Flux from non-pDep and P-dep reactions
            double[] unreactedSpeciesFlux = ps.getUnreactedSpeciesFlux();// unreacted species flux includes flux from
// both p-dep and non-pdep rxns: cf. appendUnreactedSpeciesStatus in ReactionSystem
            for (Iterator iter = cerm.getUnreactedSpeciesSet().iterator(); iter
                    .hasNext();) {
                Species us = (Species) iter.next();
                if (us.getID() < unreactedSpeciesFlux.length)
                    flux[us.getID()] = unreactedSpeciesFlux[us.getID()];
                else
                    Logger.warning("Attempted to read unreacted species flux for "
                            + us.getName()
                            + "("
                            + us.getID()
                            + "), but there are only "
                            + unreactedSpeciesFlux.length + " fluxes.");
            }
            // Determine species with maximum flux and its flux
            Species maxSpecies = null;
            double maxFlux = 0;
            for (Iterator iter = cerm.getUnreactedSpeciesSet().iterator(); iter
                    .hasNext();) {
                Species us = (Species) iter.next();
                if (Math.abs(flux[us.getID()]) >= maxFlux) {
                    maxFlux = Math.abs(flux[us.getID()]);
                    maxSpecies = us;
                }
            }
            if (maxSpecies == null)
                throw new NullPointerException();
            // Determine species with maximum leak flux and its flux
            Species maxLeakSpecies = null;
            double maxLeakFlux = 0;
            double[] leakFlux = PDepNetwork.getSpeciesLeakFluxes(ps, cerm);
            for (Iterator iter = cerm.getUnreactedSpeciesSet().iterator(); iter
                    .hasNext();) {
                Species us = (Species) iter.next();
                if (Math.abs(leakFlux[us.getID()]) >= maxLeakFlux) {
                    maxLeakFlux = Math.abs(leakFlux[us.getID()]);
                    maxLeakSpecies = us;
                }
            }
            if (maxLeakSpecies == null)
                throw new NullPointerException();
            // Output results of above calculations to console
            Logger.info(String.format("Time: %10.4e s", ps.getTime().getTime()));
            Logger.info(String.format("Rmin: %10.4e mol/cm^3*s", Rmin));
            Logger.info(String
                    .format("Edge species %s has highest flux: %10.4e mol/cm^3*s (%.6f)",
                            maxSpecies.getFullName(), maxFlux, maxFlux / Rmin));
            Logger.info(String
                    .format("Edge species %s has highest leak flux: %10.4e mol/cm^3*s (%.6f)",
                            maxLeakSpecies.getFullName(), maxLeakFlux,
                            maxLeakFlux / Rmin));
            if (maxFlux > maxLeakFlux && maxFlux > Rmin) {
                // Flux is greater than leakFlux, and big enough to matter.
                Logger.info(String.format(
                        "Species %s will be moved to the core.",
                        maxSpecies.getFullName()));
                if (!coreUpdateList.contains(maxSpecies))
                    coreUpdateList.add(maxSpecies);
                else
                    coreUpdateList.add(null);
                leakUpdateList.add(null);
            } else if (maxLeakFlux > Rmin) {
                // leakFlux is greater than Flux, and big enough to matter.
                Logger.info(String
                        .format("The pressure-dependent network of %s will be explored.",
                                maxLeakSpecies.getFullName()));
                if (!leakUpdateList.contains(maxLeakSpecies))
                    leakUpdateList.add(maxLeakSpecies);
                else
                    leakUpdateList.add(null);
                coreUpdateList.add(null);
            } else {
                // neither leakFlux nor Flux is big enough to matter
                leakUpdateList.add(null);
                coreUpdateList.add(null);
            }
        }
        // Check that species don't exist in both coreUpdateList and leakUpdateList
        // If any are in both lists, then remove from leakUpdateList
        for (int i = 0; i < rxnSystemList.size(); i++) {
            if (leakUpdateList.get(i) != null
                    && coreUpdateList.contains(leakUpdateList.get(i)))
                leakUpdateList.set(i, null);
        }
        // Make sure we're about to do something to the reaction model
        boolean found = false;
        for (int i = 0; i < rxnSystemList.size(); i++) {
            if (coreUpdateList.get(i) != null || leakUpdateList.get(i) != null)
                found = true;
        }
        if (!found) {
            Logger.critical("Could not find any species to add to core or leak species to explore. Stopping to avoid infinite loop.");
            System.exit(0);
        }
        // Update the reaction model
        for (int i = 0; i < rxnSystemList.size(); i++) {
            // Add species to core
            if (coreUpdateList.get(i) != null)
                addSpeciesToCore((Species) coreUpdateList.get(i), cerm,
                        (ReactionSystem) rxnSystemList.get(i));
            // Explore species in networks
            else if (leakUpdateList.get(i) != null)
                makeSpeciesIncluded((Species) leakUpdateList.get(i), cerm,
                        (ReactionSystem) rxnSystemList.get(i));
        }
    }

    public void addSpeciesToCore(Species maxSpecies,
            CoreEdgeReactionModel cerm, ReactionSystem rxnSystem) {
        // Add a species to the core
        Logger.info("\nAdd a new species to the model core: "
                + maxSpecies.getFullName());
        Temperature temp = new Temperature(715, "K");
        double H = maxSpecies.calculateH(temp);
        double S = maxSpecies.calculateS(temp);
        double G = maxSpecies.calculateG(temp);
        double Cp = maxSpecies.calculateCp(temp);
        Logger.debug("Thermo of species at 715K (H, S, G, Cp, respectively)\t"
                + String.valueOf(H) + '\t' + String.valueOf(S) + '\t'
                + String.valueOf(G) + '\t' + String.valueOf(Cp));
        if (cerm.containsAsReactedSpecies(maxSpecies))
            Logger.info("Species " + maxSpecies.getFullName()
                    + " is already present in reaction model");
        else {
            // Move the species and appropriate reactions from the edge to the core
            cerm.moveFromUnreactedToReactedSpecies(maxSpecies);
            cerm.moveFromUnreactedToReactedReaction();
            // Adding a species to the core automatically makes it
            // included in all networks it is contained in
            // Therefore we need to merge all networks containing that
            // species as a unimolecular isomer together
            PDepNetwork network = null;
            LinkedList<PDepNetwork> networksToRemove = new LinkedList<PDepNetwork>();
            for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter
                    .hasNext();) {
                PDepNetwork pdn = (PDepNetwork) iter.next();
                if (pdn.contains(maxSpecies)) {
                    if (network == null)
                        network = pdn; // first pdn to contain maxSpecies
                    else { // second or later pdn to contain maxSpecies. merge it with network
                        for (int j = 0; j < pdn.getIsomers().size(); j++)
                            network.addIsomer(pdn.getIsomers().get(j));
                        for (int j = 0; j < pdn.getPathReactions().size(); j++)
                            network.addReaction(pdn.getPathReactions().get(j),
                                    false);
                        networksToRemove.add(pdn);
                    }
                }
            }
            if (network != null) {
                network.getIsomer(maxSpecies).setIncluded(true);
                try {
                    network.updateReactionLists(cerm);
                } catch (PDepException e) {
                    Logger.logStackTrace(e);
                    Logger.error(e.getMessage());
                    Logger.error("Attempt to update reaction list failed "
                            + "for the following network:\n"
                            + network.toString());
                    System.exit(0);
                }
            }
            // Generate new reaction set; partition into core and edge
            LinkedHashSet newReactionSet_nodup;
            if (rxnSystem.getLibraryReactionGenerator().getReactionLibrary() != null) {
                Logger.info("Checking Reaction Library "
                        + rxnSystem.getLibraryReactionGenerator()
                                .getReactionLibrary().getName()
                        + " for reactions of " + maxSpecies.getFullName()
                        + " with the core.");
                // At this point the core (cerm.getReactedSpeciesSet()) already contains maxSpecies, so we can just
// react the entire core.
                LinkedHashSet newReactionSet = rxnSystem
                        .getLibraryReactionGenerator().react(
                                cerm.getReactedSpeciesSet());
                // LinkedHashSet newReactionSet =
// rxnSystem.getLibraryReactionGenerator().react(cerm.getReactedSpeciesSet(),maxSpecies,"All");
                // Report only those that contain the new species (maxSpecies)
                Iterator ReactionIter = newReactionSet.iterator();
                while (ReactionIter.hasNext()) {
                    Reaction current_reaction = (Reaction) ReactionIter.next();
                    if (current_reaction.contains(maxSpecies)) {
                        Logger.info("Library Reaction: "
                                + current_reaction.toString());
                    }
                }
                Logger.info("Generating reactions using reaction family templates.");
                // Iterate through the reaction templates
                newReactionSet.addAll(rxnSystem.getReactionGenerator().react(
                        cerm.getReactedSpeciesSet(), maxSpecies, "All"));
                // To remove the duplicates that are found in Reaction Library and Reaction Template
                // Preference given to Reaction Library over Template Reaction
                newReactionSet_nodup = rxnSystem.getLibraryReactionGenerator()
                        .RemoveDuplicateReac(newReactionSet);
            } else {
                // When no Reaction Library is present
                Logger.info("Generating reactions using reaction family templates.");
                newReactionSet_nodup = rxnSystem.getReactionGenerator().react(
                        cerm.getReactedSpeciesSet(), maxSpecies, "All");
            }
            // shamel 6/22/2010 Suppressed output , line is only for debugging
            // System.out.println("Reaction Set For Pdep PdepRME "+newReactionSet_nodup);
            Iterator rxnIter = newReactionSet_nodup.iterator();
            while (rxnIter.hasNext()) {
                Reaction r = (Reaction) rxnIter.next();
                if (FastMasterEqn.isReactionPressureDependent(r)) {
                    cerm.categorizeReaction(r.getStructure()); // ensure all products are in model edge.
                    PDepNetwork.addReactionToNetworks(r);
                } else
                    cerm.addReaction(r);
            }
        }
        Logger.info("");
    }

    public void makeSpeciesIncluded(Species species,
            CoreEdgeReactionModel cerm, ReactionSystem rxnSystem) {
        Logger.info("\nAdd a new included species: " + species.getFullName());
        LinkedList<PDepNetwork> networksToRemove = new LinkedList<PDepNetwork>();
        // Find the networks that this species occurs in as a nonincluded isomer
        // There should be at least one
        for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter
                .hasNext();) {
            PDepNetwork network = (PDepNetwork) iter.next();
            if (network.contains(species)) {
                PDepIsomer isomer = network.getIsomer(species);
                if (isomer.isUnimolecular() && !isomer.getIncluded()) {
                    networksToRemove.add(network);
                    if (network.getAltered()) {
                        Logger.info("\nNetwork "
                                + network.getID()
                                + " has been altered already this step, so will not be expanded until next step.");
                        return;
                    }
                } else {
                    Logger.info("Isomer "
                            + species.toString()
                            + " is in network #"
                            + network.getID()
                            + ", but is not unimolecular AND nonIncluded,\n\tso RMG will not remove this network.");
                }
            }
        }
        if (networksToRemove.size() == 0) {
            Logger.warning("Tried to make species "
                    + species.toString()
                    + " included, but it was not found as a nonincluded species in any partial network.");
            return;
        }
        // Keep the first identified network; remove the others
        PDepNetwork maxNetwork = networksToRemove.get(0);
        networksToRemove.remove(0);
        try {
            // Making a species included in one network automatically
            // makes it included in all networks it is contained in
            // Therefore we need to merge all networks containing that
            // species as a unimolecular isomer together
            for (Iterator iter = networksToRemove.iterator(); iter.hasNext();) {
                PDepNetwork pdn = (PDepNetwork) iter.next();
                PDepNetwork.getNetworks().remove(pdn);
            }
            // Make the isomer included
            // This will cause any other reactions of the form
            // isomer -> products that don't yet exist to be created
            maxNetwork.makeIsomerIncluded(maxNetwork.getIsomer(species));
            maxNetwork.updateReactionLists(cerm);
        } catch (PDepException e) {
            Logger.logStackTrace(e);
            Logger.error(e.getMessage());
            Logger.verbose(maxNetwork.toString());
            System.exit(0);
        }
    }
}

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
package jing.rxn;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.ListIterator;
import jing.chem.Species;
import jing.rxnSys.CoreEdgeReactionModel;
import jing.rxnSys.ReactionModel;
import jing.rxnSys.ReactionSystem;
import jing.rxnSys.SystemSnapshot;

/**
 * A PDepNetwork object represents a single pressure-dependent reaction network. Such a network is made up of three
 * types of reactions: isomerizations (A1 --> A2), associations (B + C --> A), and dissociations (A --> B + C). Thus
 * each network is defined by
 * <ul>
 * <li>A list of unimolecular isomers { A1, A2, ... }
 * <li>A list of multimolecular isomers { B1 + C1, B2 + C2, ... }
 * <li>A list of path reactions (isomerizations, associations, and dissociations)
 * </ul>
 * PDepNetwork can also be used in a static manner to interact with all of the existing PDepNetwork objects at once.
 * This is a rewrite of the old PDepNetwork class, which was tailored too closely to the CHEMDIS way of treating
 * pressure-dependent networks. The new class is more general and is tailored to the FAME way of treating
 * pressure-dependent networks
 * 
 * @author jwallen
 */
public class PDepNetwork {
    // ==========================================================================
    //
    // Data members
    //
    /**
     * These are used by the network to check the core/edge states of individual isomers and species and to generate
     * pathways when isomers gain included status.
     */
    public static ReactionModel reactionModel;
    public static ReactionSystem reactionSystem;
    /**
     * A count of the number of pressure-dependent networks that have been created since the inception of the current
     * instance of RMG.
     */
    private static int networkCount = 0;
    /**
     * A unique identifier integer for the network.
     */
    private int id;
    /**
     * A hash map containing all of the PDepNetworks created. Each network is identified by the chemical formula of the
     * unimolecular isomers.
     */
    protected static LinkedList<PDepNetwork> networks = new LinkedList<PDepNetwork>();
    /**
     * Set to true if RMG is ready to allow pressure-dependent networks to be created and false if not. A holdover from
     * the original PDepNetwork class.
     */
    public static boolean generateNetworks = false;
    /**
     * The list of unimolecular isomers.
     */
    private LinkedList<PDepIsomer> isomerList;
    /**
     * The list of path reactions (isomerizations, associations, and dissociations that directly connect two isomers).
     */
    private LinkedList<PDepReaction> pathReactionList;
    /**
     * The list of net reactions (allowing reactions between two isomers not directly connected by a path reaction) that
     * belong in the core or the edge of the current reaction model.
     */
    private LinkedList<PDepReaction> netReactionList;
    /**
     * The list of net reactions (allowing reactions between two isomers not directly connected by a path reaction) that
     * are neither in the core nor on the edge of the current reaction model.
     */
    private LinkedList<PDepReaction> nonincludedReactionList;
    /**
     * True if the network has been modified in such a way as to require a new pressure-dependent calculation, and false
     * otherwise. Examples include changing the number of isomers or the number of path reactions.
     */
    private boolean altered;

    // ==========================================================================
    //
    // Constructor
    //
    /**
     * Creates an empty pressure-dependent network. Does not automatically add the network to the PDepNetwork.networks
     * collection.
     */
    public PDepNetwork() {
        isomerList = new LinkedList<PDepIsomer>();
        isomerList.clear();
        pathReactionList = new LinkedList<PDepReaction>();
        pathReactionList.clear();
        netReactionList = new LinkedList<PDepReaction>();
        netReactionList.clear();
        nonincludedReactionList = new LinkedList<PDepReaction>();
        nonincludedReactionList.clear();
        altered = false;
        id = networkCount + 1;
        networkCount++;
    }

    // ==========================================================================
    //
    // Get accessor methods
    //
    /**
     * Returns the unique identifier for this network.
     * 
     * @return The unique identifier for this network.
     */
    public int getID() {
        return id;
    }

    /**
     * Returns a list of all of the species in the network.
     * 
     * @return The list of the species in the network
     */
    public LinkedList<Species> getSpeciesList() {
        LinkedList<Species> speciesList = new LinkedList<Species>();
        for (Iterator<PDepIsomer> iter = isomerList.iterator(); iter.hasNext();) {
            PDepIsomer isomer = iter.next();
            for (int i = 0; i < isomer.getNumSpecies(); i++) {
                if (!speciesList.contains(isomer.getSpecies(i)))
                    speciesList.add(isomer.getSpecies(i));
            }
        }
        return speciesList;
    }

    /**
     * Returns the list of isomers.
     * 
     * @return The list of isomers
     */
    public LinkedList<PDepIsomer> getIsomers() {
        return isomerList;
    }

    public int getNumUniIsomers() {
        int count = 0;
        for (Iterator<PDepIsomer> iter = isomerList.iterator(); iter.hasNext();) {
            PDepIsomer isomer = iter.next();
            if (isomer.isUnimolecular())
                count++;
        }
        return count;
    }

    public int getNumMultiIsomers() {
        return isomerList.size() - getNumUniIsomers();
    }

    /**
     * Returns the list of path reactions.
     * 
     * @return The list of path reactions
     */
    public LinkedList<PDepReaction> getPathReactions() {
        return pathReactionList;
    }

    public void removeFromPathReactionList(PDepReaction pdr) {
        pathReactionList.remove(pdr);
        return;
    }

    /**
     * Returns the list of net reactions (in the core or on the edge).
     * 
     * @return The list of net reactions
     */
    public LinkedList<PDepReaction> getNetReactions() {
        return netReactionList;
    }

    public void removeFromNetReactionList(PDepReaction pdr) {
        netReactionList.remove(pdr);
        return;
    }

    /**
     * Returns the list of nonincluded reactions (neither in the core nor on the edge).
     * 
     * @return The list of nonincluded reactions
     */
    public LinkedList<PDepReaction> getNonincludedReactions() {
        return nonincludedReactionList;
    }

    public void removeFromNonincludedReactionList(PDepReaction pdr) {
        nonincludedReactionList.remove(pdr);
        return;
    }

    /**
     * Returns the status of the altered flag: true if the network requires a new pressure-dependent calculation, false
     * if not.
     * 
     * @return The status of the altered flag
     */
    public boolean getAltered() {
        return altered;
    }

    /**
     * Returns the isomer that contains the same species as those in speciesList.
     * 
     * @param speciesList
     *            The species to check the isomers for
     * @return The isomer that contains the same species as those in speciesList
     */
    public PDepIsomer getIsomer(LinkedList speciesList) {
        if (speciesList.size() == 1) {
            for (ListIterator<PDepIsomer> iter = isomerList.listIterator(); iter
                    .hasNext();) {
                PDepIsomer isomer = iter.next();
                if (isomer.getSpeciesList().equals(speciesList))
                    return isomer;
            }
        }
        return null;
    }

    /**
     * Returns the unimolecular isomer that contains the indicated species.
     * 
     * @param species
     *            The species to check the isomers for
     * @return The unimolecular isomer corresponding to species
     */
    public PDepIsomer getIsomer(Species species) {
        for (ListIterator<PDepIsomer> iter = isomerList.listIterator(); iter
                .hasNext();) {
            PDepIsomer isomer = iter.next();
            if (isomer.getSpecies(0).equals(species) && isomer.isUnimolecular())
                return isomer;
        }
        return null;
    }

    public void removeFromIsomerList(PDepIsomer pdi) {
        isomerList.remove(pdi);
        return;
    }

    /**
     * Checks to see if the network contains species as a unimolecular isomer.
     * 
     * @param species
     *            The species to check for
     * @return True if species is a unimolecular isomer in the network, false if not
     */
    public boolean contains(Species species) {
        for (ListIterator<PDepIsomer> iter = isomerList.listIterator(); iter
                .hasNext();) {
            PDepIsomer isomer = iter.next();
            if (isomer.getSpecies(0).equals(species) && isomer.isUnimolecular())
                return true;
        }
        return false;
    }

    /**
     * Checks to see if the network contains species as a unimolecular isomer.
     * 
     * @param species
     *            The species to check for
     * @return True if species is a unimolecular isomer in the network, false if not
     */
    public boolean contains(Reaction reaction) {
        if (reaction == null)
            return false;
        Reaction reverse = reaction.getReverseReaction();
        for (ListIterator<PDepReaction> iter = pathReactionList.listIterator(); iter
                .hasNext();) {
            PDepReaction rxn = iter.next();
            if (reaction.equals(rxn))
                return true;
            if (reverse != null) {
                if (reverse.equals(rxn))
                    return true;
            }
        }
        return false;
    }

    // ==========================================================================
    //
    // Set accessor methods
    //
    /**
     * Adds an isomer (unimolecular or multimolecular) to the appropriate list in the network if it is not already
     * present.
     * 
     * @param isomer
     *            The isomer to add
     */
    public void addIsomer(PDepIsomer isomer) {
        // Don't add if isomer is already in network
        if (isomerList.contains(isomer))
            return;
        // Add isomer: keep unimolecular isomers first
        if (isomer.isUnimolecular()) {
            int index = -1;
            for (int i = 0; i < isomerList.size(); i++) {
                if (isomerList.get(i).isMultimolecular() && index < 0)
                    index = i;
            }
            if (index >= 0)
                isomerList.add(index, isomer);
            else
                isomerList.add(isomer);
        } else
            isomerList.add(isomer);
        // Mark network as changed so that updated rates can be determined
        altered = true;
    }

    /**
     * Adds a path reaction to the network if it is not already present.
     * 
     * @param newRxn
     *            The path reaction to add
     * @param addKinetics
     *            : In the event the reaction is already present, determine whether to add the kinetics to the already
     *            present rxn
     */
    public void addReaction(PDepReaction newRxn, boolean addKinetics) {
        // Check to ensure that reaction is not already present
        for (ListIterator<PDepReaction> iter = pathReactionList.listIterator(); iter
                .hasNext();) {
            PDepReaction rxn = iter.next();
            if (rxn.equals(newRxn)) {
                if (addKinetics)
                    rxn.addAdditionalKinetics(newRxn.getKinetics()[0], 1, false);
                return;
            }
        }
        // Add reaction
        pathReactionList.add(newRxn);
        // Mark network as changed so that updated rates can be determined
        altered = true;
    }

    /**
     * Updates the status of the altered flag: true if the network requires a new pressure-dependent calculation, false
     * if not.
     * 
     * @param alt
     *            The new status of the altered flag
     */
    public void setAltered(boolean alt) {
        altered = alt;
    }

    /**
     * Elevates the status of the designated isomer from nonincluded to included, and generates pathways for this
     * isomer. Generally a large number of pathways are generated.
     * 
     * @param isomer
     *            The isomer to make included.
     */
    public void makeIsomerIncluded(PDepIsomer isomer) {
        if (!isomer.isUnimolecular() || isomer.getIncluded())
            return;
        isomer.setIncluded(true);
        LinkedHashSet reactionSet = isomer.generatePaths(reactionSystem);
        for (Iterator iter = reactionSet.iterator(); iter.hasNext();) {
            Reaction rxn = (Reaction) iter.next();
            if (!contains(rxn)) {
                addReactionToNetworks(rxn);
                /*
                 * The reactions (of form A --> B or A --> C + D) could form a species that is not otherwise in the edge
                 * of the CoreEdgeReactionModel. We want to leave all of the reactions alone (i.e. not add them to the
                 * core OR edge) but need to check whether any of the new species should be added to the edge.
                 */
                LinkedList rxnProducts = rxn.getProductList();
                for (int numProds = 0; numProds < rxnProducts.size(); numProds++) {
                    if (!((CoreEdgeReactionModel) reactionModel)
                            .containsAsReactedSpecies((Species) rxnProducts
                                    .get(numProds)))
                        if (!((CoreEdgeReactionModel) reactionModel)
                                .containsAsUnreactedSpecies((Species) rxnProducts
                                        .get(numProds)))
                            ((CoreEdgeReactionModel) reactionModel)
                                    .addUnreactedSpecies((Species) rxnProducts
                                            .get(numProds));
                }
            }
        }
    }

    // ==========================================================================
    //
    // Other methods
    //
    /**
     * Redistributes the net reactions based on the current core and edge reaction models. Especially useful when one or
     * more species has been moved from the edge to the core since the last pressure-dependent calculation.
     * 
     * @param cerm
     *            The current core/edge reaction model
     */
    public void updateReactionLists(CoreEdgeReactionModel cerm)
            throws PDepException {
        // Merge the net reaction and nonincluded reactinon lists together
        // We will recalculate how to distribute them
        LinkedList<PDepReaction> reactionList = new LinkedList<PDepReaction>();
        for (int i = 0; i < netReactionList.size(); i++) {
            PDepReaction rxn = netReactionList.get(i);
            reactionList.add(rxn);
        }
        for (int i = 0; i < nonincludedReactionList.size(); i++) {
            PDepReaction rxn = nonincludedReactionList.get(i);
            reactionList.add(rxn);
        }
        netReactionList.clear();
        nonincludedReactionList.clear();
        for (int i = 0; i < reactionList.size(); i++) {
            PDepReaction forward = reactionList.get(i);
            PDepReaction reverse = (PDepReaction) forward.getReverseReaction();
            if (reverse != null) {
                if (forward.isCoreReaction(cerm)
                        || reverse.isCoreReaction(cerm))
                    netReactionList.add(forward);
                else if (forward.getReactant().getIncluded()
                        && forward.getProduct().getIncluded())
                    netReactionList.add(forward);
                else
                    nonincludedReactionList.add(forward);
            } else {
                if (forward.isCoreReaction(cerm))
                    netReactionList.add(forward);
                else if (forward.getReactant().getIncluded()
                        && forward.getProduct().getIncluded())
                    netReactionList.add(forward);
                else
                    nonincludedReactionList.add(forward);
            }
        }
    }

    /**
     * Calculates the leak flux for this network. The leak flux is defined as the maximum possible flux to all
     * nonincluded species. The maximum modifier implies that only the forward reaction to the nonincluded species is
     * used to generate the flux, rather than the combined forward and backward reaction.
     * 
     * @param ss
     *            A system snapshot (T, P, concentrations, etc.) to use to calculate the flux.
     * @return The leak flux for this network
     */
    public double getLeakFlux(SystemSnapshot ss) {
        double rLeak = 0.0;
        // If there is only one path reaction (and thus only one nonincluded
        // reaction), use the high-pressure limit rate as the flux rather than
        // the k(T,P) value to ensure that we are considering the maximum
        // possible flux entering the network
        if (pathReactionList.size() == 1 && netReactionList.size() == 0) {
            PDepReaction rxn = pathReactionList.get(0);
            if (!rxn.getProduct().getIncluded())
                rLeak += rxn.calculateForwardFlux(ss);
            else if (!rxn.getReactant().getIncluded())
                rLeak += rxn.calculateReverseFlux(ss);
            else
                // If both the reactant and the product are included, then the
                // leak flux is zero
                rLeak = 0.0;
        }
        // Otherwise use the set of k(T,P) values
        else {
            for (ListIterator<PDepReaction> iter = nonincludedReactionList
                    .listIterator(); iter.hasNext();) {
                PDepReaction rxn = iter.next();
                if (rxn.getReactant().getIncluded()
                        && !rxn.getProduct().getIncluded())
                    rLeak += rxn.calculateForwardFlux(ss);
                else if (!rxn.getReactant().getIncluded()
                        && rxn.getProduct().getIncluded())
                    rLeak += rxn.calculateReverseFlux(ss);
            }
        }
        return rLeak;
    }

    public static double[] getSpeciesLeakFluxes(SystemSnapshot ss,
            CoreEdgeReactionModel cerm) {
        int len = cerm.getMaxSpeciesID() + 1;
        double[] leakFlux = new double[len];
        for (int n = 0; n < len; n++)
            leakFlux[n] = 0.0;
        for (ListIterator<PDepNetwork> iter0 = networks.listIterator(); iter0
                .hasNext();) {
            PDepNetwork pdn = iter0.next();
            // If there is only one path reaction (and thus only one nonincluded
            // reaction), use the high-pressure limit rate as the flux rather than
            // the k(T,P) value to ensure that we are considering the maximum
            // possible flux entering the network
            if (pdn.getPathReactions().size() == 1
                    && pdn.getNetReactions().size() == 0) {
                PDepReaction rxn = pdn.getPathReactions().get(0);
                if (!rxn.getProduct().getIncluded())
                    leakFlux[rxn.getProduct().getSpecies(0).getID()] += rxn
                            .calculateForwardFlux(ss);
                else if (!rxn.getReactant().getIncluded())
                    leakFlux[rxn.getReactant().getSpecies(0).getID()] += rxn
                            .calculateReverseFlux(ss);
            }
            // Otherwise use the set of k(T,P) values
            else {
                for (ListIterator<PDepReaction> iter = pdn
                        .getNonincludedReactions().listIterator(); iter
                        .hasNext();) {
                    PDepReaction rxn = iter.next();
                    if (rxn.getReactant().getIncluded()
                            && !rxn.getProduct().getIncluded())
                        leakFlux[rxn.getProduct().getSpecies(0).getID()] += rxn
                                .calculateForwardFlux(ss);
                    else if (!rxn.getReactant().getIncluded()
                            && rxn.getProduct().getIncluded())
                        leakFlux[rxn.getReactant().getSpecies(0).getID()] += rxn
                                .calculateReverseFlux(ss);
                }
            }
        }
        return leakFlux;
    }

    /**
     * Calculates the isomer with the largest leak flux for this network. The reaction with the maximum flux is used to
     * select the isomer. This isomer is the candidate for elevating to included status.
     * 
     * @param ss
     *            A system snapshot (T, P, concentrations, etc.) to use to calculate the flux.
     * @return The isomer with the largest leak flux
     */
    public PDepIsomer getMaxLeakIsomer(SystemSnapshot ss) throws PDepException {
        if (nonincludedReactionList.size() == 0) {
            if (pathReactionList.size() == 1 && isomerList.size() == 2) {
                PDepIsomer isomer1 = isomerList.get(0);
                PDepIsomer isomer2 = isomerList.get(1);
                if (isomer1.isUnimolecular() && !isomer1.getIncluded())
                    return isomer1;
                else if (isomer2.isUnimolecular() && !isomer2.getIncluded())
                    return isomer2;
            }
            throw new PDepException(
                    "Tried to determine nonincluded isomer with maximum leak flux, but there are no nonincluded reactions, so no isomer can be identified.");
        }
        PDepReaction maxReaction = null;
        double maxLeak = 0.0;
        for (ListIterator<PDepReaction> iter = nonincludedReactionList
                .listIterator(); iter.hasNext();) {
            PDepReaction rxn = iter.next();
            if (!rxn.getReactant().getIncluded()
                    || !rxn.getProduct().getIncluded()) {
                if (Math.abs(rxn.calculateFlux(ss)) > maxLeak) {
                    maxReaction = rxn;
                    maxLeak = rxn.calculateFlux(ss);
                }
            }
        }
        if (maxReaction == null)
            throw new PDepException(
                    "Tried to determine nonincluded isomer with maximum leak flux, but no suitable nonincluded reaction has been found.");
        else if (!maxReaction.getReactant().getIncluded())
            return maxReaction.getReactant();
        else if (!maxReaction.getProduct().getIncluded())
            return maxReaction.getProduct();
        else
            throw new PDepException(
                    "Tried to determine nonincluded isomer with maximum leak flux, but nonincluded reaction with maximum leak flux has no nonincluded isomers.");
    }

    public String getSpeciesType() {
        for (Iterator<PDepIsomer> iter = isomerList.iterator(); iter.hasNext();) {
            PDepIsomer isomer = iter.next();
            if (isomer.isUnimolecular())
                return isomer.getSpecies(0).getName();
        }
        return "";
    }

    public String toString() {
        String str = "PDepNetwork #" + Integer.toString(id) + ":\n";
        str += "\tIsomers:\n";
        for (ListIterator<PDepIsomer> iter = isomerList.listIterator(); iter
                .hasNext();) {
            PDepIsomer isomer = iter.next();
            str += "\t\t" + isomer.toString() + "\n";
        }
        str += "\tPath reactions:\n";
        for (ListIterator<PDepReaction> iter = pathReactionList.listIterator(); iter
                .hasNext();) {
            PDepReaction rxn = iter.next();
            str += "\t\t" + rxn.toString() + "\n";
        }
        str += "\tNet reactions:\n";
        for (ListIterator<PDepReaction> iter = netReactionList.listIterator(); iter
                .hasNext();) {
            PDepReaction rxn = iter.next();
            str += "\t\t" + rxn.toString() + "\n";
        }
        str += "\tNonincluded reactions:\n";
        for (ListIterator<PDepReaction> iter = nonincludedReactionList
                .listIterator(); iter.hasNext();) {
            PDepReaction rxn = iter.next();
            str += "\t\t" + rxn.toString() + "\n";
        }
        return str;
    }

    // ==========================================================================
    //
    // Static methods (for access to PDepNetwork.networks)
    //
    /**
     * Returns the linked list containing the currently-existing pressure- dependent networks
     * 
     * @return The currently-existing pressure-dependent networks
     */
    public static LinkedList<PDepNetwork> getNetworks() {
        return networks;
    }

    /**
     * Used to add a reaction to the appropriate pressure-dependent network. If no such network exists, a new network is
     * created. For isomerization reactions connecting two existing networks, the networks are merged. This function is
     * to be called whenever a new reaction is added to the edge.
     * 
     * @param reaction
     *            The reaction to add
     * @return The network the reaction was added to
     */
    public static PDepNetwork addReactionToNetworks(Reaction reaction) {
        // Expect that most reactions passed to this function will be already
        // present in a network
        Reaction reaction0 = reaction;
        // Fail if neither reactant nor product are unimolecular
        Species species = null;
        if (reaction.getReactantNumber() == 1)
            species = (Species) reaction.getReactantList().get(0);
        if (reaction.getProductNumber() == 1)
            species = (Species) reaction.getProductList().get(0);
        if (species == null)
            return null;
        if (reaction.getReactantNumber() > 1)
            reaction = reaction.getReverseReaction();
        	if(reaction == null)
        		System.out.println("Irreversible reaction sent to pdep module, most likely it is comming from one of the libraries "+reaction0.toString());
        PDepNetwork pdn = null;
        if (reaction.getProductNumber() == 1) {
            // Isomerization reactions should cause networks to be merged together
            // This means that each unimolecular isomer should only appear in one network
            // Get the appropriate pressure-dependent network(s)
            PDepNetwork reac_pdn = null;
            PDepNetwork prod_pdn = null;
            Species reactant = (Species) reaction.getReactantList().get(0);
            Species product = (Species) reaction.getProductList().get(0);
            for (ListIterator<PDepNetwork> iter = networks.listIterator(); iter
                    .hasNext();) {
                PDepNetwork n = iter.next();
                if (n.contains(reactant)) {
                    if (n.getIsomer(reactant).getIncluded()) {
                        reac_pdn = n;
                        if (prod_pdn != null)
                            break; // have now found both prod_pdn and reac_pdn.
                    }
                }
                if (n.contains(product)) {
                    if (n.getIsomer(product).getIncluded()) {
                        prod_pdn = n;
                        if (reac_pdn != null)
                            break; // have now found both reac_pdn and prod_pdn.
                    }
                }
            }
            if (reac_pdn != null && prod_pdn != null && reac_pdn != prod_pdn) {
                // Two distinct networks found; must join them together
                pdn = reac_pdn;
                for (int i = 0; i < prod_pdn.getIsomers().size(); i++)
                    pdn.addIsomer(prod_pdn.getIsomers().get(i));
                for (int i = 0; i < prod_pdn.getPathReactions().size(); i++)
                    pdn.addReaction(prod_pdn.getPathReactions().get(i), false);
                // Also remove the second network from the list of networks
                networks.remove(prod_pdn);
            } else if (reac_pdn != null && prod_pdn != null
                    && reac_pdn == prod_pdn) {
                // Both species already present as unimolecular isomers in the same network, so use that network
                pdn = reac_pdn;
            } else if (reac_pdn != null) {
                // Only reactant species found in a network, so use that network
                pdn = reac_pdn;
            } else if (prod_pdn != null) {
                // Only product species found in a network, so use that network
                pdn = reac_pdn;
            } else {
                // No networks found for either species; will create a new network
                pdn = null;
            }
        } else if (reaction.getProductNumber() > 1) {
            // Dissociation reactions are added to the network containing that unimolecular isomer
            // Since each unimolecular isomer should only appear in one network, there should only be one such addition
            // If no existing network is found, a new one may be created
            // Get the appropriate pressure-dependent network
            Species reactant = (Species) reaction.getReactantList().get(0);
            for (ListIterator<PDepNetwork> iter = networks.listIterator(); iter
                    .hasNext();) {
                PDepNetwork n = iter.next();
                if (n.contains(reactant)) {
                    if (n.getIsomer(reactant).getIncluded())
                        pdn = n;
                }
            }
        }
        // The above check may have caused the reaction to be reversed
        // We want to add the reaction to the network in the direction in which
        // kinetics is known, so make sure we are using the reaction as
        // originally passed to this method
        reaction = reaction0;
        // If network not found, create a new network
        if (pdn == null) {
            pdn = new PDepNetwork();
            PDepIsomer isomer = new PDepIsomer(species);
            pdn.addIsomer(isomer);
            networks.add(pdn);
        }
        // Add the reaction to the network
        PDepIsomer reactantIsomer = pdn.getIsomer(reaction.getReactantList());
        if (reactantIsomer == null) {
            reactantIsomer = new PDepIsomer(reaction.getReactantList());
            pdn.addIsomer(reactantIsomer);
        }
        PDepIsomer productIsomer = pdn.getIsomer(reaction.getProductList());
        if (productIsomer == null) {
            productIsomer = new PDepIsomer(reaction.getProductList());
            pdn.addIsomer(productIsomer);
        }
        // Always add the reaction in the direction for which we have the kinetics
        PDepReaction rxn = null;
        if (reaction.isForward()) {
            rxn = new PDepReaction(reactantIsomer, productIsomer, reaction);
            rxn.setRedundancy(reaction.getRedundancy());
        } else {
            rxn = new PDepReaction(productIsomer, reactantIsomer,
                    reaction.getReverseReaction());
            rxn.setRedundancy(reaction.getReverseReaction().getRedundancy());
        }
        pdn.addReaction(rxn, false);
        // Fill in partial network if necessary
        if (reactantIsomer.isCore((CoreEdgeReactionModel) reactionModel)
                && reactantIsomer.isUnimolecular()
                && !reactantIsomer.getIncluded())
            pdn.makeIsomerIncluded(reactantIsomer);
        if (productIsomer.isCore((CoreEdgeReactionModel) reactionModel)
                && productIsomer.isUnimolecular()
                && !productIsomer.getIncluded())
            pdn.makeIsomerIncluded(productIsomer);
        // Return the created network
        return pdn;
    }

    /**
     * Useful for debugging, this function prints the isomers of each network to the console window.
     */
    public static void printNetworks() {
        int index = 0;
        for (ListIterator<PDepNetwork> iter0 = networks.listIterator(); iter0
                .hasNext();) {
            PDepNetwork pdn = iter0.next();
            index++;
            System.out.print("Network #" + Integer.toString(index) + ": ");
            for (ListIterator<PDepIsomer> iter = pdn.getIsomers()
                    .listIterator(); iter.hasNext();) {
                PDepIsomer isomer = iter.next();
                System.out.print(isomer.toString());
                if (iter.hasNext())
                    System.out.print(", ");
            }
            System.out.print("\n");
        }
    }

    /**
     * Checks to see if there are any core reactions hidden amongst those net reactions which are found in the
     * pressure-dependent networks. This is particularly useful in the initialization of the reaction model, in which
     * the core must have at least one reaction in it before the dynamic simulator can be executed.
     * 
     * @param cerm
     *            The current core/edge reaction model
     * @return True if core reactions are found, false if not
     */
    public static boolean hasCoreReactions(CoreEdgeReactionModel cerm) {
        return (getCoreReactions(cerm).size() > 0);
    }

    /**
     * Counts the number of core reactions that are hidden amongst those net reactions which are found in the
     * pressure-dependent networks. This is particularly useful in the initialization of the reaction model, in which
     * the core must have at least one reaction in it before the dynamic simulator can be executed.
     * 
     * @param cerm
     *            The current core/edge reaction model
     * @return The number of core reactions found
     */
    public static int getNumCoreReactions(CoreEdgeReactionModel cerm) {
        return getCoreReactions(cerm).size();
    }

    /**
     * Returns the core reactions that are hidden amongst those net reactions which are found in the pressure-dependent
     * networks. This is particularly useful in the initialization of the reaction model, in which the core must have at
     * least one reaction in it before the dynamic simulator can be executed.
     * 
     * @param cerm
     *            The current core/edge reaction model
     * @return The number of core reactions found
     */
    public static LinkedList<PDepReaction> getCoreReactions(
            CoreEdgeReactionModel cerm) {
        LinkedList<PDepReaction> coreReactions = new LinkedList<PDepReaction>();
        for (ListIterator<PDepNetwork> iter0 = networks.listIterator(); iter0
                .hasNext();) {
            PDepNetwork pdn = iter0.next();
            for (ListIterator<PDepReaction> iter = pdn.getNetReactions()
                    .listIterator(); iter.hasNext();) {
                PDepReaction rxn = iter.next();
                if (rxn.isCoreReaction(cerm) && !coreReactions.contains(rxn))
                    coreReactions.add(rxn);
            }
        }
        return coreReactions;
    }

    /**
     * Counts the number of edge reactions that are hidden amongst those net reactions which are found in the
     * pressure-dependent networks.
     * 
     * @param cerm
     *            The current core/edge reaction model
     * @return The number of edge reactions found
     */
    public static int getNumEdgeReactions(CoreEdgeReactionModel cerm) {
        return getEdgeReactions(cerm).size();
    }

    /**
     * Returns the edge reactions that are hidden amongst those net reactions which are found in the pressure-dependent
     * networks.
     * 
     * @param cerm
     *            The current core/edge reaction model
     * @return The list of edge reactions found
     */
    public static LinkedList<PDepReaction> getEdgeReactions(
            CoreEdgeReactionModel cerm) {
        LinkedList<PDepReaction> edgeReactions = new LinkedList<PDepReaction>();
        for (ListIterator<PDepNetwork> iter0 = networks.listIterator(); iter0
                .hasNext();) {
            PDepNetwork pdn = iter0.next();
            for (ListIterator<PDepReaction> iter = pdn.getNetReactions()
                    .listIterator(); iter.hasNext();) {
                PDepReaction rxn = iter.next();
                if (rxn.getReactant().getIncluded()
                        && rxn.getProduct().getIncluded()) {
                    if (rxn.isEdgeReaction(cerm)
                            && !edgeReactions.contains(rxn))
                        edgeReactions.add(rxn);
                }
            }
        }
        return edgeReactions;
    }

    /**
     * Counts the number of total path reactions in the pressure-dependent networks.
     * 
     * @param cerm
     *            The current core/edge reaction model
     * @return The number of path reactions found
     */
    public static int getNumPathReactions(CoreEdgeReactionModel cerm) {
        int count = 0;
        for (ListIterator<PDepNetwork> iter0 = networks.listIterator(); iter0
                .hasNext();) {
            PDepNetwork pdn = iter0.next();
            count += pdn.getPathReactions().size();
        }
        return count;
    }

    /**
     * Counts the number of total net reactions in the pressure-dependent networks, including all core-to-core ("core"),
     * core-to-edge ("edge"), and edge-to-edge reactions.
     * 
     * @param cerm
     *            The current core/edge reaction model
     * @return The number of net reactions found
     */
    public static int getNumNetReactions(CoreEdgeReactionModel cerm) {
        int count = 0;
        for (ListIterator<PDepNetwork> iter0 = networks.listIterator(); iter0
                .hasNext();) {
            PDepNetwork pdn = iter0.next();
            count += pdn.getNetReactions().size();
        }
        return count;
    }

    /**
     * Check whether or not a given species is an included (fully explored) unimolecular isomer in any
     * currently-existing network.
     * 
     * @param species
     *            The species to check for included status
     * @return true if the species is included in any existing network, false if not
     */
    public static boolean isSpeciesIncludedInAnyNetwork(Species species) {
        for (Iterator iter = networks.iterator(); iter.hasNext();) {
            PDepNetwork network = (PDepNetwork) iter.next();
            if (network.contains(species)) {
                PDepIsomer isomer = network.getIsomer(species);
                if (isomer.isUnimolecular() && isomer.getIncluded())
                    // We've identified a network wherein the species exists as
                    // a unimolecular isomer, and that its path reactions have
                    // been fully explored
                    // This satisfies all of the conditions, so we return true
                    return true;
            }
        }
        // No suitable match for all conditions was found, so we return false
        return false;
    }
}

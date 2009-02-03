////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

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
 * A PDepNetwork object represents a single pressure-dependent reaction network.
 * Such a network is made up of three types of reactions: isomerizations
 * (A1 --> A2), associations (B + C --> A), and dissociations (A --> B + C).
 * Thus each network is defined by
 * <ul>
 * <li>A list of unimolecular isomers { A1, A2, ... }
 * <li>A list of multimolecular isomers { B1 + C1, B2 + C2, ... }
 * <li>A list of path reactions (isomerizations, associations, and dissociations)
 * </ul>
 * PDepNetwork can also be used in a static manner to interact with all of the
 * existing PDepNetwork objects at once.
 * 
 * This is a rewrite of the old PDepNetwork class, which was tailored too 
 * closely to the CHEMDIS way of treating pressure-dependent networks. The new
 * class is more general and is tailored to the FAME way of treating 
 * pressure-dependent networks
 * @author jwallen
 */
public class PDepNetwork {

	//==========================================================================
	//
	//	Data members
	//
	
	public static ReactionModel reactionModel;
	public static ReactionSystem reactionSystem;

	private static int networkCount = 0;

	private int id;

	/**
	 * A hash map containing all of the PDepNetworks created. Each network is
	 * identified by the chemical formula of the unimolecular isomers.
	 */
	protected static LinkedList<PDepNetwork> networks = new LinkedList<PDepNetwork>();
    
	/**
	 * Set to true if RMG is ready to allow pressure-dependent networks to be 
	 * created and false if not. A holdover from the original PDepNetwork class. 
	 */
	public static boolean generateNetworks = false;
    
	/**
	 * The list of unimolecular isomers.
	 */
	private LinkedList<PDepIsomer> uniIsomerList;
	
	/**
	 * The list of multimolecular isomers.
	 */
	private LinkedList<PDepIsomer> multiIsomerList;
	
	/**
	 * The list of path reactions (isomerizations, associations, and 
	 * dissociations that directly connect two isomers).
	 */
	private LinkedList<PDepReaction> pathReactionList;
	
	/**
	 * The list of net reactions (allowing reactions between two isomers not
	 * directly connected by a path reaction) that belong in the core or the
	 * edge of the current reaction model.
	 */
	private LinkedList<PDepReaction> netReactionList;
	
	/**
	 * The list of net reactions (allowing reactions between two isomers not
	 * directly connected by a path reaction) that are neither in the core nor
	 * on the edge of the current reaction model.
	 */
	private LinkedList<PDepReaction> nonincludedReactionList;
	
	/**
	 * True if the network has been modified in such a way as to require a
	 * new pressure-dependent calculation, and false otherwise. Examples include
	 * changing the number of isomers or the number of path reactions.
	 */
	private boolean altered;
	
	//==========================================================================
	//
	//	Constructor
	//
	
	/**
	 * Creates an empty pressure-dependent network. Does not automatically add
	 * the network to the PDepNetwork.networks collection.
	 */
	public PDepNetwork() {
		uniIsomerList = new LinkedList<PDepIsomer>();
		uniIsomerList.clear();
		multiIsomerList = new LinkedList<PDepIsomer>();
		multiIsomerList.clear();
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
	
	//==========================================================================
	//
	//	Get accessor methods
	//
	
	public int getID() {
		return id;
	}

	/**
	 * Returns the list of unimolecular isomers.
	 * @return The list of unimolecular isomers
	 */
	public LinkedList<PDepIsomer> getUniIsomers() {
		return uniIsomerList;
	}
	
	/**
	 * Returns the list of multimolecular isomers.
	 * @return The list of multimolecular isomers
	 */
	public LinkedList<PDepIsomer> getMultiIsomers() {
		return multiIsomerList;
	}
	
	/**
	 * Returns the list of path reactions.
	 * @return The list of path reactions
	 */
	public LinkedList<PDepReaction> getPathReactions() {
		return pathReactionList;
	}
	
	/**
	 * Returns the list of net reactions (in the core or on the edge).
	 * @return The list of net reactions
	 */
	public LinkedList<PDepReaction> getNetReactions() {
		return netReactionList;
	}
	
	/**
	 * Returns the list of nonincluded reactions (neither in the core nor on 
	 * the edge).
	 * @return The list of nonincluded reactions
	 */
	public LinkedList<PDepReaction> getNonincludedReactions() {
		return nonincludedReactionList;
	}
	
	/**
	 * Returns the status of the altered flag: true if the network requires a
	 * new pressure-dependent calculation, false if not.
	 * @return The status of the altered flag
	 */
	public boolean getAltered() {
		return altered;
	}
	
	/**
	 * Returns the isomer that contains the same species as those in 
	 * speciesList.
	 * @param speciesList The species to check the isomers for
	 * @return The isomer that contains the same species as those in 
	 * speciesList
	 */
	public PDepIsomer getIsomer(LinkedList speciesList) {
		if (speciesList.size() == 1) {
			for (ListIterator<PDepIsomer> iter = uniIsomerList.listIterator(); iter.hasNext(); ) {
				PDepIsomer isomer = iter.next();
				if (isomer.getSpeciesList().equals(speciesList))
					return isomer;
			}
		}
		else if (speciesList.size() > 1) {
			for (ListIterator<PDepIsomer> iter = multiIsomerList.listIterator(); iter.hasNext(); ) {
				PDepIsomer isomer = iter.next();
				if (isomer.getSpeciesList().equals(speciesList))
					return isomer;
			}
		}
		return null;
	}

	public PDepIsomer getIsomer(Species species) {
		for (ListIterator<PDepIsomer> iter = uniIsomerList.listIterator(); iter.hasNext(); ) {
			PDepIsomer isomer = iter.next();
			if (isomer.getSpecies(0).equals(species))
				return isomer;
		}
		return null;
	}

	/**
	 * Checks to see if the network contains species as a unimolecular isomer.
	 * @param species The species to check for
	 * @return True if species is a unimolecular isomer in the network, false if
	 * not
	 */
	public boolean contains(Species species) {
		for (ListIterator<PDepIsomer> iter = uniIsomerList.listIterator(); iter.hasNext(); ) {
			PDepIsomer isomer = iter.next();
			if (isomer.getSpecies(0).equals(species))
				return true;
		}
		return false;
	}

	/**
	 * Checks to see if the network contains species as a unimolecular isomer.
	 * @param species The species to check for
	 * @return True if species is a unimolecular isomer in the network, false if
	 * not
	 */
	public boolean contains(Reaction reaction) {
		for (ListIterator<PDepReaction> iter = pathReactionList.listIterator(); iter.hasNext(); ) {
			PDepReaction rxn = iter.next();
			if (reaction.equals(rxn) || reaction.getReverseReaction().equals(rxn))
				return true;
		}
		return false;
	}
	
	//==========================================================================
	//
	//	Set accessor methods
	//
	
	
	/**
	 * Adds an isomer (unimolecular or multimolecular) to the appropriate 
	 * list in the network if it is not already present.
	 * @param isomer The isomer to add
	 */
	public void addIsomer(PDepIsomer isomer) {
		
		// Don't add if isomer is already in network
		if (uniIsomerList.contains(isomer) || multiIsomerList.contains(isomer))
			return;
		
		// Add isomer
		if (isomer.isUnimolecular())
			uniIsomerList.add(isomer);
		else if (isomer.isMultimolecular())
			multiIsomerList.add(isomer);
		
		// Mark network as changed so that updated rates can be determined
		altered = true;
	}
	
	/**
	 * Adds a path reaction to the network if it is not already present.
	 * @param newRxn The path reaction to add
	 */
	public void addReaction(PDepReaction newRxn) {
		
		// Check to ensure that reaction is not already present
		for (ListIterator<PDepReaction> iter = pathReactionList.listIterator(); iter.hasNext(); ) {
			PDepReaction rxn = iter.next();
			if (rxn.equals(newRxn))
				return;
		}
		
		// Add reaction
		pathReactionList.add(newRxn);
	}
	
	/**
	 * Updates the status of the altered flag: true if the network requires a
	 * new pressure-dependent calculation, false if not.
	 * @param alt The new status of the altered flag
	 */
	public void setAltered(boolean alt) {
		altered = alt;
	}

	public void makeIsomerIncluded(PDepIsomer isomer) {

		if (!isomer.isUnimolecular() || isomer.getIncluded())
			return;

		isomer.setIncluded(true);
		
		LinkedHashSet reactionSet = isomer.generatePaths(reactionSystem);

		for (Iterator iter = reactionSet.iterator(); iter.hasNext(); ) {
			Reaction rxn = (Reaction) iter.next();
			if (!contains(rxn))
				addReactionToNetworks(rxn);
		}
		
	}

	//==========================================================================
	//
	//	Other methods
	//

	/**
	 * Redistributes the net reactions based on the current core and edge
	 * reaction models. Especially useful when one or more species has been
	 * moved from the edge to the core since the last pressure-dependent
	 * calculation.
	 * @param cerm The current core/edge reaction model
	 */
	public void updateReactionLists(CoreEdgeReactionModel cerm) {
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
			if (forward == null || reverse == null)
				return;
			if (forward.isCoreReaction(cerm) || reverse.isCoreReaction(cerm))
				netReactionList.add(forward);
			else if (forward.getReactant().getIncluded() && forward.getProduct().getIncluded())
				netReactionList.add(forward);
			else
				nonincludedReactionList.add(forward);
		}

		System.out.println(netReactionList.size() + " included and " +
				nonincludedReactionList.size() + " nonincluded reactions in network.");

	}


	public double getLeakFlux(SystemSnapshot ss) {
		double rLeak = 0.0;
		for (ListIterator<PDepReaction> iter = nonincludedReactionList.listIterator(); iter.hasNext(); ) {
			PDepReaction rxn = iter.next();
			if (rxn.getReactant().getIncluded() && !rxn.getProduct().getIncluded())
				rLeak += rxn.calculateForwardFlux(ss);
			else if (!rxn.getReactant().getIncluded() && rxn.getProduct().getIncluded())
				rLeak += rxn.calculateReverseFlux(ss);
		}
		return rLeak;
	}

	public PDepIsomer getMaxLeakIsomer(SystemSnapshot ss) {
		PDepReaction maxReaction = null;
		double maxLeak = 0.0;
		for (ListIterator<PDepReaction> iter = nonincludedReactionList.listIterator(); iter.hasNext(); ) {
			PDepReaction rxn = iter.next();
			if (!rxn.getReactant().getIncluded() || !rxn.getProduct().getIncluded()) {
				if (Math.abs(rxn.calculateFlux(ss)) > maxLeak) {
					maxReaction = rxn;
					maxLeak = rxn.calculateFlux(ss);
				}
			}
		}
		if (maxReaction == null)
			return null;
		else if (!maxReaction.getReactant().getIncluded())
			return maxReaction.getReactant();
		else if (!maxReaction.getProduct().getIncluded())
			return maxReaction.getProduct();
		else
			return null;
	}

	//==========================================================================
	//
	//	Static methods (for access to PDepNetwork.networks)
	//
	
	/**
	 * Returns the linked list containing the currently-existing pressure-
	 * dependent networks
	 * @return The currently-existing pressure-dependent networks
	 */
	public static LinkedList<PDepNetwork> getNetworks() {
		return networks;
	}
	
	/**
	 * Used to add a reaction to the appropriate pressure-dependent network. If
	 * no such network exists, a new network is created. For isomerization
	 * reactions connecting two existing networks, the networks are merged. This
	 * function is to be called whenever a new reaction is added to the edge.
	 * @param reaction The reaction to add
	 * @return The network the reaction was added to
	 */
	public static PDepNetwork addReactionToNetworks(Reaction reaction) {
		
		// Expect that most reactions passed to this function will be already
		// present in a network

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

		PDepNetwork pdn = null;
		if (reaction.getProductNumber() == 1) {
			// Isomerization reactions should cause networks to be merged together
			// This means that each unimolecular isomer should only appear in one network

			// Get the appropriate pressure-dependent network(s)
			PDepNetwork reac_pdn = null;
			PDepNetwork prod_pdn = null;
			Species reactant = (Species) reaction.getReactantList().get(0);
			Species product = (Species) reaction.getProductList().get(0);
			for (ListIterator<PDepNetwork> iter = networks.listIterator(); iter.hasNext(); ) {
				PDepNetwork n = iter.next();
				if (n.contains(reactant))
					reac_pdn = n;
				if (n.contains(product))
					prod_pdn = n;
			}

			if (reac_pdn != null && prod_pdn != null && reac_pdn != prod_pdn) {
				// Two distinct networks found; must join them together
				pdn = reac_pdn;
				for (int i = 0; i < prod_pdn.getUniIsomers().size(); i++)
					pdn.addIsomer(prod_pdn.getUniIsomers().get(i));
				for (int i = 0; i < prod_pdn.getMultiIsomers().size(); i++)
					pdn.addIsomer(prod_pdn.getMultiIsomers().get(i));
				for (int i = 0; i < prod_pdn.getPathReactions().size(); i++)
					pdn.addReaction(prod_pdn.getPathReactions().get(i));
				// Also remove the second network from the list of networks
				networks.remove(prod_pdn);
			}
			else if (reac_pdn != null && prod_pdn != null && reac_pdn == prod_pdn) {
				// Both species already present as unimolecular isomers in the same network, so use that network
				pdn = reac_pdn;
			}
			else if (reac_pdn != null) {
				// Only reactant species found in a network, so use that network
				pdn = reac_pdn;
			}
			else if (prod_pdn != null) {
				// Only product species found in a network, so use that network
				pdn = reac_pdn;
			}
			else {
				// No networks found for either species; will create a new network
				pdn = null;
			}

		}
		else if (reaction.getProductNumber() > 1) {
			// Dissociation reactions are added to the network containing that unimolecular isomer
			// Since each unimolecular isomer should only appear in one network, there should only be one such addition
			// If no existing network is found, a new one may be created

			// Get the appropriate pressure-dependent network
			Species reactant = (Species) reaction.getReactantList().get(0);
			for (ListIterator<PDepNetwork> iter = networks.listIterator(); iter.hasNext(); ) {
				PDepNetwork n = iter.next();
				if (n.contains(reactant))
					pdn = n;
			}
		}

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
		PDepReaction rxn = new PDepReaction(reactantIsomer, productIsomer, reaction);
		pdn.addReaction(rxn);

		// Fill in partial network if necessary
		if (reactantIsomer.isCore((CoreEdgeReactionModel) reactionModel) && reactantIsomer.isUnimolecular())
			pdn.makeIsomerIncluded(reactantIsomer);
		if (productIsomer.isCore((CoreEdgeReactionModel) reactionModel) && productIsomer.isUnimolecular())
			pdn.makeIsomerIncluded(productIsomer);

		// Return the created network
		return pdn;

	}
	
	/**
	 * Useful for debugging, this function prints the isomers of each network
	 * to the console window.
	 */
	public static void printNetworks() {
		int index = 0;
		for (ListIterator<PDepNetwork> iter0 = networks.listIterator(); iter0.hasNext(); ) {
			PDepNetwork pdn = iter0.next();
			index++;
			System.out.print("Network #" + Integer.toString(index) + ": ");
			for (ListIterator<PDepIsomer> iter = pdn.getUniIsomers().listIterator(); iter.hasNext(); ) {
				PDepIsomer isomer = iter.next();
				System.out.print(isomer.toString());
				if (iter.hasNext())
					System.out.print(", ");
			}
			System.out.print("; ");
			for (ListIterator<PDepIsomer> iter = pdn.getMultiIsomers().listIterator(); iter.hasNext(); ) {
				PDepIsomer isomer = iter.next();
				System.out.print(isomer.toString());
				if (iter.hasNext())
					System.out.print(", ");
			}
			System.out.print("\n");
		}
	}
	
	/**
	 * Checks to see if there are any core reactions hidden amongst those
	 * net reactions which are found in the pressure-dependent networks.
	 * This is particularly useful in the initialization of the reaction model,
	 * in which the core must have at least one reaction in it before the
	 * dynamic simulator can be executed.
	 * @param cerm The current core/edge reaction model
	 * @return True if core reactions are found, false if not
	 */
	public static boolean hasCoreReactions(CoreEdgeReactionModel cerm) {
		return (getCoreReactions(cerm).size() > 0);
	}
	
	/**
	 * Counts the number of core reactions that are hidden amongst those
	 * net reactions which are found in the pressure-dependent networks.
	 * This is particularly useful in the initialization of the reaction model,
	 * in which the core must have at least one reaction in it before the
	 * dynamic simulator can be executed.
	 * @param cerm The current core/edge reaction model
	 * @return The number of core reactions found
	 */
	public static int getNumCoreReactions(CoreEdgeReactionModel cerm) {
		return getCoreReactions(cerm).size();
	}
	
	/**
	 * Returns the core reactions that are hidden amongst those
	 * net reactions which are found in the pressure-dependent networks.
	 * This is particularly useful in the initialization of the reaction model,
	 * in which the core must have at least one reaction in it before the
	 * dynamic simulator can be executed.
	 * @param cerm The current core/edge reaction model
	 * @return The number of core reactions found
	 */
	public static LinkedList<PDepReaction> getCoreReactions(CoreEdgeReactionModel cerm) {
		LinkedList<PDepReaction> coreReactions = new LinkedList<PDepReaction>();
		for (ListIterator<PDepNetwork> iter0 = networks.listIterator(); iter0.hasNext(); ) {
			PDepNetwork pdn = iter0.next();
			for (ListIterator<PDepReaction> iter = pdn.getNetReactions().listIterator(); iter.hasNext(); ) {
				PDepReaction rxn = iter.next();
				if (rxn.isCoreReaction(cerm) && !coreReactions.contains(rxn))
					coreReactions.add(rxn);
			}
		}
		return coreReactions;
	}
	
}

////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

package jing.rxn;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.StringTokenizer;

import jing.chem.ChemGraph;
import jing.chem.Species;
import jing.chem.SpeciesDictionary;
import jing.param.Temperature;
import jing.rxnSys.CoreEdgeReactionModel;
import jing.rxnSys.ReactionSystem;

/**
 * Represents a single isomer (unimolecular or multimolecular) of a pressure-
 * dependent network. In this context an isomer refers to one or more
 * species wherein the total number of each element in all species matches
 * the total number of that species in every other isomer of the network.
 * @author jwallen
 */
public class PDepIsomer {

	//==========================================================================
	//
	//	Data members
	//
	
	/**
	 * A linked list containing the species that make up the isomer. 
	 */
	private LinkedList<Species> speciesList;

	/**
	 * Set to true if the pathways from this isomer is included in the network.
	 */
	private boolean included;

	//==========================================================================
	//
	//	Constructors
	//
	
	/**
	 * Constructor for unimolecular isomers.
	 * @param species The species the isomer represents.
	 */
	public PDepIsomer(Species species) {
		speciesList = new LinkedList<Species>();
		speciesList.add(species);
		included = false;
	}
	
	/**
	 * Constructor for bimolecular isomers.
	 * @param species1 The first species the isomer represents.
	 * @param species2 The second species the isomer represents.
	 */
	public PDepIsomer(Species species1, Species species2) {
		speciesList = new LinkedList<Species>();
		speciesList.add(species1);
		speciesList.add(species2);
		included = true;
	}
	
	/**
	 * Constructor for generic isomers.
	 * @param speList A list of species the isomer represents.
	 */
	public PDepIsomer(LinkedList speList) {
		speciesList = new LinkedList<Species>();
		for (int i = 0; i < speList.size(); i++) {
			Object obj = speList.get(i);
			if (obj instanceof Species)
				speciesList.add((Species) obj);
			else if (obj instanceof ChemGraph)
				speciesList.add(SpeciesDictionary.getSpecies((ChemGraph) obj));
		}
		included = (speList.size() > 1);
	}
	
	/**
	 * Constructor for unimolecular isomers, read in from Restart files.
	 * @param species The species the isomer represents.
	 * @param included Whether the isomer is included or notIncluded 
	 */
	public PDepIsomer(Species species, boolean p_included) {
		speciesList = new LinkedList<Species>();
		speciesList.add(species);
		included = p_included;
	}
	
	/**
	 * Constructor for bimolecular isomers, read in from Restart file.
	 * @param species1 The first species the isomer represents.
	 * @param species2 The second species the isomer represents.
	 * @param included Whether the isomer is included or notIncluded
	 */
	public PDepIsomer(Species species1, Species species2, boolean p_included) {
		speciesList = new LinkedList<Species>();
		speciesList.add(species1);
		speciesList.add(species2);
		included = p_included;
	}
	
	//==========================================================================
	//
	//	Accessors
	//
	
	/**
	 * Returns the ith species of the isomer.
	 * @param i The index of the species to return
	 * @return The species corresponding to index i
	 */
	public Species getSpecies(int i) {
		return speciesList.get(i);
	}
	
	/**
	 * Returns the entire list of species in the isomer.
	 * @return The list of species in the isomer
	 */
	public LinkedList<Species> getSpeciesList() {
		return speciesList;
	}
	
	/**
	 * Returns an iterator to the first species in the isomer
	 * @return An iterator to the first species in the isomer
	 */
	public ListIterator<Species> getSpeciesListIterator() {
		return speciesList.listIterator();
	}
	
	/**
	 * Returns the number of species in the isomer.
	 * @return The number of species in the isomer.
	 */
	public int getNumSpecies() {
		return speciesList.size();
	}
	
	/**
	 * Returns a string containing the names and IDs of each species in the 
	 * isomer in the form "A(1) + B(2) + ...".
	 * @return The names and IDs of each species in the isomer
	 */
	public String getSpeciesNames() {
		String str = speciesList.get(0).getName() + "(" + 
				Integer.toString(speciesList.get(0).getID()) + ")";
		for (int i = 1; i < speciesList.size(); i++) {
			str += " + " + ((Species) speciesList.get(i)).getName() + "(" + 
					Integer.toString(speciesList.get(i).getID()) + ")";
		}
		return str;
	}
	
	/**
	 * Returns a string containing the names and IDs of each species in the 
	 * isomer in the form "A(1) + B(2) + ...".
	 * @return The names and IDs of each species in the isomer
	 */
	@Override
	public String toString() {
		return getSpeciesNames() + " (included =" + Boolean.toString(getIncluded()) + ")";
	}

	/**
	 * Returns the included status the isomer.
	 * @return The included status the isomer.
	 */
	public boolean getIncluded() {
		return included;
	}

	/**
	 * Sets the included status the isomer.
	 * @param incl The new included status the isomer.
	 */
	public void setIncluded(boolean incl) {
		included = incl;
	}

	//==========================================================================
	//
	//	Other methods
	//
	
	/**
	 * Returns true if the isomer is unimolecular and false otherwise.
	 * @return True if the isomer is unimolecular and false otherwise
	 */
	public boolean isUnimolecular() {
		if (speciesList == null) return false;
		return speciesList.size() == 1;
	}
	
	/**
	 * Returns true if the isomer is multimolecular and false otherwise.
	 * @return True if the isomer is multimolecular and false otherwise
	 */
	public boolean isMultimolecular() {
		if (speciesList == null) return false;
		return speciesList.size() > 1;
	}
	
	/**
	 * Calculates the total Gibbs free energy of the isomer at the specified
	 * temperature.
	 * @param temperature The temperature at which to do the calculation
	 * @return The total Gibbs free energy of the isomer
	 */
	public double calculateG(Temperature temperature) {
		double G = 0;
		for (int i = 0; i < speciesList.size(); i++) {
			G += speciesList.get(i).calculateG(temperature);
		}
		return G;
	}

	/**
	 * Calculates the Gibbs free energy of the ith species of the isomer at the 
	 * specified temperature.
	 * @param species The index of the species
	 * @param temperature The temperature at which to do the calculation
	 * @return The Gibbs free energy of the ith species
	 */
	public double calculateG(int species, Temperature temperature) {
		return speciesList.get(species).calculateG(temperature);
	}

	/**
	 * Calculates the total enthalpy of the isomer at the specified
	 * temperature.
	 * @param temperature The temperature at which to do the calculation
	 * @return The total enthalpy of the isomer
	 */
	public double calculateH(Temperature temperature) {
		double H = 0;
		for (int i = 0; i < speciesList.size(); i++) {
			H += speciesList.get(i).calculateH(temperature);
		}
		return H;
	}

	/**
	 * Calculates the enthalpy of the ith species of the isomer at the 
	 * specified temperature.
	 * @param species The index of the species
	 * @param temperature The temperature at which to do the calculation
	 * @return The enthalpy of the ith species
	 */
	public double calculateH(int species, Temperature temperature) {
		return speciesList.get(species).calculateH(temperature);
	}
	
	/**
	 * Returns the total molecular weight of the isomer.
	 * @return The total molecular weight of the isomer
	 */
	public double getMolecularWeight() {
		double molWt = 0;
		for (int i = 0; i < speciesList.size(); i++) {
			molWt += speciesList.get(i).getMolecularWeight();
		}
		return molWt;
	}

	/**
	 * Returns the molecular weight of the ith species of the isomer.
	 * @param species The index of the species
	 * @return The molecular weight of the ith species
	 */
	public double getMolecularWeight(int species) {
		return speciesList.get(species).getMolecularWeight();
	}
	
	/**
	 * Implements the equals method common to all Java objects. This method
	 * calls other versions of equals() based on the type of the object passed.
	 * @param object An object to compare the isomer to.
	 * @return True if the isomer and the object are equal, false if not
	 */
	@Override
	public boolean equals(Object object) {
		if (object instanceof PDepIsomer)
			return equals((PDepIsomer) object);
		else
			return false;
	}
	
	/**
	 * Returns true if the current object and isomer have all of the same
	 * species, no matter what order the species are stored in each object's
	 * species list.
	 * @param isomer The isomer to compare the current isomer to
	 * @return True if the two isomers have all of the same species, false if not
	 */
	public boolean equals(PDepIsomer isomer) {
		if (getNumSpecies() != isomer.getNumSpecies())
			return false;

		// We make a copy of the species list of the second isomer so that we
		// can ensure that A + B does not match B + B by removing each match
		// so that it cannot be reused
		LinkedList<Species> speciesList2 = new LinkedList<Species>(isomer.speciesList);

		for (int i = 0; i < getNumSpecies(); i++) {
			Species speciesI = speciesList.get(i);
			int index2 = -1;
			for (int j = 0; j < speciesList2.size(); j++) {
				Species speciesJ = speciesList2.get(j);
				if (speciesI.equals(speciesJ))
					index2 = j;
			}
			if (index2 == -1)
				return false;
			else
				speciesList2.remove(index2);
		}

		return true;
		
	}

	/**
	 * Returns all of the pressure-dependent pathways involving the current isomer.
	 * @param rxnSystem A reaction system to use to access the template and library reaction generators.
	 * @return A hash set of the pressure-dependent pathways involving the current isomer
	 */
	public LinkedHashSet generatePaths(ReactionSystem rxnSystem) {
		if (!isUnimolecular())
			return new LinkedHashSet();
		
		if(rxnSystem.getLibraryReactionGenerator()!= null){
			// First iterate through the Reaction Library and find all reactions which include the species being considered
		LinkedHashSet reactionSet = ((LibraryReactionGenerator) rxnSystem.getLibraryReactionGenerator()).generatePdepReactions(getSpecies(0));
		if(!reactionSet.isEmpty()){
			System.out.println("Reaction Set Found from Reaction Library "+reactionSet);
		}
				
		// Iterate through the reaction template
		reactionSet.addAll(((TemplateReactionGenerator) rxnSystem.getReactionGenerator()).generatePdepReactions(getSpecies(0)));
		
		// To remove the duplicates that are found in Reaction Library and Reaction Template
		// Preference given to Reaction Library over Template Reaction 
		
		LinkedHashSet newReactionSet_nodup = rxnSystem.getLibraryReactionGenerator().RemoveDuplicateReac(reactionSet);
		
		// shamel 6/22/2010 Suppressed output , line is only for debugging
		//System.out.println("Reaction Set For PdepIsomer "+newReactionSet_nodup);
		
		return newReactionSet_nodup;
	    }
		else{
			LinkedHashSet reactionSet = ((TemplateReactionGenerator) rxnSystem.getReactionGenerator()).generatePdepReactions(getSpecies(0));
			
			// shamel 6/22/2010 Suppressed output , line is only for debugging
			//System.out.println("Reaction Set For PdepIsomer "+reactionSet);
			
			return reactionSet;
		}
	}
	
	
   
   

	/**
	 * Checks to see if the isomer is in the model core (i.e. all of the 
	 * species in the isomer are core species).
	 * @param cerm The current core-edge reaction model
	 * @return True if all species are in the model core, false if not
	 */
	public boolean isCore(CoreEdgeReactionModel cerm) {
		for (int i = 0; i < speciesList.size(); i++) {
			if (!cerm.containsAsReactedSpecies(speciesList.get(i)))
				return false;
		}
		return true;
	}

}

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

import java.util.ListIterator;
import java.util.StringTokenizer;
import jing.chem.Species;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxnSys.CoreEdgeReactionModel;
import jing.rxnSys.NegativeConcentrationException;
import jing.rxnSys.ReactionModelGenerator;
import jing.rxnSys.SystemSnapshot;

/**
 * Represents either a pressure-dependent path reaction (a reaction connecting
 * two isomers directly) or a pressure-dependent net reaction (a reaction not
 * necessarily connecting two isomers directly). The former is provided to
 * a pressure-dependent kinetics estimator in order to generate the latter.
 * The latter are then treated as core/edge reactions based on the status
 * of the isomers they connect.
 * 
 * There are three types of path reactions: isomerizations
 * (A1 --> A2), associations (B + C --> A), and dissociations (A --> B + C).
 *
 * @author jwallen
 */
public class PDepReaction extends Reaction {

	/**
	 * An enumeraction of pressure-dependent reaction types:
	 * <ul>
	 * <li>NONE - The reaction type could not be assessed.
	 * <li>ISOMERIZATION - The reaction has the form A --> B.
	 * <li>ASSOCIATION - The reaction has the form B + C [+ ...] --> A.
	 * <li>DISSOCIATION - The reaction has the form A --> B + C [+ ...].
	 * <li>OTHER - The reaction has the form A + B [+ ...] --> P + Q [+ ...].
	 * </ul>
	 */
	public enum Type { NONE, ISOMERIZATION, ASSOCIATION, DISSOCIATION, OTHER };
	
	/**
	 * A reference to the reactant well.
	 */
	private PDepIsomer reactant;
	
	/**
	 * A reference to the product well.
	 */
	private PDepIsomer product;
	
	/**
	 * The pressure-dependent rate coefficient.
	 */
	private PDepRateConstant pDepRate;
	
	/**
	 * The reverse PDepReaction. Is used instead of Reaction.reverse when
	 * dealing with PDepReaction objects.
	 */
	private PDepReaction pDepReverse;
			
	//==========================================================================
	//
	//	Constructors
	//
	
	/**
	 * Create a pressure-dependent path reaction connecting isomers reac and
	 * prod and having high-pressure Arrhenius kinetics as found in rxn.
	 * @param reac The reactant PDepIsomer
	 * @param prod The product PDepIsomer
	 * @param rxn A Reaction containing the appropriate high-pressure kinetics
	 */
	public PDepReaction(PDepIsomer reac, PDepIsomer prod, Reaction rxn) {
		super();
		structure = rxn.structure;
		kinetics = rxn.getKinetics();
		reverseReaction = rxn.reverseReaction;
		if (structure == null)
			structure = new Structure(reac.getSpeciesList(), prod.getSpeciesList(), 1);
		setReactant(reac);
		setProduct(prod);
		pDepRate = null;
		kineticsFromPrimaryKineticLibrary = rxn.isFromPrimaryKineticLibrary();
	}
	
	/**
	 * Create a pressure-dependent path reaction connecting isomers reac and
	 * prod and having high-pressure Arrhenius kinetics as found in rxn.
	 * @param reac The reactant PDepIsomer
	 * @param prod The product PDepIsomer
	 * @param kin The high-pressure kinetics for the forward reaction
	 */
	public PDepReaction(PDepIsomer reac, PDepIsomer prod, Kinetics[] kin) {
		super();
		structure = null;
		kinetics = kin;
		reverseReaction = null;
		if (structure == null)
			structure = new Structure(reac.getSpeciesList(), prod.getSpeciesList(), 1);
		setReactant(reac);
		setProduct(prod);
		pDepRate = null;
		kineticsFromPrimaryKineticLibrary = kin[0].getFromPrimaryKineticLibrary();
	}
	
	/**
	 * Create a pressure-dependent net reaction connecting isomers reac and
	 * prod and having k(T, P) kinetics as approximated by the Chebyshev 
	 * polynomials of cheb.
	 * @param reac The reactant PDepIsomer
	 * @param prod The product PDepIsomer
	 * @param k The k(T, P) rate constant for the reaction
	 */
	public PDepReaction(PDepIsomer reac, PDepIsomer prod, PDepRateConstant k) {
		super();
		structure = new Structure(reac.getSpeciesList(), prod.getSpeciesList(), 1);
		setReactant(reac);
		setProduct(prod);
		setHighPKinetics(null);
		pDepRate = k;
	}
	
	//==========================================================================
	//
	//	Accessors
	//
	
	/**
	 * Returns the reaction type as one of the options given in the Type
	 * enumeration.
	 * @return The reaction type
	 */
	public Type getType() {
		if (reactant == null || product == null)
			return Type.NONE;
		else if (reactant.isUnimolecular() && product.isUnimolecular())
			return Type.ISOMERIZATION;
		else if (reactant.isMultimolecular() && product.isUnimolecular())
			return Type.ASSOCIATION;
		else if (reactant.isUnimolecular() && product.isMultimolecular())
			return Type.DISSOCIATION;
		else if (reactant.isMultimolecular() && product.isMultimolecular())
			return Type.OTHER;
		else
			return Type.NONE;
	}
	
	/**
	 * Returns the reactant isomer
	 * @return The reactant isomer
	 */
	public PDepIsomer getReactant() {
		return reactant;
	}
	
	/**
	 * Returns the product isomer
	 * @return The reactant isomer
	 */
	public PDepIsomer getProduct() {
		return product;
	}
	
	/**
	 * Sets the reactant isomer to r.
	 * @param r The new reactant isomer
	 */
	public void setReactant(PDepIsomer r) {
		reactant = r;
	}
	
	/**
	 * Sets the product isomer to r.
	 * @param r The new product isomer
	 */
	public void setProduct(PDepIsomer p) {
		product = p;
	}
	
	/**
	 * An alias for Reaction.getKinetics() that emphasizes that the returned
	 * kinetics are only valid in the high-pressure limit.
	 * @return The high-pressure Arrhenius kinetics for the reaction
	 */
	public Kinetics[] getHighPKinetics() {
		return getKinetics();
	}
	
	/**
	 * An alias for Reaction.setKinetics() that emphasizes that the 
	 * kinetics are only valid in the high-pressure limit.
	 * @param kin The new high-pressure Arrhenius kinetics for the reaction
	 */
	public void setHighPKinetics(Kinetics kin) {
		setKinetics(kin,-1);
	}
	
	/** 
	 * Returns the Chebyshev polynomial fit to k(T, P) for this net reaction.
	 * @return The Chebyshev polynomial fit to k(T, P) for this net reaction
	 */
	public PDepRateConstant setPDepRateConstant() {
		return pDepRate;
	}
	
	/** 
	 * Sets the Chebyshev polynomial fit to k(T, P) for this net reaction.
	 * @param cheb The new Chebyshev polynomial fit to k(T, P) for this net reaction
	 */
	public void setPDepRateConstant(PDepRateConstant r) {
		pDepRate = r;
	}
	
	/**
	 * Gets the current reverse reaction, using PDepReaction.pDepReverse rather
	 * than Reaction.reverse.
	 * @return The current pressure-dependent reverse reaction
	 */
	@Override
	public PDepReaction getReverseReaction() {
		return pDepReverse;
	}
	
	/**
	 * If the provided reaction is pressure-dependent, sets the reverse reaction
	 * to that reaction.
	 * @param rxn The new pressure-dependent reverse reaction
	 */
	@Override
	public void setReverseReaction(Reaction rxn) {
		if (rxn instanceof PDepReaction)
			pDepReverse = (PDepReaction) rxn;
	}
	
	//==========================================================================
	//
	//	Other methods
	//
	
	/**
	 * Returns true if this reaction represents an edge reaction (that is,
	 * all species in the reactant isomer are in the model core and all species 
	 * in the product isomer are in the current model edge), and false
	 * otherwise.
	 * @param cerm The current core/edge reaction model
	 * @return True if an edge reaction, false if not
	 */
	public boolean isEdgeReaction(CoreEdgeReactionModel cerm) {
		// All reactant species must be in the core
		for (int i = 0; i < reactant.getNumSpecies(); i++) 
			if (cerm.containsAsReactedSpecies(reactant.getSpecies(i)) == false)
				return false;
		// At least one product species must be in the edge
		for (int i = 0; i < product.getNumSpecies(); i++) 
			if (cerm.containsAsUnreactedSpecies(product.getSpecies(i)) == true)
				return true;
		// If here, then reaction is not on the edge
		return false;
	}
	
	/**
	 * Returns true if this reaction represents a core reaction (that is,
	 * all species in the reactant and product isomers are in the current
	 * model core), and false otherwise.
	 * @param cerm The current core/edge reaction model
	 * @return True if a core reaction, false if not
	 */
	public boolean isCoreReaction(CoreEdgeReactionModel cerm) {
		// All reactant species must be in the core
		for (int i = 0; i < reactant.getNumSpecies(); i++) 
			if (cerm.containsAsReactedSpecies(reactant.getSpecies(i)) == false)
				return false;
		// All product species must be in the edge
		for (int i = 0; i < product.getNumSpecies(); i++) 
			if (cerm.containsAsReactedSpecies(product.getSpecies(i)) == false)
				return false;
		// If here, then reaction is in the core
		return true;
	}
	
	/**
	 * Returns true if the reaction is a net reaction, signified by it having
	 * a non-null Chebyshev polynomial fit for k(T, P).
	 * @return True if the reaction is a net reaction, false otherwise
	 */
	public boolean isNetReaction() {
		return (pDepRate != null);
	}
	
	/**
	 * Returns the reaction as an ASCII string.
	 * @return A string representing the reaction equation in ASCII test.
	 */
	@Override
	public String toString() {
		if (reactant == null || product == null)
			return "";
		else {
			/*
			 * Distinguish between reversible and irreversible pdep reactions
			 * 	In particular, this is important when writing/reading
			 * 	Restart/pdepnetworks.txt file
			 */
			if (this.getReverseReaction() == null)
				return (reactant.toString() + " --> " + product.toString());
			else
				return (reactant.toString() + " <=> " + product.toString());
		}
	}
	
	/**
	 * Calculates the rate coefficient for the forward reaction at the 
	 * specified temperature and pressure. Uses the Chebyshev polynomial fit
	 * if present, then reverts to the high-pressure Arrhenius fit if present.
	 * @param temperature The temperature to determine the rate coefficient at
	 * @param pressure The pressure to determine the rate coefficient at
	 * @return The calculated rate coefficient for the forward reaction
	 */
	public double calculateRate(Temperature temperature, Pressure pressure) {
		double k = 0.0;
		try {
			if (pDepRate != null)
				k = pDepRate.calculateRate(temperature, pressure);
			else if (kinetics != null) {
				for (int numKinetics=0; numKinetics<kinetics.length; ++numKinetics) {
					k += kinetics[numKinetics].calculateRate(temperature);
				}
			}
		}
		catch (Exception e) {
			System.err.println(e.getMessage());
			System.err.println("Reaction: "+this.toChemkinString(temperature, pressure));
			System.exit(0);
		}

		return k;
	}
	
	/**
	 * Calculates the flux of this reaction given the provided system snapshot.
	 * The system snapshot contains the temperature, pressure, and 
	 * concentrations of each core species.
	 * @param ss The system snapshot at which to determine the reaction flux
	 * @return The determined reaction flux
	 */
	public double calculateFlux(SystemSnapshot ss) {
		return calculateForwardFlux(ss) - calculateReverseFlux(ss);
	}

	/**
	 * Calculates the forward flux of this reaction given the provided system snapshot.
	 * The system snapshot contains the temperature, pressure, and
	 * concentrations of each core species.
	 * @param ss The system snapshot at which to determine the reaction flux
	 * @return The determined reaction flux
	 */
	public double calculateForwardFlux(SystemSnapshot ss) {
		Temperature T = ss.getTemperature();
		Pressure P = ss.getPressure();
		double forwardFlux = calculateRate(T, P);
		for (ListIterator<Species> iter = reactant.getSpeciesListIterator(); iter.hasNext(); ) {
			Species spe = iter.next();
			double conc = 0.0;
			if (ss.getSpeciesStatus(spe) != null)
				conc = ss.getSpeciesStatus(spe).getConcentration();
			if (conc < 0) {
				double aTol = ReactionModelGenerator.getAtol();
				//if (Math.abs(conc) < aTol) conc = 0;
				//else throw new NegativeConcentrationException(spe.getName() + ": " + String.valueOf(conc));
				if (conc < -100.0 * aTol)
					throw new NegativeConcentrationException("Species " + spe.getName() + " has negative concentration: " + String.valueOf(conc));
			}
			forwardFlux *= conc;
		}
		return forwardFlux;
	}

	/**
	 * Calculates the flux of this reaction given the provided system snapshot.
	 * The system snapshot contains the temperature, pressure, and
	 * concentrations of each core species.
	 * @param ss The system snapshot at which to determine the reaction flux
	 * @return The determined reaction flux
	 */
	public double calculateReverseFlux(SystemSnapshot ss) {
		if (pDepReverse != null)
			return pDepReverse.calculateForwardFlux(ss);
		else
			return 0.0;
	}

	/**
	 * Returns true if either the forward or reverse reaction matches the
	 * provided reaction.
	 * @param rxn The reaction to compare the current reaction to
	 * @return True if the reactions are the same, false if not
	 */
	public boolean equals(PDepReaction rxn) {
		if (rxn.reactant.equals(reactant) && rxn.product.equals(product))
			return true;
		else if (rxn.reactant.equals(product) && rxn.product.equals(reactant))
			return true;
		else
			return false;
	}

	/**
	 * Returns true if either the forward or reverse reaction matches the
	 * provided reaction.
	 * @param rxn The reaction to compare the current reaction to
	 * @return True if the reactions are the same, false if not
	 */
	public boolean equals(Reaction rxn) {
		return super.equals(rxn);
	}
	
	/**
	 * Generates the reverse PDepReaction, overriding Reaction.generateReverseReaction().
	 */
	@Override
	public void generateReverseReaction() {
        if (pDepRate != null) {
			PDepReaction r = new PDepReaction(product, reactant, pDepRate);
			setReverseReaction(r);
			r.setReverseReaction(this);
		}
		else {
			super.generateReverseReaction();
			PDepReaction r = new PDepReaction(product, reactant, super.getReverseReaction()); 
			setReverseReaction(r); 
			r.setReverseReaction(this);
		}
    }
	
	/**
	 * Returns true if the PDepReaction has a previously-defined reverse reaction.
	 * @return
	 */
	@Override
	public boolean hasReverseReaction() {
		return (pDepReverse != null);
	}
	
	/**
	 * A holdover from the old PDepNetReaction class, used by
	 * PDepReaction.toChemkinString().
	 * @param p_string A Chemkin string from a reaction structure
	 * @return The parsed version of the string
	 */
	public String formPDepSign(String p_string) {
        StringTokenizer st = new StringTokenizer(p_string, "=");
        String s1 = st.nextToken();
        s1 += "(+m)=";
        String s2 = st.nextToken();
        s2 += "(+m)";
        return (s1+s2);
    }
	
	/**
	 * A holdover from the old PDepNetReaction class, used to generate a 
	 * Chemkin string. For a path reaction Reaction.toChemkinString() is used,
	 * while for a net reaction a different string is constructed.
	 * @param t A temperature at which to evaluate the Chemkin string at
	 * @return The resulting Chemkin string
	 */
	@Override
	public String toChemkinString(Temperature t) {
        if (pDepRate != null) {

			String result = getStructure().toChemkinString(true).toString();
			if (PDepRateConstant.getMode() == PDepRateConstant.Mode.CHEBYSHEV)
				result = formPDepSign(result);
			result += '\t' + "1.0E0 0.0 0.0" ;
			result += "\t!" + getComments().toString()  + '\n';
			if (PDepRateConstant.getMode() == PDepRateConstant.Mode.CHEBYSHEV)
				result += pDepRate.getChebyshev().toChemkinString() + '\n';
			else if (PDepRateConstant.getMode() == PDepRateConstant.Mode.PDEPARRHENIUS)
				result += pDepRate.getPDepArrheniusKinetics().toChemkinString();

			return result;
		}
		else if (kinetics != null)
			//return super.toChemkinString(t);
			// MRH 18Jan2010:
			//	Changed from toChemkinString(t) to toRestartString(t) to avoid bug in 
			//		reporting chem.inp file (issue of reporting A(single event) vs.
			//		A(single event) * (# events)
			return super.toRestartString(t,true);
		else
			return "";
    
    }
	
	public String toRestartString(Temperature t) {
        if (pDepRate != null) {

			String result = getStructure().toRestartString(true).toString();
			if (PDepRateConstant.getMode() == PDepRateConstant.Mode.CHEBYSHEV)
				result = formPDepSign(result);
			result += '\t' + "1.0E0 0.0 0.0" ;
			result += "\t!" + getComments().toString() +
				"\tdeltaHrxn(T=298K) = " + 
				calculateHrxn(new Temperature(298.0,"K")) + " kcal/mol\n";
			if (PDepRateConstant.getMode() == PDepRateConstant.Mode.CHEBYSHEV)
				result += pDepRate.getChebyshev().toChemkinString() + '\n';
			else if (PDepRateConstant.getMode() == PDepRateConstant.Mode.PDEPARRHENIUS)
				result += pDepRate.getPDepArrheniusKinetics().toChemkinString();

			return result;
		}
		else if (kinetics != null)
			//return super.toChemkinString(t);
			// MRH 18Jan2010:
			//	Changed from toChemkinString(t) to toRestartString(t) to avoid bug in 
			//		reporting chem.inp file (issue of reporting A(single event) vs.
			//		A(single event) * (# events)
			return super.toRestartString(t,true);
		else
			return "";
    
    }
	
	// 6Jul2009-MRH:
	//	This toChemkinString function is identical to the toChemkinString(Temperature)
	//		function, with the exception of the if(PDep.getMode() == Mode.RATE) statement
	//		The passed pressure p is only read if the mode == RATE
	@Override
	public String toChemkinString(Temperature t, Pressure p) {
        if (pDepRate != null) {

			String result = getStructure().toChemkinString(true).toString();
			if (PDepRateConstant.getMode() == PDepRateConstant.Mode.CHEBYSHEV)
				result = formPDepSign(result);
			if (PDepRateConstant.getMode() == PDepRateConstant.Mode.RATE) {
				result = formPDepSign(result) + "\t" + calculateRate(t,p) + " 0.0 0.0";
				result += "\t!" + getComments().toString()  + '\n';
				return result;
			}
			result += '\t' + "1.0E0 0.0 0.0";
			result += "\t!" + getComments().toString()  + '\n';

			if (PDepRateConstant.getMode() == PDepRateConstant.Mode.CHEBYSHEV)
				result += pDepRate.getChebyshev().toChemkinString() + '\n';
			else if (PDepRateConstant.getMode() == PDepRateConstant.Mode.PDEPARRHENIUS)
				result += pDepRate.getPDepArrheniusKinetics().toChemkinString();
			
			return result;
		}
		else if (kinetics != null)
			return super.toChemkinString(t);
		else
			return "";
    
    }
	
	public PDepRateConstant getPDepRate() {
		return pDepRate;
	}
}

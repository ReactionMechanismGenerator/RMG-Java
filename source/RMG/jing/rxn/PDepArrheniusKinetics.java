/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jing.rxn;

import jing.param.Pressure;
import jing.param.Temperature;

/**
 * A model of pressure-dependent reaction kinetics k(T, P) of the form
 * 
 * k(T, P) = A(P) * T ^ n(P) * exp( -Ea(P) / R * T )
 *
 * The values of A(P), n(P), and Ea(P) are stored for multiple pressures and
 * interpolation on a log P scale is used.
 *
 * @author jwallen
 */
public class PDepArrheniusKinetics implements PDepKinetics {

	/**
	 * The list of pressures at which we have Arrhenius parameters.
	 */
	public static Pressure[] pressures;

	/**
	 * The list of Arrhenius kinetics fitted at each pressure.
	 */
	private ArrheniusKinetics[] kinetics;
	
	protected static int numPressures = 0;

	public PDepArrheniusKinetics(int numP) {
		pressures = new Pressure[numP];
		kinetics = new ArrheniusKinetics[numP];
		setNumPressures(numP);
	}
	
	public void setKinetics(int index, Pressure P, ArrheniusKinetics kin) {
		if (index < 0 || index >= pressures.length)
			return;
		pressures[index] = P;
		kinetics[index] = kin;
	}

	/**
	 * Calculate the rate cofficient at the specified conditions.
	 * @param T The temperature of interest
	 * @param P The pressure of interest
	 * @return The rate coefficient evaluated at T and P
	 */
	public double calculateRate(Temperature T, Pressure P) {
		int index1 = -1; int index2 = -1;

		for (int i = 0; i < pressures.length - 1; i++) {
			if (pressures[i].getBar() <= P.getBar() && P.getBar() <= pressures[i+1].getBar()) {
				index1 = i; index2 = i + 1;
			}
		}

		if (index1 < 0 || index2 < 0)
			return 0.0;

		double logk1 = Math.log10(kinetics[index1].calculateRate(T));
		double logk2 = Math.log10(kinetics[index2].calculateRate(T));
		double logP0 = Math.log10(P.getBar());
		double logP1 = Math.log10(pressures[index1].getBar());
		double logP2 = Math.log10(pressures[index2].getBar());

		double logk0 = logk1 + (logk2 - logk1) / (logP2 - logP1) * (logP0 - logP1);

		return Math.pow(10, logk0);
	}

    public String toChemkinString() {
        String result = "";
		for (int i = 0; i < pressures.length; i++) {
			result += "PLOG / " + Double.toString(pressures[i].getAtm());
//			result += " / " + Double.toString(kinetics[i].getAValue());
//			result += " / " + Double.toString(kinetics[i].getNValue());
//			result += " / " + Double.toString(kinetics[i].getEValue());
			// 6Jul2009-MRH:
			//	PLOG format does not need "/" between parameters
			result += " " + Double.toString(kinetics[i].getAValue());
			result += " " + Double.toString(kinetics[i].getNValue());
			double Ea_in_kcalmol = kinetics[i].getEValue();
			double Ea = 0.0;
			if (ArrheniusKinetics.getEaUnits().equals("kcal/mol"))		Ea = Ea_in_kcalmol;
			else if (ArrheniusKinetics.getEaUnits().equals("cal/mol"))	Ea = Ea_in_kcalmol * 1000.0;
			else if (ArrheniusKinetics.getEaUnits().equals("kJ/mol"))	Ea = Ea_in_kcalmol * 4.184;
			else if (ArrheniusKinetics.getEaUnits().equals("J/mol"))	Ea = Ea_in_kcalmol * 4184.0;
			else if (ArrheniusKinetics.getEaUnits().equals("Kelvins"))	Ea = Ea_in_kcalmol / 1.987e-3;
			result += " " + Double.toString(Ea);
			result += " /\n";
		}
		return result;
    }
    
    public static void setNumPressures(int numP) {
    	if (numP > getNumPressures()) numPressures = numP;
    }
    
    public static int getNumPressures() {
    	return numPressures;
    }
    
    public ArrheniusKinetics getKinetics(int i) {
    	return kinetics[i];
    }
    
    public static void setPressures(Pressure[] ListOfPressures) {
    	pressures = ListOfPressures;
    }
    
    public void setRateCoefficients(ArrheniusKinetics[] ListOfKinetics) {
    	kinetics = ListOfKinetics;
    }
    
    public Pressure getPressures(int i) {
    	return pressures[i];
    }

}

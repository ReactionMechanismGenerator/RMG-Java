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

import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxnSys.Logger;

/**
 * A model of pressure-dependent reaction kinetics k(T, P) of the form k(T, P) = A(P) * T ^ n(P) * exp( -Ea(P) / R * T )
 * The values of A(P), n(P), and Ea(P) are stored for multiple pressures and interpolation on a log P scale is used.
 * 
 * @author jwallen
 */
public class PDepArrheniusKinetics implements PDepKinetics {
    /**
     * The list of pressures at which we have Arrhenius parameters.
     */
    public Pressure[] pressures;
    /**
     * The list of Arrhenius kinetics fitted at each pressure.
     */
    private ArrheniusKinetics[] kinetics;
    protected int numPressures = 0;

    public PDepArrheniusKinetics(int numP) {
        pressures = new Pressure[numP];
        kinetics = new ArrheniusKinetics[numP];
        setNumPressures(numP);
    }

    public void setKinetics(int index, Pressure P, ArrheniusKinetics kin) {
        if (index < 0 || index >= pressures.length)
            throw new RuntimeException(
                    String.format(
                            "Cannot set kinetics with index %s because array is only of size %s",
                            index, pressures.length));
        pressures[index] = P;
        kinetics[index] = kin;
    }

    /**
     * Calculate the rate cofficient at the specified conditions.
     * 
     * @param T
     *            The temperature of interest
     * @param P
     *            The pressure of interest
     * @return The rate coefficient evaluated at T and P
     */
    public double calculateRate(Temperature T, Pressure P) {
        int index1 = -1;
        int index2 = -1;
        for (int i = 0; i < pressures.length - 1; i++) {
            if (pressures[i].getBar() <= P.getBar()
                    && P.getBar() <= pressures[i + 1].getBar()) {
                index1 = i;
                index2 = i + 1;
                break;
            }
        }
        /*
         * Chemkin 4 theory manual specifies: "If the rate of the reaction is desired for a pressure lower than any of
         * those provided, the rate parameters provided for the lowest pressure are used. Likewise, if rate of the
         * reaction is desired for a pressure higher than any of those provided, the rate parameters provided for the
         * highest pressure are used." We take the same approach here, but warn the user (so they can fix their input
         * file).
         */
        if (P.getPa() < pressures[0].getPa()) {
            Logger.warning(String
                    .format("Tried to evaluate rate coefficient at P=%.3g Atm, which is below minimum for this PLOG rate.",
                            P.getAtm()));
            Logger.warning(String.format(
                    "Using rate for minimum %s Atm instead",
                    pressures[0].getAtm()));
            return kinetics[0].calculateRate(T);
        }
        if (P.getPa() > pressures[pressures.length - 1].getPa()) {
            Logger.warning(String
                    .format("Tried to evaluate rate coefficient at P=%.3g Atm, which is above maximum for this PLOG rate.",
                            P.getAtm()));
            Logger.warning(String.format(
                    "Using rate for maximum %s Atm instead",
                    pressures[pressures.length - 1].getAtm()));
            return kinetics[pressures.length - 1].calculateRate(T);
        }
        double logk1 = Math.log10(kinetics[index1].calculateRate(T));
        double logk2 = Math.log10(kinetics[index2].calculateRate(T));
        double logP0 = Math.log10(P.getBar());
        double logP1 = Math.log10(pressures[index1].getBar());
        double logP2 = Math.log10(pressures[index2].getBar());
        // We can't take logarithms of k=0 and get meaningful interpolation, so we have to do something weird.
        // The approach used here is arbitrary, but at least it gives a continuous k(P) function.
        //
        // If interpolating between k1=0 and k2=0, return k=0
        if (logk1 == Double.NEGATIVE_INFINITY
                && logk2 == Double.NEGATIVE_INFINITY)
            return 0.0;
        // if interpolating between k1=0 and k2>0, set k1 to something small but nonzero.
        else if (logk1 == Double.NEGATIVE_INFINITY)
            logk1 = Math.min(0, logk2 - 1); // k1 is a small 1 cm3/mol/sec, or k2/10 if that's even smaller.
        // if interpolating between k1>0 and k2=0, set k2 to something small but nonzero.
        else if (logk2 == Double.NEGATIVE_INFINITY)
            logk2 = Math.min(0, logk1 - 1); // k2 is a small 1 cm3/mol/sec, or k1/10 if that's even smaller.
        double logk0 = logk1 + (logk2 - logk1) / (logP2 - logP1)
                * (logP0 - logP1);
        return Math.pow(10, logk0);
    }

    public String toChemkinString(int numReac) {
        String result = "";
        for (int i = 0; i < pressures.length; i++) {
            double Ea_in_kcalmol = kinetics[i].getEValue();
            // ***note: PLOG uses the same units for Ea and A as Arrhenius expressions; this has been a persistent
// source of confusion; see
// https://github.com/GreenGroup/RMG-Java/commit/2947e7b8d5b1e3e19543f2489990fa42e43ecad2#commitcomment-844009
            double A = kinetics[i].getAValue();
            double A_multiplier = 1.0;// for "moles", multiplier is 1
            if (ArrheniusKinetics.getAUnits().equals("molecules")) {
                A_multiplier = 1 / 6.022e23;
            }
            // convert the units (cf. similar code in ChemParser)
            if (numReac == 1) {
                // do nothing, no conversion needed
            } else if (numReac == 2) {
                A = A * A_multiplier;
            } else if (numReac == 3) {
                A = A * A_multiplier * A_multiplier;
            } else {
                Logger.error("Unsupported number of reactants:" + numReac);
                System.exit(0);
            }
            double Ea = 0.0;
            if (ArrheniusKinetics.getEaUnits().equals("kcal/mol"))
                Ea = Ea_in_kcalmol;
            else if (ArrheniusKinetics.getEaUnits().equals("cal/mol"))
                Ea = Ea_in_kcalmol * 1000.0;
            else if (ArrheniusKinetics.getEaUnits().equals("kJ/mol"))
                Ea = Ea_in_kcalmol * 4.184;
            else if (ArrheniusKinetics.getEaUnits().equals("J/mol"))
                Ea = Ea_in_kcalmol * 4184.0;
            else if (ArrheniusKinetics.getEaUnits().equals("Kelvins"))
                Ea = Ea_in_kcalmol / 1.987e-3;
            result += String.format("PLOG / %10s    %10.2e  %10s  %10s /\n",
                    pressures[i].getAtm(), A, kinetics[i].getNValue(), Ea);
        }
        return result;
    }

    public void setNumPressures(int numP) {
        numPressures = numP;
    }

    public int getNumPressures() {
        return numPressures;
    }

    public ArrheniusKinetics getKinetics(int i) {
        return kinetics[i];
    }

    public void setPressures(Pressure[] ListOfPressures) {
        pressures = ListOfPressures;
    }

    public void setRateCoefficients(ArrheniusKinetics[] ListOfKinetics) {
        kinetics = ListOfKinetics;
    }

    public Pressure getPressure(int i) {
        return pressures[i];
    }
}

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

import java.util.*;
import jing.param.*;
import jing.mathTool.UncertainDouble;
import jing.rxnSys.Logger;

// ## package jing::rxn
// ----------------------------------------------------------------------------
// jing\rxn\ArrheniusEPKinetics.java
// ----------------------------------------------------------------------------
// ## class ArrheniusEPKinetics
public class ArrheniusEPKinetics extends ArrheniusKinetics {
    protected UncertainDouble alpha = new UncertainDouble(0, 0, "Adder"); // ## attribute alpha

    // Constructors
    // ## operation
// ArrheniusEPKinetics(UncertainDouble,UncertainDouble,UncertainDouble,UncertainDouble,String,int,String,String)
    public ArrheniusEPKinetics(UncertainDouble p_A, UncertainDouble p_n,
            UncertainDouble p_alpha, UncertainDouble p_E, String p_TRange,
            int p_rank, String p_source, String p_comment) {
        super(p_A, p_n, p_E, p_TRange, p_rank, p_source, p_comment);
        alpha = p_alpha;
    }

    public ArrheniusEPKinetics() {
    }

    // ## operation toChemkinString(double, Temperature, boolean)
    public String toChemkinString(double p_Hrxn, Temperature p_temperature,
            boolean includeComments) {
        ArrheniusKinetics arrhenius = fixBarrier(p_Hrxn);
        return arrhenius
                .toChemkinString(p_Hrxn, p_temperature, includeComments);
    }

    // ## operation calculateRate(Temperature,double)
    public double calculateRate(Temperature p_temperature, double p_Hrxn) {
        ArrheniusKinetics arrhenius = fixBarrier(p_Hrxn);
        return arrhenius.calculateRate(p_temperature);
    }

    // ## operation getAlphaUncertainty()
    public double getAlphaUncertainty() {
        // #[ operation getAlphaUncertainty()
        return alpha.getUncertainty();
        // #]
    }

    // ## operation getAlphaValue()
    public double getAlphaValue() {
        // #[ operation getAlphaValue()
        return alpha.getValue();
        // #]
    }

    // ## operation multiply(double)
    public Kinetics multiply(double p_multiple) {
        // #[ operation multiply(double)
        UncertainDouble newA = getA().multiply(p_multiple);
        Kinetics newK = new ArrheniusEPKinetics(newA, getN(), getAlpha(),
                getE(), getTRange(), getRank(), getSource(), getComment());
        return newK;
        // #]
    }

    // ## operation toString()
    public String toString() {
        // #[ operation toString()
        String string = "A = " + A.toString() + "; n = " + n.toString()
                + "; E = " + E.toString() + ";" + " alpha = "
                + alpha.toString() + '\n';
        string = string + "Source: " + source + '\n';
        string = string + "Comments: " + comment + '\n';
        return string;
        // #]
    }

    public UncertainDouble getAlpha() {
        return alpha;
    }

    public double getEValue() {
        throw new InvalidKineticsTypeException(
                "Convert ArrheniusEPKinetics to ArrheniusKinetics with fixBarrier(Hrxn).");
    }

    public double getEaValue(double p_Hrxn) {
        ArrheniusKinetics arrhenius = fixBarrier(p_Hrxn);
        return arrhenius.getEValue();
    }

    public ArrheniusKinetics fixBarrier(double p_Hrxn) {
        // create a new ArrheniusKinetics object with a corrected barrier and return it.
        double al = alpha.getValue();
        double Eo = E.getValue();
        double Ea = Eo + al * p_Hrxn;
        UncertainDouble newEa = getE();
        String newComment = getComment();
        String warning = "";
        if (al != 0.0) {
            warning = String
                    .format("Ea computed using Evans-Polanyi dHrxn(298K)=%.1f kcal/mol and alpha=%.2f.",
                            p_Hrxn, al);
            newComment += " " + warning;
            Logger.info(warning);
            newEa = newEa.plus((UncertainDouble) getAlpha().multiply(p_Hrxn));
        }
        if (Eo >= 0 && Ea < 0) {
            // Negative barrier estimated by Evans-Polanyi, despite non-negative intrinsic barrier.
            warning = String.format("Ea raised from %.1f kcal/mol to 0.0.", Ea);
            newComment += " Warning: " + warning;
            Logger.info(warning);
            newEa = newEa.plus((-Ea));
            Ea = 0.0;
        }
        if (Eo < 0 && Ea < Eo) {
            // Negative barrier estimated by Evans-Polanyi is even more negative than negative intrinsic barrier.
            warning = String.format("Ea raised from %.1f to Eo=%.1f kcal/mol",
                    Ea, Eo);
            newComment += " Warning: " + warning;
            Logger.info(warning);
            newEa = E;
            Ea = Eo;
        }
        if (p_Hrxn > 0 && Ea < p_Hrxn) {
            // Reaction is endothermic and the barrier is less than the endothermicity.
            warning = String
                    .format("Ea raised by %.1f from %.1f to dHrxn(298K)=%.1f kcal/mol.",
                            p_Hrxn - Ea, Ea, p_Hrxn);
            newComment += " Warning: " + warning;
            Logger.info(warning);
            newEa = newEa.plus((p_Hrxn - Ea));
            Ea = p_Hrxn;
        }
        assert (Ea == newEa.getValue()) : String.format(
                "Something wrong with Ea calculation. %e != %e", Ea,
                newEa.getValue());
        ArrheniusKinetics newK = new ArrheniusKinetics(getA(), getN(), newEa,
                getTRange(), getRank(), getSource(), newComment);
        return newK;
    }
}
/*********************************************************************
 * File Path : RMG\RMG\jing\rxn\ArrheniusEPKinetics.java
 *********************************************************************/

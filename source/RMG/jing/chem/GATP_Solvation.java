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
package jing.chem;

import java.util.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.mathTool.*;
import jing.param.Temperature;
import jing.rxnSys.*;

public class GATP_Solvation implements GeneralSolvationGAPP {
    private static GATP_Solvation INSTANCE = new GATP_Solvation();

    // Constructors
    private GATP_Solvation() {
    }

    public ThermoData generateSolvThermoData(ChemGraph p_chemGraph) {
        /*
         * AJ 16JULY2010 The Pierotti method has now been replaced with the method of MIntz et al. based on some recent
         * comparisons between the 2 methods
         */
        double R = 8.314; // Gas constant units J/mol K
        double T = 298; // Standard state temperature
        // Generation of Abraham Solute Parameters
        AbramData result_Abraham = new AbramData();
        result_Abraham = p_chemGraph.getAbramData();
        // Solute descriptors from the Abraham Model
        double S = result_Abraham.S;
        double B = result_Abraham.B;
        double E = result_Abraham.E;
        double L = result_Abraham.L;
        double A = result_Abraham.A;
        double V = result_Abraham.V;
        SolventData solvent = ReactionModelGenerator.getSolvent();
        String solventname = solvent.name;
        double c_g = solvent.c_g;
        double e_g = solvent.e_g;
        double s_g = solvent.s_g;
        double a_g = solvent.a_g;
        double b_g = solvent.b_g;
        double l_g = solvent.l_g;
        double logK = c_g + s_g * S + b_g * B + e_g * E + l_g * L + a_g * A; // Implementation of Abraham Model for
// calculation of partition coefficient
        double deltaG0 = -8.314 * 298 * 2.303 * logK; // J/mol
        deltaG0 = deltaG0 / 4180; // conversion from kJ/mol to kcal/mol
        // System.out.println("The free energy of solvation in decane at 298K w/o reference state corrections  = " +
// deltaG0_decane +" J/mol for " );
        /*
         * AJ 16JULY 2010 (Mintz method for calculation of solution phase enthalpy
         */
        double c_h = solvent.c_h;
        double e_h = solvent.e_h;
        double s_h = solvent.s_h;
        double a_h = solvent.a_h;
        double b_h = solvent.b_h;
        double l_h = solvent.l_h;
        double deltaH0 = c_g + s_g * S + b_g * B + e_g * E + l_g * L + a_g * A; // Implementation of Mintz model for
// calculation of solution phase enthalpy (kJ/mol)
        deltaH0 = deltaH0 / 4.18; // Conversion from kJ/mol to kcal/mol
        double deltaS0 = (deltaH0 - deltaG0) / T; // kcal/mol/K
        deltaS0 = deltaS0 * 1000; // conversion from kcal/mol/K to cal/mol/K
        // Generation of Gas Phase data to add to the solution phase quantities
        ThermoData solvationCorrection = new ThermoData(deltaH0, deltaS0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                "Solvation correction");
        // Now, solvationCorrection contains solution phase estimates of CORRECTION TO H298, S298 and all the gas phase
// heat capacities.
        // Assuming the solution phase heat capcities to be the same as that in the gas phase we wouls now want to pass
// on this
        // modified version of result to the kinetics codes. This might require reading in a keyword from the
// condition.txt file.
        // Exactly how this will be done is yet to be figured out.
        return solvationCorrection;
    }

    protected static GATP_Solvation getINSTANCE() {
        return INSTANCE;
    }
}

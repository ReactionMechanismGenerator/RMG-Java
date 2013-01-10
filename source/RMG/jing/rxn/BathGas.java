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

import java.util.LinkedHashMap;
import java.util.Iterator;
import jing.chem.Species;
import jing.rxnSys.Logger;
import jing.rxnSys.ReactionSystem;
import jing.rxn.DeltaEDown;

/**
 * Represents the bath gas. Calculates parameters of the bath gas as a weighted average of the components in the bath
 * gas, where the weights are the mole fractions.
 * 
 * @author jwallen
 */
public class BathGas {
    private DeltaEDown dEdown = new DeltaEDown(0, 0, 0);
    private double ljSigma = 0.0;
    private double ljEpsilon = 0.0;
    private double molWt = 0.0;
    private LinkedHashMap colliders;

    public BathGas() {
        colliders = null;
    }

    public BathGas(String inertGas) {
        LinkedHashMap bathGas4PopRxns = new LinkedHashMap();
        bathGas4PopRxns.put(inertGas, new Double(1e-6));
        setColliders(bathGas4PopRxns);
        update();
    }

    public BathGas(ReactionSystem rxnSystem) {
        if (colliders == null)
            colliders = rxnSystem.identifyColliders();
        update();
    }

    public DeltaEDown getDeltaEdown() {
        return dEdown;
    }

    public void setDeltaEdown(DeltaEDown dEdown) {
        this.dEdown = dEdown;
    }

    public double getLJSigma() {
        return ljSigma;
    }

    public double getLJEpsilon() {
        return ljEpsilon;
    }

    public double getMolecularWeight() {
        return molWt;
    }

    public LinkedHashMap getColliders() {
        return colliders;
    }

    public void setColliders(ReactionSystem rxnSystem) {
        colliders = rxnSystem.identifyColliders();
        update();
    }

    public void setColliders(LinkedHashMap Colliders) {
        colliders = Colliders;
    }

    /**
     * Updates the bath gas parameters via a weighted average. What type of average *should* it be? Geometric?
     * Arithmetic?...
     */
    private void update() {
        // Check for null pointer (i.e. colliders not set)
        if (colliders == null)
            throw new NullPointerException();
        // Clear parameters
        dEdown = new DeltaEDown(0, 0, 0);
        ljSigma = 0.0;
        ljEpsilon = 0.0;
        molWt = 0.0;
        // Determine bath gas concentration (i.e. total concentration of colliders)
        double totalConc = 0;
        for (Iterator iter = colliders.values().iterator(); iter.hasNext();) {
            totalConc += ((Double) iter.next()).doubleValue();
        }
        // Calculate values as weighted average
        for (Iterator iter = colliders.keySet().iterator(); iter.hasNext();) {
            Object key = iter.next();
            if (key instanceof Species) {
                Species spe = (Species) key;
                double conc = ((Double) colliders.get(spe)).doubleValue();
                double mf = conc / totalConc;
                molWt += mf * spe.getMolecularWeight();
                ljSigma += mf * spe.getChemkinTransportData().getSigma();
                ljEpsilon += mf * spe.getChemkinTransportData().getEpsilon();
            } else if (key instanceof String) {
                String name = (String) key;
                double conc = ((Double) colliders.get(key)).doubleValue();
                double mf = conc / totalConc;
                if (name.equals("Ar") || name.equals("AR")) {
                    // Taken from West et al 2009 (doi:10.1016/j.combustflame.2009.04.011)
                    dEdown.setParameters(150.0, 300, 0.85);
                    molWt += mf * 39.95;
                    /*
                     * Numbers from Table K.2 of "Fundamentals of Momentum, Heat, and Mass Transfer", Welty, Wicks,
                     * Wilson, Rorrer 4th ed. Numbers come from R.C. Reid and T.K. Sherwood, "The Properties of Gases
                     * and Liquids", McGraw-Hill Book Company, New York, 1958. MRH 16-Jun-2009
                     */
                    ljSigma += mf * 3.418; // Units of Angstroms
                    ljEpsilon += mf * 124; // Units of Kelvin (actually epsilon/boltzmann constant)
                } else if (name.equals("N2")) {
                    // Taken from Zhang et al 2011 (doi:doi:10.1016/j.proci.2010.05.010)
                    dEdown.setParameters(200.0, 300, 0.85);
                    molWt += mf * 28.01;
// ljEpsilon += mf * 97.5;
// ljSigma += mf * 3.62;
                    /*
                     * Numbers from Table K.2 of "Fundamentals of Momentum, Heat, and Mass Transfer", Welty, Wicks,
                     * Wilson, Rorrer 4th ed. Numbers come from R.C. Reid and T.K. Sherwood, "The Properties of Gases
                     * and Liquids", McGraw-Hill Book Company, New York, 1958. Changed the numbers already stored in RMG
                     * for consistency MRH 16-Jun-2009
                     */
                    ljSigma += mf * 3.681; // Units of Angstroms
                    ljEpsilon += mf * 91.5; // Units of Kelvin (actually epsilon/boltzmann constant)
                } else if (name.equals("He") || name.equals("HE")) {
                    // Taken from Jasper and Miller 2009 (doi:10.1021/jp900802f)
                    dEdown.setParameters(110.0, 300, 0.81);
                    molWt += mf * 4.00;
                    /*
                     * Numbers from Table K.2 of "Fundamentals of Momentum, Heat, and Mass Transfer", Welty, Wicks,
                     * Wilson, Rorrer 4th ed. Numbers come from R.C. Reid and T.K. Sherwood, "The Properties of Gases
                     * and Liquids", McGraw-Hill Book Company, New York, 1958. MRH 16-Jun-2009
                     */
                    ljSigma += mf * 2.576; // Units of Angstroms
                    ljEpsilon += mf * 10.22; // Units of Kelvin (actually epsilon/boltzmann constant)
                } else if (name.equals("Ne") || name.equals("NE")) {
                    // Taken from Giri et al. 2008 (doi:10.1039/B808168A)
                    dEdown.setParameters(150.0, 300, 0.4);
                    molWt += mf * 20.18;
                    /*
                     * Numbers from Table K.2 of "Fundamentals of Momentum, Heat, and Mass Transfer", Welty, Wicks,
                     * Wilson, Rorrer 4th ed. Numbers come from R.C. Reid and T.K. Sherwood, "The Properties of Gases
                     * and Liquids", McGraw-Hill Book Company, New York, 1958. MRH 30-Jul-2009
                     */
                    ljSigma += mf * 2.789; // Units of Angstroms
                    ljEpsilon += mf * 35.7; // Units of Kelvin (actually epsilon/boltzmann constant)
                } else {
                    Logger.critical("Unknown colliders: " + name);
                    System.exit(0);
                }
            } else {
                Logger.critical("Unknown colliders: " + key.toString());
                System.exit(0);
            }
        }
        // Convert to units used by FAME
        dEdown.setAlpha(dEdown.getAlpha() * 2.9979e10 * 6.626e-34 * 6.022e23
                / 1000); // cm^-1 --> kJ/mol
        ljSigma *= 1e-10; // A --> m
        ljEpsilon *= 1.381e-23; // K --> J
    }
}

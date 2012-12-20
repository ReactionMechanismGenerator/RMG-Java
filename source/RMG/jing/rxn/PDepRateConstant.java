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

/**
 * A representation of a pressure-dependent rate coefficient k(T, P).
 * 
 * @author jwallen
 */
public class PDepRateConstant {
    // ==========================================================================
    //
    // Data members
    //
    /**
     * An enumeraction of pressure-dependent rate coefficient evaluation methods:
     * <ul>
     * <li>NONE - The method could not be assessed.
     * <li>INTERPOLATE - Bilinear interpolation on T^-1, log P, and log k(T, P) axes.
     * <li>CHEBYSHEV - Evaluation of Chebyshev polynomials.
     * <li>PDEPARRHENIUS - Pressure-dependent Arrhenius equations.
     * </ul>
     */
    public enum Mode {
        NONE, INTERPOLATE, CHEBYSHEV, PDEPARRHENIUS, RATE
    };

    /**
     * The default mode to use for evaluation the pressure-dependent rate coefficients. Default is interpolate.
     */
    private static Mode mode = Mode.INTERPOLATE;
    // The mode for a specific rate coefficient (so it can be non-default)
    private Mode thisMode;
    /**
     * A list of the temperatures at which the rate coefficient has been explicitly calculated. For now this is the same
     * for all pressure- dependent rates, and so is static to conserve memory.
     */
    private static Temperature[] temperatures;
    /**
     * A list of the pressures at which the rate coefficient has been explicitly calculated. For now this is the same
     * for all pressure- dependent rates, and so is static to conserve memory.
     */
    private static Pressure[] pressures;
    /**
     * A matrix of explicitly-evaluated rate coefficient values in s^-1, cm^3 mol^-1 s^-1, etc. at each temperature
     * (rows) and pressure (columns) in the above lists. These values will be used for interpolation.
     */
    private double[][] rateConstants;
    /**
     * The rate coefficient as represented by Chebyshev polynomials.
     */
    private ChebyshevPolynomials chebyshev;
    /**
     * The rate coefficient as represented by pressure-dependent Arrhenius kinetics.
     */
    private PDepArrheniusKinetics[] pDepArrhenius;
    private static Temperature TMIN;
    private static Temperature TMAX;
    private static Pressure PMIN;
    private static Pressure PMAX;

    // ==========================================================================
    //
    // Constructors
    //
    /**
     * Constructor.
     * 
     * @param rates
     * @param cheb
     */
    public PDepRateConstant(double[][] rates) {
        rateConstants = rates;
        chebyshev = null;
        pDepArrhenius = null;
    }

    public PDepRateConstant() {
        rateConstants = null;
        chebyshev = null;
        pDepArrhenius = null;
    }

    public PDepRateConstant(double[][] rates, ChebyshevPolynomials chebyPols) {
        rateConstants = rates;
        chebyshev = chebyPols;
        pDepArrhenius = null;
        setThisMode(Mode.CHEBYSHEV);
    }

    public PDepRateConstant(double[][] rates, PDepArrheniusKinetics plogKinetics) {
        rateConstants = rates;
        chebyshev = null;
        pDepArrhenius = new PDepArrheniusKinetics[1];
        pDepArrhenius[0] = plogKinetics;
        setThisMode(Mode.PDEPARRHENIUS);
    }

    public PDepRateConstant(PDepArrheniusKinetics plogKinetics) {
        rateConstants = null;
        chebyshev = null;
        pDepArrhenius = new PDepArrheniusKinetics[1];
        pDepArrhenius[0] = plogKinetics;
        setThisMode(Mode.PDEPARRHENIUS);
    }

    public PDepRateConstant(ChebyshevPolynomials chebyPols) {
        rateConstants = null;
        chebyshev = chebyPols;
        pDepArrhenius = null;
        setThisMode(Mode.CHEBYSHEV);
    }

    // ==========================================================================
    //
    // Static methods
    //
    public Mode getMode() {
        if (thisMode != null)
            return thisMode;
        else
            return mode;
    }

    public static Mode getDefaultMode() {
        return mode;
    }

    public static void setDefaultMode(Mode m) {
        mode = m;
    }

    public void setThisMode(Mode m) {
        thisMode = m;
    }

    public static Temperature[] getTemperatures() {
        return temperatures;
    }

    public static void setTemperatures(Temperature[] temps) {
        temperatures = temps;
    }

    public static Pressure[] getPressures() {
        return pressures;
    }

    public static void setPressures(Pressure[] press) {
        pressures = press;
    }

    // ==========================================================================
    //
    // Other methods
    //
    public double calculateRate(Temperature temperature, Pressure pressure) {
        double rate = 0.0;
        if ((TMIN != null && temperature.getK() < TMIN.getK())
                || (TMAX != null && temperature.getK() > TMAX.getK()))
            throw new IllegalArgumentException(
                    String.format(
                            "Tried to evaluate P-dep rate coefficient at T=%.1f K but Pdep calculations only valid from %.1f to %.1f K",
                            temperature.getK(), TMIN.getK(), TMAX.getK()));
        if ((PMIN != null && pressure.getBar() < PMIN.getBar())
                || (PMAX != null && pressure.getBar() > PMAX.getBar()))
            throw new IllegalArgumentException(
                    String.format(
                            "Tried to evaluate P-dep rate coefficient at P=%.2g bar but Pdep calculations only valid from %s to %s bar",
                            pressure.getBar(), PMIN.getBar(), PMAX.getBar()));
        if (getMode() == Mode.INTERPOLATE || getMode() == Mode.RATE) {
            /*
             * MRH 10Feb2010 I am initializing the t1, t2, p1, p2 indices to be zero. In the case of the temperature of
             * interest being equal to the lowest temperature, t1 would not be re-defined (ditto for the pressure)
             */
            int t1 = 0, t2 = 0, p1 = 0, p2 = 0;
            // int t1 = -1, t2 = -1, p1 = -1, p2 = -1;
            double x = 0.0, x1 = 0.0, x2 = 0.0, y = 0.0, y1 = 0.0, y2 = 0.0;
            double z11 = 0.0, z12 = 0.0, z21 = 0.0, z22 = 0.0;
            for (int t = 0; t < temperatures.length - 1; t++) {
                if (temperatures[t].getK() < temperature.getK()) {
                    t1 = t;
                }
            }
            for (int p = 0; p < pressures.length - 1; p++) {
                if (pressures[p].getBar() < pressure.getBar()) {
                    p1 = p;
                }
            }
            x1 = 1.0 / temperatures[t1].getK();
            t2 = t1 + 1;
            x2 = 1.0 / temperatures[t2].getK();
            if (pressures.length == 1) {
                x = 1.0 / temperature.getK();
                double z1 = Math.log10(rateConstants[t1][0]);
                double z2 = Math.log10(rateConstants[t2][0]);
                rate = (z2 - z1) * (x - x1) / (x2 - x1) + z1;
                rate = Math.pow(10, rate);
            } else {
                y1 = Math.log10(pressures[p1].getPa());
                p2 = p1 + 1;
                y2 = Math.log10(pressures[p2].getPa());
                x = 1.0 / temperature.getK();
                y = Math.log10(pressure.getPa());
                z11 = Math.log10(rateConstants[t1][p1]);
                z12 = Math.log10(rateConstants[t1][p2]);
                z21 = Math.log10(rateConstants[t2][p1]);
                z22 = Math.log10(rateConstants[t2][p2]);
                rate = (z11 * (x2 - x) * (y2 - y) + z21 * (x - x1) * (y2 - y)
                        + z12 * (x2 - x) * (y - y1) + z22 * (x - x1) * (y - y1))
                        / ((x2 - x1) * (y2 - y1));
                rate = Math.pow(10, rate);
            }
        } else if (getMode() == Mode.CHEBYSHEV && chebyshev != null) {
            rate = chebyshev.calculateRate(temperature, pressure);
        } else if (getMode() == Mode.PDEPARRHENIUS && pDepArrhenius != null) {
            rate = 0;
            for (int i = 0; i < pDepArrhenius.length; i++) {
                rate += pDepArrhenius[i].calculateRate(temperature, pressure);
            }
        } else {
            throw new RuntimeException(
                    "Failed to evaluate P-dep rate coefficient with type "
                            + getMode());
        }
        return rate;
    }

    public ChebyshevPolynomials getChebyshev() {
        return chebyshev;
    }

    public void setChebyshev(ChebyshevPolynomials cheb) {
        chebyshev = cheb;
    }

    public void setChebyshev(double[][] cheb) {
        if (TMIN == null) {
            chebyshev = new ChebyshevPolynomials(cheb.length, temperatures[0],
                    temperatures[temperatures.length - 1], cheb[0].length,
                    pressures[0], pressures[pressures.length - 1], cheb);
        } else {
            chebyshev = new ChebyshevPolynomials(cheb.length, TMIN, TMAX,
                    cheb[0].length, PMIN, PMAX, cheb);
        }
    }

    public PDepArrheniusKinetics[] getPDepArrheniusKinetics() {
        return pDepArrhenius;
    }

    public void setPDepArrheniusKinetics(PDepArrheniusKinetics[] kin) {
        // Store the passed in list of kinetics
        pDepArrhenius = kin;
    }

    public void setPDepArrheniusKinetics(PDepArrheniusKinetics kin) {
        // Make a new list of size 1, containing only the passed in kinetics
        pDepArrhenius = new PDepArrheniusKinetics[1];
        pDepArrhenius[0] = kin;
    }

    public void addPDepArrheniusKinetics(PDepArrheniusKinetics p_kinetics) {
        // Add the passed in kinetics to the list (create the list if none exists)
        if (p_kinetics == null)
            return;
        if (getMode() != PDepRateConstant.Mode.PDEPARRHENIUS)
            throw new RuntimeException(
                    String.format(
                            "Cannot add PDepArrheniusKinetics Rate constants to %s rate.",
                            getMode()));
        if (pDepArrhenius == null) {
            pDepArrhenius = new PDepArrheniusKinetics[1];
            pDepArrhenius[0] = p_kinetics;
        } else {
            PDepArrheniusKinetics[] tempKinetics = pDepArrhenius;
            pDepArrhenius = new PDepArrheniusKinetics[tempKinetics.length + 1];
            for (int i = 0; i < tempKinetics.length; i++) {
                pDepArrhenius[i] = tempKinetics[i];
            }
            pDepArrhenius[pDepArrhenius.length - 1] = p_kinetics;
        }
    }

    public static void setPMin(Pressure p) {
        PMIN = p;
    }

    public static void setPMax(Pressure p) {
        PMAX = p;
    }

    public static void setTMin(Temperature t) {
        TMIN = t;
    }

    public static void setTMax(Temperature t) {
        TMAX = t;
    }

    public static Pressure getPMin() {
        return PMIN;
    }

    public static Pressure getPMax() {
        return PMAX;
    }

    public static Temperature getTMin() {
        return TMIN;
    }

    public static Temperature getTMax() {
        return TMAX;
    }

    public double[][] getRateConstants() {
        return rateConstants;
    }
}

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

import jing.param.Temperature;

/**
 * A class for representing the exponential down parameter (used in pressure dependence calculations) as a power law of
 * the form dEdown = alpha * (T / T0)^n where the parameters are alpha, T0 (in K), and n. The units of dEdown are the
 * same as the units of alpha; this is left to the end user to choose and maintain.
 */
public class DeltaEDown {
    private double alpha;
    private double T0;
    private double n;

    public DeltaEDown(double alpha, double T0, double n) {
        setParameters(alpha, T0, n);
    }

    public double getAlpha() {
        return alpha;
    }

    public double getT0() {
        return T0;
    }

    public double getN() {
        return n;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public void setT0(double T0) {
        this.T0 = T0;
    }

    public void setN(double n) {
        this.n = n;
    }

    public void setParameters(double alpha, double T0, double n) {
        setAlpha(alpha);
        setT0(T0);
        setN(n);
    }

    public double evaluate(Temperature temperature) {
        double T = temperature.getK();
        if (T0 == 0)
            return alpha;
        else
            return alpha * Math.pow(T / T0, n);
    }
}

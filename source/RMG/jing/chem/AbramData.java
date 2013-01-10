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

public class AbramData extends AbrahamGAValue {
    /*
     * Contains AbramGAValue parameters (S,B,E,L,A) and also V, the McGowan's Volume. If created from an AbrahamGAValue
     * instance, V=0
     */
    public double V = 0;
    String comments = "";
    String source = "";
    String name = "";
    boolean fromPrimaryAbrahamLibrary = false;

    public AbramData() {
        super();
        V = 0;
    }

    public AbramData(double S, double B, double E, double L, double A) {
        super(S, B, E, L, A);
        V = 0;
    }

    public AbramData(double S, double B, double E, double L, double A,
            double p_V) {
        super(S, B, E, L, A);
        V = p_V;
    }

    public AbramData(String p_name, AbramData p_abramdata, String p_comments,
            String p_source, boolean truefalse) {
        name = p_name;
        comments = p_comments;
        source = p_source;
        S = p_abramdata.S;
        E = p_abramdata.E;
        L = p_abramdata.L;
        B = p_abramdata.B;
        A = p_abramdata.A;
        V = p_abramdata.V;
        fromPrimaryAbrahamLibrary = truefalse;
    }

    public AbramData(double p_S, double p_B, double p_E, double p_L,
            double p_A, double p_V, boolean truefalse) {
        S = p_S;
        E = p_E;
        A = p_A;
        B = p_B;
        L = p_L;
        V = p_V;
        fromPrimaryAbrahamLibrary = truefalse;
    }

    public AbramData(AbrahamGAValue p_ga) {
        super(p_ga);
        V = 0;
    }

    public void plus(AbrahamGAValue p_ga) {
        if (p_ga == null)
            return;
        S += p_ga.S;
        B += p_ga.B;
        E += p_ga.E;
        L += p_ga.L;
        A += p_ga.A;
    }

    public String getName() {
        return name;
    }

    public String getSource() {
        return source;
    }

    public String toString() {
        // Get the parent class's toString() and append the V value.
        String s = super.toString();
        s = s + '\t' + String.valueOf(V);
        return s;
    }
}

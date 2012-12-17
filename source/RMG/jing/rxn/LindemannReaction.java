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
import java.util.StringTokenizer;
import jing.mathTool.MathTool;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxnSys.SystemSnapshot;

public class LindemannReaction extends ThirdBodyReaction {
    protected ArrheniusKinetics low;

    private LindemannReaction() {
    }

    public double calculateRate(SystemSnapshot p_presentStatus) {
        Temperature temp = p_presentStatus.getTemperature();
        double rate = super.calculateTotalRate(temp);
        rate *= calculateLindemannFallOff(p_presentStatus);
        return rate;
    }

    public double calculateLindemannFallOff(SystemSnapshot p_presentStatus) {
        Temperature temp = p_presentStatus.getTemperature();
        double M = calculateThirdBodyCoefficient(p_presentStatus);
        double kZero = low.calculateRate(temp, -1);
        double kInf = 0.0;
        Kinetics[] k_array = getKinetics();
        for (int i = 0; i < k_array.length; i++) {
            kInf += k_array[i].calculateRate(temp, -1);
        }
        double Pr = kZero * M / kInf;
        double fallOffFactor = (Pr / (1.0 + Pr));
        return fallOffFactor;
    }

    public String formPDepSign(String p_string) {
        StringTokenizer st = new StringTokenizer(p_string, "=");
        String s1 = st.nextToken();
        s1 += "(+m)=";
        String s2 = st.nextToken();
        s2 += "(+m)";
        return (s1 + s2);
    }

    public void generateReverseReaction() {
        LindemannReaction r = new LindemannReaction();
        r.structure = getStructure().generateReverseStructure();
        r.kinetics = getKinetics();
        for (int i = 0; i < r.kinetics.length; i++) {
            r.comments = "Reverse reaction";
        }
        r.weightMap = weightMap;
        r.low = low;
        r.setReverseReaction(this);
        this.setReverseReaction(r);
        return;
    }

    public static LindemannReaction make(Reaction p_reaction,
            LinkedHashMap p_weightMap, final ArrheniusKinetics p_low) {
        LindemannReaction lr = new LindemannReaction();
        lr.structure = p_reaction.getStructure();
        lr.kinetics = p_reaction.getKinetics();
        lr.comments = p_reaction.getComments();
        lr.weightMap = p_weightMap;
        lr.low = p_low;
        lr.generateReverseReaction();
        return lr;
    }

    public String toChemkinString(Temperature p_temperature) {
        String s = super.toChemkinString(p_temperature) + '\n';
        s += "LOW/ "
                + low.toChemkinString(calculateHrxn(p_temperature),
                        p_temperature, false) + " /\n";
        return s;
    }

    public String toRestartString(Temperature p_temperature) {
        String s = super.toRestartString(p_temperature) + '\n';
        s += "LOW/ "
                + low.toChemkinString(calculateHrxn(p_temperature),
                        p_temperature, false) + " /\n";
        return s;
    }

    public String toChemkinString(Temperature p_temperature, Pressure p_pressure) {
        String s = super.toChemkinString(p_temperature) + "\n";
        s += "LOW/ "
                + low.toChemkinString(calculateHrxn(p_temperature),
                        p_temperature, false) + "/\n";
        return s;
    }

    public String toString(Temperature p_temperature) {
        String s = getStructure().toChemkinString(true).toString() + '\n';
        for (int i = 0; i < getKinetics().length; i++) {
            s += "kInf = "
                    + getKinetics()[i].toChemkinString(
                            calculateHrxn(p_temperature), p_temperature, false)
                    + '\n';
            s += "kZero = "
                    + low.toChemkinString(calculateHrxn(p_temperature),
                            p_temperature, false) + '\n';
        }
        return s;
    }

    public ArrheniusKinetics getLow() {
        return low;
    }

    public void setLow(ArrheniusKinetics p_low) {
        low = p_low;
    }
}

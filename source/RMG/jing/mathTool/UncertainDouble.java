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
package jing.mathTool;

import java.util.*;

// ## package jing::mathTool
// ----------------------------------------------------------------------------
// jing\mathTool\UncertainDouble.java
// ----------------------------------------------------------------------------
// ## class UncertainDouble
public class UncertainDouble {
    protected String type; // ## attribute type
    protected double uncertainty; // ## attribute uncertainty
    protected double value; // ## attribute value

    // Constructors
    // ## operation UncertainDouble(double,double,String)
    public UncertainDouble(double p_value, double p_uncertainty, String p_type) {
        // #[ operation UncertainDouble(double,double,String)
        value = p_value;
        uncertainty = Math.abs(p_uncertainty);
        if (p_type.compareToIgnoreCase("Multiplying") == 0
                || p_type.compareToIgnoreCase("Multiplier") == 0
                || p_type.compareToIgnoreCase("M") == 0) {
            type = "Multiplier";
        } else if (p_type.compareToIgnoreCase("Adding") == 0
                || p_type.compareToIgnoreCase("Adder") == 0
                || p_type.compareToIgnoreCase("A") == 0) {
            type = "Adder";
        } else
            throw new InvalidUncertaintyTypeException(p_type);
        // #]
    }

    public UncertainDouble() {
    }

    // ## operation getLowerBound()
    public double getLowerBound() {
        // #[ operation getLowerBound()
        if (isAddingUncertainty()) {
            return getValue() - getUncertainty();
        } else if (isMultiplyingUncertainty()) {
            if (getUncertainty() < 1.0)
                throw new InvalidUncertaintyException(
                        "Uncertainty mutiplier is less than one");
            return getValue() / getUncertainty();
        } else
            throw new InvalidUncertaintyTypeException();
        // #]
    }

    // ## operation getUpperBound()
    public double getUpperBound() {
        // #[ operation getUpperBound()
        if (isAddingUncertainty()) {
            return getValue() + getUncertainty();
        } else if (isMultiplyingUncertainty()) {
            if (getUncertainty() < 1.0)
                throw new InvalidUncertaintyException(
                        "Uncertainty mutiplier is less than one");
            return getValue() * getUncertainty();
        } else
            throw new InvalidUncertaintyTypeException();
        // #]
    }

    // ## operation getValue()
    public double getValue() {
        // #[ operation getValue()
        return value;
        // #]
    }

    // ## operation isAddingUncertainty()
    public boolean isAddingUncertainty() {
        // #[ operation isAddingUncertainty()
        return type.equals("Adder");
        // #]
    }

    // ## operation isMultiplyingUncertainty()
    public boolean isMultiplyingUncertainty() {
        // #[ operation isMultiplyingUncertainty()
        return type.equals("Multiplier");
        // #]
    }

    // ## operation multiply(double)
    public UncertainDouble multiply(double p_multiplier) {
        // #[ operation multiply(double)
        if (isAddingUncertainty())
            return new UncertainDouble(value * p_multiplier, uncertainty
                    * p_multiplier, type);
        if (isMultiplyingUncertainty())
            return new UncertainDouble(value * p_multiplier, uncertainty, type);
        throw new RuntimeException("Can't multiply UncertainDouble with type "
                + type);
        // #]
    }

    // ## operation plus(double)
    public UncertainDouble plus(double p_adder) {
        // #[ operation plus(double)
        return new UncertainDouble(value + p_adder, uncertainty, type);
        // #]
    }

    // ## operation plus(UncertainDouble)
    public UncertainDouble plus(UncertainDouble p_adder) {
        double newValue = value + p_adder.getValue();
        double newUncertainty;
        if (isMultiplyingUncertainty() && p_adder.isMultiplyingUncertainty()) {
            newUncertainty = Math.max(uncertainty, p_adder.getUncertainty()); // may overestimate error slightly?
            return new UncertainDouble(newValue, uncertainty, "Multiplier");
        } else {
            newUncertainty = getAddingUncertainty()
                    + p_adder.getAddingUncertainty();
            return new UncertainDouble(newValue, newUncertainty, "Adder");
        }
    }

    // ## operation getAddingUncertainty()
    public double getAddingUncertainty() {
        // return the uncertainty in an adding type
        // i.e. if it's a multiplying type, subtract one and multiply it by 'value'.
        // This exaggerates how low the lower bound will be.
        // (e.g. 10*/2 will become 10+-10.
        if (isAddingUncertainty()) {
            return uncertainty;
        } else if (isMultiplyingUncertainty()) {
            return value * (uncertainty - 1.0);
        } else
            throw new InvalidUncertaintyTypeException();
    }

    // ## operation toString()
    public String toString() {
        // #[ operation toString()
        if (isAddingUncertainty()) {
            return String.valueOf(value) + "+|-" + String.valueOf(uncertainty);
        } else if (isMultiplyingUncertainty()) {
            return String.valueOf(value) + "*|/" + String.valueOf(uncertainty);
        } else
            throw new InvalidUncertaintyTypeException(type);
        // #]
    }

    public String getType() {
        return type;
    }

    public double getUncertainty() {
        return uncertainty;
    }
}
/*********************************************************************
 * File Path : RMG\RMG\jing\mathTool\UncertainDouble.java
 *********************************************************************/

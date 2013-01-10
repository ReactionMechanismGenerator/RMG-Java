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
package jing.rxnSys;

import java.util.*;
import jing.chem.Species;

// ## package jing::rxnSys
// ----------------------------------------------------------------------------
// jing\rxnSys\SensitivityStatus.java
// ----------------------------------------------------------------------------
// ## class SensitivityStatus
// svp
public class SensitivityStatus {
    protected double sensitivity;
    protected double sflux;
    protected int species_number;
    protected int reaction_number;

    // ## operation SensitivityStatus(double,double,int,int)
    public SensitivityStatus(double p_sensitivity, double p_sflux, int p_snum,
            int p_rnum) {
        // #[ operation SensitivityStatus(double,double,int,int)
        sensitivity = p_sensitivity;
        sflux = p_sflux;
        species_number = p_snum;
        reaction_number = p_rnum;
        // #]
    }

    public double getSensitivity() {
        return sensitivity;
    }

    public double getSFlux() {
        return sflux;
    }

    public int getSpeciesNumber() {
        return species_number;
    }

    public int getReactionNumber() {
        return reaction_number;
    }
}
/*********************************************************************
 * File Path : RMG\RMG\jing\rxnSys\SensitivityStatus.java
 *********************************************************************/

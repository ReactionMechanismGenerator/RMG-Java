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

public class LindemannODEReaction extends ODEReaction {
    protected int[] colliders = null;
    protected double[] efficiency = null;
    protected int numCollider;
    protected double inertColliderEfficiency;
    protected double lowRate;
    protected double highRate;

    public LindemannODEReaction(int p_rNum, int p_pNum, int[] p_rID,
            int[] p_pID, int p_direction, double p_Keq, int[] p_collider,
            double[] p_efficiency, int p_numCollider,
            double p_inertColliderEfficiency, double p_highRate,
            double p_lowRate) {
        rNum = p_rNum;
        pNum = p_pNum;
        rID = p_rID;
        pID = p_pID;
        direction = p_direction;
        Keq = p_Keq;
        setRate = false;
        highRate = p_highRate;
        lowRate = p_lowRate;
        colliders = p_collider;
        efficiency = p_efficiency;
        numCollider = p_numCollider;
        inertColliderEfficiency = p_inertColliderEfficiency;
    }
}

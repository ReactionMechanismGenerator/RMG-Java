////////////////////////////////////////////////////////////////////////////////
//
//  RMG - Reaction Mechanism Generator
//
//  Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
//  RMG Team (rmg_dev@mit.edu)
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

package jing.chem;

import java.util.HashMap;
import java.util.Map;

/**
 * A set of Benson-like ring corrections, along with a flag for indicating if
 * the corrections are approximate. This allows us to both apply the best
 * available ring correction (BensonTDGenerator) while still falling back to
 * QMTP if the corrections are approximate (HybridTDGenerator). *
 */
public class BensonRingCorrections {

    /**
     * Default constructor.
     */
    public BensonRingCorrections() {
        corrections = new HashMap<ThermoGAValue, Integer>();
        clear();
    }
    
    /**
     * Clear the ring corrections and the approximate match flag.
     */
    public void clear() {
        corrections.clear();
        imperfectMatch = false;
    }
    
    /**
     * Return the set of ring corrections, a mapping of group additivity values
     * with associated integer degeneracies.
     */
    public Map<ThermoGAValue, Integer> getCorrections() {
        return corrections;
    }
    
    /**
     * Add a ring correction to the set of corrections.
     * @param ga The Benson group additivity values for the ring correction
     * @param value The degeneracy of the correction
     */
    public void addCorrection(ThermoGAValue ga, int value) {
        if (corrections.containsKey(ga)) {
            corrections.put(ga, corrections.get(ga) + value);
        }
        else{
            corrections.put(ga,value);
        }
    }
    
    /**
     * Return true if there are no ring corrections, or false otherwise.
     */
    public boolean isEmpty() {
        return corrections.isEmpty();
    }
    
    /**
     * Return true if the ring corrections are approximate, or false otherwise.
     */
    public boolean isImperfectMatch() {
        return imperfectMatch;
    }
    
    /**
     * Set the flag indicating whether or not the corrections are approximate.
     */
    public void setImperfectMatch(boolean imperfect) {
        imperfectMatch = imperfect;
    }
    
    /** The set of group additivity ring corrections, mapped to integer degeneracies. */
    private Map<ThermoGAValue, Integer> corrections;
    
    /** A flag indicating the quality of the matches. */
    private boolean imperfectMatch;
            
}

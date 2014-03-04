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

// ## package jing::chem
// ----------------------------------------------------------------------------
// jing\chem\FunctionalGroupCollection.java
// ----------------------------------------------------------------------------
// ## class FunctionalGroupCollection
public class FunctionalGroupCollection implements Matchable {
    protected String name; // ## attribute name
    protected LinkedHashSet functionalGroups;

    // Constructors
    // ## operation FunctionalGroupCollection(String)
    public FunctionalGroupCollection(String p_name) {
        {
            functionalGroups = new LinkedHashSet();
        }
        // #[ operation FunctionalGroupCollection(String)
        name = p_name;
        // #]
    }

    public FunctionalGroupCollection() {
        {
            functionalGroups = new LinkedHashSet();
        }
    }

    public String toString() {
        Iterator iter = ((FunctionalGroupCollection) this)
                 .getFunctionalGroups();
        String t = new String();
        while (iter.hasNext()) {
             FunctionalGroup fg = (FunctionalGroup) iter.next();
             t = t + fg.getGraph().toString();
        }
        return t;
    }

    // ## operation addFunctionalGroups(FunctionalGroup)
    public void addFunctionalGroups(FunctionalGroup p_FunctionalGroup) {
        // #[ operation addFunctionalGroups(FunctionalGroup)
        functionalGroups.add(p_FunctionalGroup);
        // #]
    }

    // ## operation addFunctionalGroups(Matchable)
    public void addFunctionalGroups(Matchable p_functionalGroup)
            throws InvalidFunctionalGroupException {
        // #[ operation addFunctionalGroups(Matchable)
        if (p_functionalGroup instanceof FunctionalGroup) {
            FunctionalGroup fg = (FunctionalGroup) p_functionalGroup;
            addFunctionalGroups(fg);
        } else if (p_functionalGroup instanceof FunctionalGroupCollection) {
            FunctionalGroupCollection fgc = (FunctionalGroupCollection) p_functionalGroup;
            addFunctionalGroups(fgc);
        } else
            throw new InvalidFunctionalGroupException();
        // #]
    }

    // ## operation addFunctionalGroups(FunctionalGroupCollection)
    public void addFunctionalGroups(FunctionalGroupCollection p_fgc) {
        // #[ operation addFunctionalGroups(FunctionalGroupCollection)
        Iterator iter = p_fgc.getFunctionalGroups();
        while (iter.hasNext()) {
            FunctionalGroup fg = (FunctionalGroup) iter.next();
            addFunctionalGroups(fg);
        }
        // #]
    }

    // ## operation contains(FunctionalGroup)
    public boolean contains(FunctionalGroup p_fg) {
        // #[ operation contains(FunctionalGroup)
        if (functionalGroups == null)
            return false;
        else
            return functionalGroups.contains(p_fg);
        // #]
    }

    // ## operation isSubAtCentralNodes(Matchable)
    public boolean isSubAtCentralNodes(Matchable p_functional) {
        // #[ operation isSubAtCentralNodes(Matchable)
        if (this == p_functional)
            return false;
        if (!(p_functional instanceof FunctionalGroupCollection)) {
            boolean found = false;
            Iterator iter1 = functionalGroups.iterator();
            while (iter1.hasNext()) {
                FunctionalGroup fg1 = (FunctionalGroup) iter1.next();
                if (fg1.equals(p_functional)
                        || fg1.isSubAtCentralNodes(p_functional)) {
                    found = true;
                } else {
                    found = false;
                    break;
                }
            }
            if (!found)
                return false;
            return true;
        }
        Collection c1 = functionalGroups;
        Collection c2 = ((FunctionalGroupCollection) p_functional).functionalGroups;
        if (c2.size() == c1.size() && c2.containsAll(c1))
            return false;
        boolean found = true;
        Iterator iter1 = c1.iterator();
        while (iter1.hasNext()) {
            found = false;
            FunctionalGroup fg1 = (FunctionalGroup) iter1.next();
            Iterator iter2 = c2.iterator();
            while (iter2.hasNext()) {
                FunctionalGroup fg2 = (FunctionalGroup) iter2.next();
                if (fg1.equals(fg2) || fg1.isSubAtCentralNodes(fg2)) {
                    found = true;
                    break;
                }
            }
            if (!found)
                return false;
        }
        return true;
        // #]
    }

    public String getName() {
        return name;
    }

    public void setName(String p_name) {
        name = p_name;
    }

    public Iterator getFunctionalGroups() {
        Iterator iter = functionalGroups.iterator();
        return iter;
    }

    public void removeFunctionalGroups(FunctionalGroup p_FunctionalGroup) {
        functionalGroups.remove(p_FunctionalGroup);
    }

    public void clearFunctionalGroups() {
        functionalGroups.clear();
    }
}
/*********************************************************************
 * File Path : RMG\RMG\jing\chem\FunctionalGroupCollection.java
 *********************************************************************/

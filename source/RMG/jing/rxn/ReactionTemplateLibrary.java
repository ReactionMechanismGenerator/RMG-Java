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

import java.io.*;
import java.util.*;
import jing.chemParser.ChemParser;
import jing.param.Temperature;
import jing.rxnSys.Logger;

// ## package jing::rxn
// ----------------------------------------------------------------------------
// jing\rxn\ReactionTemplateLibrary.java
// ----------------------------------------------------------------------------
// ## class ReactionTemplateLibrary
public class ReactionTemplateLibrary {
    protected boolean ForceTheremoConsistence = false; // ## attribute ForceTheremoConsistence
    protected static ReactionTemplateLibrary INSTANCE = new ReactionTemplateLibrary(); // ## attribute INSTANCE
    protected LinkedHashSet reactionTemplate;

    // Constructors
    // ## operation ReactionTemplateLibrary()
    public ReactionTemplateLibrary() {
        reactionTemplate = new LinkedHashSet();
        // Get the directory which contains the reaction families
        String kineticsDirectory = System
                .getProperty("jing.rxn.ReactionTemplateLibrary.pathName");
        if (kineticsDirectory == null) {
            Logger.critical("Undefined system property: jing.rxn.ReactionTemplateLibrary.pathName!");
            System.exit(0);
        }
        String separator = System.getProperty("file.separator");
        if (!kineticsDirectory.endsWith(separator))
            kineticsDirectory = kineticsDirectory + separator;
        Logger.info("Reading kinetics database from " + kineticsDirectory);
        // Read the file families.txt
        try {
            Logger.info("\nReading kinetics groups database from: "
                    + kineticsDirectory + "\n");
            String familiesPath = kineticsDirectory + "families.txt";
            FileReader in = new FileReader(familiesPath);
            BufferedReader data = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(data, true);
            while (line != null) {
                // Tokenize string into its component tokens
                StringTokenizer token = new StringTokenizer(line);
                String index = token.nextToken();
                String status = token.nextToken().toLowerCase();
                String forward = token.nextToken();
                // Check that family directory exists
                String forwardPath = kineticsDirectory + forward + separator;
                File forwardDir = new File(forwardPath);
                if (!forwardDir.exists())
                    throw new IOException("Reaction family \"" + forward
                            + "\" specified in families.txt not found!");
                if (status.equals("on")) {
                    ReactionTemplate rt = new ReactionTemplate();
                    rt.read(forward, kineticsDirectory + forward + separator);
                    addReactionTemplate(rt);
                    ReactionTemplate reverse_rt = rt
                            .getReverseReactionTemplate();
                    if (reverse_rt != null)
                        addReactionTemplate(reverse_rt);
                } else {
                    Logger.info("Skipping reaction family: " + forward);
                }
                line = ChemParser.readMeaningfulLine(data, true);
            }
            in.close();
        } catch (IOException e) {
            Logger.critical(e.getMessage());
            System.exit(0);
        }
        Logger.info("");
    }

    // ## operation isEmpty()
    public boolean isEmpty() {
        // #[ operation isEmpty()
        return (size() == 0);
        // #]
    }

    /**
     * Requires: Effects: read in all the defined reaction templates in the pass-in directory Modifies:
     */
    // ## operation read(String)
    public void read(String p_directoryName) {
        // #[ operation read(String)
        File f = new File(p_directoryName);
        String[] fileNames = f.list();
        Arrays.sort(fileNames);
        if (f == null) {
            throw new RuntimeException(String.format(
                    "Empty reaction template directory %s!", p_directoryName));
        }
        for (int i = 0; i < fileNames.length; i++) {
            String fullName = p_directoryName + "/" + fileNames[i];
            File d = new java.io.File(fullName);
            if (d.isDirectory() && !d.getName().toUpperCase().equals("CVS")) {
                ReactionTemplate rt = new ReactionTemplate();
                rt.read(fileNames[i], fullName);
                addReactionTemplate(rt);
                ReactionTemplate reverse_rt = rt.getReverseReactionTemplate();
                if (reverse_rt != null)
                    addReactionTemplate(reverse_rt);
            }
        }
        return;
        // #]
    }

    // ## operation size()
    public int size() {
        // #[ operation size()
        return reactionTemplate.size();
        // #]
    }

    public boolean getForceTheremoConsistence() {
        return ForceTheremoConsistence;
    }

    public void setForceTheremoConsistence(boolean p_ForceTheremoConsistence) {
        ForceTheremoConsistence = p_ForceTheremoConsistence;
    }

    public static ReactionTemplateLibrary getINSTANCE() {
        return INSTANCE;
    }

    public static void setINSTANCE(ReactionTemplateLibrary p_INSTANCE) {
        INSTANCE = p_INSTANCE;
    }

    public Iterator getReactionTemplate() {
        Iterator iter = reactionTemplate.iterator();
        return iter;
    }

    public void addReactionTemplate(ReactionTemplate p_ReactionTemplate) {
        reactionTemplate.add(p_ReactionTemplate);
    }

    public void removeReactionTemplate(ReactionTemplate p_ReactionTemplate) {
        reactionTemplate.remove(p_ReactionTemplate);
    }

    public void clearReactionTemplate() {
        reactionTemplate.clear();
    }

    public int removeFromAllReactionDictionariesByStructure(Structure structure) {
        // remove the given reaction structure from all the ReactionTemplate.reactionDictionaryByStructure sets
        // Returns the number of times it was found and removed.
        int number_of_times_removed = 0;
        Iterator<ReactionTemplate> iterRT = getReactionTemplate();
        while (iterRT.hasNext()) {
            ReactionTemplate rt = (ReactionTemplate) iterRT.next();
            if (rt.removeFromReactionDictionaryByStructure(structure)) {
                number_of_times_removed += 1;
                // If reaction structures were in only one template dictionary, we could break now
                // but they may be in more than one so we need to check them all.
            }
        }
        return number_of_times_removed;
    }
}
/*********************************************************************
 * File Path : RMG\RMG\jing\rxn\ReactionTemplateLibrary.java
 *********************************************************************/

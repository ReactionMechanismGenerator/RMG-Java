////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////



package jing.rxnSys;


import java.io.*;

import jing.mathTool.UncertainDouble;
import jing.rxn.*;
import jing.chem.*;
import java.util.*;

import jing.chemUtil.*;
import jing.chemParser.*;



/**
This is the primary kinetics set that any reaction system has to include into its model.  For example, in combustion system, we build a primary small molecule reaction set, and every combustion/oxydation system should include such a primary kinetic library.  The reaction / rates are basically from Leeds methane oxidation mechanism.
*/
public class PrimaryKineticLibrary {
    
    protected String name; 
    
    protected static LinkedHashSet reactionSet = new LinkedHashSet(); 
    
    protected HashMap speciesSet = new HashMap(); 
    
    
    // Constructors
    
    public  PrimaryKineticLibrary(String p_libraryName, String p_directoryPath) throws IOException {
        name = p_libraryName;

        if ( p_directoryPath == null) throw new NullPointerException("PrimaryKineticLibrary directory path");
        try {
        	read(p_directoryPath,p_libraryName);
        }
        catch (IOException e) {
        	throw new IOException("Error reading Primary Kinetic Library: " + name + '\n' + e.getMessage());
        }
    }
    public  PrimaryKineticLibrary() {
    }
    
    public void appendPrimaryKineticLibrary(String new_p_libraryName, String p_directoryPath) throws IOException {
    	// Appends the current PKLib with an additional one, allowing the user
    	// to combine separate PKLibs easily.  GJB 10/05.
     	setName(name+"/"+new_p_libraryName);
    	try {
    		read(p_directoryPath,new_p_libraryName);	
    	}
        catch (IOException e) {
        	throw new IOException("Error reading Primary Kinetic Library: " + new_p_libraryName + '\n' + e.getMessage());
        }
    }
    
    
    
    public LinkedHashSet getSpeciesSet() {
        
        return new LinkedHashSet(speciesSet.values());
        
    }
    
    public void read(String p_directoryName, String p_name) throws IOException {
        try {
        	if (!p_directoryName.endsWith("/")) p_directoryName = p_directoryName + "/";
			System.out.println("Reading Primary Kinetic Library from: "+p_directoryName);
        	
            String speciesFile = p_directoryName + "species.txt";
            String reactionFile = p_directoryName + "reactions.txt";
            String pdepreactionFile = p_directoryName + "pdepreactions.txt";
            
            SeedMechanism sm = new SeedMechanism();
        	speciesSet.putAll(sm.readSpecies(speciesFile,p_name,"Primary Kinetic Library: "));
        	reactionSet.addAll(sm.readReactions(reactionFile,p_name,speciesSet,"Primary Kinetic Library: ",true));
        	reactionSet.addAll(sm.readPdepReactions(pdepreactionFile,p_name,speciesSet,"Primary Kinetic Library: ",true));
        	return;
        }
        catch (Exception e) {
        	throw new IOException("Can't read primary kinetic library.\n" + e.getMessage());
        }
    }
    
    public static int size() {
        return reactionSet.size();
    }
    
    public String getName() {
        return name;
    }
    
    public void setName(String p_name) {
        name = p_name;
    }
    
    public static LinkedHashSet getReactionSet() {
        return reactionSet;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\PrimaryKineticLibrary.java
*********************************************************************/


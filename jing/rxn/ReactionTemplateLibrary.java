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



package jing.rxn;


import java.io.*;
import java.util.*;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ReactionTemplateLibrary.java                                                                  
//----------------------------------------------------------------------------

//## class ReactionTemplateLibrary 
public class ReactionTemplateLibrary {
    
    protected boolean ForceTheremoConsistence = false;		//## attribute ForceTheremoConsistence 
    
    protected static ReactionTemplateLibrary INSTANCE = new ReactionTemplateLibrary();		//## attribute INSTANCE 
    
    protected LinkedHashSet reactionTemplate;
    
    // Constructors
    
    //## operation ReactionTemplateLibrary() 
    private  ReactionTemplateLibrary() {
        {
            reactionTemplate=new LinkedHashSet();
        }
        //#[ operation ReactionTemplateLibrary() 
        String kineticsDirectory = System.getProperty("jing.rxn.ReactionTemplateLibrary.pathName");
        if (kineticsDirectory == null) {
        	System.out.println("undefined system property: jing.rxn.ReactionTemplateLibrary.pathName, exit!");
        	System.exit(0);
        }
        
        String separator = System.getProperty("file.separator");
        if (!kineticsDirectory.endsWith(separator)) kineticsDirectory = kineticsDirectory + separator;
	      System.out.println("\nReading kinetics database from "+kineticsDirectory);        
        read(kineticsDirectory);
        //#]
    }
    
    //## operation isEmpty() 
    public boolean isEmpty() {
        //#[ operation isEmpty() 
        return (size() == 0);
        //#]
    }
    
    /**
    Requires:
    Effects: read in all the defined reaction templates in the pass-in directory
    Modifies:
    */
    //## operation read(String) 
    public void read(String p_directoryName) {
        //#[ operation read(String) 
        File f = new File(p_directoryName);
        String[] fileNames = f.list();
        Arrays.sort(fileNames);
        if (f==null) {
        	System.err.println("Empty reaction template directory!");
        	System.exit(0);
        }
        for (int i=0; i<fileNames.length;i++) {
        	String fullName = p_directoryName + "/" + fileNames[i];
        	File d = new java.io.File(fullName);
        	if (d.isDirectory() && !d.getName().toUpperCase().equals("CVS")) {
        		ReactionTemplate rt = new ReactionTemplate();
        		rt.read(fileNames[i],fullName);
        		addReactionTemplate(rt);
        		ReactionTemplate reverse_rt = rt.getReverseReactionTemplate();
        		if (reverse_rt != null) addReactionTemplate(reverse_rt);
        	}
        }
        return;
        
        
        //#]
    }
    
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return reactionTemplate.size();
        //#]
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
        Iterator iter=reactionTemplate.iterator();
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
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ReactionTemplateLibrary.java
*********************************************************************/


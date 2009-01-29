//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator                                            
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!     
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!  
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED 
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR  
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;  
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//! 
//!******************************************************************************



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


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


import jing.chem.*;
import java.util.*;
import jing.mathTool.*;
import jing.chem.Species;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\PDepWell.java                                                                  
//----------------------------------------------------------------------------

//## class PDepWell 
public class PDepWell {
    
    protected Species isomer;		//## attribute isomer 
    
    protected HashSet paths;
    
    // Constructors
    
    //## operation PDepWell(Species) 
    protected  PDepWell(Species p_species) {
        {
            paths=new HashSet();
        }
        
        isomer = p_species;
        paths = new HashSet();
        boolean added = false;
        ReactionTemplateLibrary rtl = ReactionTemplateLibrary.getINSTANCE();
        for (Iterator iter = rtl.getReactionTemplate(); iter.hasNext(); ) {
        	ReactionTemplate rt = (ReactionTemplate)iter.next();
        	if (rt.ableToGeneratePDepWellPaths()) {
        		HashSet rxnSet = rt.reactOneReactant(p_species);
        		for (Iterator pathIter = rxnSet.iterator();pathIter.hasNext();) {
        			TemplateReaction tr = (TemplateReaction)pathIter.next();
        			PDepPathReaction pdpr = new PDepPathReaction(tr);
        			paths.add(pdpr);
        			added = true;
        		}
        	}
        }
        if (!added) {
        	System.out.println("no path for well found!");
        	System.out.println(p_species);
        	System.exit(0);
        }
        return;
        //#]
    }
    
    
    
    protected  PDepWell(Species p_species, Reaction tr) {
    	{
            paths=new HashSet();
        }
        
        isomer = p_species;
        PDepPathReaction pdpr = new PDepPathReaction(tr);
        paths.add(pdpr);
        return;
    }
    
    protected PDepWell(Species p_species, HashSet tsSet){
    	{
    		paths = new HashSet();
    	}
    	
    	isomer = p_species;
    	Iterator iter = tsSet.iterator();
    	while (iter.hasNext()){
    		Reaction tr = (Reaction) iter.next();
    		PDepPathReaction pdpr = new PDepPathReaction(tr);
            paths.add(pdpr);
    	}
    }
    
    protected void addPath(Reaction tr){
    	PDepPathReaction pdpr = new PDepPathReaction(tr);
		paths.add(pdpr);
    }
    
    public  PDepWell() {
        {
            paths=new HashSet();
        }
    }
    
    //## operation toChemDisString() 
    public String toChemDisString() {
        //#[ operation toChemDisString() 
        Species spe = getIsomer();
        String s = spe.toChemDisString() + '\n';
        s += "FREQ\n"; 
        ThreeFrequencyModel tfm = spe.getThreeFrequencyModel();
        int fnumber = tfm.getFrequencyNumber();
        s += " " + String.valueOf(fnumber) + " ";
        double [] freq = tfm.getFrequency();
        double [] deg = tfm.getDegeneracy();
        for (int i=0; i<fnumber; i++) {
        	double f = freq[i];
        	double d = deg[i];
        	s += MathTool.formatDouble(f, 7, 1) + " ";
        	s += MathTool.formatDouble(d, 7, 1) + " ";
        }
        s += '\n';
        for (Iterator iter = getPaths(); iter.hasNext(); ) {
        	PDepPathReaction pdpr = (PDepPathReaction)iter.next();
        	s += pdpr.toChemDisString();
        }
        return s;  
        
        	
        
        
        
        
        
        
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        return toChemDisString();
        //#]
    }
    
    public Species getIsomer() {
        return isomer;
    }
    
    public Iterator getPaths() {
        Iterator iter=paths.iterator();
        return iter;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\PDepWell.java
*********************************************************************/


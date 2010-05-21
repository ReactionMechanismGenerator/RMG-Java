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

package jing.chem;

public class TransportData {
    
    /**
    N2's value as default
    */
    protected double epsilon = 97.5;		//## attribute epsilon 
    
    /**
    N2's value as default
    */
    protected double sigma = 3.62;		//## attribute sigma 
    
    /**
     * Assumed values by MRH
     */
    double dipoleMoment = 0.0;
    double polarizability = 0.0;
    double rotrelaxcollnum = 1.0;
    int shapeIndex = 0;
    String comments = "";
    String source = "";
    String name = "";
    LJData gaData;
    boolean fromPrimaryTransportLibrary = false;
    
    // Constructors
    public TransportData(LJData lj_data, boolean truefalse) {
    	epsilon = lj_data.calculateEpsilon();
    	sigma = lj_data.calculateSigma();
    	// No group additivity scheme for dipoleMoment, polarizability, or rotrelaxcollnum
    	shapeIndex = lj_data.shapeIndex;
    	comments = lj_data.name + ")";
    	source = lj_data.toString();
    	gaData = lj_data;
    	fromPrimaryTransportLibrary = truefalse;
    }
        
    public TransportData(double p_sigma, double p_epsilon) {
        sigma = p_sigma;
        epsilon = p_epsilon;
    }
    
    public TransportData(int p_index, double p_epsilon, double p_sigma, double p_dipole, double p_polar, double p_rotrelax, boolean truefalse) {
    	shapeIndex = p_index;
    	epsilon = p_epsilon;
    	sigma = p_sigma;
    	dipoleMoment = p_dipole;
    	polarizability = p_polar;
    	rotrelaxcollnum = p_rotrelax;
    	fromPrimaryTransportLibrary = truefalse;
    }
    
    public TransportData(String p_name, TransportData p_transdata, String p_comments, String p_source, boolean truefalse) {
    	name = p_name;
    	comments = p_comments;
    	source = p_source;
    	sigma = p_transdata.sigma;
    	epsilon = p_transdata.epsilon;
    	dipoleMoment = p_transdata.dipoleMoment;
    	polarizability = p_transdata.polarizability;
    	rotrelaxcollnum = p_transdata.rotrelaxcollnum;
    	shapeIndex = p_transdata.shapeIndex;
    	fromPrimaryTransportLibrary = truefalse;
    }
    
    public TransportData() {
	}

	// Lennard-Jones potential well depth [=] Kelvin
    public double getEpsilon() {
        return epsilon;
    }
    
    // Lennard-Jones collision diameter [=] angstroms
    public double getSigma() {
        return sigma;
    }
    
    // dipole moment [=] Debye
    public double getDipoleMoment() {
    	return dipoleMoment;
    }
    
    // polarizability [=] (angstroms)^3
    public double getPolarizability() {
    	return polarizability;
    }
    
    // rotational relaxation collision number at 298 K
    public double getRotationalRelaxationCollisionNumber() {
    	return rotrelaxcollnum;
    }
    
    // shape index of molecule, as defined by CHEMKIN
    public int getShapeIndex() {
    	// if the molecule comes from a primary transport library, return the value as is
    	if (fromPrimaryTransportLibrary) return shapeIndex;
    	// else, the data was estimated using group additivity
    	if (shapeIndex == 0)
    		return 1;	// linear
    	else
    		return 2;	// nonlinear
    }
    
    public String getName() {
    	return name;
    }
    
    public String toString() {
    	/*
    	 *  The Lennard-Jones sigma/epsilon parameters are estimated by RMG
    	 *  The others are assumed values (MRH 17-May-2010)
    	 */
    	return  getShapeIndex() + " " +
    			String.format("%9.3f", getEpsilon()) + " " +
    			String.format("%9.3f", getSigma()) + " " +    					
    			String.format("%9.3f", getDipoleMoment()) + " " +
    			String.format("%9.3f", getPolarizability()) + " " +
    			String.format("%9.3f", getRotationalRelaxationCollisionNumber());
    }

	public String getSource() {
		return source;
	}
	
	public String getComment() {
		return comments;
	}
	
	public LJData getGAData() {
		return gaData;
	}
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\LennardJones.java
*********************************************************************/


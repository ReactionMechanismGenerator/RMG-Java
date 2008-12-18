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

package jing.chem;

/**
 * Provides functionality for storing and manipulating the spectroscopic data
 * for a given molecule, including the vibrational, external rotation,
 * interal rotation modes, and the symmetry number. This information is used,
 * among other things, to estimate the partition function and density of states 
 * for a molecule.
 * @author jwallen
 */
public class SpectroscopicData {
	
	//==========================================================================
	//
	//	Data members
	//
	
	/**
	 * A flag that specifies whether or not to utilize the three-frequency
	 * model for the spectroscopic data. Setting this to true makes the code
	 * dependent on the THERFIT module; setting this to false makes the code
	 * use Frankie instead.
	 */
	public static boolean useThreeFrequencyModel = false;
	
	/**
	 * An array of vibrational frequencies, each representing one degree of
	 * freedom modeled as a simple harmonic oscillator. Units are cm^-1.
	 */
	private double[] vibrations;
	
	/**
	 * An array of external rotational frequencies, each representing one degree
	 * of freedom modeled as a rigid rotor. This array should contain either
	 * one frequency (linear molecule) or three frequencies (nonlinear 
	 * molecule). Units are cm^-1.
	 */
	private double[] rotations;
	
	/**
	 * The symmetry number of the molecule, used to correct for multiple
	 * counting of rotational states caused by rotational symmetry.
	 */
	private int symmetry;
	
	/**
	 * An array of internal rotational frequencies, each representing one degree
	 * of freedom modeled as a hindered rotor. Units are cm^-1.
	 */
	private double[] hinderedFrequencies;
	
	/**
	 * An array of internal rotational barrier energies, each corresponding to
	 * one of the hindered frequencies above. Units are cm^-1.
	 */
	private double[] hinderedBarriers;
	
	//==========================================================================
	//
	//	Constructors
	//
	
	/**
	 * Default constructor. Initializes all data to empty.
	 */
	public SpectroscopicData() {
		clear();
	}
        
        /**
	 * gmagoon 12/11/08: alternative constructor
	 */
	public SpectroscopicData(double[] p_vibrations, double[] p_hinderedFrequencies, double[] p_hinderedBarriers) {
                vibrations = p_vibrations;
		rotations = null;
		symmetry = 1;
		hinderedFrequencies = p_hinderedFrequencies;
		hinderedBarriers = p_hinderedBarriers;
	}
	
	//==========================================================================
	//
	//	Accessors
	//
	
	/**
	 * Returns the number of vibrational frequencies.
	 * @return The number of vibrational frequencies.
	 */
	public int getVibrationCount() {
		if (vibrations == null)
			return 0;
		else
			return vibrations.length;
	}
	
	/**
	 * Returns the vibrational frequency at the specified index in the array.
	 * @return The desired vibrational frequency.
	 */
	public double getVibration(int index) {
		return vibrations[index];
	}
	
	/**
	 * Returns the vibrational frequencies.
	 * @return The vibrational frequencies.
	 */
	public double[] getVibrations() {
		return vibrations;
	}
	
	/**
	 * Returns the number of (external) rotational frequencies.
	 * @return The number of rotational frequencies.
	 */
	public int getRotationCount() {
		if (rotations == null)
			return 0;
		else
			return rotations.length;
	}
	
	/**
	 * Returns the (external) rotational frequency at the specified index in
	 * the array.
	 * @return The desired rotational frequency.
	 */
	public double getRotation(int index) {
		return rotations[index];
	}
	
	/**
	 * Returns the (external) rotational frequencies.
	 * @return The rotational frequencies.
	 */
	public double[] getRotations() {
		return rotations;
	}
	
	/**
	 * Returns the symmetry number.
	 * @return The symmetry number.
	 */
	public int getSymmetryNumber() {
		return symmetry;
	}
	
	/**
	 * Returns the number of hindered rotor frequency-barrier pairs.
	 * @return The number of hindered rotor frequency-barrier pairs.
	 */
	public int getHinderedCount() {
		if (hinderedFrequencies == null)
			return 0;
		else
			return hinderedFrequencies.length;
	}
	
	/**
	 * Returns the (external) rotational frequency at the specified index in
	 * the array.
	 * @return The desired rotational frequency.
	 */
	public double getHinderedFrequency(int index) {
		return hinderedFrequencies[index];
	}
	
	/**
	 * Returns the (external) rotational frequencies.
	 * @return The rotational frequencies.
	 */
	public double[] getHinderedFrequencies() {
		return hinderedFrequencies;
	}
	
	/**
	 * Returns the (external) rotational frequency at the specified index in
	 * the array.
	 * @return The desired rotational frequency.
	 */
	public double getHinderedBarrier(int index) {
		return hinderedBarriers[index];
	}
	
	/**
	 * Returns the (external) rotational frequencies.
	 * @return The rotational frequencies.
	 */
	public double[] getHinderedBarriers() {
		return hinderedBarriers;
	}
	
	//==========================================================================
	//
	//	Other methods
	//
	
	/**
	 * Clears any existing data from the object.
	 */
	public void clear() {
		vibrations = null;
		rotations = null;
		symmetry = 1;
		hinderedFrequencies = null;
		hinderedBarriers = null;
	}
	
	/**
	 * Converts a three frequency model (i.e. from a THERFIT calculation) to
	 * a list of vibrational frequencies and stores them in the current object.
	 * @param model The three frequency model to take the data from.
	 */
	public void loadThreeFrequencyModel(ThreeFrequencyModel model) {
		
		clear();
		
		// Get data from three frequency model
		double[] freq = model.getFrequency();
		double[] deg = model.getDegeneracy();
		
		// Determine total number of degrees of freedom represented
		int nVib = 0;
		for (int i = 0; i < model.getFrequencyNumber(); i++)
			nVib += (int) deg[i];
		
		// Transfer frequencies
		vibrations = new double[nVib];
		int index = 0;
		for (int i = 0; i < model.getFrequencyNumber(); i++) {
			for (int j = 0; j < deg[i]; j++) {
				vibrations[index] = freq[i];
				index++;
			}
		}
	}
	
}

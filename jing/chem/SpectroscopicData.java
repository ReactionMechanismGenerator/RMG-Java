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
	
	public enum Mode { OFF, THREEFREQUENCY, FREQUENCYGROUPS };

	/**
	 * A flag that specifies the method to use to estimate the spectroscopic
	 * data. Setting this to THREEFREQUENCY makes the code dependent on the
	 * THERFIT module; setting this to FREQUENCYGROUPS makes the code use
	 * Frankie instead. If spectroscopic data is not needed (i.e. non-pressure-
	 * dependent cases) then OFF is used.
	 */
	public static Mode mode = Mode.OFF;

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
	
	/**
	 * Default constructor. Initializes all data to empty.
	 */
	public SpectroscopicData(ThreeFrequencyModel tfm) {
		loadThreeFrequencyModel(tfm);
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

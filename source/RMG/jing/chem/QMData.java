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

/**
 * @author User1
 */
public class QMData {
    public int natoms; // number of atoms from Mopac file; in principle, this should agree with number of chemGraph
// atoms
    public ArrayList atomicNumber; // vector of atomic numbers (integers) (apparently Vector is thread-safe; cf.
// http://answers.yahoo.com/question/index?qid=20081214065127AArZDT3; ...should I be using this instead?)
    public ArrayList x_coor; // vectors of x-, y-, and z-coordinates (doubles) (Angstroms) (in order corresponding to
// above atomic numbers)
    public ArrayList y_coor;
    public ArrayList z_coor;
    public double energy; // energy (Hf298) in Hartree
    public double stericEnergy;// steric energy in Hartree
    public double molmass; // molecular mass in amu
    public ArrayList freqs; // list of frequencies in units of cm^-1
    public double rotCons_1;// rotational constants in (1/s)
    public double rotCons_2;
    public double rotCons_3;

    // constructor class
    public QMData(int p_natoms, ArrayList p_atomicNumber, ArrayList p_x_coor,
            ArrayList p_y_coor, ArrayList p_z_coor, double p_energy,
            double p_stericEnergy, double p_molmass, ArrayList p_freqs,
            double p_rotCons_1, double p_rotCons_2, double p_rotCons_3) {
        natoms = p_natoms;
        atomicNumber = p_atomicNumber;
        x_coor = p_x_coor;
        y_coor = p_y_coor;
        z_coor = p_z_coor;
        energy = p_energy;
        stericEnergy = p_stericEnergy;
        molmass = p_molmass;
        freqs = p_freqs;
        rotCons_1 = p_rotCons_1;
        rotCons_2 = p_rotCons_2;
        rotCons_3 = p_rotCons_3;
    }

    public String getSYMMETRYinput() {
        String geom = natoms + "\n";
        for (int i = 0; i < natoms; i++) {
            geom += atomicNumber.get(i) + " " + x_coor.get(i) + " "
                    + y_coor.get(i) + " " + z_coor.get(i) + "\n";
        }
        return geom;
    }
}

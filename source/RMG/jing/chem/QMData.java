/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jing.chem;

import java.util.*;
/**
 *
 * @author User1
 */
public class QMData {
    public int natoms; //number of atoms from Mopac file; in principle, this should agree with number of chemGraph atoms
    public ArrayList atomicNumber; //vector of atomic numbers (integers) (apparently Vector is thread-safe; cf. http://answers.yahoo.com/question/index?qid=20081214065127AArZDT3; ...should I be using this instead?)
    public ArrayList x_coor; //vectors of x-, y-, and z-coordinates (doubles) (Angstroms) (in order corresponding to above atomic numbers)
    public ArrayList y_coor;
    public ArrayList z_coor;
    public double energy; // energy (Hf298) in Hartree
    public double stericEnergy;//steric energy in Hartree
    public double molmass; //molecular mass in amu
    public ArrayList freqs; //list of frequencies in units of cm^-1
    public double rotCons_1;//rotational constants in (1/s)
    public double rotCons_2;
    public double rotCons_3;

    //constructor class
    public QMData(int p_natoms, ArrayList p_atomicNumber, ArrayList p_x_coor, ArrayList p_y_coor, ArrayList p_z_coor, double p_energy, double p_stericEnergy, double p_molmass, ArrayList p_freqs, double p_rotCons_1, double p_rotCons_2, double p_rotCons_3){
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

    public String getSYMMETRYinput(){
	String geom = natoms + "\n";
	for(int i=0; i < natoms; i++){
	    geom += atomicNumber.get(i) + " "+ x_coor.get(i) + " " + y_coor.get(i) + " " +z_coor.get(i) + "\n";
	}
	return geom;
    }
}


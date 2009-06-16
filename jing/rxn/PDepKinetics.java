/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jing.rxn;

import jing.param.Pressure;
import jing.param.Temperature;

/**
 * An interface for representing pressure-dependent kinetics k(T, P).
 * 
 * @author jwallen
 */
public interface PDepKinetics {

	public double calculateRate(Temperature p_temperature, Pressure p_pressure);

}

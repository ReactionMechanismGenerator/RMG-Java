package jing.chem;

/**
 * The abstract class for every strategy to estimate ideal gas thermochemistry for species
 * @author nmvdewie
 *
 */
public abstract class TDGenerator {
	
	GeneralGAPP thermoGAPP;
	
	GeneralGAPP thermoQM;

	public abstract ThermoData generateThermo(ChemGraph chemGraph) ;
}

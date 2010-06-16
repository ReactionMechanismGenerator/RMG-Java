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
import jing.chem.*;
import java.util.*;

import jing.param.*;
import jing.mathTool.*;
import jing.chemParser.*;
import jing.chem.Species;
import jing.param.Temperature;
import jing.rxnSys.NegativeConcentrationException;
import jing.rxnSys.ReactionModelGenerator;
import jing.rxnSys.SystemSnapshot;

//## package jing::rxn

//----------------------------------------------------------------------------
//jing\rxn\Reaction.java
//----------------------------------------------------------------------------

/**
Immutable objects.
*/
//## class Reaction
public class Reaction {

  protected static double TRIMOLECULAR_RATE_UPPER = 1.0E100;
  protected static double BIMOLECULAR_RATE_UPPER = 1.0E100;		//## attribute BIMOLECULAR_RATE_UPPER

  protected static double UNIMOLECULAR_RATE_UPPER = 1.0E100;		//## attribute UNIMOLECULAR_RATE_UPPER

  protected String comments = "No comment";		//## attribute comments

  protected Kinetics[] fittedReverseKinetics = null;		//## attribute fittedReverseKinetics

  protected double rateConstant; 
  
  protected Reaction reverseReaction = null;		//## attribute reverseReaction

  protected Kinetics[] kinetics;
  protected Structure structure;
  protected double UpperBoundRate;//svp
  protected double LowerBoundRate;//svp
  //protected Kinetics additionalKinetics = null;  //This is incase a reaction has two completely different transition states.
  protected boolean finalized = false;
  protected String ChemkinString = null;
  protected boolean ratesForKineticsAndAdditionalKineticsCross = false; //10/29/07 gmagoon: added variable to keep track of whether both rate constants are maximum for some temperature in the temperature range
  
  protected boolean kineticsFromPrimaryKineticLibrary = false;
  protected ReactionTemplate rxnTemplate;
  // Constructors

  //## operation Reaction()
  public  Reaction() {
      //#[ operation Reaction()
      //#]
  }
  //## operation Reaction(Structure,RateConstant)
  private  Reaction(Structure p_structure, Kinetics[] p_kinetics) {
      //#[ operation Reaction(Structure,RateConstant)
      structure = p_structure;
      kinetics = p_kinetics;
	  //rateConstant = calculateTotalRate(Global.temperature);
      //#]
  }
  /*public Reaction(Reaction rxn) {
		structure = rxn.structure;
		kinetics = rxn.kinetics;

		comments = rxn.comments;

		fittedReverseKinetics = rxn.fittedReverseKinetics;

		rateConstant = rxn.rateConstant;

		reverseReaction = rxn.reverseReaction;

		UpperBoundRate = rxn.UpperBoundRate;
		LowerBoundRate = rxn.LowerBoundRate;
		additionalKinetics = rxn.additionalKinetics;
		finalized = rxn.finalized;
		ChemkinString  = rxn.ChemkinString;
		ratesForKineticsAndAdditionalKineticsCross = rxn.ratesForKineticsAndAdditionalKineticsCross;

  }*/
  
  
  //## operation allProductsIncluded(HashSet)
  public boolean allProductsIncluded(HashSet p_speciesSet) {
      //#[ operation allProductsIncluded(HashSet)
      Iterator iter = getProducts();
      while (iter.hasNext()) {
      	Species spe = ((Species)iter.next());
      	if (!p_speciesSet.contains(spe)) return false;
      }
      return true;
      //#]
  }

  //## operation allReactantsIncluded(HashSet)
  public boolean allReactantsIncluded(HashSet p_speciesSet) {
      //#[ operation allReactantsIncluded(HashSet)
      if (p_speciesSet == null) throw new NullPointerException();
      Iterator iter = getReactants();
      while (iter.hasNext()) {
      	Species spe = ((Species)iter.next());
      	if (!p_speciesSet.contains(spe)) return false;
      }
      return true;
      //#]
  }

  /**
  Calculate this reaction's thermo parameter.  Basically, make addition of the thermo parameters of all the reactants and products.
  */
  //## operation calculateHrxn(Temperature)
  public double calculateHrxn(Temperature p_temperature) {
      //#[ operation calculateHrxn(Temperature)
      return structure.calculateHrxn(p_temperature);
      //#]
  }

  //## operation calculateKeq(Temperature)
  public double calculateKeq(Temperature p_temperature) {
      //#[ operation calculateKeq(Temperature)
      return structure.calculateKeq(p_temperature);
      //#]
  }

  //## operation calculateKeqUpperBound(Temperature)
  //svp
    public double calculateKeqUpperBound(Temperature p_temperature) {
      //#[ operation calculateKeqUpperBound(Temperature)
      return structure.calculateKeqUpperBound(p_temperature);
      //#]
    }

  //## operation calculateKeqLowerBound(Temperature)
  //svp
    public double calculateKeqLowerBound(Temperature p_temperature) {
      //#[ operation calculateKeqLowerBound(Temperature)
      return structure.calculateKeqLowerBound(p_temperature);
      //#]
    }

  
  public double calculateTotalRate(Temperature p_temperature){
  	double rate =0;
	Temperature stdtemp = new Temperature(298,"K");
	double Hrxn = calculateHrxn(stdtemp);
	/*
	 * 29Jun2009-MRH: Added a kinetics from PRL check
	 * 	If the kinetics for this reaction is from a PRL, use those numbers
	 * 		to compute the rate.  Else, proceed as before.
	 */
	
	/*
	 * MRH 18MAR2010:
	 * Changing the structure of a reaction's kinetics
	 * 	If the kinetics are from a primary kinetic library, we assume the user
	 * 		has supplied the total pre-exponential factor for the reaction (and
	 * 		not the per-event pre-exponential facor).
	 *  If the kinetics were estimated by RMG, the pre-exponential factor must
	 *  	be multiplied by the "redundancy" (# of events)
	 */
	if (kineticsFromPrimaryKineticLibrary) {
		Kinetics[] k_All = kinetics;
		for (int numKinetics=0; numKinetics<kinetics.length; numKinetics++) {
			Kinetics k = k_All[numKinetics];
			if (k instanceof ArrheniusEPKinetics)
				rate += k.calculateRate(p_temperature,Hrxn);
			else
				rate += k.calculateRate(p_temperature);
		}
		return rate;
	}
	else if (isForward()){
		Kinetics[] k_All = kinetics;
		for (int numKinetics=0; numKinetics<kinetics.length; numKinetics++) {
			Kinetics k = k_All[numKinetics].multiply(structure.redundancy);
			if (k instanceof ArrheniusEPKinetics)
				rate += k.calculateRate(p_temperature,Hrxn);
			else
				rate += k.calculateRate(p_temperature);
		}
		return rate;
//  		Iterator kineticIter = getAllKinetics().iterator();
//      	while (kineticIter.hasNext()){
//      		Kinetics k = (Kinetics)kineticIter.next();
//			if (k instanceof ArrheniusEPKinetics)
//				rate = rate + k.calculateRate(p_temperature,Hrxn);
//			else 
//				rate = rate + k.calculateRate(p_temperature);
//      	}
//		
//      	return rate;
  	}
  	else if (isBackward()){
  		Reaction r = getReverseReaction();
  		rate = r.calculateTotalRate(p_temperature);
  		return rate*calculateKeq(p_temperature);
  	}
  	else {
  		throw new InvalidReactionDirectionException();
  	}
  	
  }

  //## operation calculateUpperBoundRate(Temperature)
  //svp
    public double calculateUpperBoundRate(Temperature p_temperature){
      //#[ operation calculateUpperBoundRate(Temperature)
      if (isForward()){
        double A;
        double E;
        double n;
        
        for (int numKinetics=0; numKinetics<kinetics.length; numKinetics++) {
	        A = kinetics[numKinetics].getA().getUpperBound();
	        E = kinetics[numKinetics].getE().getLowerBound();      
	        n = kinetics[numKinetics].getN().getUpperBound();
	              
	        if (A > 1E300) {
	        	A = kinetics[numKinetics].getA().getValue()*1.2;
	        }
	        //Kinetics kinetics = getRateConstant().getKinetics();
	        if (kinetics[numKinetics] instanceof ArrheniusEPKinetics){
	        	ArrheniusEPKinetics arrhenius = (ArrheniusEPKinetics)kinetics[numKinetics];
	        	double H = calculateHrxn(p_temperature);
	        	if (H < 0) {
	        		if (arrhenius.getAlpha().getValue() > 0){
	        			double alpha = arrhenius.getAlpha().getUpperBound();
	        			E = E + alpha * H;
	        		}
	        		else{
	        			double alpha = arrhenius.getAlpha().getLowerBound();
	        			E = E + alpha*H;
	        		}
	        	}
	        	else {
	        		if (arrhenius.getAlpha().getValue() > 0){
	        			double alpha = arrhenius.getAlpha().getLowerBound();
	        			E = E + alpha * H;
	        		}
	        		else{
	        			double alpha = arrhenius.getAlpha().getUpperBound();
	        			E = E + alpha*H;
	        		}
	        	}
	        }
	        if (E < 0){
	        	E = 0;
	        }
	        double indiv_k = 0.0;
	        indiv_k = A*Math.pow(p_temperature.getK(),n)*Math.exp(-E/GasConstant.getKcalMolK()/p_temperature.getK());
	        indiv_k *= getStructure().getRedundancy();
	        UpperBoundRate += indiv_k;
        }
        return UpperBoundRate;

      }
      else if (isBackward()) {
        Reaction r = getReverseReaction();
        if (r == null) throw new NullPointerException("Reverse reaction is null.\n" + structure.toString());
        if (!r.isForward()) throw new InvalidReactionDirectionException();
        
        
        for (int numKinetics=0; numKinetics<kinetics.length; numKinetics++) {
	        double A = kinetics[numKinetics].getA().getUpperBound();
	        double E = kinetics[numKinetics].getE().getLowerBound();
	        double n = kinetics[numKinetics].getN().getUpperBound();
	        if (A > 1E300) {
	        	A = kinetics[numKinetics].getA().getValue()*1.2;
	        }
	        //Kinetics kinetics = getRateConstant().getKinetics();
	        if (kinetics[numKinetics] instanceof ArrheniusEPKinetics){
	        	ArrheniusEPKinetics arrhenius = (ArrheniusEPKinetics)kinetics[numKinetics];
	        	double H = calculateHrxn(p_temperature);
	        	if (H < 0) {
	        		if (arrhenius.getAlpha().getValue() > 0){
	        			double alpha = arrhenius.getAlpha().getUpperBound();
	        			E = E + alpha * H;
	        		}
	        		else{
	        			double alpha = arrhenius.getAlpha().getLowerBound();
	        			E = E + alpha*H;
	        		}
	        	}
	        	else {
	        		if (arrhenius.getAlpha().getValue() > 0){
	        			double alpha = arrhenius.getAlpha().getLowerBound();
	        			E = E + alpha * H;
	        		}
	        		else{
	        			double alpha = arrhenius.getAlpha().getUpperBound();
	        			E = E + alpha*H;
	        		}
	        	}
	        }
	        if (E < 0){
	        	E = 0;
	        }
	        double indiv_k = 0.0;
	        indiv_k = A*Math.pow(p_temperature.getK(),n)*Math.exp(-E/GasConstant.getKcalMolK()/p_temperature.getK());
	        indiv_k *= getStructure().getRedundancy();
	        UpperBoundRate += indiv_k*calculateKeqUpperBound(p_temperature);
        }
        return UpperBoundRate;
      }
      else{
        throw new InvalidReactionDirectionException();
      }
      //#]
    }

  //## operation calculateLowerBoundRate(Temperature)
  //svp
    public double calculateLowerBoundRate(Temperature p_temperature){
      //#[ operation calculateLowerBoundRate(Temperature)
      if (isForward()){
    	  for (int numKinetics=0; numKinetics<kinetics.length; ++numKinetics) {
	        double A = kinetics[numKinetics].getA().getLowerBound();
	        double E = kinetics[numKinetics].getE().getUpperBound();
	        double n = kinetics[numKinetics].getN().getLowerBound();
	        if (A > 1E300 || A <= 0) {
	          A = kinetics[numKinetics].getA().getValue()/1.2;
	        }
	        //Kinetics kinetics = getRateConstant().getKinetics();
	        if (kinetics[numKinetics] instanceof ArrheniusEPKinetics){
	          ArrheniusEPKinetics arrhenius = (ArrheniusEPKinetics)kinetics[numKinetics];
	          double H = calculateHrxn(p_temperature);
	          if (H < 0) {
	            if (arrhenius.getAlpha().getValue()>0){
	              double alpha = arrhenius.getAlpha().getLowerBound();
	              E = E + alpha * H;
	            }
	            else{
	              double alpha = arrhenius.getAlpha().getUpperBound();
	              E = E + alpha*H;
	            }
	          }
	          else {
	            if (arrhenius.getAlpha().getValue()>0){
	              double alpha = arrhenius.getAlpha().getUpperBound();
	              E = E + alpha * H;
	            }
	            else{
	              double alpha = arrhenius.getAlpha().getLowerBound();
	              E = E + alpha*H;
	            }
	          }
	        }

	        double indiv_k = 0.0;
	        indiv_k = A*Math.pow(p_temperature.getK(),n)*Math.exp(-E/GasConstant.getKcalMolK()/p_temperature.getK());
	        indiv_k *= getStructure().getRedundancy();
	        LowerBoundRate += indiv_k;
    	  }
        return LowerBoundRate;
      }
      else if (isBackward()) {
        Reaction r = getReverseReaction();
        if (r == null) throw new NullPointerException("Reverse reaction is null.\n" + structure.toString());
        if (!r.isForward()) throw new InvalidReactionDirectionException();
        for (int numKinetics=0; numKinetics<kinetics.length; ++numKinetics) {
	        double A = kinetics[numKinetics].getA().getLowerBound();
	        double E = kinetics[numKinetics].getE().getUpperBound();
	        double n = kinetics[numKinetics].getN().getLowerBound();
	        if (A > 1E300) {
	           A = kinetics[numKinetics].getA().getValue()/1.2;
	        }
	        if (kinetics[numKinetics] instanceof ArrheniusEPKinetics){
	          ArrheniusEPKinetics arrhenius = (ArrheniusEPKinetics)kinetics[numKinetics];
	          double H = calculateHrxn(p_temperature);
	          if (H < 0) {
	            if (arrhenius.getAlpha().getValue() > 0){
	              double alpha = arrhenius.getAlpha().getLowerBound();
	              E = E + alpha * H;
	            }
	            else{
	              double alpha = arrhenius.getAlpha().getUpperBound();
	              E = E + alpha*H;
	            }
	          }
	          else {
	            if (arrhenius.getAlpha().getValue() > 0){
	              double alpha = arrhenius.getAlpha().getUpperBound();
	              E = E + alpha * H;
	            }
	            else{
	              double alpha = arrhenius.getAlpha().getLowerBound();
	              E = E + alpha*H;
	            }
	          }
	        }
	
	        double indiv_k = 0.0;
	        indiv_k = A*Math.pow(p_temperature.getK(),n)*Math.exp(-E/GasConstant.getKcalMolK()/p_temperature.getK());
	        indiv_k *= getStructure().getRedundancy();
	        LowerBoundRate += indiv_k*calculateKeqLowerBound(p_temperature);
        }
        return LowerBoundRate;
      }

      else{
        throw new InvalidReactionDirectionException();
      }
      //#]
    }


  //## operation calculateSrxn(Temperature)
  public double calculateSrxn(Temperature p_temperature) {
      //#[ operation calculateSrxn(Temperature)
      return structure.calculateSrxn(p_temperature);
      //#]
  }

  //## operation calculateThirdBodyCoefficient(SystemSnapshot)
  public double calculateThirdBodyCoefficient(SystemSnapshot p_presentStatus) {
      //#[ operation calculateThirdBodyCoefficient(SystemSnapshot)
      if (!(this instanceof ThirdBodyReaction)) return 1;
      else {
      	return ((ThirdBodyReaction)this).calculateThirdBodyCoefficient(p_presentStatus);
      }
      //#]
  }

  //## operation checkRateRange()
  public boolean checkRateRange() {
      //#[ operation checkRateRange()
      Temperature t = new Temperature(1500,"K");

      double rate = calculateTotalRate(t);

      if (getReactantNumber() == 2) {
      	if (rate > BIMOLECULAR_RATE_UPPER) return false;
      }
      else if (getReactantNumber() == 1) {
      	if (rate > UNIMOLECULAR_RATE_UPPER) return false;
      }
	  else if (getReactantNumber() == 3) {
		  if (rate > TRIMOLECULAR_RATE_UPPER) return false;
	  }
      else throw new InvalidReactantNumberException();

      return true;
      //#]
  }

  //## operation contains(Species)
  public boolean contains(Species p_species) {
      //#[ operation contains(Species)
      if (containsAsReactant(p_species) || containsAsProduct(p_species)) return true;
      else return false;


      //#]
  }

  //## operation containsAsProduct(Species)
  public boolean containsAsProduct(Species p_species) {
      //#[ operation containsAsProduct(Species)
      Iterator iter = getProducts();
      while (iter.hasNext()) {
      	//ChemGraph cg = (ChemGraph)iter.next();
      	Species spe = (Species)iter.next();
      	if (spe.equals(p_species)) return true;
      }

      return false;
      //#]
  }

  //## operation containsAsReactant(Species)
  public boolean containsAsReactant(Species p_species) {
      //#[ operation containsAsReactant(Species)
      Iterator iter = getReactants();
      while (iter.hasNext()) {
      	//ChemGraph cg = (ChemGraph)iter.next();
      	Species spe = (Species)iter.next();
      	if (spe.equals(p_species)) return true;
      }

      return false;
      //#]
  }

  /**
   * Checks if the structure of the reaction is the same. Does not check the rate constant.
   * Two reactions with the same structure but different rate constants will be equal.
   */
  //## operation equals(Object)
  public boolean equals(Object p_reaction) {
      //#[ operation equals(Object)
      if (this == p_reaction) return true;

      if (!(p_reaction instanceof Reaction)) return false;

      Reaction r = (Reaction)p_reaction;

      if (!getStructure().equals(r.getStructure())) return false;

      return true;


      //#]
  }

/*  
 // fitReverseKineticsPrecisely and getPreciseReverseKinetics 
 // are not used, not maintained, and we have clue what they do,
 // so we're commenting them out so we don't keep looking at them.
 // (oh, and they look pretty similar to each other!)
 //            - RWest & JWAllen, June 2009
 
  //## operation fitReverseKineticsPrecisely()
  public void fitReverseKineticsPrecisely() {
      //#[ operation fitReverseKineticsPrecisely()
      if (isForward()) {
      	fittedReverseKinetics = null;
      }
      else {

      	String result = "";
      	for (double t = 300.0; t<1500.0; t+=50.0) {
      		double rate = calculateTotalRate(new Temperature(t,"K"));
      		result += String.valueOf(t) + '\t' + String.valueOf(rate) + '\n';
      	}

          // run fit3p
      	String dir = System.getProperty("RMG.workingDirectory");

      	File fit3p_input;

      	try {
      	        // prepare fit3p input file, "input.dat" is the input file name
      	        fit3p_input = new File("fit3p/input.dat");
      	        FileWriter fw = new FileWriter(fit3p_input);
      	        fw.write(result);
      	        fw.close();
      	}
      	catch (IOException e) {
      	        System.out.println("Wrong input file for fit3p!");
      	        System.out.println(e.getMessage());
      	        System.exit(0);
      	}

      	try {
      	    // system call for fit3p
      		String[] command = {dir+ "/bin/fit3pbnd.exe"};
      		File runningDir = new File("fit3p");
      	    Process fit = Runtime.getRuntime().exec(command, null, runningDir);
      	    int exitValue = fit.waitFor();
      	}
      	catch (Exception e) {
      	    System.out.println("Error in run fit3p!");
      	    System.out.println(e.getMessage());
      	    System.exit(0);
      	}

      	// parse the output file from chemdis
      	try {
      	        String fit3p_output = "fit3p/output.dat";

      	        FileReader in = new FileReader(fit3p_output);
      	        BufferedReader data = new BufferedReader(in);

      	        String line = ChemParser.readMeaningfulLine(data);
      	        line = line.trim();
      	    	StringTokenizer st = new StringTokenizer(line);
      	        String A = st.nextToken();
      	        String temp = st.nextToken();
      	        temp = st.nextToken();
      	        temp = st.nextToken();
      	        double Ar = Double.parseDouble(temp);

      			line = ChemParser.readMeaningfulLine(data);
      	        line = line.trim();
      	        st = new StringTokenizer(line);
      	        String n = st.nextToken();
      	        temp = st.nextToken();
      	        temp = st.nextToken();
      	        double nr = Double.parseDouble(temp);

      			line = ChemParser.readMeaningfulLine(data);
      	        line = line.trim();
      	        st = new StringTokenizer(line);
      	        String E = st.nextToken();
      	        temp = st.nextToken();
      	        temp = st.nextToken();
      	        temp = st.nextToken();
      	        double Er = Double.parseDouble(temp);
      			if (Er < 0) {
      				System.err.println(getStructure().toString());
      				System.err.println("fitted Er < 0: "+Double.toString(Er));
      				double increase = Math.exp(-Er/GasConstant.getKcalMolK()/715.0);
      				double deltan = Math.log(increase)/Math.log(715.0);
      				System.err.println("n enlarged by factor of: " + Double.toString(deltan));
      				nr += deltan;
      				Er = 0;

      			}


      	    	UncertainDouble udAr = new UncertainDouble(Ar, 0, "Adder");
      	        UncertainDouble udnr = new UncertainDouble(nr, 0, "Adder");
      	        UncertainDouble udEr = new UncertainDouble(Er, 0, "Adder");

      	        fittedReverseKinetics = new ArrheniusKinetics(udAr, udnr , udEr, "300-1500", 1, "fitting from forward and thermal",null);
      	        in.close();
      	}
      	catch (Exception e) {
      	        System.out.println("Error in read output.dat from fit3p!");
      	        System.out.println(e.getMessage());
      	        System.exit(0);
      	}
      }

      return;
      //#]
  }

  //## operation fitReverseKineticsPrecisely()
  public Kinetics getPreciseReverseKinetics() {
      //#[ operation fitReverseKineticsPrecisely()
      
      Kinetics fittedReverseKinetics =null;
      

      	String result = "";
      	for (double t = 300.0; t<1500.0; t+=50.0) {
      		double rate = calculateTotalRate(new Temperature(t,"K"));
      		result += String.valueOf(t) + '\t' + String.valueOf(rate) + '\n';
      	}

          // run fit3p
      	String dir = System.getProperty("RMG.workingDirectory");

      	File fit3p_input;

      	try {
      	        // prepare fit3p input file, "input.dat" is the input file name
      	        fit3p_input = new File("fit3p/input.dat");
      	        FileWriter fw = new FileWriter(fit3p_input);
      	        fw.write(result);
      	        fw.close();
      	}
      	catch (IOException e) {
      	        System.out.println("Wrong input file for fit3p!");
      	        System.out.println(e.getMessage());
      	        System.exit(0);
      	}

      	try {
      	    // system call for fit3p
      		String[] command = {dir+ "/bin/fit3pbnd.exe"};
      		File runningDir = new File("fit3p");
      	    Process fit = Runtime.getRuntime().exec(command, null, runningDir);
      	    int exitValue = fit.waitFor();
      	}
      	catch (Exception e) {
      	    System.out.println("Error in run fit3p!");
      	    System.out.println(e.getMessage());
      	    System.exit(0);
      	}

      	// parse the output file from chemdis
      	try {
      	        String fit3p_output = "fit3p/output.dat";

      	        FileReader in = new FileReader(fit3p_output);
      	        BufferedReader data = new BufferedReader(in);

      	        String line = ChemParser.readMeaningfulLine(data);
      	        line = line.trim();
      	    	StringTokenizer st = new StringTokenizer(line);
      	        String A = st.nextToken();
      	        String temp = st.nextToken();
      	        temp = st.nextToken();
      	        temp = st.nextToken();
      	        double Ar = Double.parseDouble(temp);

      			line = ChemParser.readMeaningfulLine(data);
      	        line = line.trim();
      	        st = new StringTokenizer(line);
      	        String n = st.nextToken();
      	        temp = st.nextToken();
      	        temp = st.nextToken();
      	        double nr = Double.parseDouble(temp);

      			line = ChemParser.readMeaningfulLine(data);
      	        line = line.trim();
      	        st = new StringTokenizer(line);
      	        String E = st.nextToken();
      	        temp = st.nextToken();
      	        temp = st.nextToken();
      	        temp = st.nextToken();
      	        double Er = Double.parseDouble(temp);
      			if (Er < 0) {
      				System.err.println(getStructure().toString());
      				System.err.println("fitted Er < 0: "+Double.toString(Er));
      				double increase = Math.exp(-Er/GasConstant.getKcalMolK()/715.0);
      				double deltan = Math.log(increase)/Math.log(715.0);
      				System.err.println("n enlarged by factor of: " + Double.toString(deltan));
      				nr += deltan;
      				Er = 0;

      			}


      	    	UncertainDouble udAr = new UncertainDouble(Ar, 0, "Adder");
      	        UncertainDouble udnr = new UncertainDouble(nr, 0, "Adder");
      	        UncertainDouble udEr = new UncertainDouble(Er, 0, "Adder");

      	        fittedReverseKinetics = new ArrheniusKinetics(udAr, udnr , udEr, "300-1500", 1, "fitting from forward and thermal",null);
      	        in.close();
      	}
      	catch (Exception e) {
      	        System.out.println("Error in read output.dat from fit3p!");
      	        System.out.println(e.getMessage());
      	        System.exit(0);
      	}

      return fittedReverseKinetics;
      //#]
  }
// */	
	
  //## operation fitReverseKineticsRoughly()
  public void fitReverseKineticsRoughly() {
      //#[ operation fitReverseKineticsRoughly()
      // now is a rough fitting
      if (isForward()) {
      	fittedReverseKinetics = null;
      }
      else {
    //double temp = 715;
    //    double temp = 298.15; //10/29/07 gmagoon: Sandeep made change to temp = 298 on his computer locally
    //    double temp = 1350; //11/6/07 gmagoon:**** changed to actual temperature in my condition file to create agreement with old version; apparently, choice of temp has large effect; //11/9/07 gmagoon: commented out
          double temp = 298.15; //11/9/07 gmagoon: restored use of 298.15 per discussion with Sandeep
          //double temp = Global.temperature.getK();
    	Kinetics[] k = getKinetics();
    	fittedReverseKinetics = new Kinetics[k.length];
      	double doubleAlpha;
      	for (int numKinetics=0; numKinetics<k.length; numKinetics++) {
	      	if (k[numKinetics] instanceof ArrheniusEPKinetics) doubleAlpha = ((ArrheniusEPKinetics)k[numKinetics]).getAlphaValue();
	      	else doubleAlpha = 0;
	      	double Hrxn = calculateHrxn(new Temperature(temp,"K"));
	      	double Srxn = calculateSrxn(new Temperature(temp, "K"));
			  // for EvansPolyani kinetics (Ea = Eo + alpha * Hrxn) remember that  k.getEValue() gets Eo not Ea
			  // this Hrxn is for the reverse reaction (ie. -Hrxn_forward)
	      	double doubleEr = k[numKinetics].getEValue() - (doubleAlpha-1)*Hrxn;
	      	if (doubleEr < 0) {
	      		System.err.println("fitted Er < 0: "+Double.toString(doubleEr));
	      		System.err.println(getStructure().toString());
	      		//doubleEr = 0;
	      	}
	
	      	UncertainDouble Er = new UncertainDouble(doubleEr, k[numKinetics].getE().getUncertainty(), k[numKinetics].getE().getType());
	      	UncertainDouble n = new UncertainDouble(0,0, "Adder");
	
	      	double doubleA = k[numKinetics].getAValue()* Math.pow(temp, k[numKinetics].getNValue())* Math.exp(Srxn/GasConstant.getCalMolK());
	      	doubleA *= Math.pow(GasConstant.getCCAtmMolK()*temp, -getStructure().getDeltaN());  // assumes Ideal gas law concentration and 1 Atm reference state
	      	fittedReverseKinetics[numKinetics] = new ArrheniusKinetics(new UncertainDouble(doubleA, 0, "Adder"), n , Er, "300-1500", 1, "fitting from forward and thermal",null);
      	}
      }
      return;
      //#]
  }

  /**
   * Generates a reaction whose structure is opposite to that of the present reaction.
   * Just appends the rate constant of this reaction to the reverse reaction.
   *
   */
  //## operation generateReverseReaction()
  public void generateReverseReaction() {
      //#[ operation generateReverseReaction()
      Structure s = getStructure();

      //Kinetics k = getKinetics();
	  Kinetics[] k = kinetics;
	  if (kinetics == null)
		throw new NullPointerException();
		Structure newS = s.generateReverseStructure();
	  newS.setRedundancy(s.getRedundancy());
      Reaction r = new Reaction(newS, k);

//	  if (hasAdditionalKinetics()){
//		  r.addAdditionalKinetics(additionalKinetics,1);
//	  }
	  
      r.setReverseReaction(this);
      this.setReverseReaction(r);

      return;
      //#]
  }

  

 
  
  public String getComments() {
     
      return comments;
     
  }

  //## operation getDirection()
  public int getDirection() {
      //#[ operation getDirection()
      return getStructure().getDirection();
      //#]
  }

  //## operation getFittedReverseKinetics()
  public Kinetics[] getFittedReverseKinetics() {
      //#[ operation getFittedReverseKinetics()
      if (fittedReverseKinetics == null) fitReverseKineticsRoughly();
      return fittedReverseKinetics;
      //#]
  }

  /*//## operation getForwardRateConstant()
  public Kinetics getForwardRateConstant() {
      //#[ operation getForwardRateConstant()
      if (isForward()) return kinetics;
      else return null;
      //#]
  }
*/
  //## operation getKinetics()
  public Kinetics[] getKinetics() {
      //#[ operation getKinetics()
	  // returns the kinetics OF THE FORWARD REACTION
	  // ie. if THIS is reverse, it calls this.getReverseReaction().getKinetics()
	  
	  /*
	   * 29Jun2009-MRH:
	   * 	When getting kinetics, check whether it comes from a PRL or not.
	   * 	If so, return the kinetics.  We are not worried about the redundancy
	   * 	because I assume the user inputs the Arrhenius kinetics for the overall
	   * 	reaction A = B + C
	   * 
	   * 	E.g. CH4 = CH3 + H		A n E
	   * 		The Arrhenius parameters would be for the overall decomposition of CH4,
	   * 		not for each carbon-hydrogen bond fission
	   */
	  if (isFromPrimaryKineticLibrary()) {
		  return kinetics;
	  }
      if (isForward()) {
      	int red = structure.getRedundancy();
      	Kinetics[] kinetics2return = new Kinetics[kinetics.length];
      	for (int numKinetics=0; numKinetics<kinetics.length; ++numKinetics) {
      		kinetics2return[numKinetics] = kinetics[numKinetics].multiply(red);
      	}
		return kinetics2return;
      }
      else if (isBackward()) {
      	Reaction rr = getReverseReaction();
      	//	Added by MRH on 7/Sept/2009
      	//	Required when reading in the restart files
      	if (rr == null) {
      		generateReverseReaction();
      		rr = getReverseReaction();
      	}
      	if (rr == null) throw new NullPointerException("Reverse reaction is null.\n" + structure.toString());
      	if (!rr.isForward()) throw new InvalidReactionDirectionException(structure.toString());
      	return rr.getKinetics();
      }
      else 
		  throw new InvalidReactionDirectionException(structure.toString());



      //#]
  }

  public void setKineticsComments(String p_string, int num_k){
	  kinetics[num_k].setComments(p_string);
  }
  
  
  // shamel: Added this function 6/10/2010, to get Kinetics Source to identify duplicates 
  // in Reaction Library, Seed Mech and Template Reaction with Library Reaction feature
  
  public String getKineticsSource(int num_k){
	  // The num_k is the number for different kinetics stored for one type of reaction but formed due to different families
	  // Check if the kinetics exits
	  if (kinetics != null) {
		  // Check if the "source" string is not null
		  if(kinetics[num_k].getSource() != null){
			  return kinetics[num_k].getSource();  
		  }
		  else{
			  // This is mostly done for case of H Abstraction where forward kinetic source is null
			  //we might need to check if this "source" is also null (Can be source of Bug)
			  return this.reverseReaction.kinetics[num_k].getSource();  
		  }
	  } else
		  // Returns Source as null when there are no Kinetics at all!
		  return null;
	  
  }
  
  
  public void setKineticsSource(String p_string, int num_k){
	  kinetics[num_k].setSource(p_string);
  }
  
  //## operation getUpperBoundRate(Temperature)
  public double getUpperBoundRate(Temperature p_temperature){//svp
    //#[ operation getUpperBoundRate(Temperature)
    if (UpperBoundRate == 0.0){
      calculateUpperBoundRate(p_temperature);
    }
    return UpperBoundRate;
    //#]
  }

  //## operation getLowerBoundRate(Temperature)
  public double getLowerBoundRate(Temperature p_temperature){//svp
    //#[ operation getLowerBoundRate(Temperature)
    if (LowerBoundRate == 0.0){
      calculateLowerBoundRate(p_temperature);
    }
    return LowerBoundRate;
    //#]
  }


  //## operation getProductList()
  public LinkedList getProductList() {
      //#[ operation getProductList()
      return structure.getProductList();
      //#]
  }

  //10/26/07 gmagoon: changed to have temperature and pressure passed as parameters (part of eliminating use of Global.temperature)
  public double getRateConstant(Temperature p_temperature){
  //public double getRateConstant(){
	  if (rateConstant == 0)
                  rateConstant = calculateTotalRate(p_temperature);
	//	  rateConstant = calculateTotalRate(Global.temperature);
	  return rateConstant;
  }
  
  //## operation getProductNumber()
  public int getProductNumber() {
      //#[ operation getProductNumber()
      return getStructure().getProductNumber();
      //#]
  }

  //## operation getProducts()
  public ListIterator getProducts() {
      //#[ operation getProducts()
      return structure.getProducts();
      //#]
  }

  /*//## operation getRateConstant()
  public Kinetics getRateConstant() {
      //#[ operation getRateConstant()
      if (isForward()) {
      	return rateConstant;
      }
      else if (isBackward()) {
      	Reaction rr = getReverseReaction();
      	if (rr == null) throw new NullPointerException("Reverse reaction is null.\n" + structure.toString());
      	if (!rr.isForward()) throw new InvalidReactionDirectionException(structure.toString());
      	return rr.getRateConstant();
      }
      else throw new InvalidReactionDirectionException(structure.toString());
      //#]
  }*/

  //## operation getReactantList()
  public LinkedList getReactantList() {
      //#[ operation getReactantList()
      return structure.getReactantList();
      //#]
  }

  //## operation getReactantNumber()
  public int getReactantNumber() {
      //#[ operation getReactantNumber()
      return getStructure().getReactantNumber();
      //#]
  }

  //## operation getReactants()
  public ListIterator getReactants() {
      //#[ operation getReactants()
      return structure.getReactants();
      //#]
  }

  //## operation getRedundancy()
  public int getRedundancy() {
      //#[ operation getRedundancy()
      return getStructure().getRedundancy();
      //#]
  }

	 public boolean hasResonanceIsomer() {
	        //#[ operation hasResonanceIsomer()
	        return (hasResonanceIsomerAsReactant() || hasResonanceIsomerAsProduct());
	        //#]
	    }

	    //## operation hasResonanceIsomerAsProduct()
	    public boolean hasResonanceIsomerAsProduct() {
	        //#[ operation hasResonanceIsomerAsProduct()
	        for (Iterator iter = getProducts(); iter.hasNext();) {
	        	Species spe = ((Species)iter.next());
	        	if (spe.hasResonanceIsomers()) return true;
	        }
	        return false;
	        //#]
	    }

	    //## operation hasResonanceIsomerAsReactant()
	    public boolean hasResonanceIsomerAsReactant() {
	        //#[ operation hasResonanceIsomerAsReactant()
	        for (Iterator iter = getReactants(); iter.hasNext();) {
				Species spe = ((Species)iter.next());
	        	if (spe.hasResonanceIsomers()) return true;
	        }
	        return false;
	        //#]
	    }

	    //## operation hasReverseReaction()
	    public boolean hasReverseReaction() {
	        //#[ operation hasReverseReaction()
	        return reverseReaction != null;
	        //#]
	    }


  //## operation hashCode()
  public int hashCode() {
      //#[ operation hashCode()
      // just use the structure's hashcode
      return structure.hashCode();


      //#]
  }

	   //## operation isDuplicated(Reaction)
  /*public boolean isDuplicated(Reaction p_reaction) {
      //#[ operation isDuplicated(Reaction)
      // the same structure, return true
      Structure str1 = getStructure();
      Structure str2 = p_reaction.getStructure();

      //if (str1.isDuplicate(str2)) return true;

      // if not the same structure, check the resonance isomers
      if (!hasResonanceIsomer()) return false;

      if (str1.equals(str2)) return true;
      else return false;


      //#]
  }*/

  //## operation isBackward()
  public boolean isBackward() {
      //#[ operation isBackward()
      return structure.isBackward();
      //#]
  }

  //## operation isForward()
  public boolean isForward() {
      //#[ operation isForward()
      return structure.isForward();
      //#]
  }

  //## operation isIncluded(HashSet)
  public boolean isIncluded(HashSet p_speciesSet) {
      //#[ operation isIncluded(HashSet)
      return (allReactantsIncluded(p_speciesSet) && allProductsIncluded(p_speciesSet));
      //#]
  }

  //## operation makeReaction(Structure,Kinetics,boolean)
  public static Reaction makeReaction(Structure p_structure, Kinetics[] p_kinetics, boolean p_generateReverse) {
      //#[ operation makeReaction(Structure,Kinetics,boolean)
      if (!p_structure.repOk()) throw new InvalidStructureException(p_structure.toChemkinString(false).toString());
      for (int numKinetics=0; numKinetics<p_kinetics.length; numKinetics++) {
    	  if (!p_kinetics[numKinetics].repOk()) throw new InvalidKineticsException(p_kinetics[numKinetics].toString());
      }


      Reaction r = new Reaction(p_structure, p_kinetics);

      if (p_generateReverse) {
      	r.generateReverseReaction();
      }
      else {
      	r.setReverseReaction(null);
      }

      return r;
      //#]
  }

  //## operation reactantEqualsProduct()
  public boolean reactantEqualsProduct() {
      //#[ operation reactantEqualsProduct()
      return getStructure().reactantEqualsProduct();
      //#]
  }

  //## operation repOk()
  public boolean repOk() {
      //#[ operation repOk()
      if (!structure.repOk()) {
      	System.out.println("Invalid Reaction Structure:" + structure.toString());
      	return false;
      }

      if (!isForward() && !isBackward()) {
      	System.out.println("Invalid Reaction Direction: " + String.valueOf(getDirection()));
      	return false;
      }
      if (isBackward() && reverseReaction == null) {
      	System.out.println("Backward Reaction without a reversed reaction defined!");
      	return false;
      }

      /*if (!getRateConstant().repOk()) {
      	System.out.println("Invalid Rate Constant: " + getRateConstant().toString());
      	return false;
      }*/

      Kinetics[] allKinetics = getKinetics();
      for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
	      if (!allKinetics[numKinetics].repOk()) {
	      	System.out.println("Invalid Kinetics: " + allKinetics[numKinetics].toString());
	      	return false;
	      }
      }

      if (!checkRateRange()) {
      	System.out.println("reaction rate is higher than the upper rate limit!");
      	System.out.println(getStructure().toString());
      	Temperature tup = new Temperature(1500,"K");
      	if (isForward()) {
      		System.out.println("k(T=1500) = " + String.valueOf(calculateTotalRate(tup)));
      	}
      	else {
      		System.out.println("k(T=1500) = " + String.valueOf(calculateTotalRate(tup)));
      		System.out.println("Keq(T=1500) = " + String.valueOf(calculateKeq(tup)));
      		System.out.println("krev(T=1500) = " + String.valueOf(getReverseReaction().calculateTotalRate(tup)));
      	}
      	System.out.println(getKinetics());
      	return false;
      }
      return true;
      //#]
  }

  //## operation setReverseReaction(Reaction)
  public void setReverseReaction(Reaction p_reverseReaction) {
      //#[ operation setReverseReaction(Reaction)
      reverseReaction = p_reverseReaction;
      if (p_reverseReaction != null) reverseReaction.reverseReaction = this;
      //#]
  }


	  //## operation toChemkinString()
  public String toChemkinString(Temperature p_temperature) {
      //#[ operation toChemkinString()
	  if (ChemkinString != null)
		  return ChemkinString;
	  	StringBuilder result = new StringBuilder();
      	StringBuilder strucString = getStructure().toChemkinString(hasReverseReaction());
		Temperature stdtemp = new Temperature(298,"K");
		double Hrxn = calculateHrxn(stdtemp);
		Kinetics[] allKinetics = getKinetics();
		for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
			String k = allKinetics[numKinetics].toChemkinString(Hrxn,p_temperature,true);
			if (allKinetics.length == 1)
				result.append(strucString + " " + k);
			else
				result.append(strucString + " " + k + "\nDUP\n");
		}		
		ChemkinString = result.toString();
		return result.toString();
  }
  
  
  public String toChemkinString(Temperature p_temperature, Pressure p_pressure) {
	  // For certain PDep cases it's helpful to be able to call this with a temperature and pressure
	  // but usually (and in this case) the pressure is irrelevant, so we just call the above toChemkinString(Temperature) method:
	  return toChemkinString(p_temperature);
  }

	public String toRestartString(Temperature p_temperature, boolean pathReaction) {
		/*
		 * Edited by MRH on 18Jan2010
		 * 
		 * Writing restart files was causing a bug in the RMG-generated chem.inp file
		 * 	For example, H+CH4=CH3+H2 in input file w/HXD13, CH4, and H2
		 * 	RMG would correctly multiply the A factor by the structure's redundancy
		 * 		when calculating the rate to place in the ODEsolver input file.  However,
		 * 		the A reported in the chem.inp file would be the "per event" A.  This was
		 * 		due to the reaction.toChemkinString() method being called when writing the
		 * 		Restart coreReactions.txt and edgeReactions.txt files.  At the first point of
		 * 		writing the chemkinString for this reaction (when it is still an edge reaction),
		 * 		RMG had not yet computed the redundancy of the structure (as H was not a core
		 * 		species at the time, but CH3 and H2 were).  When RMG tried to write the chemkinString
		 * 		for the above reaction, using the correct redundancy, the chemkinString already existed
		 * 		and thus the toChemkinString() method was exited immediately.
		 * 	MRH is replacing the reaction.toChemkinString() call with reaction.toRestartString()
		 * 		when writing the Restart files, to account for this bug.
		 */
		
		String result = getStructure().toRestartString(hasReverseReaction()).toString(); //+ " "+getStructure().direction + " "+getStructure().redundancy;
		// MRH 18Jan2010: Restart files do not look for direction/redundancy
		/*
		 * MRH 14Feb2010: Handle reactions with multiple kinetics
		 */
		String totalResult = "";
		Kinetics[] allKinetics = getKinetics();
		for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
			totalResult += result + " " + allKinetics[numKinetics].toChemkinString(calculateHrxn(p_temperature),p_temperature,true);
			if (allKinetics.length != 1) {
				if (pathReaction) {
					totalResult += "\n";
					if (numKinetics != allKinetics.length-1) totalResult += getStructure().direction + "\t"; 
				}
				else totalResult += "\n\tDUP\n";
			}
		}
		return totalResult;
	}

	/*
	 * MRH 23MAR2010:
	 * 	Method not used in RMG
	 */
//  //## operation toFullString()
//  public String toFullString() {
//      //#[ operation toFullString()
//      return getStructure().toString() + getKinetics().toString() + getComments().toString();
//
//
//
//      //#]
//  }

 
  //## operation toString()
  public String toString(Temperature p_temperature) {
      //#[ operation toString()
	  String string2return = "";
      Kinetics[] k = getKinetics();
      for (int numKinetics=0; numKinetics<k.length; ++numKinetics) {
    	  string2return += getStructure().toString() + "\t";
    	  string2return += k[numKinetics].toChemkinString(calculateHrxn(p_temperature),p_temperature,false);
    	  if (k.length > 1) string2return += "\n";
      }
      return string2return;
  }

  /*
   * MRH 23MAR2010:
   * 	This method is redundant to toString()
   */
  //10/26/07 gmagoon: changed to take temperature as parameter (required changing function name from toString to reactionToString
//  public String reactionToString(Temperature p_temperature) {
// // public String toString() {
//      //#[ operation toString()
//	//  Temperature p_temperature = Global.temperature;
//      Kinetics k = getKinetics();
//      String kString = k.toChemkinString(calculateHrxn(p_temperature),p_temperature,false);
//      
//      return getStructure().toString() + '\t' + kString;
//      //#]
//  }
  
  public static double getBIMOLECULAR_RATE_UPPER() {
      return BIMOLECULAR_RATE_UPPER;
  }

  public static double getUNIMOLECULAR_RATE_UPPER() {
      return UNIMOLECULAR_RATE_UPPER;
  }

  public void setComments(String p_comments) {
      comments = p_comments;
  }
  
  

  /**
   * Returns the reverse reaction of this reaction. If there is no reverse reaction present
   * then a null object is returned.
   * @return
   */
  public Reaction getReverseReaction() {
      return reverseReaction;
  }

  public void setKinetics(Kinetics p_kinetics, int k_index) {
	  if (p_kinetics == null) {
		  kinetics = null;
	  }
	  else {
		  kinetics[k_index] = p_kinetics;
	  }
  }
	
	public void addAdditionalKinetics(Kinetics p_kinetics, int red) {
		if (finalized)
			return;
		if (p_kinetics == null)
			return;
		if (kinetics == null){
			kinetics = new Kinetics[1];
			kinetics[0] = p_kinetics;
			structure.redundancy = 1;
		}
		else {
			boolean kineticsAlreadyPresent = false;
			for (int numKinetics=0; numKinetics<kinetics.length; ++numKinetics) {
				if (kinetics[numKinetics].equals(p_kinetics)) {
					structure.increaseRedundancy(red);
					kineticsAlreadyPresent = true;
				}
			}
			if (!kineticsAlreadyPresent) {
				Kinetics[] tempKinetics = kinetics;
				kinetics = new Kinetics[tempKinetics.length+1];
				for (int i=0; i<tempKinetics.length; i++) {
					kinetics[i] = tempKinetics[i];
				}
				kinetics[kinetics.length-1] = p_kinetics;
				structure.redundancy = 1;
			}
		}
		
		/*
		 * MRH 24MAR2010:
		 * 	Commented out.  As RMG will be able to handle more than 2 Kinetics
		 * 		per reaction, the code below is no longer necessary
		 */
		//10/29/07 gmagoon: changed to use Global.highTemperature, Global.lowTemperature (versus Global.temperature); apparently this function chooses the top two rates when there are multiple reactions with same reactants and products; the reactions with the top two rates are used; use of high and low temperatures would be less than ideal in cases where temperature of system changes over the course of reaction
		//10/31/07 gmagoon: it is assumed that two rate constants vs. temperature cross each other at most one time over the temperature range of interest
		//if there at least three different reactions/rates and rate crossings/intersections occur, the two rates used are based on the lowest simulation temperature and a warning is displayed
//		else if (additionalKinetics == null){
//			if (p_kinetics.calculateRate(Global.lowTemperature) > kinetics.calculateRate(Global.lowTemperature)){
//				if (p_kinetics.calculateRate(Global.highTemperature) < kinetics.calculateRate(Global.highTemperature))
//					ratesForKineticsAndAdditionalKineticsCross = true;
//				additionalKinetics = kinetics;
//				kinetics = p_kinetics;
//				structure.redundancy = 1;
//			}
//			else{
//				if(p_kinetics.calculateRate(Global.highTemperature) > kinetics.calculateRate(Global.highTemperature))
//					ratesForKineticsAndAdditionalKineticsCross = true;
//				additionalKinetics = p_kinetics;
//				}
//			}
//		else if (additionalKinetics.equals(p_kinetics))
//			return;
//		else {
//			if(ratesForKineticsAndAdditionalKineticsCross){
//				if(p_kinetics.calculateRate(Global.lowTemperature) > kinetics.calculateRate(Global.lowTemperature)){
//					if(p_kinetics.calculateRate(Global.highTemperature) > kinetics.calculateRate(Global.highTemperature))
//						ratesForKineticsAndAdditionalKineticsCross = false;
//					System.out.println("WARNING: reaction that may be significant at higher temperatures within the provided range is being neglected; see following for details:");
//					System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//					System.out.println("Ignoring the rate constant " + additionalKinetics.toString() );				
//					additionalKinetics = kinetics;
//					kinetics = p_kinetics;
//				}
//				else if (p_kinetics.calculateRate(Global.lowTemperature) < additionalKinetics.calculateRate(Global.lowTemperature)){
//					if(p_kinetics.calculateRate(Global.highTemperature) > kinetics.calculateRate(Global.highTemperature))
//						System.out.println("WARNING: reaction that may be significant at higher temperatures within the provided range is being neglected; see following for details:");
//					System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//					System.out.println("Ignoring the rate constant " + p_kinetics.toString() );
//				}
//				else{//else p_kinetics @ low temperature is between kinetics and additional kinetics at low temperature
//					if(p_kinetics.calculateRate(Global.highTemperature) < kinetics.calculateRate(Global.highTemperature))
//						ratesForKineticsAndAdditionalKineticsCross = false;
//					System.out.println("WARNING: reaction that may be significant at higher temperatures within the provided range is being neglected; see following for details:");
//					System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//					System.out.println("Ignoring the rate constant " + additionalKinetics.toString() );
//					additionalKinetics = p_kinetics;
//				}                        
//			}
//			else{
//				if ((p_kinetics.calculateRate(Global.lowTemperature) > kinetics.calculateRate(Global.lowTemperature)) && (p_kinetics.calculateRate(Global.highTemperature) > kinetics.calculateRate(Global.highTemperature))){
//					System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//					System.out.println("Ignoring the rate constant " + additionalKinetics.toString() ); //10/29/07 gmagoon: note that I have moved this before reassignment of variables; I think this was a minor bug in original code
//					additionalKinetics = kinetics;
//					kinetics = p_kinetics;
//				}
//				else if ((p_kinetics.calculateRate(Global.lowTemperature) < additionalKinetics.calculateRate(Global.lowTemperature))&&(p_kinetics.calculateRate(Global.highTemperature) < additionalKinetics.calculateRate(Global.highTemperature))){
//					System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//					System.out.println("Ignoring the rate constant " + p_kinetics.toString() );
//				}
//				else if ((p_kinetics.calculateRate(Global.lowTemperature) > additionalKinetics.calculateRate(Global.lowTemperature))&&(p_kinetics.calculateRate(Global.highTemperature) > additionalKinetics.calculateRate(Global.highTemperature))&&(p_kinetics.calculateRate(Global.lowTemperature) < kinetics.calculateRate(Global.lowTemperature))&&(p_kinetics.calculateRate(Global.highTemperature) < kinetics.calculateRate(Global.highTemperature))){
//					System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//					System.out.println("Ignoring the rate constant " + additionalKinetics.toString() );
//					additionalKinetics = p_kinetics;
//				}
//				else //else there is at least one crossing in the temperature range of interest between p_kinetics and either kinetics or additionalKinetics; base which reaction is kept on the lowest temperature
//                {
//					if(p_kinetics.calculateRate(Global.lowTemperature) > kinetics.calculateRate(Global.lowTemperature)){
//						if(p_kinetics.calculateRate(Global.highTemperature) < additionalKinetics.calculateRate(Global.highTemperature))
//							System.out.println("WARNING: reaction that may be significant at higher temperatures within the provided range is being neglected; see following for details:");
//						System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//						System.out.println("Ignoring the rate constant " + additionalKinetics.toString() );
//						additionalKinetics = kinetics;
//						kinetics = p_kinetics;
//						ratesForKineticsAndAdditionalKineticsCross = true;
//					}
//					else if(p_kinetics.calculateRate(Global.lowTemperature) < additionalKinetics.calculateRate(Global.lowTemperature)){
//						System.out.println("WARNING: reaction that may be significant at higher temperatures within the provided range is being neglected; see following for details:");
//						System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//						System.out.println("Ignoring the rate constant " + p_kinetics.toString() );
//					}
//					else //else p_kinetics at low temperature is between kinetics and additional kinetics at low temperature
//                    {
//						if(p_kinetics.calculateRate(Global.highTemperature) > kinetics.calculateRate(Global.highTemperature))
//							ratesForKineticsAndAdditionalKineticsCross = true;
//						else//else p_kinetics crosses additional kinetics
//							System.out.println("WARNING: reaction that may be significant at higher temperatures within the provided range is being neglected; see following for details:");
//						System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//						System.out.println("Ignoring the rate constant " + additionalKinetics.toString() );
//						additionalKinetics = p_kinetics;
//					}
//				}
//			}
//		}                
//                else if (additionalKinetics == null){
//			if (p_kinetics.calculateRate(Global.temperature) > kinetics.calculateRate(Global.temperature)){
//				additionalKinetics = kinetics;
//				kinetics = p_kinetics;
//				structure.redundancy = 1;
//			}
//			else additionalKinetics = p_kinetics;
//              }
//		else if (additionalKinetics.equals(p_kinetics))
//			return;
//		else {
//			if (p_kinetics.calculateRate(Global.temperature) > kinetics.calculateRate(Global.temperature)){
//				additionalKinetics = kinetics;
//				kinetics = p_kinetics;
//				System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//				System.out.println("Ignoring the rate constant " + additionalKinetics.toString() );
//			
//			}
//			else if (p_kinetics.calculateRate(Global.temperature) < additionalKinetics.calculateRate(Global.temperature)){
//				System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//				System.out.println("Ignoring the rate constant " + p_kinetics.toString() );
//			}
//			else {
//				System.out.println("More than 2 kinetics provided for reaction " + structure.toChemkinString(true));
//				System.out.println("Ignoring the rate constant " + additionalKinetics.toString() );
//				additionalKinetics = p_kinetics;
//			}
//		}
	}

	/*
	 * MRH 18MAR2010:
	 * A reaction's kinetics is now an array.  The additionalKinetics is now
	 * 	obsolete and this method will be removed (to ensure nothing calls it)
	 */
//	public boolean hasAdditionalKinetics(){
//		return (additionalKinetics != null);
//	}
	
	//public int totalNumberOfKinetics(){ //7/26/09 gmagoon: this is not used, and appears to incorrectly assume that there are a maximum of two kinetics...I think this was old and I have changed it since
	//	if (hasAdditionalKinetics())
	//		return 2;
	//	else
	//		return 1;
	//}
	
	/*
	 * MRH 18MAR2010:
	 * Restructuring a reaction's kinetics
	 * 	With a reaction's kinetics now being defined as an array of Kinetics,
	 * 		instead of "kinetics" and "additionalKinetics", the getAllKinetics()
	 * 		method is now obsolete.  
	 */
//	public HashSet getAllKinetics(){
//		HashSet allKinetics = new HashSet();
//		allKinetics.add(kinetics.multiply(structure.redundancy));
//		if ( hasAdditionalKinetics()){
//			allKinetics.add(additionalKinetics.multiply(structure.redundancy));
//			
//		}
//		
//		return allKinetics;	
//	}
	
	

	public void setFinalized(boolean p_finalized) {
		finalized = p_finalized;
		return;
	}
	
  public Structure getStructure() {
      return structure;
  }

  public void setStructure(Structure p_Structure) {
      structure = p_Structure;
  }
  
  /**
	 * Returns the reaction as an ASCII string.
	 * @return A string representing the reaction equation in ASCII test.
	 */
	@Override
	public String toString() {
		if (getReactantNumber() == 0 || getProductNumber() == 0)
			return "";
		
		String rxn = "";
		Species species = (Species) structure.getReactantList().get(0);
		rxn = rxn + species.getName() + "(" + Integer.toString(species.getID()) + ")";
		for (int i = 1; i < getReactantNumber(); i++) {
			species = (Species) structure.getReactantList().get(i);
			rxn += " + " + species.getName() + "(" + Integer.toString(species.getID()) + ")";
		}
		rxn += " --> ";
		species = (Species) structure.getProductList().get(0);
		rxn = rxn + species.getName() + "(" + Integer.toString(species.getID()) + ")";
		for (int i = 1; i < getProductNumber(); i++) {
			species = (Species) structure.getProductList().get(i);
			rxn += " + " + species.getName() + "(" + Integer.toString(species.getID()) + ")";
		}
		
		return rxn;
	}
	
	public String toInChIString() {
		if (getReactantNumber() == 0 || getProductNumber() == 0)
			return "";
		
		String rxn = "";
		Species species = (Species) structure.getReactantList().get(0);
		rxn = rxn + species.getInChI();
		for (int i = 1; i < getReactantNumber(); i++) {
			species = (Species) structure.getReactantList().get(i);
			rxn += " + " + species.getInChI();
		}
		rxn += " --> ";
		species = (Species) structure.getProductList().get(0);
		rxn = rxn + species.getInChI();
		for (int i = 1; i < getProductNumber(); i++) {
			species = (Species) structure.getProductList().get(i);
			rxn += " + " + species.getInChI();
		}
		
		return rxn;
	}

	/**
	 * Calculates the flux of this reaction given the provided system snapshot.
	 * The system snapshot contains the temperature, pressure, and
	 * concentrations of each core species.
	 * @param ss The system snapshot at which to determine the reaction flux
	 * @return The determined reaction flux
	 */
	public double calculateFlux(SystemSnapshot ss) {
		return calculateForwardFlux(ss) - calculateReverseFlux(ss);
	}

	/**
	 * Calculates the forward flux of this reaction given the provided system snapshot.
	 * The system snapshot contains the temperature, pressure, and
	 * concentrations of each core species.
	 * @param ss The system snapshot at which to determine the reaction flux
	 * @return The determined reaction flux
	 */
	public double calculateForwardFlux(SystemSnapshot ss) {
		Temperature T = ss.getTemperature();
		double forwardFlux = getRateConstant(T);
		for (ListIterator<Species> iter = getReactants(); iter.hasNext(); ) {
			Species spe = iter.next();
			double conc = 0.0;
			if (ss.getSpeciesStatus(spe) != null)
				conc = ss.getSpeciesStatus(spe).getConcentration();
			if (conc < 0) {
				double aTol = ReactionModelGenerator.getAtol();
				//if (Math.abs(conc) < aTol) conc = 0;
				//else throw new NegativeConcentrationException(spe.getName() + ": " + String.valueOf(conc));
				if (conc < -100.0 * aTol)
					throw new NegativeConcentrationException("Species " + spe.getName() + " has negative concentration: " + String.valueOf(conc));
			}
			forwardFlux *= conc;
		}
		return forwardFlux;
	}

	/**
	 * Calculates the flux of this reaction given the provided system snapshot.
	 * The system snapshot contains the temperature, pressure, and
	 * concentrations of each core species.
	 * @param ss The system snapshot at which to determine the reaction flux
	 * @return The determined reaction flux
	 */
	public double calculateReverseFlux(SystemSnapshot ss) {
		if (hasReverseReaction())
			return reverseReaction.calculateForwardFlux(ss);
		else
			return 0.0;
	}
	
	public boolean isFromPrimaryKineticLibrary() {
		return kineticsFromPrimaryKineticLibrary;
	}
	
	public void setIsFromPrimaryKineticLibrary(boolean p_boolean) {
		kineticsFromPrimaryKineticLibrary = p_boolean;
	}
	
	public ReactionTemplate getReactionTemplate() {
		return rxnTemplate;
	}

	public void setReactionTemplate(ReactionTemplate rt) {
		rxnTemplate = rt;
	}
	
	public boolean hasMultipleKinetics() {
		if (getKinetics().length > 1) return true;
		else return false;
	}

}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\Reaction.java
*********************************************************************/


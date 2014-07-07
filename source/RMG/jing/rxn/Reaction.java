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
package jing.rxn;

import java.io.*;
import jing.chem.*;
import java.util.*;
import jing.param.*;
import jing.mathTool.*;
import jing.chemParser.*;
import jing.chem.Species;
import jing.param.Temperature;
import jing.rxnSys.Logger;
import jing.rxnSys.NegativeConcentrationException;
import jing.rxnSys.ReactionModelGenerator;
import jing.rxnSys.SystemSnapshot;

// ## package jing::rxn
// ----------------------------------------------------------------------------
// jing\rxn\Reaction.java
// ----------------------------------------------------------------------------
/**
 * Immutable objects.
 */
// ## class Reaction
public class Reaction {
    protected static double TRIMOLECULAR_RATE_UPPER = 1.0E100;
    protected static double BIMOLECULAR_RATE_UPPER = 1.0E100; // ## attribute BIMOLECULAR_RATE_UPPER
    protected static double UNIMOLECULAR_RATE_UPPER = 1.0E100; // ## attribute UNIMOLECULAR_RATE_UPPER
    protected String comments = "No comment"; // ## attribute comments
    protected Kinetics[] fittedReverseKinetics = null; // ## attribute fittedReverseKinetics
    protected double rateConstant;
    protected Reaction reverseReaction = null; // ## attribute reverseReaction
    protected Kinetics[] kinetics;
    protected Structure structure;
    protected double UpperBoundRate;// svp
    protected double LowerBoundRate;// svp
    protected boolean finalized = false;
    protected String ChemkinString = null;
    protected boolean kineticsFromPrimaryKineticLibrary = false;
    protected boolean expectDuplicate = false;

    // Constructors
    // ## operation Reaction()
    public Reaction() {
        // #[ operation Reaction()
        // #]
    }

    // ## operation Reaction(Structure,RateConstant)
    private Reaction(Structure p_structure, Kinetics[] p_kinetics) {
        // #[ operation Reaction(Structure,RateConstant)
        structure = p_structure;
        kinetics = p_kinetics;
        // rateConstant = calculateTotalRate(Global.temperature);
        // #]
    }

    // ## operation allProductsIncluded(LinkedHashSet)
    public boolean allProductsIncluded(LinkedHashSet p_speciesSet) {
        // #[ operation allProductsIncluded(LinkedHashSet)
        Iterator iter = getProducts();
        while (iter.hasNext()) {
            Species spe = ((Species) iter.next());
            if (!p_speciesSet.contains(spe))
                return false;
        }
        return true;
        // #]
    }

    // ## operation allReactantsIncluded(LinkedHashSet)
    public boolean allReactantsIncluded(LinkedHashSet p_speciesSet) {
        // #[ operation allReactantsIncluded(LinkedHashSet)
        if (p_speciesSet == null)
            throw new NullPointerException();
        Iterator iter = getReactants();
        while (iter.hasNext()) {
            Species spe = ((Species) iter.next());
            if (!p_speciesSet.contains(spe))
                return false;
        }
        return true;
        // #]
    }

    /**
     * Calculate this reaction's thermo parameter. Basically, make addition of the thermo parameters of all the
     * reactants and products.
     */
    // ## operation calculateHrxn(Temperature)
    public double calculateHrxn(Temperature p_temperature) {
        // #[ operation calculateHrxn(Temperature)
        return structure.calculateHrxn(p_temperature);
        // #]
    }

    // ## operation calculateKeq(Temperature)
    public double calculateKeq(Temperature p_temperature) {
        // #[ operation calculateKeq(Temperature)
        return structure.calculateKeq(p_temperature);
        // #]
    }

    // ## operation calculateKeqUpperBound(Temperature)
    // svp
    public double calculateKeqUpperBound(Temperature p_temperature) {
        // #[ operation calculateKeqUpperBound(Temperature)
        return structure.calculateKeqUpperBound(p_temperature);
        // #]
    }

    // ## operation calculateKeqLowerBound(Temperature)
    // svp
    public double calculateKeqLowerBound(Temperature p_temperature) {
        // #[ operation calculateKeqLowerBound(Temperature)
        return structure.calculateKeqLowerBound(p_temperature);
        // #]
    }

    public double calculateTotalRate(Temperature p_temperature) {
        double rate = 0;
        Temperature stdtemp = new Temperature(298, "K");
        double Hrxn = calculateHrxn(stdtemp);
        Temperature sys_temp = ReactionModelGenerator.getTemp4BestKinetics();
        /*
         * AJ 12JULY2010: Added diffusive limits from previous RMG version by replacing function calculateTotalRate
         * Checks the exothermicity and molecularity of reaction to determine the diffusive rate limit
         */
        /*
         * 29Jun2009-MRH: Added a kinetics from PRL check If the kinetics for this reaction is from a PRL, use those
         * numbers to compute the rate. Else, proceed as before.
         */
        /*
         * MRH 18MAR2010: Changing the structure of a reaction's kinetics If the kinetics are from a primary kinetic
         * library, we assume the user has supplied the total pre-exponential factor for the reaction (and not the
         * per-event pre-exponential facor). If the kinetics were estimated by RMG, the pre-exponential factor must be
         * multiplied by the "redundancy" (# of events)
         */
        if (kineticsFromPrimaryKineticLibrary) {
            Kinetics[] k_All = kinetics;
            for (int numKinetics = 0; numKinetics < kinetics.length; numKinetics++) {
                Kinetics k = k_All[numKinetics];
                rate += k.calculateRate(p_temperature, Hrxn);
            }            
            if(isBackward()){
            	rate = rate * calculateKeq(p_temperature);	
            }
            return rate;
        } else if (isForward()) {
            Kinetics[] k_All = kinetics;
            for (int numKinetics = 0; numKinetics < kinetics.length; numKinetics++) {
                Kinetics k = k_All[numKinetics];
                if ((int) structure.redundancy != 1) {
                    k = k.multiply(structure.redundancy);
                }
                rate += k.calculateRate(p_temperature, Hrxn);
            }
            /*
             * Diffusion limits added by AJ on July 12, 2010 Requires correction in the forward direction only, reverse
             * reaction corrects itself If ReactionModelGenerator.useDiffusion is true (solvation is on) compute kd and
             * return keff
             */
            if (ReactionModelGenerator.getUseDiffusion()) {
                int numReacts = structure.getReactantNumber();
                int numProds = structure.getProductNumber();
                double keff = 0.0;
                double DiffFactor = 0.0;
                if (numReacts == 1 && numProds == 1) {
                    keff = rate;
                    rate = keff;
                    DiffFactor = 1;
                } else if (numReacts == 1 && numProds == 2) {
                    double k_back = rate / calculateKeq(p_temperature);
                    LinkedList reactantsInBackRxn = structure.products;
                    double k_back_diff = calculatediff(reactantsInBackRxn);
                    double k_back_eff = k_back * k_back_diff
                            / (k_back + k_back_diff);
                    keff = k_back_eff * calculateKeq(p_temperature);
                    DiffFactor = keff / rate;
                    rate = keff;
                } else if (numReacts == 2 && numProds == 1) {
                    double k_forw = rate;
                    LinkedList reactantsInForwRxn = structure.reactants;
                    double k_forw_diff = calculatediff(reactantsInForwRxn);
                    double k_forw_eff = k_forw * k_forw_diff
                            / (k_forw + k_forw_diff);
                    keff = k_forw_eff;
                    DiffFactor = keff / rate;
                    rate = keff;
                } else if (numReacts == 2 && numProds == 2) {
                    double rxn_Keq = structure.calculateKeq(p_temperature);
                    double deltaHrxn = structure.calculateHrxn(p_temperature);
                    if (rxn_Keq > 1) { // Forward reaction is exothermic hence the corresponding diffusion limit applies
                        double k_forw = rate;
                        LinkedList reactantsInForwRxn = structure.reactants;
                        double k_forw_diff = calculatediff(reactantsInForwRxn);
                        double k_forw_eff = k_forw * k_forw_diff
                                / (k_forw + k_forw_diff);
                        keff = k_forw_eff;
                        DiffFactor = keff / rate;
                        rate = keff;
                    } else if (rxn_Keq < 1) { // Reverse reaction is exothermic and the corresponding diffusion limit
// should be used
                        double k_back = rate / calculateKeq(p_temperature);
                        LinkedList reactantsInBackRxn = structure.products;
                        double k_back_diff = calculatediff(reactantsInBackRxn);
                        double k_back_eff = k_back * k_back_diff
                                / (k_back + k_back_diff);
                        keff = k_back_eff * calculateKeq(p_temperature);
                        DiffFactor = keff / rate;
                        rate = keff;
                    }
                } else if (numReacts == 2 && numProds == 3) {
                    double k_forw = rate;
                    LinkedList reactantsInForwRxn = structure.reactants;
                    double k_forw_diff = calculatediff(reactantsInForwRxn);
                    double k_forw_eff = k_forw * k_forw_diff
                            / (k_forw + k_forw_diff);
                    keff = k_forw_eff;
                    DiffFactor = keff / rate;
                    rate = keff;
                } else if (numReacts == 2 && numProds == 4) {
                    double k_forw = rate;
                    LinkedList reactantsInForwRxn = structure.reactants;
                    double k_forw_diff = calculatediff(reactantsInForwRxn);
                    double k_forw_eff = k_forw * k_forw_diff
                            / (k_forw + k_forw_diff);
                    keff = k_forw_eff;
                    DiffFactor = keff / rate;
                    rate = keff;
                }
                // Add comments only if the temperature at which the function has been called corresponds to the system
// temperature specified in the condition file
                // And if they haven't already been added
                if (p_temperature.getK() == sys_temp.getK()) {
                    String oldComments = k_All[0].getComment();
                    if (!oldComments.contains("Diffusion")) {
                        if (numReacts == 1 && numProds == 1)
                            oldComments += " Unimolecular: no diffusion limitation.";
                        String newComments = oldComments
                                + String.format(
                                        " Diffusion factor=%.3g, rate=%.3g, Keq=%.3g, at T=%.0fK.",
                                        DiffFactor, rate,
                                        calculateKeq(p_temperature),
                                        p_temperature.getK());
                        setKineticsComments(newComments, 0);
                    }
                }
            }
            return rate;
        } else if (isBackward()) {
            Reaction r = getReverseReaction();
            rate = r.calculateTotalRate(p_temperature);
            return rate * calculateKeq(p_temperature);
        } else {
            throw new InvalidReactionDirectionException();
        }
    }

    // ## operation calculatediff(LinkedList)
    public double calculatediff(LinkedList p_struct) {
        if (p_struct.size() != 2) {
            Logger.warning("Cannot compute diffusive limit if number of reactants is not equal to 2");
        }
        // Array containing the radii of the two species passed in p_struct
        double[] r;
        double[] d;
        r = new double[2];
        d = new double[2];
        int i = 0;
        for (Iterator iter = p_struct.iterator(); iter.hasNext();) {
            Species sp = (Species) iter.next();
            ChemGraph cg = sp.getChemGraph();
            r[i] = cg.getRadius();
            d[i] = cg.getDiffusivity();
            i = i + 1;
        }
        double kdiff;
        kdiff = (88 / 7) * (d[0] + d[1]) * (r[0] + r[1]) * 6.023e29; // units of r[i]=m; d[1]=m2/sec; kdiff=cm3/mole sec
        return kdiff;
    }

    // ## operation calculateUpperBoundRate(Temperature)
    // svp
    public double calculateUpperBoundRate(Temperature p_temperature) {
        // #[ operation calculateUpperBoundRate(Temperature)
        if (isForward()) {
            double A;
            double E;
            double n;
            for (int numKinetics = 0; numKinetics < kinetics.length; numKinetics++) {
                A = kinetics[numKinetics].getA().getUpperBound();
                E = kinetics[numKinetics].getE().getLowerBound();
                n = kinetics[numKinetics].getN().getUpperBound();
                if (A > 1E300) {
                    A = kinetics[numKinetics].getA().getValue() * 1.2;
                }
                // Kinetics kinetics = getRateConstant().getKinetics();
                if (kinetics[numKinetics] instanceof ArrheniusEPKinetics) {
                    ArrheniusEPKinetics arrhenius = (ArrheniusEPKinetics) kinetics[numKinetics];
                    double H = calculateHrxn(p_temperature);
                    if (H < 0) {
                        if (arrhenius.getAlpha().getValue() > 0) {
                            double alpha = arrhenius.getAlpha().getUpperBound();
                            E = E + alpha * H;
                        } else {
                            double alpha = arrhenius.getAlpha().getLowerBound();
                            E = E + alpha * H;
                        }
                    } else {
                        if (arrhenius.getAlpha().getValue() > 0) {
                            double alpha = arrhenius.getAlpha().getLowerBound();
                            E = E + alpha * H;
                        } else {
                            double alpha = arrhenius.getAlpha().getUpperBound();
                            E = E + alpha * H;
                        }
                    }
                }
                if (E < 0) {
                    E = 0;
                }
                double indiv_k = 0.0;
                indiv_k = A
                        * Math.pow(p_temperature.getK(), n)
                        * Math.exp(-E / GasConstant.getKcalMolK()
                                / p_temperature.getK());
                indiv_k *= getStructure().getRedundancy();
                UpperBoundRate += indiv_k;
            }
            return UpperBoundRate;
        } else if (isBackward()) {
            Reaction r = getReverseReaction();
            if (r == null)
                throw new NullPointerException("Reverse reaction is null.\n"
                        + structure.toString());
            if (!r.isForward())
                throw new InvalidReactionDirectionException();
            for (int numKinetics = 0; numKinetics < kinetics.length; numKinetics++) {
                double A = kinetics[numKinetics].getA().getUpperBound();
                double E = kinetics[numKinetics].getE().getLowerBound();
                double n = kinetics[numKinetics].getN().getUpperBound();
                if (A > 1E300) {
                    A = kinetics[numKinetics].getA().getValue() * 1.2;
                }
                // Kinetics kinetics = getRateConstant().getKinetics();
                if (kinetics[numKinetics] instanceof ArrheniusEPKinetics) {
                    ArrheniusEPKinetics arrhenius = (ArrheniusEPKinetics) kinetics[numKinetics];
                    double H = calculateHrxn(p_temperature);
                    if (H < 0) {
                        if (arrhenius.getAlpha().getValue() > 0) {
                            double alpha = arrhenius.getAlpha().getUpperBound();
                            E = E + alpha * H;
                        } else {
                            double alpha = arrhenius.getAlpha().getLowerBound();
                            E = E + alpha * H;
                        }
                    } else {
                        if (arrhenius.getAlpha().getValue() > 0) {
                            double alpha = arrhenius.getAlpha().getLowerBound();
                            E = E + alpha * H;
                        } else {
                            double alpha = arrhenius.getAlpha().getUpperBound();
                            E = E + alpha * H;
                        }
                    }
                }
                if (E < 0) {
                    E = 0;
                }
                double indiv_k = 0.0;
                indiv_k = A
                        * Math.pow(p_temperature.getK(), n)
                        * Math.exp(-E / GasConstant.getKcalMolK()
                                / p_temperature.getK());
                indiv_k *= getStructure().getRedundancy();
                UpperBoundRate += indiv_k
                        * calculateKeqUpperBound(p_temperature);
            }
            return UpperBoundRate;
        } else {
            throw new InvalidReactionDirectionException();
        }
        // #]
    }

    // ## operation calculateLowerBoundRate(Temperature)
    // svp
    public double calculateLowerBoundRate(Temperature p_temperature) {
        // #[ operation calculateLowerBoundRate(Temperature)
        if (isForward()) {
            for (int numKinetics = 0; numKinetics < kinetics.length; ++numKinetics) {
                double A = kinetics[numKinetics].getA().getLowerBound();
                double E = kinetics[numKinetics].getE().getUpperBound();
                double n = kinetics[numKinetics].getN().getLowerBound();
                if (A > 1E300 || A <= 0) {
                    A = kinetics[numKinetics].getA().getValue() / 1.2;
                }
                // Kinetics kinetics = getRateConstant().getKinetics();
                if (kinetics[numKinetics] instanceof ArrheniusEPKinetics) {
                    ArrheniusEPKinetics arrhenius = (ArrheniusEPKinetics) kinetics[numKinetics];
                    double H = calculateHrxn(p_temperature);
                    if (H < 0) {
                        if (arrhenius.getAlpha().getValue() > 0) {
                            double alpha = arrhenius.getAlpha().getLowerBound();
                            E = E + alpha * H;
                        } else {
                            double alpha = arrhenius.getAlpha().getUpperBound();
                            E = E + alpha * H;
                        }
                    } else {
                        if (arrhenius.getAlpha().getValue() > 0) {
                            double alpha = arrhenius.getAlpha().getUpperBound();
                            E = E + alpha * H;
                        } else {
                            double alpha = arrhenius.getAlpha().getLowerBound();
                            E = E + alpha * H;
                        }
                    }
                }
                double indiv_k = 0.0;
                indiv_k = A
                        * Math.pow(p_temperature.getK(), n)
                        * Math.exp(-E / GasConstant.getKcalMolK()
                                / p_temperature.getK());
                indiv_k *= getStructure().getRedundancy();
                LowerBoundRate += indiv_k;
            }
            return LowerBoundRate;
        } else if (isBackward()) {
            Reaction r = getReverseReaction();
            if (r == null)
                throw new NullPointerException("Reverse reaction is null.\n"
                        + structure.toString());
            if (!r.isForward())
                throw new InvalidReactionDirectionException();
            for (int numKinetics = 0; numKinetics < kinetics.length; ++numKinetics) {
                double A = kinetics[numKinetics].getA().getLowerBound();
                double E = kinetics[numKinetics].getE().getUpperBound();
                double n = kinetics[numKinetics].getN().getLowerBound();
                if (A > 1E300) {
                    A = kinetics[numKinetics].getA().getValue() / 1.2;
                }
                if (kinetics[numKinetics] instanceof ArrheniusEPKinetics) {
                    ArrheniusEPKinetics arrhenius = (ArrheniusEPKinetics) kinetics[numKinetics];
                    double H = calculateHrxn(p_temperature);
                    if (H < 0) {
                        if (arrhenius.getAlpha().getValue() > 0) {
                            double alpha = arrhenius.getAlpha().getLowerBound();
                            E = E + alpha * H;
                        } else {
                            double alpha = arrhenius.getAlpha().getUpperBound();
                            E = E + alpha * H;
                        }
                    } else {
                        if (arrhenius.getAlpha().getValue() > 0) {
                            double alpha = arrhenius.getAlpha().getUpperBound();
                            E = E + alpha * H;
                        } else {
                            double alpha = arrhenius.getAlpha().getLowerBound();
                            E = E + alpha * H;
                        }
                    }
                }
                double indiv_k = 0.0;
                indiv_k = A
                        * Math.pow(p_temperature.getK(), n)
                        * Math.exp(-E / GasConstant.getKcalMolK()
                                / p_temperature.getK());
                indiv_k *= getStructure().getRedundancy();
                LowerBoundRate += indiv_k
                        * calculateKeqLowerBound(p_temperature);
            }
            return LowerBoundRate;
        } else {
            throw new InvalidReactionDirectionException();
        }
        // #]
    }

    // ## operation calculateSrxn(Temperature)
    public double calculateSrxn(Temperature p_temperature) {
        // #[ operation calculateSrxn(Temperature)
        return structure.calculateSrxn(p_temperature);
        // #]
    }

    // ## operation calculateThirdBodyCoefficient(SystemSnapshot)
    public double calculateThirdBodyCoefficient(SystemSnapshot p_presentStatus) {
        // #[ operation calculateThirdBodyCoefficient(SystemSnapshot)
        if (!(this instanceof ThirdBodyReaction))
            return 1;
        else {
            return ((ThirdBodyReaction) this)
                    .calculateThirdBodyCoefficient(p_presentStatus);
        }
        // #]
    }

    // ## operation checkRateRange()
    public boolean checkRateRange() {
        // #[ operation checkRateRange()
        Temperature t = new Temperature(1500, "K");
        double rate = calculateTotalRate(t);
        if (getReactantNumber() == 2) {
            if (rate > BIMOLECULAR_RATE_UPPER)
                return false;
        } else if (getReactantNumber() == 1) {
            if (rate > UNIMOLECULAR_RATE_UPPER)
                return false;
        } else if (getReactantNumber() == 3) {
            if (rate > TRIMOLECULAR_RATE_UPPER)
                return false;
        } else
            throw new InvalidReactantNumberException();
        return true;
        // #]
    }

    // ## operation contains(Species)
    public boolean contains(Species p_species) {
        // #[ operation contains(Species)
        if (containsAsReactant(p_species) || containsAsProduct(p_species))
            return true;
        else
            return false;
        // #]
    }

    // ## operation containsAsProduct(Species)
    public boolean containsAsProduct(Species p_species) {
        // #[ operation containsAsProduct(Species)
        Iterator iter = getProducts();
        while (iter.hasNext()) {
            // ChemGraph cg = (ChemGraph)iter.next();
            Species spe = (Species) iter.next();
            if (spe.equals(p_species))
                return true;
        }
        return false;
        // #]
    }

    // ## operation containsAsReactant(Species)
    public boolean containsAsReactant(Species p_species) {
        // #[ operation containsAsReactant(Species)
        Iterator iter = getReactants();
        while (iter.hasNext()) {
            // ChemGraph cg = (ChemGraph)iter.next();
            Species spe = (Species) iter.next();
            if (spe.equals(p_species))
                return true;
        }
        return false;
        // #]
    }

    /**
     * Checks if the structure of the reaction is the same. Does not check the rate constant. Two reactions with the
     * same structure but different rate constants will be equal.
     */
    // ## operation equals(Object)
    public boolean equals(Object p_reaction) {
        // #[ operation equals(Object)
        if (this == p_reaction)
            return true;
        if (!(p_reaction instanceof Reaction))
            return false;
        Reaction r = (Reaction) p_reaction;
        if (!getStructure().equals(r.getStructure()))
            return false;
        return true;
    }

/*
 * // fitReverseKineticsPrecisely and getPreciseReverseKinetics // are not used, not maintained, and we have clue what
 * they do, // so we're commenting them out so we don't keep looking at them. // (oh, and they look pretty similar to
 * each other!) // - RWest & JWAllen, June 2009 //## operation fitReverseKineticsPrecisely() public void
 * fitReverseKineticsPrecisely() { //#[ operation fitReverseKineticsPrecisely() if (isForward()) { fittedReverseKinetics
 * = null; } else { String result = ""; for (double t = 300.0; t<1500.0; t+=50.0) { double rate = calculateTotalRate(new
 * Temperature(t,"K")); result += String.valueOf(t) + '\t' + String.valueOf(rate) + '\n'; } // run fit3p String dir =
 * System.getProperty("RMG.workingDirectory"); File fit3p_input; try { // prepare fit3p input file, "input.dat" is the
 * input file name fit3p_input = new File("fit3p/input.dat"); FileWriter fw = new FileWriter(fit3p_input);
 * fw.write(result); fw.close(); } catch (IOException e) { System.out.println("Wrong input file for fit3p!");
 * System.out.println(e.getMessage()); System.exit(0); } try { // system call for fit3p String[] command = {dir+
 * "/bin/fit3pbnd.exe"}; File runningDir = new File("fit3p"); Process fit = Runtime.getRuntime().exec(command, null,
 * runningDir); int exitValue = fit.waitFor(); } catch (Exception e) { System.out.println("Error in run fit3p!");
 * System.out.println(e.getMessage()); System.exit(0); } // parse the output file from chemdis try { String fit3p_output
 * = "fit3p/output.dat"; FileReader in = new FileReader(fit3p_output); BufferedReader data = new BufferedReader(in);
 * String line = ChemParser.readMeaningfulLine(data); line = line.trim(); StringTokenizer st = new
 * StringTokenizer(line); String A = st.nextToken(); String temp = st.nextToken(); temp = st.nextToken(); temp =
 * st.nextToken(); double Ar = Double.parseDouble(temp); line = ChemParser.readMeaningfulLine(data); line = line.trim();
 * st = new StringTokenizer(line); String n = st.nextToken(); temp = st.nextToken(); temp = st.nextToken(); double nr =
 * Double.parseDouble(temp); line = ChemParser.readMeaningfulLine(data); line = line.trim(); st = new
 * StringTokenizer(line); String E = st.nextToken(); temp = st.nextToken(); temp = st.nextToken(); temp =
 * st.nextToken(); double Er = Double.parseDouble(temp); if (Er < 0) { System.err.println(getStructure().toString());
 * System.err.println("fitted Er < 0: "+Double.toString(Er)); double increase =
 * Math.exp(-Er/GasConstant.getKcalMolK()/715.0); double deltan = Math.log(increase)/Math.log(715.0);
 * System.err.println("n enlarged by factor of: " + Double.toString(deltan)); nr += deltan; Er = 0; } UncertainDouble
 * udAr = new UncertainDouble(Ar, 0, "Adder"); UncertainDouble udnr = new UncertainDouble(nr, 0, "Adder");
 * UncertainDouble udEr = new UncertainDouble(Er, 0, "Adder"); fittedReverseKinetics = new ArrheniusKinetics(udAr, udnr
 * , udEr, "300-1500", 1, "fitting from forward and thermal",null); in.close(); } catch (Exception e) {
 * System.out.println("Error in read output.dat from fit3p!"); System.out.println(e.getMessage()); System.exit(0); } }
 * return; //#] } //## operation fitReverseKineticsPrecisely() public Kinetics getPreciseReverseKinetics() { //#[
 * operation fitReverseKineticsPrecisely() Kinetics fittedReverseKinetics =null; String result = ""; for (double t =
 * 300.0; t<1500.0; t+=50.0) { double rate = calculateTotalRate(new Temperature(t,"K")); result += String.valueOf(t) +
 * '\t' + String.valueOf(rate) + '\n'; } // run fit3p String dir = System.getProperty("RMG.workingDirectory"); File
 * fit3p_input; try { // prepare fit3p input file, "input.dat" is the input file name fit3p_input = new
 * File("fit3p/input.dat"); FileWriter fw = new FileWriter(fit3p_input); fw.write(result); fw.close(); } catch
 * (IOException e) { System.out.println("Wrong input file for fit3p!"); System.out.println(e.getMessage());
 * System.exit(0); } try { // system call for fit3p String[] command = {dir+ "/bin/fit3pbnd.exe"}; File runningDir = new
 * File("fit3p"); Process fit = Runtime.getRuntime().exec(command, null, runningDir); int exitValue = fit.waitFor(); }
 * catch (Exception e) { System.out.println("Error in run fit3p!"); System.out.println(e.getMessage()); System.exit(0);
 * } // parse the output file from chemdis try { String fit3p_output = "fit3p/output.dat"; FileReader in = new
 * FileReader(fit3p_output); BufferedReader data = new BufferedReader(in); String line =
 * ChemParser.readMeaningfulLine(data); line = line.trim(); StringTokenizer st = new StringTokenizer(line); String A =
 * st.nextToken(); String temp = st.nextToken(); temp = st.nextToken(); temp = st.nextToken(); double Ar =
 * Double.parseDouble(temp); line = ChemParser.readMeaningfulLine(data); line = line.trim(); st = new
 * StringTokenizer(line); String n = st.nextToken(); temp = st.nextToken(); temp = st.nextToken(); double nr =
 * Double.parseDouble(temp); line = ChemParser.readMeaningfulLine(data); line = line.trim(); st = new
 * StringTokenizer(line); String E = st.nextToken(); temp = st.nextToken(); temp = st.nextToken(); temp =
 * st.nextToken(); double Er = Double.parseDouble(temp); if (Er < 0) { System.err.println(getStructure().toString());
 * System.err.println("fitted Er < 0: "+Double.toString(Er)); double increase =
 * Math.exp(-Er/GasConstant.getKcalMolK()/715.0); double deltan = Math.log(increase)/Math.log(715.0);
 * System.err.println("n enlarged by factor of: " + Double.toString(deltan)); nr += deltan; Er = 0; } UncertainDouble
 * udAr = new UncertainDouble(Ar, 0, "Adder"); UncertainDouble udnr = new UncertainDouble(nr, 0, "Adder");
 * UncertainDouble udEr = new UncertainDouble(Er, 0, "Adder"); fittedReverseKinetics = new ArrheniusKinetics(udAr, udnr
 * , udEr, "300-1500", 1, "fitting from forward and thermal",null); in.close(); } catch (Exception e) {
 * System.out.println("Error in read output.dat from fit3p!"); System.out.println(e.getMessage()); System.exit(0); }
 * return fittedReverseKinetics; //#] } //
 */
    // ## operation fitReverseKineticsRoughly()
    public void fitReverseKineticsRoughly() {
        // #[ operation fitReverseKineticsRoughly()
        // now is a rough fitting
        if (isForward()) {
            fittedReverseKinetics = null;
        } else {
            // double temp = 715;
            double temp = 298.15; // 11/9/07 gmagoon: restored use of 298.15 per discussion with Sandeep
            // double temp = Global.temperature.getK();
            Kinetics[] k = getKinetics();
            fittedReverseKinetics = new Kinetics[k.length];
            double doubleAlpha;
            for (int numKinetics = 0; numKinetics < k.length; numKinetics++) {
                if (k[numKinetics] instanceof ArrheniusEPKinetics)
                    doubleAlpha = ((ArrheniusEPKinetics) k[numKinetics])
                            .getAlphaValue();
                else
                    doubleAlpha = 0;
                double Hrxn = calculateHrxn(new Temperature(temp, "K"));
                double Srxn = calculateSrxn(new Temperature(temp, "K"));
                // for EvansPolyani kinetics (Ea = Eo + alpha * Hrxn) remember that k.getEValue() gets Eo not Ea
                // this Hrxn is for the reverse reaction (ie. -Hrxn_forward)
                double doubleEr = k[numKinetics].getEValue()
                        - (doubleAlpha - 1) * Hrxn;
                if (doubleEr < 0) {
                    Logger.warning("fitted Er < 0: "
                            + Double.toString(doubleEr));
                    Logger.warning(getStructure().toString());
                    // doubleEr = 0;
                }
                UncertainDouble Er = new UncertainDouble(doubleEr,
                        k[numKinetics].getE().getUncertainty(), k[numKinetics]
                                .getE().getType());
                UncertainDouble n = new UncertainDouble(0, 0, "Adder");
                double doubleA = k[numKinetics].getAValue()
                        * Math.pow(temp, k[numKinetics].getNValue())
                        * Math.exp(Srxn / GasConstant.getCalMolK());
                doubleA *= Math.pow(GasConstant.getCCAtmMolK() * temp,
                        -getStructure().getDeltaN()); // assumes Ideal gas law concentration and 1 Atm reference state
                fittedReverseKinetics[numKinetics] = new ArrheniusKinetics(
                        new UncertainDouble(doubleA, 0, "Adder"), n, Er,
                        "300-1500", 1, "fitting from forward and thermal", null);
            }
        }
        return;
        // #]
    }

    /**
     * Generates a reaction whose structure is opposite to that of the present reaction. Just appends the rate constant
     * of this reaction to the reverse reaction.
     */
    // ## operation generateReverseReaction()
    public void generateReverseReaction() {
        // #[ operation generateReverseReaction()
        Structure s = getStructure();
        // Kinetics k = getKinetics();
        Kinetics[] k = kinetics;
        if (kinetics == null)
            throw new NullPointerException();
        Structure newS = s.generateReverseStructure();
        newS.setRedundancy(s.getRedundancy());
        Reaction r = new Reaction(newS, k);
// if (hasAdditionalKinetics()){
// r.addAdditionalKinetics(additionalKinetics,1);
// }
        r.setReverseReaction(this);
        this.setReverseReaction(r);
        return;
        // #]
    }

    public String getComments() {
        return comments;
    }

    // ## operation getDirection()
    public int getDirection() {
        // #[ operation getDirection()
        return getStructure().getDirection();
        // #]
    }

    // ## operation getFittedReverseKinetics()
    public Kinetics[] getFittedReverseKinetics() {
        // #[ operation getFittedReverseKinetics()
        if (fittedReverseKinetics == null)
            fitReverseKineticsRoughly();
        return fittedReverseKinetics;
        // #]
    }

    /*
     * //## operation getForwardRateConstant() public Kinetics getForwardRateConstant() { //#[ operation
     * getForwardRateConstant() if (isForward()) return kinetics; else return null; //#] }
     */
    // ## operation getKinetics()
    public Kinetics[] getKinetics() {
        // #[ operation getKinetics()
        // returns the kinetics OF THE FORWARD REACTION
        // ie. if THIS is reverse, it calls this.getReverseReaction().getKinetics()
        /*
         * 29Jun2009-MRH: When getting kinetics, check whether it comes from a PRL or not. If so, return the kinetics.
         * We are not worried about the redundancy because I assume the user inputs the Arrhenius kinetics for the
         * overall reaction A = B + C E.g. CH4 = CH3 + H A n E The Arrhenius parameters would be for the overall
         * decomposition of CH4, not for each carbon-hydrogen bond fission
         */
        if (isFromPrimaryKineticLibrary()) {
            return kinetics;
        }
        if (isForward()) {
            int red = structure.getRedundancy();
            if (red == 1)
                return kinetics; // Don't waste time multiplying by 1
            Kinetics[] kinetics2return = new Kinetics[kinetics.length];
            for (int numKinetics = 0; numKinetics < kinetics.length; ++numKinetics) {
                kinetics2return[numKinetics] = kinetics[numKinetics]
                        .multiply(red);
            }
            return kinetics2return;
        } else if (isBackward()) {
            Reaction rr = getReverseReaction();
            // Added by MRH on 7/Sept/2009
            // Required when reading in the restart files
            if (rr == null) {
                generateReverseReaction();
                rr = getReverseReaction();
            }
            if (rr == null)
                throw new NullPointerException("Reverse reaction is null.\n"
                        + structure.toString());
            if (!rr.isForward())
                throw new InvalidReactionDirectionException(
                        structure.toString());
            return rr.getKinetics();
        } else
            throw new InvalidReactionDirectionException(structure.toString());
    }

    public void setKineticsComments(String p_string, int num_k) {
        kinetics[num_k].setComments(p_string);
    }

    // shamel: Added this function 6/10/2010, to get Kinetics Source to identify duplicates
    // in Reaction Library, Seed Mech and Template Reaction with Library Reaction feature
    public String getKineticsSource(int num_k) {
        // The num_k is the number for different kinetics stored for one type of reaction but formed due to different
// families
        // Returns Source as null when there are no Kinetics at all!
        if (kinetics == null) {
            return null;
        }
        String source = null;
        if (kinetics[num_k] != null) {
            source = kinetics[num_k].getSource();
        }
        // This is mostly done for case of H Abstraction where forward kinetic source may be null
        if (source == null) {
            try {
                source = this.reverseReaction.kinetics[num_k].getSource();
            } catch (NullPointerException e) {
                // this catches several possibilities:
                // this.reverseReaction == null
                // this.reverseReaction.kinetics == null
                // this.reverseReaction.kinetics[num_k] == null
                source = null;
            }
        }
        return source;
    }

    public void setKineticsSource(String p_string, int num_k) {
        kinetics[num_k].setSource(p_string);
    }

    // ## operation getUpperBoundRate(Temperature)
    public double getUpperBoundRate(Temperature p_temperature) {// svp
        // #[ operation getUpperBoundRate(Temperature)
        if (UpperBoundRate == 0.0) {
            calculateUpperBoundRate(p_temperature);
        }
        return UpperBoundRate;
        // #]
    }

    // ## operation getLowerBoundRate(Temperature)
    public double getLowerBoundRate(Temperature p_temperature) {// svp
        // #[ operation getLowerBoundRate(Temperature)
        if (LowerBoundRate == 0.0) {
            calculateLowerBoundRate(p_temperature);
        }
        return LowerBoundRate;
        // #]
    }

    // ## operation getProductList()
    public LinkedList getProductList() {
        // #[ operation getProductList()
        return structure.getProductList();
        // #]
    }

    // ## operation getProductNumber()
    public int getProductNumber() {
        // #[ operation getProductNumber()
        return getStructure().getProductNumber();
        // #]
    }

    // ## operation getProducts()
    public ListIterator getProducts() {
        // #[ operation getProducts()
        return structure.getProducts();
        // #]
    }

    // ## operation getReactantList()
    public LinkedList getReactantList() {
        // #[ operation getReactantList()
        return structure.getReactantList();
        // #]
    }

    // ## operation getReactantNumber()
    public int getReactantNumber() {
        // #[ operation getReactantNumber()
        return getStructure().getReactantNumber();
        // #]
    }

    // ## operation getReactants()
    public ListIterator getReactants() {
        // #[ operation getReactants()
        try {
            return structure.getReactants();
        } catch (RuntimeException e) {
            Logger.critical("******DEBUGGING LINES FOLLOW******");
            Logger.logStackTrace(e);
            Logger.critical("This.toString:" + this.toString());
            Logger.critical("Comments:" + comments);
            Logger.critical("ChemkinString:" + ChemkinString);
            Logger.critical("Rate constant:" + rateConstant);
            Logger.critical("Finalized:" + finalized);
            Logger.critical("KineticsFromPrimaryKineticLibrary:"
                    + kineticsFromPrimaryKineticLibrary);
            Logger.critical("ExpectDuplicate:" + expectDuplicate);
            if (kinetics != null) {
                Logger.critical("Kinetics size:" + kinetics.length);
                Logger.critical("First kinetics info:" + kinetics[0].toString());
                Logger.critical("First kinetics source:"
                        + kinetics[0].getSource());
                Logger.critical("First kinetics comment:"
                        + kinetics[0].getComment());
            }
            if (fittedReverseKinetics != null) {
                Logger.critical("Reverse kinetics size:"
                        + fittedReverseKinetics.length);
                Logger.critical("First reverse kinetics info:"
                        + fittedReverseKinetics[0].toString());
                Logger.critical("First reverse kinetics source:"
                        + fittedReverseKinetics[0].getSource());
                Logger.critical("First reverse kinetics comment:"
                        + fittedReverseKinetics[0].getComment());
            }
            if (reverseReaction != null) {
                Logger.critical("***Reverse reaction info below***");
                Logger.critical("reverseReaction.toString:"
                        + reverseReaction.toString());
                Logger.critical("Comments:" + reverseReaction.comments);
                Logger.critical("ChemkinString:"
                        + reverseReaction.ChemkinString);
                Logger.critical("Rate constant:" + reverseReaction.rateConstant);
                Logger.critical("Finalized:" + reverseReaction.finalized);
                Logger.critical("KineticsFromPrimaryKineticLibrary:"
                        + reverseReaction.kineticsFromPrimaryKineticLibrary);
                Logger.critical("ExpectDuplicate:"
                        + reverseReaction.expectDuplicate);
                if (reverseReaction.structure != null) {
                    Logger.critical("(for reverse reaction): structure.toString()"
                            + reverseReaction.structure.toString());
                }
                if (reverseReaction.kinetics != null) {
                    Logger.critical("Kinetics size:"
                            + reverseReaction.kinetics.length);
                    Logger.critical("First kinetics info:"
                            + reverseReaction.kinetics[0].toString());
                    Logger.critical("First kinetics source:"
                            + reverseReaction.kinetics[0].getSource());
                    Logger.critical("First kinetics comment:"
                            + reverseReaction.kinetics[0].getComment());
                }
                if (reverseReaction.fittedReverseKinetics != null) {
                    Logger.critical("Reverse kinetics size:"
                            + reverseReaction.fittedReverseKinetics.length);
                    Logger.critical("First reverse kinetics info:"
                            + reverseReaction.fittedReverseKinetics[0]
                                    .toString());
                    Logger.critical("First reverse kinetics source:"
                            + reverseReaction.fittedReverseKinetics[0]
                                    .getSource());
                    Logger.critical("First reverse kinetics comment:"
                            + reverseReaction.fittedReverseKinetics[0]
                                    .getComment());
                }
            }
            Logger.flush();
            throw e;
        }
        // #]
    }

    // ## operation getRedundancy()
    public int getRedundancy() {
        // #[ operation getRedundancy()
        return getStructure().getRedundancy();
        // #]
    }

    public void setRedundancy(int redundancy) {
        getStructure().setRedundancy(redundancy);
    }

    public boolean hasResonanceIsomer() {
        // #[ operation hasResonanceIsomer()
        return (hasResonanceIsomerAsReactant() || hasResonanceIsomerAsProduct());
        // #]
    }

    // ## operation hasResonanceIsomerAsProduct()
    public boolean hasResonanceIsomerAsProduct() {
        // #[ operation hasResonanceIsomerAsProduct()
        for (Iterator iter = getProducts(); iter.hasNext();) {
            Species spe = ((Species) iter.next());
            if (spe.hasResonanceIsomers())
                return true;
        }
        return false;
        // #]
    }

    // ## operation hasResonanceIsomerAsReactant()
    public boolean hasResonanceIsomerAsReactant() {
        // #[ operation hasResonanceIsomerAsReactant()
        for (Iterator iter = getReactants(); iter.hasNext();) {
            Species spe = ((Species) iter.next());
            if (spe.hasResonanceIsomers())
                return true;
        }
        return false;
        // #]
    }

    // ## operation hasReverseReaction()
    public boolean hasReverseReaction() {
        // #[ operation hasReverseReaction()
        return reverseReaction != null;
        // #]
    }

    // ## operation hashCode()
    public int hashCode() {
        // #[ operation hashCode()
        // just use the structure's hashcode
        return structure.hashCode();
        // #]
    }

    // ## operation isDuplicated(Reaction)
    /*
     * public boolean isDuplicated(Reaction p_reaction) { //#[ operation isDuplicated(Reaction) // the same structure,
     * return true Structure str1 = getStructure(); Structure str2 = p_reaction.getStructure(); //if
     * (str1.isDuplicate(str2)) return true; // if not the same structure, check the resonance isomers if
     * (!hasResonanceIsomer()) return false; if (str1.equals(str2)) return true; else return false; //#] }
     */
    // ## operation isBackward()
    public boolean isBackward() {
        // #[ operation isBackward()
        return structure.isBackward();
        // #]
    }

    // ## operation isForward()
    public boolean isForward() {
        // #[ operation isForward()
        return structure.isForward();
        // #]
    }

    // ## operation isIncluded(LinkedHashSet)
    public boolean isIncluded(LinkedHashSet p_speciesSet) {
        // #[ operation isIncluded(LinkedHashSet)
        return (allReactantsIncluded(p_speciesSet) && allProductsIncluded(p_speciesSet));
        // #]
    }

    // ## operation makeReaction(Structure,Kinetics,boolean)
    public static Reaction makeReaction(Structure p_structure,
            Kinetics[] p_kinetics, boolean p_generateReverse) {
        // #[ operation makeReaction(Structure,Kinetics,boolean)
        if (!p_structure.repOk())
            throw new InvalidStructureException(p_structure.toChemkinString(
                    false).toString());
        for (int numKinetics = 0; numKinetics < p_kinetics.length; numKinetics++) {
            if (!p_kinetics[numKinetics].repOk())
                throw new InvalidKineticsException(
                        p_kinetics[numKinetics].toString());
        }
        Reaction r = new Reaction(p_structure, p_kinetics);
        if (p_generateReverse) {
            r.generateReverseReaction();
        } else {
            r.setReverseReaction(null);
        }
        return r;
        // #]
    }

    // ## operation reactantEqualsProduct()
    public boolean reactantEqualsProduct() {
        // #[ operation reactantEqualsProduct()
        return getStructure().reactantEqualsProduct();
        // #]
    }

    // ## operation repOk()
    public boolean repOk() {
        // #[ operation repOk()
        if (!structure.repOk()) {
            Logger.error("Invalid Reaction Structure:" + structure.toString());
            return false;
        }
        if (!isForward() && !isBackward()) {
            Logger.error("Invalid Reaction Direction: "
                    + String.valueOf(getDirection()));
            return false;
        }
        if (isBackward() && reverseReaction == null) {
            Logger.error("Backward Reaction without a reversed reaction defined!");
            return false;
        }
        /*
         * if (!getRateConstant().repOk()) { Logger.error("Invalid Rate Constant: " + getRateConstant().toString());
         * return false; }
         */
        Kinetics[] allKinetics = getKinetics();
        if (allKinetics == null) return false;
        
        for (int numKinetics = 0; numKinetics < allKinetics.length; ++numKinetics) {
            if (!allKinetics[numKinetics].repOk()) {
                Logger.error("Invalid Kinetics: "
                        + allKinetics[numKinetics].toString());
                return false;
            }
        }
        if (!checkRateRange()) {
            Logger.error("reaction rate is higher than the upper rate limit!");
            Logger.info(getStructure().toString());
            Temperature tup = new Temperature(1500, "K");
            if (isForward()) {
                Logger.info("k(T=1500) = "
                        + String.valueOf(calculateTotalRate(tup)));
            } else {
                Logger.info("k(T=1500) = "
                        + String.valueOf(calculateTotalRate(tup)));
                Logger.info("Keq(T=1500) = "
                        + String.valueOf(calculateKeq(tup)));
                Logger.info("krev(T=1500) = "
                        + String.valueOf(getReverseReaction()
                                .calculateTotalRate(tup)));
            }
            Logger.info(getKinetics().toString());
            return false;
        }
        return true;
        // #]
    }

    // ## operation setReverseReaction(Reaction)
    public void setReverseReaction(Reaction p_reverseReaction) {
        // #[ operation setReverseReaction(Reaction)
        reverseReaction = p_reverseReaction;
        if (p_reverseReaction != null)
            reverseReaction.reverseReaction = this;
        // #]
    }

    // ## operation toChemkinString()
    public String toChemkinString(Temperature p_temperature) {
        // #[ operation toChemkinString()
        if (ChemkinString != null)
            return ChemkinString;
        StringBuilder result = new StringBuilder();
        String strucString = String.format("%-52s", getStructure()
                .toChemkinString(hasReverseReaction()));
        Temperature stdtemp = new Temperature(298, "K");
        double Hrxn = calculateHrxn(stdtemp);
        Kinetics[] allKinetics = getKinetics();
        for (int numKinetics = 0; numKinetics < allKinetics.length; ++numKinetics) {
            String k = allKinetics[numKinetics].toChemkinString(Hrxn,
                    p_temperature, true);
            if (allKinetics.length == 1)
                result.append(strucString + " " + k);
            else
                result.append(strucString + " " + k + "\nDUP\n");
        }
        if (result.charAt(result.length() - 1) == '\n')
            result.deleteCharAt(result.length() - 1);
        ChemkinString = result.toString();
        return result.toString();
    }

    public String toChemkinString(Temperature p_temperature, Pressure p_pressure) {
        // For certain PDep cases it's helpful to be able to call this with a temperature and pressure
        // but usually (and in this case) the pressure is irrelevant, so we just call the above
// toChemkinString(Temperature) method:
        return toChemkinString(p_temperature);
    }

    public String toRestartString(Temperature p_temperature,
            boolean pathReaction) {
        /*
         * Edited by MRH on 18Jan2010 Writing restart files was causing a bug in the RMG-generated chem.inp file For
         * example, H+CH4=CH3+H2 in input file w/HXD13, CH4, and H2 RMG would correctly multiply the A factor by the
         * structure's redundancy when calculating the rate to place in the ODEsolver input file. However, the A
         * reported in the chem.inp file would be the "per event" A. This was due to the reaction.toChemkinString()
         * method being called when writing the Restart coreReactions.txt and edgeReactions.txt files. At the first
         * point of writing the chemkinString for this reaction (when it is still an edge reaction), RMG had not yet
         * computed the redundancy of the structure (as H was not a core species at the time, but CH3 and H2 were). When
         * RMG tried to write the chemkinString for the above reaction, using the correct redundancy, the chemkinString
         * already existed and thus the toChemkinString() method was exited immediately. MRH is replacing the
         * reaction.toChemkinString() call with reaction.toRestartString() when writing the Restart files, to account
         * for this bug.
         */
        String result = String.format("%-52s",
                getStructure().toRestartString(hasReverseReaction())); // + " "+getStructure().direction +
// " "+getStructure().redundancy;
        // MRH 18Jan2010: Restart files do not look for direction/redundancy
        /*
         * MRH 14Feb2010: Handle reactions with multiple kinetics
         */
        String totalResult = "";
        Kinetics[] allKinetics = getKinetics();
        for (int numKinetics = 0; numKinetics < allKinetics.length; ++numKinetics) {
            totalResult += result
                    + " "
                    + allKinetics[numKinetics].toChemkinString(
                            calculateHrxn(p_temperature), p_temperature, true);
            if (allKinetics.length != 1) {
                if (pathReaction) {
                    totalResult += "\n";
                    if (numKinetics != allKinetics.length - 1)
                        totalResult += getStructure().direction + "\t";
                } else
                    totalResult += "\n\tDUP\n";
            }
        }
        return totalResult;
    }

    /*
     * MRH 23MAR2010: Method not used in RMG
     */
// //## operation toFullString()
// public String toFullString() {
// //#[ operation toFullString()
// return getStructure().toString() + getKinetics().toString() + getComments().toString();
//
//
//
// //#]
// }
    // ## operation toString()
    public String toString(Temperature p_temperature) {
        // #[ operation toString()
        String string2return = "";
        Kinetics[] k = getKinetics();
        for (int numKinetics = 0; numKinetics < k.length; ++numKinetics) {
            string2return += getStructure().toString() + "\t";
            string2return += k[numKinetics].toChemkinString(
                    calculateHrxn(p_temperature), p_temperature, false);
            if (k.length > 1)
                string2return += "\n";
        }
        return string2return;
    }

    /*
     * MRH 23MAR2010: This method is redundant to toString()
     */
    // 10/26/07 gmagoon: changed to take temperature as parameter (required changing function name from toString to
// reactionToString
// public String reactionToString(Temperature p_temperature) {
// // public String toString() {
// //#[ operation toString()
// // Temperature p_temperature = Global.temperature;
// Kinetics k = getKinetics();
// String kString = k.toChemkinString(calculateHrxn(p_temperature),p_temperature,false);
//
// return getStructure().toString() + '\t' + kString;
// //#]
// }
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
     * Returns the reverse reaction of this reaction. If there is no reverse reaction present then a null object is
     * returned.
     * 
     * @return
     */
    public Reaction getReverseReaction() {
        return reverseReaction;
    }

    public void setKinetics(Kinetics p_kinetics, int k_index) {
        if (p_kinetics == null) {
            kinetics = null;
        } else {
            kinetics[k_index] = p_kinetics;
        }
    }

    public void addAdditionalKinetics(Kinetics p_kinetics, int red,
            boolean readingFromUserLibrary) {
        if (finalized)
            return;
        if (p_kinetics == null)
            return;
        if (kinetics == null) {
            kinetics = new Kinetics[1];
            kinetics[0] = p_kinetics;
            structure.redundancy = 1;
        } else if (readingFromUserLibrary
                || (!readingFromUserLibrary && !p_kinetics
                        .isFromPrimaryKineticLibrary())) {
            boolean kineticsAlreadyPresent = false;
            for (int i = 0; i < kinetics.length; i++) {
                Kinetics old_kinetics = kinetics[i];
                if (!(old_kinetics instanceof ArrheniusKinetics))
                    continue; // we can't compare or increment
                if (((ArrheniusKinetics) old_kinetics)
                        .equalNESource(p_kinetics)) {
                    ((ArrheniusKinetics) old_kinetics).addToA(p_kinetics.getA()
                            .multiply((double) red));
                    kineticsAlreadyPresent = true;
                }
            }
            if (!kineticsAlreadyPresent) {
                Kinetics[] tempKinetics = kinetics;
                kinetics = new Kinetics[tempKinetics.length + 1];
                for (int i = 0; i < tempKinetics.length; i++) {
                    kinetics[i] = tempKinetics[i];
                }
                kinetics[kinetics.length - 1] = p_kinetics;
                structure.redundancy = 1;
            }
        }
    }

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
     * 
     * @return A string representing the reaction equation in ASCII test.
     */
    @Override
    public String toString() {
        if (getReactantNumber() == 0 || getProductNumber() == 0)
            return "";
        String rxn = "";
        Species species = (Species) structure.getReactantList().get(0);
        rxn = rxn + species.getFullName();
        for (int i = 1; i < getReactantNumber(); i++) {
            species = (Species) structure.getReactantList().get(i);
            rxn += " + " + species.getFullName();
        }
        rxn += " --> ";
        species = (Species) structure.getProductList().get(0);
        rxn = rxn + species.getFullName();
        for (int i = 1; i < getProductNumber(); i++) {
            species = (Species) structure.getProductList().get(i);
            rxn += " + " + species.getFullName();
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
     * Calculates the flux of this reaction given the provided system snapshot. The system snapshot contains the
     * temperature, pressure, and concentrations of each core species.
     * 
     * @param ss
     *            The system snapshot at which to determine the reaction flux
     * @return The determined reaction flux
     */
    public double calculateFlux(SystemSnapshot ss) {
        return calculateForwardFlux(ss) - calculateReverseFlux(ss);
    }

    /**
     * Calculates the forward flux of this reaction given the provided system snapshot. The system snapshot contains the
     * temperature, pressure, and concentrations of each core species.
     * 
     * @param ss
     *            The system snapshot at which to determine the reaction flux
     * @return The determined reaction flux
     */
    public double calculateForwardFlux(SystemSnapshot ss) {
        Temperature T = ss.getTemperature();
        double forwardFlux = calculateTotalRate(T);
        for (ListIterator<Species> iter = getReactants(); iter.hasNext();) {
            Species spe = iter.next();
            double conc = 0.0;
            if (ss.getSpeciesStatus(spe) != null)
                conc = ss.getSpeciesStatus(spe).getConcentration();
            if (conc < 0) {
                double aTol = ReactionModelGenerator.getAtol();
                // if (Math.abs(conc) < aTol) conc = 0;
                // else throw new NegativeConcentrationException(spe.getFullName() + ": " + String.valueOf(conc));
                if (conc < -100.0 * aTol)
                    throw new NegativeConcentrationException("Species "
                            + spe.getFullName()
                            + " has negative concentration: "
                            + String.valueOf(conc));
            }
            forwardFlux *= conc;
        }
        return forwardFlux;
    }

    /**
     * Calculates the flux of this reaction given the provided system snapshot. The system snapshot contains the
     * temperature, pressure, and concentrations of each core species.
     * 
     * @param ss
     *            The system snapshot at which to determine the reaction flux
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
        if (reverseReaction != null) {
            reverseReaction.kineticsFromPrimaryKineticLibrary = p_boolean;
        }
    }

    public boolean hasMultipleKinetics() {
        if (getKinetics().length > 1)
            return true;
        else
            return false;
    }

    public void setExpectDuplicate(boolean b) {
        expectDuplicate = b;
    }

    public boolean getExpectDuplicate() {
        return expectDuplicate;
    }

    public void prune() {
        // Do what's necessary to prune the reaction
        // Also clears the reverse, so don't expect to be able to get it back
        // Do santy check on the isFromPrimaryKineticLibrary() method we rely on.
        // (This shouldn't be necessary, but we have no unit test framework so I'm building one in here!)
        if (!isFromPrimaryKineticLibrary()) {
            if (isForward())
                if (getKineticsSource(0) != null)
                    if (getKineticsSource(0).contains("Library"))
                        throw new RuntimeException(
                                String.format(
                                        "Reaction %s kinetics source contains 'Library' but isFromPrimaryKineticLibrary() returned false.",
                                        this));
            if (isBackward())
                if (getReverseReaction().getKineticsSource(0) != null)
                    if (getReverseReaction().getKineticsSource(0).contains(
                            "Library"))
                        throw new RuntimeException(
                                String.format(
                                        "Reverse of reaction %s kinetics source contains 'Library' but isFromPrimaryKineticLibrary() returned false.",
                                        this));
            // I'm not sure why we can't clear the reverse anyway - as long as the direction WITH kinetics still has a
// Structure we should be ok.
        }
        // Use isFromPrimaryKineticLibrary() to decide if it's safe to clear the reaction structure
        if (!isFromPrimaryKineticLibrary()) {
            setStructure(null);
            setReverseReaction(null);
        }
    }

    /**
     * Return the total number of atoms in the reactants (and products).
     */
    public int getAtomNumber() {
        int atoms = 0;
        for (ListIterator<Species> iter = getReactants(); iter.hasNext();) {
            Species spe = iter.next();
            atoms += spe.getChemGraph().getAtomNumber();
        }
        return atoms;
    }
    
    /**
     * Return the total number of Carbon atoms in the reactants (and products).
     */
    public int getCarbonAtomNumber() {
        int atoms = 0;
        for (ListIterator<Species> iter = getReactants(); iter.hasNext();) {
            Species spe = iter.next();
            atoms += spe.getChemGraph().getCarbonNumber();
        }
        return atoms;
    }
}
/*********************************************************************
 * File Path : RMG\RMG\jing\rxn\Reaction.java
 *********************************************************************/

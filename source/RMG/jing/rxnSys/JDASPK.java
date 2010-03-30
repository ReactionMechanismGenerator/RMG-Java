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

package jing.rxnSys;


import jing.rxn.*;
import jing.chem.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

import jing.chem.Species;
import jing.rxn.Reaction;
import jing.rxn.Structure;
import jing.rxn.TROEReaction;
import jing.rxn.ThirdBodyReaction;
import jing.param.Global;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.param.ParameterInfor;

//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\JDASPK.java
//----------------------------------------------------------------------------

//## class JDASPK
public class JDASPK extends JDAS {
	
	private JDASPK() {
        super();
    }
    
	public JDASPK(double p_rtol, double p_atol, int p_parameterInfor, 
			InitialStatus p_initialStatus, int p_index, ValidityTester p_vt, 
			boolean p_autoflag, Double p_termTol, Double p_coreTol) {
        super(p_rtol, p_atol, p_parameterInfor, p_initialStatus, p_index, p_vt, 
			p_autoflag, p_termTol, p_coreTol);
    }
    
    //6/25/08 gmagoon: defined alternate constructor for use with sensitivity analysis (lacks autoflag and validityTester parameters)
    //6/25/08 gmagoon: set autoflag to be false with this constructor (not used for sensitivity analysis)
    //3/22/10 gmagoon: set termTol and coreTol =null with this constructor; these variables should not be used
	public JDASPK(double p_rtol, double p_atol, int p_parameterInfor, InitialStatus p_initialStatus, int p_index) {
        super(p_rtol, p_atol, p_parameterInfor, p_initialStatus, p_index, null, 
			false, null, null);
    }

    //## operation generateSensitivityStatus(ReactionModel,double [],double [],int)
    private double [] generateSensitivityStatus(ReactionModel p_reactionModel, double [] p_y, double [] p_yprime, int p_paraNum) {
    	//#[ operation generateSensitivityStatus(ReactionModel,double [],double [],int)
    	int nequ = p_reactionModel.getSpeciesNumber()*(p_paraNum+1);
    	if (p_y.length != nequ) throw new DynamicSimulatorException();
    	if (p_yprime.length != nequ) throw new DynamicSimulatorException();

    	double [] senStatus = new double[nParameter*nState]; //gmagoon 12/21/09: this does not include volume and could be shrunk to nParameter*(nState-1)

    	for (int i = p_reactionModel.getSpeciesNumber();i<nequ;i++){
    		//double sens = p_y[i]; gmagoon 12/21/09: this doesn't seem to be used anywhere
    		int ind = i-p_reactionModel.getSpeciesNumber();
    		senStatus[ind] = p_y[i];
    	}
    	return senStatus;
    	//#]
    }

    //## operation solve(boolean,ReactionModel,boolean,SystemSnapshot,ReactionTime,ReactionTime,Temperature,Pressure,boolean)
    public SystemSnapshot solve(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, ReactionTime p_beginTime, ReactionTime p_endTime, Temperature p_temperature, Pressure p_pressure, boolean p_conditionChanged,TerminationTester tt, int p_iterationNum) {
        // set up the input file
        setupInputFile();
    	//outputString = new StringBuilder();
    	//first generate an id for all the species
    	Iterator spe_iter = p_reactionModel.getSpecies();
    	while (spe_iter.hasNext()){
    		Species spe = (Species)spe_iter.next();
    		int id = getRealID(spe);
    	}
    	double startTime = System.currentTimeMillis();
		
		ReactionTime rt = p_beginStatus.getTime();
        if (!rt.equals(p_beginTime)) throw new InvalidBeginStatusException();

        double tBegin = p_beginTime.getStandardTime();
        double tEnd = p_endTime.getStandardTime();

		double T = p_temperature.getK();
		double P = p_pressure.getAtm();
		
		LinkedList initialSpecies = new LinkedList();
		
        // set reaction set
        if (p_initialization || p_reactionChanged || p_conditionChanged) {
			
        	nState = p_reactionModel.getSpeciesNumber();
			
//        	troeString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10, alpha, Tstar, T2star, T3star, lowRate  (21 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0), ncollider, c1, c2,..c10, troe(0=T or 1=F) (21 elements)
        	troeString = generateTROEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
        	
//        	tbrString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10  (16 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0), ncollider, c1, c2,..c10 (20 elements)
			tbrString = generateThirdBodyReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			lindemannString = generateLindemannReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			//rString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0)
			rString = generatePDepODEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			nParameter = 0;
			if (parameterInfor != 0) {
				nParameter = rList.size() + thirdBodyList.size() + troeList.size() + lindemannList.size() + p_reactionModel.getSpeciesNumber();				
			}
			neq = nState*(nParameter + 1);
			initializeWorkSpace();
			initializeConcentrations(p_beginStatus, p_reactionModel, p_beginTime, p_endTime, initialSpecies);
			
			
        }
        //6/25/08 gmagoon: (next two lines) for autoflag, get binary 0 or 1 corresponding to boolean false/true
        int af = 0;
        if (autoflag) af = 1;
        try{
            if (tt instanceof ConversionTT){
                    SpeciesConversion sc = (SpeciesConversion)((ConversionTT)tt).speciesGoalConversionSet.get(0);
                bw.write(nState + "\t" + neq + "\t" +  getRealID(sc.species) + "\t 1" +"\t"+ af + "\t0\n"); //6/25/08 gmagoon: added autoflag, needed when using daspkAUTO.exe; 080509 gmagoon: added sensitivity flag = 0
                bw.write(conversionSet[p_iterationNum]+"\n");

            }
            else{
                    bw.write(nState + "\t" + neq + "\t" +  -1 + "\t" +1+"\t"+ af + "\t0\n");//6/25/08 gmagoon: added autoflag, needed when using daspkAUTO.exe; 080509 gmagoon: added sensitivity flag = 0
                    bw.write(0+"\n");

            }
            bw.write( tBegin+" "+tEnd+"\n" );
                    for (int i=0; i<nState; i++)
                            bw.write(y[i]+" ");
                    bw.write("\n");
                    for (int i=0; i<nState; i++)
                            bw.write(yprime[i]+" ");
                    bw.write("\n");
                    for (int i=0; i<30; i++)
                            bw.write(info[i]+" ");
                    bw.write("\n"+ rtol + " "+atol);

                    bw.write("\n" + thermoString.toString() + "\n" + p_temperature.getK() + " " + p_pressure.getPa() + "\n" + rList.size() + "\n" + rString.toString() + "\n" + thirdBodyList.size() + "\n"+tbrString.toString() + "\n" + troeList.size() + "\n" + troeString.toString()+"\n" + lindemannList.size() + "\n" + lindemannString.toString() + "\n");
        }
        catch (IOException e) {
            System.err.println("Problem writing Solver Input File!");
            e.printStackTrace();
        }
		///4/30/08 gmagoon: code for providing edge reaction info to DASPK in cases if the automatic time stepping flag is set to true
		if (autoflag)
			getAutoEdgeReactionInfo((CoreEdgeReactionModel) p_reactionModel, p_temperature, p_pressure);
	
               // Add flags that specify whether the concentrations are constant or not
		getConcentrationFlags(p_reactionModel);
                //this should be the end of the input file
                try{
                    bw.flush();
                    bw.close();
                    fw.close();
                }
                catch (IOException e) {
                    System.err.println("Problem closing Solver Input File!");
                    e.printStackTrace();
		}
        int idid=0;
        LinkedHashMap speStatus = new LinkedHashMap();
        double [] senStatus = new double[nParameter*nState];
        
		int temp = 1;
        Global.solverPrepossesor = Global.solverPrepossesor + (System.currentTimeMillis() - startTime)/1000/60;
		
        startTime = System.currentTimeMillis();
        	//idid = solveDAE(p_initialization, reactionList, p_reactionChanged, thirdBodyReactionList, troeReactionList, nState, y, yprime, tBegin, tEnd, this.rtol, this.atol, T, P);
        idid = solveDAE();
        if (idid !=1 && idid != 2 && idid != 3)	{
        	System.out.println("The idid from DASPK was "+idid );
        	throw new DynamicSimulatorException("DASPK: SA off.");
        }
        System.out.println("After ODE: from " + String.valueOf(tBegin) + " SEC to " + String.valueOf(endTime) + "SEC");
        Global.solvertime = Global.solvertime + (System.currentTimeMillis() - startTime)/1000/60;
        startTime = System.currentTimeMillis();
        speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, 0);
        Global.speciesStatusGenerator = Global.speciesStatusGenerator + (System.currentTimeMillis() - startTime)/1000/60;
        
        
	SystemSnapshot sss = new SystemSnapshot(new ReactionTime(endTime, "sec"), speStatus, p_beginStatus.getTemperature(), p_beginStatus.getPressure());
	//sss.inertGas = p_beginStatus.inertGas; //gmagoon 6/23/09: copy inertGas information from initialStatus
        //total the concentrations of non-inert species
        double totalNonInertConc = totalNonInertConcentrations();
        //calculate the scale factor needed to account for volume change 
        double inertScaleFactor = 1;
        if(p_beginStatus.inertGas != null){//8/5/09 gmagoon: make sure inert gas is defined; otherwise, we will have issues with getTotalInertGas() when we start from non-zero time, as jdmo encountered (null-pointer exception)
            if(p_beginStatus.getTotalInertGas() > 0){//this check will ensure we don't try to divide by zero; otherwise, we will end up setting new concentration to be 0.0*Infinity=NaN
                inertScaleFactor = (p_beginStatus.getTotalMole()- totalNonInertConc)/p_beginStatus.getTotalInertGas();
            }
        
            //scale the initial concentrations of the inertGas to account for volume change
            for (Iterator iter = p_beginStatus.getInertGas(); iter.hasNext(); ) {
                String inertName = (String)iter.next();
                double originalInertConc = p_beginStatus.getInertGas(inertName);
                sss.putInertGas(inertName, originalInertConc*inertScaleFactor);       
            }
        }
        LinkedList reactionList = new LinkedList();
        reactionList.addAll(rList);
        reactionList.addAll(duplicates);
        reactionList.addAll(thirdBodyList);
        reactionList.addAll(troeList);
        reactionList.addAll(lindemannList);
        sss.setReactionList(reactionList);
        sss.setReactionFlux(reactionFlux);
		
        return sss;
        //#]
    }

	private int solveDAE() {
                String workingDirectory = System.getProperty("RMG.workingDirectory");
		
		// write the input file
//		File SolverInput = new File("ODESolver/SolverInput.dat");
//		try {
//			FileWriter fw = new FileWriter(SolverInput);
//			fw.write(outputString.toString());
//			fw.close();
//		} catch (IOException e) {
//			System.err.println("Problem writing Solver Input File!");
//			e.printStackTrace();
//		}
		
		// Rename RWORK and IWORK files if they exist
                renameIntermediateFilesBeforeRun();
		
		//run the solver on the input file
		boolean silentError = true;//start off assuming there is an error

                try {

                        String[] command = {workingDirectory +  "/bin/daspkAUTO.exe"};//5/5/08 gmagoon: changed to call daspkAUTO.exe
                                File runningDir = new File("ODESolver");

                                Process solver = Runtime.getRuntime().exec(command, null, runningDir);
                                InputStream is = solver.getInputStream();
                                InputStreamReader isr = new InputStreamReader(is);
                                BufferedReader br = new BufferedReader(isr);
                                String line=null;
                                while ( (line = br.readLine()) != null) {
                                        line = line.trim();
					silentError = false; //there is actual output from the ODE solver
                                        if (!(line.contains("ODESOLVER SUCCESSFUL"))) {
                                            System.err.println("Error running the ODESolver: "+line);
                                        }
                                }
				if(silentError){
				    System.err.println("Error: No stdout output from DASPK");
				}
                        int exitValue = solver.waitFor();
                }
                catch (Exception e) {
                        String err = "Error in running ODESolver \n";
                        err += e.toString();
                        e.printStackTrace();
                        System.exit(0);
                }

                //11/1/07 gmagoon: renaming RWORK and IWORK files
                renameIntermediateFilesAfterRun();
		return readOutputFile("ODESolver/SolverOutput.dat");
	}
        
        private void renameIntermediateFilesBeforeRun(){
                File f = new File("ODESolver/RWORK_"+index+".DAT");
		File newFile = new File("ODESolver/RWORK.DAT");
                boolean renameSuccess = false;
                if(f.exists()){
			if(newFile.exists())
				newFile.delete();
			renameSuccess = f.renameTo(newFile);
                        if (!renameSuccess)
                        {
                            System.out.println("Renaming of RWORK file(s) failed. renameIntermediateFilesBeforeRun()");
                            System.exit(0);
                        }
                }
                
                f = new File("ODESolver/IWORK_"+index+".DAT");
                newFile = new File("ODESolver/IWORK.DAT");
                if(f.exists()){
                    if(newFile.exists())
                            newFile.delete();
                    renameSuccess = f.renameTo(newFile);
                    if (!renameSuccess)
                    {
                        System.out.println("Renaming of IWORK file(s) failed. renameIntermediateFilesBeforeRun()");
                        System.exit(0);
                    }
                }
		
                f = new File("ODESolver/variables_"+index+".dat");
                newFile = new File("ODESolver/variables.dat");
                if(f.exists()){
                    if(newFile.exists())
                            newFile.delete();
                    renameSuccess = f.renameTo(newFile);
                    if (!renameSuccess)
                    {
                        System.out.println("Renaming of variables.dat file(s) failed. renameIntermediateFilesBeforeRun()");
                        System.exit(0);
                    }
                }
        }
        
	private void renameIntermediateFilesAfterRun() {
            File f = new File("ODESolver/RWORK.DAT");
            File newFile = new File("ODESolver/RWORK_"+index+".DAT");
            if(newFile.exists())
                newFile.delete();
            boolean renameSuccess = f.renameTo(newFile);
            if (!renameSuccess)
            {
                System.out.println("Renaming of RWORK file(s) failed. (renameIntermediateFilesAfterRun())");
                System.exit(0);
            }
            
            f = new File("ODESolver/IWORK.DAT");
            newFile = new File("ODESolver/IWORK_"+index+".DAT");
            if(newFile.exists())
                newFile.delete();
            renameSuccess = f.renameTo(newFile);
            if (!renameSuccess)
            {
                System.out.println("Renaming of IWORK file(s) failed. (renameIntermediateFilesAfterRun())");
                System.exit(0);
            }
            
            f = new File("ODESolver/variables.dat");
            newFile = new File("ODESolver/variables_"+index+".dat");
            if(newFile.exists())
                newFile.delete();
            renameSuccess = f.renameTo(newFile);
            if (!renameSuccess)
            {
                System.out.println("Renaming of variables.dat file(s) failed. (renameIntermediateFilesAfterRun())");
                System.exit(0);
            }
	}
        
	public int readOutputFile(String path) {
		
		//read the result
        File SolverOutput = new File(path);
        try {
        	FileReader fr = new FileReader(SolverOutput);
        	BufferedReader br = new BufferedReader(fr);
        	String line = br.readLine();
        	//StringTokenizer st = new StringTokenizer(line);
        	Global.solverIterations = Integer.parseInt(line.trim());
        	line = br.readLine();
        	if (Double.parseDouble(line.trim()) != neq) {
        		System.out.println("ODESolver didnt generate all species result");
        		System.exit(0);
        	}
        	endTime = Double.parseDouble(br.readLine().trim());
        	for (int i=0; i<nParameter+1; i++){
        		for (int j=0; j<nState; j++) {
        			line = br.readLine();
            		y[i*nState + j] = Double.parseDouble(line.trim());
        		}
        		line = br.readLine();
        	}
        	
        	for (int i=0; i<nParameter+1; i++){
        		for (int j=0; j<nState; j++) {
        			line = br.readLine();
            		yprime[i*nState + j] = Double.parseDouble(line.trim());
        		}
        		line = br.readLine();
        	}
        	reactionFlux = new double[rList.size()+thirdBodyList.size()+troeList.size()+lindemannList.size()];
        	for (int i=0; i<rList.size()+thirdBodyList.size()+troeList.size()+lindemannList.size(); i++){
        		line = br.readLine();
        		reactionFlux[i] = Double.parseDouble(line.trim());
        	}
        	for (int i=0; i<30; i++){
        		line = br.readLine();
        		info[i] = Integer.parseInt(line.trim());
        	}
		//for autoflag cases, there will be additional information which may be used for pruning
		if (autoflag){
		    prunableSpecies = new boolean[edgeID.size()];
		    maxEdgeFluxRatio = new double[edgeID.size()];
		    line=br.readLine();//read volume; (this is actually in the output even if AUTO is off, but is not used)
		    line=br.readLine();//read the edgeflag
		    Integer edgeflag = Integer.parseInt(line.trim());
		    if (edgeflag < 0){//if the edgeflag is negative, the ODE solver terminated by reaching the target time/concentration
			targetReached = true;
		    }
		    else{
			targetReached = false;
		    }
		    line=br.readLine();//read the time integrated to
		    double finalTime = Double.parseDouble(line.trim());
		    System.out.println("ODE solver integrated to "+ finalTime+" sec.");
		    for (int i=0; i<edgeID.size(); i++){//read the "prunability index" (0 or 1) and maximum ratio (edge flux/Rchar) for each edge species; note that edgeID only contains species, not P-dep networks, so we will not be reading in all the output from DASSL...only the flux ratio to actual edge species (vs. P-dep network pseudospecies)
			line = br.readLine().trim();//read the prunability index
			int q = Integer.parseInt(line);
			if(q > 0) prunableSpecies[i]=true; //q should be 1 or 0
			else prunableSpecies[i] = false;
			line = br.readLine().trim();//read the max edge flux ratio
			if(line.startsWith("+Inf")) maxEdgeFluxRatio[i]=Double.POSITIVE_INFINITY;
			else maxEdgeFluxRatio[i] = Double.parseDouble(line);
		    }
		}
        	
        }
        catch (IOException e) {
        	String err = "Error in reading Solver Output File! \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        
		return 1;
	}
	
	public LinkedList solveSEN(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, ReactionTime p_beginTime, ReactionTime p_endTime, Temperature p_temperature, Pressure p_pressure, boolean p_conditionChanged,TerminationTester tt) {
                setupInputFile();
	//	outputString = new StringBuilder();
		Iterator spe_iter = p_reactionModel.getSpecies();
                while (spe_iter.hasNext()){
                        Species spe = (Species)spe_iter.next();
                        int id = getRealID(spe);
                }
                double startTime = System.currentTimeMillis();

                        ReactionTime rt = p_beginStatus.getTime();
                if (!rt.equals(p_beginTime)) throw new InvalidBeginStatusException();

                double tBegin = p_beginTime.getStandardTime();
                double tEnd = p_endTime.getStandardTime();

		double T = p_temperature.getK();
		double P = p_pressure.getAtm();
		
		LinkedList initialSpecies = new LinkedList();
		
        // set reaction set
        //if (p_initialization || p_reactionChanged || p_conditionChanged) {
			
        	nState = p_reactionModel.getSpeciesNumber();
//        	troeString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10, alpha, Tstar, T2star, T3star, lowRate  (21 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0), ncollider, c1, c2,..c10, troe(0=T or 1=F) (21 elements)
        	troeString = generateTROEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
        	
			
//        	tbrString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10  (16 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0), ncollider, c1, c2,..c10 (20 elements)
			tbrString = generateThirdBodyReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			lindemannString = generateLindemannReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			//rString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, HASrev(T=1 or F=0)
			rString = generatePDepODEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			nParameter = 0;
			if (parameterInfor != 0) {
				nParameter = rList.size() + thirdBodyList.size() + troeList.size() + lindemannList.size() + p_reactionModel.getSpeciesNumber();				
			}
			neq = nState*(nParameter + 1);
			initializeWorkSpace();
			initializeConcentrations(p_beginStatus, p_reactionModel, p_beginTime, p_endTime, initialSpecies);
			
			
        //}
		int iterNum = 0;
                try{
                    if (tt instanceof ConversionTT){
                        SpeciesConversion sc = (SpeciesConversion)((ConversionTT)tt).speciesGoalConversionSet.get(0);
                        iterNum = conversionSet.length;
                        bw.write(nState + "\t" + neq + "\t" +  getRealID(sc.species) + "\t" +conversionSet.length+ "\t0\t1\n");//gmagoon 080509: added autoflag=0; later: added sensflag=1
                        for (int i=0; i<conversionSet.length; i++){
                                bw.write(conversionSet[i] + " ");
                        }
                        bw.write("\n");
                    }
                    else{
                        LinkedList timeSteps = ((ReactionTimeTT)tt).timeStep;
                        iterNum = timeSteps.size();
                        bw.write(nState + "\t" + neq + "\t" +  -1 + "\t" +timeSteps.size()+"\t0\t1\n");
                        for (int i=0; i<timeSteps.size(); i++){
                                bw.write(((ReactionTime)timeSteps.get(i)).time + " ");
                        }
                        bw.write("\n");
                     }
                    bw.write( tBegin+" "+tEnd+ "\n");
                    for (int i=0; i<nState; i++)
                            bw.write(y[i]+" ");
                    bw.write("\n");
                    for (int i=0; i<nState; i++)
                            bw.write(yprime[i]+" ");
                    bw.write("\n");
                    for (int i=0; i<30; i++)
                            bw.write(info[i]+" ");
                    bw.write("\n"+ rtol + " "+atol);

                    bw.write("\n" + thermoString.toString() + "\n" + p_temperature.getK() + " " + p_pressure.getPa() + "\n" + rList.size() + "\n" + rString.toString() + "\n" + thirdBodyList.size() + "\n"+tbrString.toString() + "\n" + troeList.size() + "\n" + troeString.toString()+"\n" + lindemannList.size() + "\n" + lindemannString.toString() + "\n");

                    // Add list of flags for constantConcentration
                    // one for each species, and a final one for the volume
                    // if 1:  will not change the number of moles of that species (or the volume)
                    // if 0:  will integrate the ODE as normal
                    // eg. liquid phase calculations with a constant concentration of O2 (the solubility limit - replenished from the gas phase)
                    // for normal use, this will be a sequence of '0 's
                    getConcentrationFlags(p_reactionModel);        
                }
                catch (IOException e) {
                    System.err.println("Problem writing Solver Input File!");
                    e.printStackTrace();
                }
       //this should be the end of the input file
        try{
            bw.close();
        }
        catch (IOException e) {
            System.err.println("Problem closing Solver Input File!");
            e.printStackTrace();
        }
        int idid=0;
        
        
		int temp = 1;
        Global.solverPrepossesor = Global.solverPrepossesor + (System.currentTimeMillis() - startTime)/1000/60;
		

        LinkedList systemSnapshotList = callSolverSEN(iterNum, p_reactionModel,  p_beginStatus);
        
        return systemSnapshotList;
        //#]
    }
	private LinkedList callSolverSEN(int p_numSteps, ReactionModel p_reactionModel, SystemSnapshot p_beginStatus) {
		double startTime = System.currentTimeMillis();
		String workingDirectory = System.getProperty("RMG.workingDirectory");
		
		LinkedList systemSnapshotList = new LinkedList();
		ReactionTime beginT = new ReactionTime(0.0, "sec");
		ReactionTime endT;
		//write the input file
	//	File SolverInput = new File("ODESolver/SolverInput.dat");
	//	try {
	//		FileWriter fw = new FileWriter(SolverInput);
	//		fw.write(outputString.toString());
	//		fw.close();
	//	} catch (IOException e) {
	//		System.err.println("Problem writing Solver Input File!");
        //			e.printStackTrace();
	//	}
		Global.writeSolverFile +=(System.currentTimeMillis()-startTime)/1000/60;
		//run the solver on the input file
		boolean error = false;
        try {
        	 // system call for therfit
        	String[] command = {workingDirectory +  "/bin/daspkAUTO.exe"};
			File runningDir = new File("ODESolver");
			
			Process ODESolver = Runtime.getRuntime().exec(command, null, runningDir);
			InputStream is = ODESolver.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				//System.out.println(line);
				line = line.trim();
				//if (!(line.contains("ODESOLVER SUCCESSFUL"))) {
					System.out.println(line);
					//error = true;
				//}
			}
        	int exitValue = 4;
        	exitValue = ODESolver.waitFor();
        	//System.out.println(br.readLine() + exitValue);
        	
        }
        catch (Exception e) {
        	String err = "Error in running ODESolver \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        
        startTime = System.currentTimeMillis();
        //read the result
        File SolverOutput = new File("ODESolver/SolverOutput.dat");
        try {
        	FileReader fr = new FileReader(SolverOutput);
        	BufferedReader br = new BufferedReader(fr);
        	String line ;
        	double presentTime = 0;
        	for (int k=0; k<p_numSteps; k++){
        		line = br.readLine();
                    if (Double.parseDouble(line.trim()) != neq) {
                            System.out.println("ODESolver didnt generate all species result");
                            System.exit(0);
                    }
                    presentTime = Double.parseDouble(br.readLine().trim());
                    endT = new ReactionTime(presentTime, "sec");
                    for (int i=0; i<nParameter+1; i++){
                            for (int j=0; j<nState; j++) {
                                    line = br.readLine();
                                    y[i*nState + j] = Double.parseDouble(line.trim());
                            }
                            line = br.readLine();//gmagoon 12/21/09: this apparently parses the line containing volume or volume sensitivity
                    }

                    for (int i=0; i<nParameter+1; i++){
                            for (int j=0; j<nState; j++) {
                                    line = br.readLine();
                                    yprime[i*nState + j] = Double.parseDouble(line.trim());
                            }
                            line = br.readLine();
                    }
                    reactionFlux = new double[rList.size()+thirdBodyList.size()+troeList.size()+lindemannList.size()];
                    for (int i=0; i<rList.size()+thirdBodyList.size()+troeList.size()+lindemannList.size(); i++){
                            line = br.readLine();
                            reactionFlux[i] = Double.parseDouble(line.trim());
                    }
                    LinkedHashMap speStatus = new LinkedHashMap();
                    double [] senStatus = new double[nParameter*nState];

                    System.out.println("After ODE: from " + String.valueOf(beginT.time) + " SEC to " + String.valueOf(endT.time) + "SEC");
                    speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, nParameter);
                    senStatus = generateSensitivityStatus(p_reactionModel,y,yprime,nParameter);
                    SystemSnapshot sss = new SystemSnapshot(endT, speStatus, senStatus, p_beginStatus.getTemperature(), p_beginStatus.getPressure());
                   // sss.inertGas = p_beginStatus.inertGas; //gmagoon 6/23/09: copy inertGas information from initialStatus
                    sss.inertGas = new LinkedHashMap();//zero out the inert gas info (in principal, this should already be null, but this is done "just-in-case")
                    //total the concentrations of non-inert species
                    double totalNonInertConc = totalNonInertConcentrations();
                    //calculate the scale factor needed to account for volume change
                    double inertScaleFactor = 1;
                    if(p_beginStatus.inertGas != null){//8/4/09 gmagoon: make sure inert gas is defined; otherwise, we will have issues with getTotalInertGas() when we start from non-zero time, as jdmo encountered (null-pointer exception)
                        if(p_beginStatus.getTotalInertGas() > 0){//this check will ensure we don't try to divide by zero; otherwise, we will end up setting new concentration to be 0.0*Infinity=NaN
                            inertScaleFactor = (p_beginStatus.getTotalMole()- totalNonInertConc)/p_beginStatus.getTotalInertGas();
                        }

                        //scale the initial concentrations of the inertGas to account for volume change
                        for (Iterator iter = p_beginStatus.getInertGas(); iter.hasNext(); ) {
                            String inertName = (String)iter.next();
                            double originalInertConc = p_beginStatus.getInertGas(inertName);
                            sss.putInertGas(inertName, originalInertConc*inertScaleFactor);       
                        }
                    }

                    sss.setIDTranslator(IDTranslator);
                    LinkedList reactionList = new LinkedList();
                    reactionList.addAll(rList);
                    reactionList.addAll(duplicates);
                    reactionList.addAll(thirdBodyList);
                    reactionList.addAll(troeList);
                    reactionList.addAll(lindemannList);
                    sss.setReactionList(reactionList);
                    systemSnapshotList.add(sss);
                    sss.setReactionFlux(reactionFlux);
                    beginT = endT;
                    //tEnd = tEnd.add(tStep);
        	}
        	
        	
        }
        catch (IOException e) {
        	String err = "Error in reading Solver Output File! \n";
        	err += e.toString();
        	e.printStackTrace();
        	System.exit(0);
        }
        Global.readSolverFile += (System.currentTimeMillis() - startTime)/1000/60;
		return systemSnapshotList;
	}
	
	
	@Override
	protected void initializeWorkSpace() {
		super.initializeWorkSpace();
		info[4] = 1; //use analytical jacobian
		if (nParameter != 0) {
			info[18] = nParameter; //the number of parameters
			info[19] = 2; //perform senstivity analysis
			info[24] = 1;//staggered corrector method is used
		}
	}
	
	@Override
	protected void initializeConcentrations(SystemSnapshot p_beginStatus, ReactionModel p_reactionModel, ReactionTime p_beginTime, ReactionTime p_endTime, LinkedList initialSpecies) {
    	
		super.initializeConcentrations(p_beginStatus, p_reactionModel,
				p_beginTime, p_endTime, initialSpecies);

		if (nParameter != 0){//svp
						
			double [] sensitivityStatus = new double[nState*nParameter];
			int speciesNumber = p_reactionModel.getSpeciesNumber();
			
			for (int i=0; i<nParameter*speciesNumber;i++){
				sensitivityStatus[i] = 0;
			}
			p_beginStatus.addSensitivity(sensitivityStatus);							
		}

	}

}

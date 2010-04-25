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


import java.awt.image.IndexColorModel;
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
import jing.rxn.PDepNetwork;
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

//## class JDASSL
public class JDASSL extends JDAS {

	private JDASSL() {
        super();
    }
    
	public JDASSL(double p_rtol, double p_atol, int p_parameterInfor, 
			InitialStatus p_initialStatus, int p_index, ValidityTester p_vt, 
			boolean p_autoflag, Double p_termTol, Double p_coreTol) {
        super(p_rtol, p_atol, p_parameterInfor, p_initialStatus, p_index, p_vt, 
			p_autoflag, p_termTol, p_coreTol);
    }

    // DASSL is not used for sensitivity analysis, so this function shouldn't be needed
	/*private LinkedList generateSensitivityStatus(ReactionModel p_reactionModel, double [] p_y, double [] p_yprime, int p_paraNum) {
    	int neq = p_reactionModel.getSpeciesNumber()*(p_paraNum+1);
    	if (p_y.length != neq) throw new DynamicSimulatorException();
    	if (p_yprime.length != neq) throw new DynamicSimulatorException();

    	LinkedList senStatus = new LinkedList();

    	if (p_paraNum > 0){
    		for (int i = p_reactionModel.getSpeciesNumber();i<neq;i++){
    			double sens = p_y[i];
    			double sflux = p_yprime[i];
    			int reaction_num = i/p_reactionModel.getSpeciesNumber();
    			int species_num = (i+1)%p_reactionModel.getSpeciesNumber();
    			if (species_num == 0){
    				species_num = p_reactionModel.getSpeciesNumber();
    			}
    			//String name = "dC"+String.valueOf(species_num)+"/dk"+String.valueOf(reaction_num);
    			//System.out.println(name + '\t' + String.valueOf(sens) + '\t' + String.valueOf(sflux));
    				SensitivityStatus ss = new SensitivityStatus(sens, sflux, species_num, reaction_num);
    				int index = i-p_reactionModel.getSpeciesNumber();
    				senStatus.add(ss);
    		}
    	}
    	return senStatus;
    }*/	 

    //## operation solve(boolean,ReactionModel,boolean,SystemSnapshot,ReactionTime,ReactionTime,Temperature,Pressure,boolean)
    public SystemSnapshot solve(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, ReactionTime p_beginTime, ReactionTime p_endTime, Temperature p_temperature, Pressure p_pressure, boolean p_conditionChanged, TerminationTester tt, int p_iterationNum) {
        
        // set up the input file
        setupInputFile();
        //outputString = new StringBuilder();
        //first generate an id for all the species
    	Iterator spe_iter = p_reactionModel.getSpecies();
    	while (spe_iter.hasNext()){ //gmagoon 021810; this section had been commented out since Jun 2009, but I am restoring it in order to have IDs generated for all the core species; later on, in transferReaction (called from, e.g. generateLindemannReactionList) we may want to check whether a species is in the core, and we do this through IDTranslator, which is the only place we can relatively easily get this information; without this fix, certain colliders may not be considered for Lindemann, Troe, and third-body reactions
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
			nParameter = 0;
			if (parameterInfor != 0) {
				nParameter = p_reactionModel.getReactionNumber(); //svp
				if (initialStatus == null) System.out.println("initialStatus = null");
				spe_iter = initialStatus.getSpeciesStatus();
				while (spe_iter.hasNext()) {
					SpeciesStatus ss = (SpeciesStatus) spe_iter.next();
					String name = ss.getSpecies().getName();
					initialSpecies.add(name);
					nParameter++;
				}
			}
			neq = nState*(nParameter + 1);
			
//			tbrString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10  (16 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, rev(=1 or -1), ncollider, c1, c2,..c10 (20 elements)
			tbrString = generateThirdBodyReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			//troeString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq, inertEfficiency, e1, e2, ..., e10, alpha, Tstar, T2star, T3star, lowRate  (21 elements)
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, rev(=1 or -1), ncollider, c1, c2,..c10, troe(0=T or 1=F) (21 elements)
        	troeString = generateTROEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
        	lindemannString = generateLindemannReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
        	
			//rString is a combination of a integer and a real array
			//real array format:  rate, A, n, Ea, Keq
			//int array format :  nReac, nProd, r1, r2, r3, p1, p2, p3, rev(=1 or -1)
			rString = generatePDepODEReactionList(p_reactionModel, p_beginStatus, p_temperature, p_pressure);
			
			initializeWorkSpace();
			initializeConcentrations(p_beginStatus, p_reactionModel, p_beginTime, p_endTime, initialSpecies);
			
			
        }
        else
        	info[0] = 1;
        //5/5/08 gmagoon: (next two lines) for autoflag, get binary 1 or 0 corresponding to boolean true/false
        int af = 0;
        if (autoflag) af = 1;
        try{
            if (tt instanceof ConversionTT){
                    SpeciesConversion sc = (SpeciesConversion)((ConversionTT)tt).speciesGoalConversionSet.get(0);    
                    bw.write(nState + "\t" + neq + "\t" +  getRealID(sc.species) + "\t" +conversionSet[p_iterationNum]+"\t"+ af+"\n");//5/5/08 gmagoon: added autoflag, needed when using dasslAUTO.exe

            }
            else{
                    bw.write(nState + "\t" + neq + "\t" +  -1 + "\t" +0+ "\t"+ af+"\n");//5/5/08 gmagoon: added autoflag, needed when using dasslAUTO.exe

            }
            for (int i=0; i<neq; i++)
                    bw.write(y[i]+" ");
            bw.write("\n");
            for (int i=0; i<neq; i++)
                    bw.write(yprime[i]+" ");
            bw.write("\n" + tBegin+" "+tEnd+"\n");
            for (int i=0; i<30; i++)
                    bw.write(info[i]+" ");
            bw.write("\n"+ rtol + " "+atol);
            bw.write("\n" + p_temperature.getK() + " " + p_pressure.getPa() + "\n" + rList.size() + "\n" + rString.toString() + "\n" + thirdBodyList.size() + "\n"+tbrString.toString() + "\n" + troeList.size() + "\n" + troeString.toString()+"\n" + lindemannList.size() + "\n" + lindemannString.toString() + "\n");
        }
        catch (IOException e) {
        	System.err.println("Problem writing Solver Input File!");
                e.printStackTrace();
        }
        //4/30/08 gmagoon: code for providing edge reaction info to DASSL in cases if the automatic time stepping flag is set to true
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
        LinkedList senStatus = new LinkedList();
        
		int temp = 1;
        Global.solverPrepossesor = Global.solverPrepossesor + (System.currentTimeMillis() - startTime)/1000/60;
		if (nParameter==0) {
			startTime = System.currentTimeMillis();
        	
			idid = solveDAE();
                        
 //                       createDotGraphs((CoreEdgeReactionModel)p_reactionModel);

                        //gmagoon 1/25/09: note that the below lines do not actually use the actual idid, but instead use the 1 returned by readOutputFile; the actual success represented by idid is effectively read in through the Fortran output "ODESOLVER SUCCESSFUL" vs. "ODESOLVER FAILED"
			if (idid !=1 && idid != 2 && idid != 3)	{
				System.out.println("The idid from DASSL was "+idid );
				throw new DynamicSimulatorException("DASSL");
        	}
            System.out.println("After ODE: from " + String.valueOf(tBegin) + " SEC to " + String.valueOf(endTime) + "SEC");
			Global.solvertime = Global.solvertime + (System.currentTimeMillis() - startTime)/1000/60;
			startTime = System.currentTimeMillis();
        	speStatus = generateSpeciesStatus(p_reactionModel, y, yprime, 0);
			Global.speciesStatusGenerator = Global.speciesStatusGenerator + (System.currentTimeMillis() - startTime)/1000/60;
        }
        

	SystemSnapshot sss = new SystemSnapshot(new ReactionTime(endTime,"sec"), speStatus, p_beginStatus.getTemperature(), p_beginStatus.getPressure());
	//sss.inertGas = p_beginStatus.inertGas; //gmagoon 6/23/09: copy inertGas information from initialStatus
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
        
        LinkedList reactionList = new LinkedList();
	reactionList.addAll(rList);
        reactionList.addAll(duplicates);
        reactionList.addAll(thirdBodyList);
        reactionList.addAll(troeList);
        reactionList.addAll(lindemannList);
        sss.setReactionList(reactionList);
	sss.setReactionFlux(reactionFlux);
                
        return sss;
    }

	private int solveDAE() {	
	//	super.solveDAE("dasslAUTO.exe");
            	String workingDirectory = System.getProperty("RMG.workingDirectory");
		
//		// write the input file
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
		boolean error = false;
                try {

                        String[] command = {workingDirectory +  "/bin/dasslAUTO.exe"};//5/5/08 gmagoon: changed to call dasslAUTO.exe
                                File runningDir = new File("ODESolver");

                                Process solver = Runtime.getRuntime().exec(command, null, runningDir);
                                InputStream is = solver.getInputStream();
                                InputStreamReader isr = new InputStreamReader(is);
                                BufferedReader br = new BufferedReader(isr);
                                String line=null;
                                while ( (line = br.readLine()) != null) {
                                        line = line.trim();
                                        if (!(line.contains("ODESOLVER SUCCESSFUL"))) {
                                                System.err.println("Error running the ODESolver: "+line);
                                                error = true;
                                        }          
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
                            System.out.println("Renaming of RWORK file(s) failed.");
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
                        System.out.println("Renaming of IWORK file(s) failed.");
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
                System.out.println("Renaming of RWORK file(s) failed. (renameIntermediateFiles())");
                System.exit(0);
            }
            
            f = new File("ODESolver/IWORK.DAT");
            newFile = new File("ODESolver/IWORK_"+index+".DAT");
            if(newFile.exists())
                newFile.delete();
            renameSuccess = f.renameTo(newFile);
            if (!renameSuccess)
            {
                System.out.println("Renaming of IWORK file(s) failed. (renameIntermediateFiles())");
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
        	Global.solverIterations = Integer.parseInt(line.trim());
        	line = br.readLine();
        	
        	if (Double.parseDouble(line.trim()) != neq) {
        		System.out.println("ODESolver didnt generate all species result");
        		System.exit(0);
        	}
        	endTime = Double.parseDouble(br.readLine().trim());
        	for (int i=0; i<neq; i++){
        		line = br.readLine();
        		y[i] = Double.parseDouble(line.trim());
        	}
        	
        	for (int i=0; i<neq; i++){
        		line = br.readLine();
        		yprime[i] = Double.parseDouble(line.trim());
        	}
        	reactionFlux = new double[rList.size()+thirdBodyList.size()+troeList.size()+lindemannList.size()];
        	for (int i=0; i<rList.size()+thirdBodyList.size()+troeList.size()+lindemannList.size(); i++){
        		line = br.readLine();
        		reactionFlux[i] = Double.parseDouble(line.trim());
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
        SolverOutput.delete();

		return 1;
	}
	
	@Override
	protected void initializeWorkSpace() {
		super.initializeWorkSpace();
		//info[4] = 1; //use analytical jacobian
	}
        
//        protected void createDotGraphs(CoreEdgeReactionModel cerm) {
//            // 1. determine the core-core connectivity (it is assumed all are "regular" non-P-dep reactions
//            
//            int n = IDTranslator.size();
//            int [][] matrix = new int[n][n];
//            for (int i = 0; i<n; i++){
//               for (int j = 0; j<n; j++){
//                    matrix[i][j]=0;
//                } 
//            }
//            
//            Iterator riter = cerm.getReaction();
//            while (riter.hasNext()){
//                //below is based on transferReaction()
//                    Reaction p_reaction = (Reaction)riter.next();
//                    int rnum = p_reaction.getReactantNumber();
//                    int pnum = p_reaction.getProductNumber();
//
//                    int [] rid = new int[rnum];
//                    int index = 0;
//                    for (Iterator r_iter = p_reaction.getReactants(); r_iter.hasNext(); ) {
//                            Species s = (Species)r_iter.next();
//                            rid[index] = getRealID(s);
//                            index++;
//                    }
//
//                    int [] pid = new int[pnum];
//                    index = 0;
//                    for (Iterator p_iter = p_reaction.getProducts(); p_iter.hasNext(); ) {
//                                    Species s = (Species)p_iter.next();
//                            pid[index] = getRealID(s);
//                            index++;
//                    }
//                    
//                   for(int i = 0; i<rid.length; i++){
//                       for(int j = 0; j < pid.length; j++){
//                           matrix[rid[i]-1][pid[j]-1]=1;
//                       }
//                   }
//            }
//            //now we want to "symmetrize" the matrix so we only have one arrow connecting each pair of species; we will only construct upper half
//           // int [][] matrixTr = new int[n][n];
//           // for (int i = 0; i<n; i++){
//           //    for (int j = 0; j<n; j++){
//           //         matrixTr[i][j]=matrix[j][i];
//           //     } 
//           // }
//            
//            int [][] matrixNew = new int[n][n];
//            for (int i = 0; i<n; i++){
//               for (int j = 0; j<n; j++){
//                    matrixNew[i][j]=0;
//                    if(j<i) matrixNew[i][j]=matrix[j][i] + matrix[i][j];
//                } 
//            }
//
//            matrix = matrixNew;
//            
//            // 2. determine the core-edge connectivity
//            int [][] matrixE = null;
//            int [][] rid = null;
//            int [][] pid = null;
//            int [] nreac = null;
//            int nEdge = 0;
//            int nEdgeReacs = 0;
//            String path = "ODESolver/SolverInput.dat";
//            File SolverInput = new File(path);
//            //String path2 = "ODESolver/EdgeReacFlux.dat";
//            //File erf = new File(path2);
//            try {
//        	FileReader fr = new FileReader(SolverInput);
//        	BufferedReader br = new BufferedReader(fr);
//              //  FileReader fr2 = new FileReader(erf);
//        //	BufferedReader br2 = new BufferedReader(fr2);
//                
//                String line = br.readLine();
//                for(int i = 0; i < 21; i++){
//                    line = br.readLine();
//                }
//                //line with nEdgeSpecies, nEdgeReactions
//                String [] split = line.split(" ");
//                nEdge = Integer.parseInt(split[0]);
//                nEdgeReacs = Integer.parseInt(split[1]);
//               
//                matrixE = new int[n][nEdge];
//                for (int i = 0; i<n; i++){
//                    for (int j = 0; j<nEdge; j++){
//                        matrixE[i][j]=0;
//                    } 
//                }
//                
//                rid = new int[nEdgeReacs][3];
//                pid = new int[nEdgeReacs][3];
//                nreac = new int[nEdgeReacs];
//                
//                line = br.readLine();
//                for(int m = 0; m < nEdgeReacs; m++){
//                    line = line.trim();
//                    split = line.split(" ");
//                    nreac[m] = Integer.parseInt(split[0]);
//                    rid[m][0] = Integer.parseInt(split[2]);
//                    rid[m][1] = Integer.parseInt(split[3]);
//                    rid[m][2] = Integer.parseInt(split[4]);
//                    pid[m][0] = Integer.parseInt(split[5]);
//                    pid[m][1] = Integer.parseInt(split[6]);
//                    pid[m][2] = Integer.parseInt(split[7]);            
//                    for(int i = 0; i < 3; i++){
//                       for(int j = 0; j < 3; j++){                           
//                           if(rid[m][i]>0 && pid[m][j]> 0) matrixE[rid[m][i]-1][pid[m][j]-1]=1;
//                       }
//                    }
//    
//                    line = br.readLine();
//                }
//                fr.close();
//                
//            }
//            catch (IOException e) {
//        	String err = "Error in reading Solver Input File! \n";
//        	err += e.toString();
//        	e.printStackTrace();
//        	System.exit(0);
//            }
//            
//            // 3. read in core conc. vs. time
//            path = "ODESolver/CoreConc.dat";
//            File f = new File(path);
//            LinkedList tvec = new LinkedList();
//            LinkedList concvec = new LinkedList();//list of arrays of concentration
//            double [] concvecElement = new double [n];
//            try {
//        	FileReader fr = new FileReader(f);
//        	BufferedReader br = new BufferedReader(fr);
//                
//                int spe = 0;
//                double t = 0;
//                double conc = 0;
//                
//                String line = br.readLine();
//
//                while(line!=null){
//                    line = line.trim();
//                    String [] split = line.split(" ");
//                    spe = Integer.parseInt(split[0]);
//                    t = Double.parseDouble(split[1]);
//                    conc = Double.parseDouble(split[2]);
//                    
//                    if (spe==1){
//                        tvec.add(t);
//                        if(tvec.size()>1) concvec.add(concvecElement.clone()); //don't add until the second time
//                    }
//                    concvecElement[spe-1]=conc;
//                    
//                    line = br.readLine();
//                }
//                concvec.add(concvecElement.clone());
//                
//                fr.close();
//            }
//            catch (IOException e) {
//        	String err = "Error in reading core conc. File! \n";
//        	err += e.toString();
//        	e.printStackTrace();
//        	System.exit(0);
//            }
//            
//            
//            // 4. read in edge fluxes vs. time
//            path = "ODESolver/EdgeFlux.dat";
//            f = new File(path);
//            LinkedList fluxvec = new LinkedList();//list of arrays of concentration
//            double [] fluxvecElement = new double [nEdge];
//            try {
//        	FileReader fr = new FileReader(f);
//        	BufferedReader br = new BufferedReader(fr);
//                
//                int spe = 0;
//                double t = 0;
//                double flux = 0;
//                
//                String line = br.readLine();
//
//                while(line!=null){
//                    line = line.trim();
//                    String [] split = line.split(" ");
//                    spe = Integer.parseInt(split[0]);
//                    t = Double.parseDouble(split[1]);
//                    flux = Double.parseDouble(split[2]);
//                    
//                    if (spe==1){
//                        if(t>0) fluxvec.add(fluxvecElement.clone()); //don't add until the second time
//                    }
//                    fluxvecElement[spe-1]=flux;
//                    
//                    line = br.readLine();
//                }
//                fluxvec.add(fluxvecElement.clone());
//                
//                fr.close();
//            }
//            catch (IOException e) {
//        	String err = "Error in reading edge flux File! \n";
//        	err += e.toString();
//        	e.printStackTrace();
//        	System.exit(0);
//            }
//            
//            //4a. read in the edge reaction fluxes
//            path = "ODESolver/EdgeReacFlux.dat";
//            f = new File(path);
//            LinkedList reacfluxvec = new LinkedList();//list of arrays of edge reaction flux
//            double [] rfvElement = new double [nEdgeReacs];
//            try {
//        	FileReader fr = new FileReader(f);
//        	BufferedReader br = new BufferedReader(fr);
//                
//                int reac = 0;
//                double t = 0;
//                double flux = 0;
//                
//                String line = br.readLine();
//
//                while(line!=null){
//                    line = line.trim();
//                    String [] split = line.split(" ");
//                    reac = Integer.parseInt(split[0]);
//                    t = Double.parseDouble(split[1]);
//                    flux = Double.parseDouble(split[2]);
//                    
//                    if (reac==1){
//                        if(t>0) reacfluxvec.add(rfvElement.clone()); //don't add until the second time
//                    }
//                    rfvElement[reac-1]=flux;
//                    
//                    line = br.readLine();
//                }
//                reacfluxvec.add(rfvElement.clone());
//                
//                fr.close();
//            }
//            catch (IOException e) {
//        	String err = "Error in reading edge reaction flux File! \n";
//        	err += e.toString();
//        	e.printStackTrace();
//        	System.exit(0);
//            }
//            
//            //4aa. make a new matrix with reaction flux weights
//            
//
//            LinkedList matrixE2vec = new LinkedList();//list of arrays of reaction flux
//            for(int m = 0; m < tvec.size(); m++){
//                double[][] matrixE2 = new double[n][nEdge];
//                for (int i = 0; i<n; i++){//zero the matrix
//                    for (int j = 0; j<nEdge; j++){
//                        matrixE2[i][j]=0;
//                    } 
//                }
//                double[] rf = (double[])reacfluxvec.get(m);
//                double[] fl = (double[])fluxvec.get(m);
//                for(int p = 0; p < nEdgeReacs; p++){
//                    for(int i = 0; i<3; i++){
//                       for(int j = 0; j < 3; j++){                           
//                           if(rid[p][i]>0 && pid[p][j]> 0 && fl[pid[p][j]-1] > 0) matrixE2[rid[p][i]-1][pid[p][j]-1]+=rf[p]/(fl[pid[p][j]-1]*nreac[p]);//includes an "optional" division by number of reactants to make arrow sizes look more "realistic"; in the end, this should give the fraction of the total edge flux from core species
//                       }
//                    }
//                }
//                matrixE2vec.add(matrixE2);//for some reason, cloning approach used before with matrix initialized outside loop doesn't seem to work for this case (since it is a 2D array?)
//            }
//            
//            // 5. write the dot file
//            //write a dot file for each time
//            for(int m = 0; m < tvec.size(); m++){
//                double t = (Double)tvec.get(m);
//                double [] fluxvecEl = (double[])fluxvec.get(m);
//                double [][] matrixE2El = (double[][])matrixE2vec.get(m);
//                double [] concvecEl = (double[])concvec.get(m); 
//                // write the input file
//                String num = "";
//                if(m<10) num = "0000"+m;
//                else if(m<100) num = "000"+m;
//                else if(m<1000) num = "00"+m;
//                else if(m<10000) num = "0"+m;
//                else if(m<100000) num = ""+m;
//                File dotFile = new File("ODESolver/dot"+num+".dot");
//                try {
//                        FileWriter fw = new FileWriter(dotFile);
//                        
//                        
//                        fw.write("digraph reaction_network {\n");
//                        
//                        //write the core connectivity
//                        for (int i = 0; i<n; i++){
//                            for (int j = 0; j<n; j++){
//                                if(matrix[i][j] > 0) fw.write((i+1) + " -> "+ (j+1)+ "[fontname=\"Arial\", style=\"setlinewidth(1.5)\", arrowsize=0.00 color=\"black\"];\n");//make label size depend on log of concentration
//                                //else fw.write((i+1) + " -> "+ (j+1)+ "[fontname=\"Arial\", style=\"setlinewidth(0.0)\", arrowsize=0.00 color=\"black\", len=5.0];\n");//use invisible lines
//                            } 
//                        }
//                        //write the core-edge connectivity
//                        for (int i = 0; i<n; i++){
//                            for (int j = 0; j<nEdge; j++){
//                                if(matrixE[i][j]==1) fw.write((i+1) + " -> E"+ (j+1)+ "a[fontname=\"Arial\", style=\"setlinewidth("+(.10+6*(Math.pow(fluxvecEl[j]*matrixE2El[i][j], 0.25)))+")\", arrowsize=0.00 color=\"red\"];\n");
//                            } 
//                        }
//                        //write the (time-dependent) arrows to the edge
//                        for (int j = 0; j<nEdge; j++){
//                            fw.write("E"+ (j+1)+ "a -> E" + (j+1)+"[fontname=\"Arial\", style=\"setlinewidth("+(.10+6*(Math.pow(fluxvecEl[j], 0.25)))+")\", arrowsize=" +Math.sqrt(.1+6*fluxvecEl[j])/2 +  " color=\"red\"];\n");
//                        }
//                        //write the core nodes
//                        for (int i = 0; i<n; i++){
//                            //find the corresponding species:
//                            Iterator iter = IDTranslator.keySet().iterator();
//                            int found = 0;
//                            Species spe = null;
//                            while (iter.hasNext()&&found==0) {
//                                spe = (Species)iter.next();
//                                if((Integer)(IDTranslator.get(spe))==i+1) found = 1;
//                            }
//                            String name = spe.getName();
//                                  
//                            double fontSize = 10;
//                            if(concvecEl[i]>1E-24) fontSize = 3*(Math.log10(concvecEl[i])+24)+10.0;//make label size depend on log of concentration
//                            fw.write((i+1)+"[fontsize="+fontSize/1.5+", width="+(.02+Math.sqrt((fontSize-10)/100))+" height="+(.02+Math.sqrt((fontSize-10)/100))+ ", shape=\"rectangle\", fontcolor = \"black\", label=\""+name.replace("JJ", ":").replace("J", "&bull;")+"\"];\n");//make label size depend on log of concentration
//                        }
//                        //make the ancillary edge nodes invisible
//                        for (int j = 0; j<nEdge; j++){
//                            //find the corresponding species:
//                            Iterator iter = edgeID.keySet().iterator();
//                            int found = 0;
//                            Species spe = null;
//                            while (iter.hasNext()&&found==0) {
//                                spe = (Species)iter.next();
//                                if((Integer)(edgeID.get(spe))==j+1) found = 1;
//                            }
//                            String name = spe.getName();
//                            
//                            String shape = "shape=\"none\",";
//                          //  if(fluxvecEl[j] > 1.0) shape="shape=\"circle\", color=\"red\",";//circle the edge node when it exceeds the threshhold
//                            fw.write("E"+(j+1)+"a[shape=\"point\", width=0.0, height=0.0, label=\"\", color=\"red\"];\n");//make the ancillary edge nodes invisible
//                            if(fluxvecEl[j] > 1.0) fw.write("E"+(j+1)+"[fontsize=8, fontcolor=\"red\", "+shape+" width=0.0, height=0.0, label=\""+name.replace("JJ", ":").replace("J", "&bull;")+"\"];\n");
//                            else fw.write("E"+(j+1)+"[fontsize=8, "+shape+" width=0.0, height=0.0, label=\""+name.replace("JJ", ":").replace("J", "&bull;")+"\"];\n");
//                        }
//                        //NumberFormat form;
//                        //form.setMinFractionDigits(2);
//                        //form.setMaxFractionDigits(2);
//                        //String str = form.format(t);
//                        //fw.write("label = \"t = " + t + " sec.\";\n");
//                        String str = String.format("label = \"t = %3.2G sec.\";\n", t);
//                        fw.write(str);
//                        //fw.write("overlap=prism;\n");
//                        //fw.write("model=circuit;\n");
//                        fw.write("}\n");
//                        fw.close();
//                        String command = "\"c:\\Program Files\\GraphViz2.24\\bin\\\\twopi.exe\" -Tgif -O " + "dot"+num+".dot";
//                        File runningDir = new File("ODESolver");
//                        Process dot = Runtime.getRuntime().exec(command, null, runningDir);
//                        int ev = dot.waitFor();               
//                } catch (Exception e) {
//                        System.err.println("Problem writing dot File!");
//                        e.printStackTrace();
//                }
//
//            }
//            
//
//            
//            
//        }
	
}

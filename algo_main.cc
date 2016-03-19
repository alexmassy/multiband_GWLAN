#include "gwlan_lib.hpp"

//Input: input_filename - model_filename - output_filename
//Output: output_file 

///Upper bound on the number of clock ticks that may be spent for solving a model.
const clock_t tlimit = 7200;

int main (int argc, char *argv[]) {

 if (argc!=4)
 {
   cout << "Wrong number of input!" << endl;
   exit(1);
 }
 else {

  string input_filename, output_filename;
  IloInt verbose_flag;

  input_filename = argv[1];
  //cout << "Input file: " << input_filename << endl;

  output_filename = argv[2];
  //cout << "Output file: " << output_filename << endl;

  verbose_flag = atoi(argv[3]);
  //cout << "Verbose flag: " << verbose_flag << endl;
 

  /********** Check if input files exist ***********/
  ifstream infile;

  infile.open (input_filename.c_str(), ifstream::in);
  if(!infile.is_open()){
    cout << "Unable to open file: "<<input_filename<<endl<<endl;		
    return 0;
  }
  infile.close();	


  /**************** Read selected input data ***************/
  IloEnv simEnv;

  IloInt K;
  IloNumArray H(simEnv);
  IloNumArray L(simEnv);
  IloNum rho;
  vector<string> cs_vect;
  vector<string> ut_vect;
  IloNumArray w(simEnv);
  IloNumArray b(simEnv);
  IloNumArray p_w(simEnv);
  IloNumArray2D r_curr(simEnv);
  IloNumArray3D r_min(simEnv);
  IloNumArray2D r_future(simEnv);

  //cout<<"Reading input.\n";
  
  readData(simEnv, input_filename, K, H, L, rho, cs_vect, ut_vect, w, b, p_w, r_curr, r_min, r_future);

  /********** Algorithm ***********/

  IloInt i, j, k;
  IloInt numAPs = cs_vect.size();
  IloInt numUTs = ut_vect.size();
  IloInt numBands = H.getSize();  
  
  try {
 

    /*----------------------------------------- Create multiband -------------------------------------*/
    IloEnv multibAlgEnv;
    IloModel multibAlgMod(multibAlgEnv, "multiband_gwlan");
    assignVarMx x(multibAlgEnv, numUTs);//no bouds specified
    IloIntVarArray y(multibAlgEnv, numAPs, 0, 1); //bounds already set (constraint 41) 

    createMultibMILP(multibAlgMod, x, y, b, p_w, rho, w, r_curr, r_min, K, H, L);
 
    IloCplex multibAlgCplex(multibAlgMod);
    
    // Silent Mode
    multibAlgCplex.setOut(multibAlgEnv.getNullStream());
    multibAlgCplex.setWarning(multibAlgEnv.getNullStream());
    
    // Set time limit
    multibAlgCplex.setParam(IloCplex::ClockType, 1); //0: automatic - 1: CPU time - 2: wall-clock
    multibAlgCplex.setParam(IloCplex::TiLim, tlimit);
    
    // Set mono-threading
    multibAlgCplex.setParam(IloCplex::Threads, 1);
	
    
    /*------------------------------------ Create nominal present---------------------------------------*/
    IloEnv algEnv;
    IloModel algMod(algEnv, "nominal_gwlan");
    assignVarMx nr_x(algEnv, numUTs);
    IloIntVarArray nr_y(algEnv, numAPs, 0, 1);

    createAlgILP(algMod, nr_x, nr_y, b, p_w, rho, w, r_curr);
    
    IloCplex algCplex(algMod);
    
    // Silent Mode
    algCplex.setOut(algEnv.getNullStream());
    algCplex.setWarning(algEnv.getNullStream());
    // Set time limit
    algCplex.setParam(IloCplex::ClockType, 1); //0: automatic - 1: CPU time - 2: wall-clock
    algCplex.setParam(IloCplex::TiLim, tlimit);
    // Set mono-threading
    algCplex.setParam(IloCplex::Threads, 1);
    
    
    /*------------------------------------- create nominal future -------------------------------------*/
    IloEnv futAlgEnv;
    IloModel futAlgMod(futAlgEnv, "future_gwlan");
    assignVarMx fut_x(futAlgEnv, numUTs);
    IloIntVarArray fut_y(futAlgEnv, numAPs, 0, 1);

    createAlgILP(futAlgMod, fut_x, fut_y, b, p_w, rho, w, r_future);
    
    IloCplex futAlgCplex(futAlgMod);
    // Silent Mode
    futAlgCplex.setOut(futAlgEnv.getNullStream());
    futAlgCplex.setWarning(futAlgEnv.getNullStream());
    // Set time limit
    futAlgCplex.setParam(IloCplex::ClockType, 1); //0: automatic - 1: CPU time - 2: wall-clock
    futAlgCplex.setParam(IloCplex::TiLim, tlimit);
    // Set mono-threading
    futAlgCplex.setParam(IloCplex::Threads, 1);
	
    
    
    // Run multiband present
    if(verbose_flag > 0) cout<<"***Solving multiband problem."<<endl;
    clock_t tstart = clock();
    multibAlgCplex.solve();     
    double multibTime = (double) (clock() - tstart) / CLOCKS_PER_SEC;
    if( multibAlgCplex.getStatus() == IloAlgorithm::Infeasible ) {
      cout<<"Algorithm ended: infeasible problem.\n";
      return 0;
    }
    else{
      if( multibAlgCplex.getStatus() == IloAlgorithm::Optimal ) {
        //cout<<"Algorithm ended: optimal solution found.\n";
      }
      else{
        cout<<"Algorithm ended: time exipred with no solution available.\n";
        return 0;
      }
    }
  
    
    // Run nominal present
    if(verbose_flag > 0) cout<<"***Solving no robust problem."<<endl;
    tstart = clock();
    algCplex.solve();
    double nrTime = (double) (clock() - tstart) / CLOCKS_PER_SEC;
    if( algCplex.getStatus() == IloAlgorithm::Infeasible ) {
      cout<<"Algorithm ended: infeasible problem.\n";
      return 0;
    }
    else{
      if( algCplex.getStatus() == IloAlgorithm::Optimal ) {
        //cout<<"Algorithm ended: optimal solution found.\n";
      }
      else{
        cout<<"Algorithm ended: time exipred with no solution available.\n";
        return 0;
      }
    }
    
 

  

    if(verbose_flag > 0) cout<<"Saving solutions."<<endl;
    // Save multiband and nominal solutions
    IloNumArray2D xSol(multibAlgEnv/*type of the elements of the matrix*/, numUTs /*cardinality of one dimension*/);
    IloNumArray ySol(multibAlgEnv, numAPs);
    IloNumArray2D nr_xSol(algEnv, numUTs);
    IloNumArray nr_ySol(algEnv, numAPs);
    
    /*putting solutions in the proper variables*/
    multibAlgCplex.getValues(y, ySol);
    algCplex.getValues(nr_y, nr_ySol);
    for (i = 0; i < numUTs; i++) {
      xSol[i] = IloNumArray(multibAlgEnv);
      multibAlgCplex.getValues( xSol[i], x[i]);
      nr_xSol[i] = IloNumArray(algEnv);
      algCplex.getValues(nr_x[i], nr_xSol[i]);
    }
 
    if(verbose_flag > 0) cout<<"Checking future feasibilities."<<endl;
    
    /*r_future reflects movements of users. If xSol(computed in the present) respects capacity constraints for all AP, 
     the solution is still good in the future :-) 
     */
    IloBool multib_feasFlag = checkFeasibility(rho, w, xSol, r_future); 
    IloBool nr_feasFlag = checkFeasibility(rho, w, nr_xSol, r_future);
    
    IloNumArray2D fut_xSol(futAlgEnv, numUTs);
    IloNumArray fut_ySol(futAlgEnv, numAPs);
    
    //Computing the nominal future solution.
    if(multib_feasFlag){
      // Run future ILP
      if(verbose_flag > 0) cout<<"***Solving future problem."<<endl;
      futAlgCplex.solve();
      if( futAlgCplex.getStatus() == IloAlgorithm::Infeasible ) {
        cout<<"Algorithm ended: infeasible problem.\n";
        return 0;
      }
      else{
        if( futAlgCplex.getStatus() == IloAlgorithm::Optimal ) {
          //cout<<"Algorithm ended: optimal solution found.\n";
        }
        else{
          cout<<"Algorithm ended: time exipred with no solution available.\n";
          return 0;
        }
      }
 
      // Save the future solutions
      futAlgCplex.getValues(fut_y, fut_ySol);
      for (i = 0; i < numUTs; i++) {
        fut_xSol[i] = IloNumArray(futAlgEnv);
        futAlgCplex.getValues(fut_x[i], fut_xSol[i]);
      }
    }
	
    
    /*At this point, we have computed the nominal(nr) and the multiband solution in the present; and the nominal solution for the future(fut), 
     * GIVEN the future positions of UTs. The latter can be computed only if at least the multiband solution is valid also for the future.
    */
    if(verbose_flag > 0) cout<<"Computing power consumptions."<<endl;
    IloNum multib_PowCons = computePowerConsumption(b, p_w, w, xSol, ySol, r_curr);//multiband present P.C.
    IloNum nr_PowCons = computePowerConsumption(b, p_w, w, nr_xSol, nr_ySol, r_curr);//nominal present P.C.
    
    IloNum mf_PowCons; 
    IloNum nrf_PowCons;
    IloNum fut_PowCons;
    if(multib_feasFlag){
      mf_PowCons = computePowerConsumption(b, p_w, w, xSol, ySol, r_future);//multiband future P.C.
      fut_PowCons = computePowerConsumption(b, p_w, w, fut_xSol, fut_ySol, r_future);// nominal future P.C. (GIVEN future positions)
    }
    if(nr_feasFlag)
      nrf_PowCons = computePowerConsumption(b, p_w, w, nr_xSol, nr_ySol, r_future);//nominal future P.C.

    //Compute maximum power consumption (as if all were working at full time)
    IloNum max_PowCons = 0;
    for (j = 0; j < numAPs; j++) {
      max_PowCons += b[j] + p_w[j];
    }

    //cout<<"Printing results.\n";
    printResults(output_filename, multib_PowCons, nr_PowCons, mf_PowCons, nrf_PowCons, fut_PowCons,
                 multib_feasFlag, nr_feasFlag, xSol, ySol, nr_xSol, nr_ySol, fut_xSol, fut_ySol, 
                 ut_vect, cs_vect, multibTime, nrTime, max_PowCons);
    
   
    //Free memory
    xSol.end();
    nr_xSol.end();
    fut_xSol.end();
    ySol.end();
    nr_ySol.end();
    fut_ySol.end();
    multibAlgEnv.end();
    algEnv.end();
    futAlgEnv.end();
 	
  }
  
  catch (const IloException& e) {
     cerr << "Exception caught: " << e << endl;
  }
  catch (...) {
     cerr << "Unknown exception caught!" << endl;
  }

  
  // Close the environment
  simEnv.end();
	
  }
  return 0;

} 



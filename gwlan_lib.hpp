#ifndef __GWLANlib_HPP__
#define __GWLANlib_HPP__

#include <cstdlib>
#include <fstream> 
#include <iostream>
#include <ctime> 
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <unistd.h>
#include <cstring>
#include <limits>
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

using namespace std;
 
///Threshold used in some computations.
/**
 * @see createRobAlgMILP createAlgILP
 */
const IloNum delta_0 = 0.000001;
 
typedef IloArray<IloIntVarArray> assignVarMx;
typedef IloArray<IloNumArray> IloNumArray2D;
typedef IloArray<IloIntArray> IloIntArray2D;
typedef IloArray<IloArray<IloNumArray> > IloNumArray3D;


/// <b>Initializes CPLEX data structures for the Multiband mathematical model.</b>
void createMultibMILP(IloModel mod, assignVarMx x, IloIntVarArray y, IloNumArray b,
                     IloNumArray  p_w, IloNum rho, IloNumArray w, IloNumArray2D r_curr, 
					 IloNumArray3D r_min, IloInt K, IloNumArray H);



/// <b>Initializes CPLEX data structures for the Nominal mathematical model.</b>
/**
 * Several settings are performed: the matrix <i>x</i> is modeled to take values from 0 to 1 (this statement will be refined later when defining constraints),
 * and, if the current rate is less than the threshold <b>delta_0</b>, the values will range from 0 to 0; 
 * the vector y is modeled to take values from 0 to 1(this statement will be refined later when defining constraints).\n
 * Then, the objective function of the nominal model (see <i>(2)</i> in the paper) is defined in the CPLEX model.\n
 * Finally, several constraints are set (see <i>(4)</i> and <i>(5)</i> in the paper).
 */
void createAlgILP(IloModel mod, assignVarMx x, IloIntVarArray y, IloNumArray b,
                  IloNumArray p_w, IloNum rho, IloNumArray w, IloNumArray2D r);




/// <b>It fetches necessary data from the file which contains an instance of a certain configuration.</b>
/**
 *Notable information fetched from the instance is: H[b]; K; traffic demand <i>w</i>; base power consuption <i>b</i>; 
 *"airtime" power consumption<i>p_w</i>; current, future, minimum and minimum-deviation rates.
 */
void readData(IloEnv env, string input_filename, IloInt& K, IloNumArray H, IloNum& rho, 
              vector<string>& cs_vect, vector<string>& ut_vect, IloNumArray w, 
              IloNumArray b, IloNumArray p_w, IloNumArray2D r_curr, 
              IloNumArray3D r_min, IloNumArray2D r_future);




/// <b> Computation of the power consuption of a given solution (given a certain matrix <i>x</i> and a certain array <i>y</i>) </b>
/**
 * According to the paper, the power consuption of an AP is modeled as:\n
 * <b> P[j] = b + p_w * a[j] </b> \n
 * where a[j] = (sum i in I) w[i] / r[i][j] 
*/
IloNum
computePowerConsumption(IloNumArray b, IloNumArray p_w, IloNumArray w, 
                        IloNumArray2D xSol, IloNumArray ySol, IloNumArray2D r);




/// <b> Checks if the capacity constraint is fulfilled by a certain solution. </b>
/**
 * Defined the airtime <i>a[j]</i> as <i>(sum i in I) w[i] / r[i][j] </i>: if it is less or equal than <b>rho</b> the solution is feasible,
 * otherwise it is not.
 */
IloBool 
checkFeasibility(IloNum rho, IloNumArray w, IloNumArray2D xSol, IloNumArray2D r);





/// <b> Computed solutions are saved into an output file (its name is specified as parameter). </b>
void 
printResults(string output_filename, IloNum multib_PowCons, IloNum nr_PowCons, 
             IloNum mf_PowCons, IloNum nrf_PowCons, IloNum fut_PowCons, 
             IloBool rob_feasFlag, IloBool nr_feasFlag,
             IloNumArray2D xSol, IloNumArray ySol,  
             IloNumArray2D nr_xSol, IloNumArray nr_ySol,
             IloNumArray2D fut_xSol, IloNumArray fut_ySol, 
             vector<string> ut_vect, vector<string> cs_vect,
             double multibTime, double nrTime, double max_PowCons);
             
#endif

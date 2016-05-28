#include "gwlan_lib.hpp"


void createMultibMILP(IloModel mod, assignVarMx x, IloIntVarArray y, IloNumArray b,
                     IloNumArray  p_w, IloNum rho, IloNumArray w, IloNumArray2D r_curr, 
					 IloNumArray3D r_min, IloInt K, IloNumArray H)
{
  
   IloEnv env = mod.getEnv();// copies the optimization environment of "mod" in "env"
   IloInt numUTs = x.getSize();
   IloInt numAPs = y.getSize();
   IloInt numBands = H.getSize();
   
   // Create variables x(i,j) forall i in I, j in J
   char varName[100];                         
   for (int i = 0; i < numUTs; ++i) {
      x[i] = IloIntVarArray(env, numAPs, 0/*lb*/, 1/*ub*/); //constraint (40)
      for (int j = 0; j < numAPs; ++j) {
         sprintf(varName, "x.%d.%d", i+1, j+1); 
         x[i][j].setName(varName);
         if(r_curr[i][j]<delta_0)
           x[i][j].setBounds(0,0);
      }
      mod.add(x[i]); //add the x array to the model
   }                                 
   
    
   // Create objective function: 
   // minimize sum(j in J) b(j) * y(j) + p_w(j) * (sum(i in I) w(i) / r_curr(i,j) * x(i,j))
   IloExpr obj(env);
   for (int j = 0; j < numAPs; ++j) { 
        obj += b[j] * y[j];
        for (int i = 0; i < numUTs; ++i) {
          if(r_curr[i][j]>delta_0)
            obj += p_w[j] * (w[i] / r_curr[i][j]) * x[i][j];
        }
   }
   mod.add(IloMinimize(env, obj)); //Add a "minimizator" for the objective function and the evironment 
   obj.end();//frees resources allocated to such variable
 
   
   // Create variables y(j) for j in J 
      for (int j = 0; j < numAPs; ++j) {
        sprintf(varName, "y.%d", j+1); 
        y[j].setName(varName);            
      }
      mod.add(y); // add the y array to the model 
      

   // Create variables pigreco(i) for i in |I| 
  IloNumVarArray pigreco(env, numUTs, 0, IloInfinity); //(constraint 37)
   for (int i = 0; i < numUTs; i++) {
     sprintf(varName, "pigreco.%d", i+1); 
     pigreco[i].setName(varName);
   }
   mod.add(pigreco);
 
   // Create variables mu(b) for b in |B|  
   IloNumVarArray mu(env, numBands, 0, IloInfinity); //(constraint 38)
   for (int b = 0; b < numBands; ++b) {
     sprintf(varName, "mu.%d", b); 
     mu[b].setName(varName);
   }
   mod.add(mu);
   
   
   //create variable delta
   IloNumVar delta(env);// constraint 39
   delta.setBounds(0,IloInfinity);
   delta.setName("delta");
   mod.add(delta);
   
   
   // Add the UT assignament constraints: #34
   // forall i in I: sum(j in J) x(i,j) = 1
   for (int i = 0; i < numUTs; ++i) {
     IloExpr expr(env);
     for (int j = 0; j < numAPs; ++j){
       expr += x[i][j];
     }
     mod.add(expr == 1);
     expr.end();
   }
 
   
   // Add multiband constraints: #35
   for(int j=0; j<numAPs; j++){
	   IloExpr expr(env); 
	   expr += ( K * delta ); 
	   for(int i=0; i<numUTs; i++)
		   if(r_curr[i][j]>delta_0)
			   expr += (w[i] * (x[i][j] / r_curr[i][j])) + pigreco[i];
	   for(int b=0; b<numBands; b++)
		   expr += (mu[b] * H[b]);
	   
	   mod.add(expr <= rho * y[j]);
	   expr.end();
   }
 
   
   // constraint #36
   for(int i=0; i<numUTs; i++){ 	//i=1,...,|I|
	   for(int b=0; b<numBands; b++){	//b=0,...,|B|
		   IloExpr l_expr(env);
		   IloExpr r_expr(env);
		   l_expr += (delta + pigreco[i] + mu[b]);
		   for(int j=0; j<numAPs; j++){
			   IloNumArray rmin_ij = r_min[i][j];
			   if(rmin_ij[b]>delta_0)
				   r_expr += w[i] * ( 1/rmin_ij[b] ) * x[i][j];
			   if(r_curr[i][j]>delta_0)
				   r_expr -= w[i] * ( 1/r_curr[i][j] ) * x[i][j];
		   }
		   mod.add(l_expr >= r_expr); //specified for each i and b
		   l_expr.end();
		   r_expr.end();
	   }
   }
}// END createMultibandAlgMILP
 

void
createAlgILP(IloModel mod, assignVarMx x, IloIntVarArray y, IloNumArray b,
             IloNumArray p_w, IloNum rho, IloNumArray w, IloNumArray2D r)
{
  
   IloInt i, j;
   IloEnv env = mod.getEnv();
   IloInt numUTs = x.getSize();
   IloInt numAPs = y.getSize();
 
   // Create variables x(i,j) forall i in I, j in J
   char varName[100];                         
   for (i = 0; i < numUTs; ++i) {
      x[i] = IloIntVarArray(env, numAPs, 0, 1);
      for (j = 0; j < numAPs; ++j) {
         sprintf(varName, "x.%d.%d", (int)i+1, (int)j+1); 
         x[i][j].setName(varName);
         if(r[i][j]<delta_0)
           x[i][j].setBounds(0,0);
      }
      mod.add(x[i]);
   }                                 
 
   // Create variables y(j) for j in J 
   for (j = 0; j < numAPs; ++j) {
     sprintf(varName, "y.%d", (int)j+1); 
     y[j].setName(varName);            
   }
   mod.add(y);
 
   
   // Create objective function: 
   // minimize sum(j in J) (b(j) * y(j) + p_w(j) * sum(i in I) w(i) / r_curr(i,j) * x(i,j))
   IloExpr obj(env);
   for (j = 0; j < numAPs; ++j) { 
     obj += b[j] *y[j];
     for (i = 0; i < numUTs; ++i) {
       if(r[i][j]>delta_0)
         obj += p_w[j] * (w[i] / r[i][j]) * x[i][j];
     }
   }
   mod.add(IloMinimize(env, obj));// add a "minimizator" :-)
   obj.end();
 
   // Add the UT assignament constraints: #4
   // forall i in I: sum(j in J) x(i,j) = 1
   for (i = 0; i < numUTs; ++i) {
     IloExpr expr(env);
     for (j = 0; j < numAPs; ++j) 
       expr += x[i][j];
     mod.add(expr == 1);
     expr.end();
   }
  
  // Add the capacity constraints: #5
  // forall j in J: sum(i in I) w(i) / r_curr(i,j) * x(i,j) <= rho *y(j)
   for (j = 0; j < numAPs; ++j) {
     IloExpr expr(env);
     for (i = 0; i < numUTs; i++) {
       if(r[i][j]>delta_0)
         expr += (w[i] / r[i][j]) * x[i][j];
     }
     mod.add(expr <= rho * y[j]);
     expr.end();
   }
  

}// END createAlgILP

 

IloNum
computePowerConsumption(IloNumArray b, IloNumArray p_w, IloNumArray w, 
                        IloNumArray2D xSol, IloNumArray ySol, IloNumArray2D r)
{
 
  IloInt i, j;
  IloInt numUTs = xSol.getSize();
  IloInt numAPs = ySol.getSize();
 
  IloNum PowCons = 0;
  for (j = 0; j < numAPs; ++j) { 
    if(ySol[j]>delta_0) 
      PowCons += b[j];
    for (i = 0; i < numUTs; ++i) {
      if(xSol[i][j]>delta_0) 
        PowCons += p_w[j] * (w[i] / r[i][j]);
    }
  }
  
  return PowCons;

}


IloBool
checkFeasibility(IloNum rho, IloNumArray w, IloNumArray2D xSol, IloNumArray2D r)
{
 
  IloInt i, j;
  IloInt numUTs = xSol.getSize();
  IloInt numAPs = xSol[0].getSize();
 

  IloBool feas_flag = IloTrue;
  for (j = 0; j < numAPs; ++j) { 
    IloNum airtime = 0;
    for (i = 0; i < numUTs; ++i) 
      if(xSol[i][j] > delta_0) airtime += (w[i] / r[i][j]); 
    if(airtime > rho) {
      feas_flag = IloFalse;
      break;
    }
  }
  

  return feas_flag;

}



void 
printResults(string output_filename, IloNum multib_PowCons, IloNum nr_PowCons, 
             IloNum mf_PowCons, IloNum nrf_PowCons, IloNum fut_PowCons, 
             IloBool multib_feasFlag, IloBool nr_feasFlag,
             IloNumArray2D xSol, IloNumArray ySol,  
             IloNumArray2D nr_xSol, IloNumArray nr_ySol,
             IloNumArray2D fut_xSol, IloNumArray fut_ySol, 
             vector<string> ut_vect, vector<string> cs_vect,
             double multibTime, double nrTime, double max_PowCons)
{
 
  IloInt i, j;
  IloInt numUTs = xSol.getSize();
  IloInt numAPs = ySol.getSize();
 
  ofstream outfile;
  outfile.open (output_filename.c_str(), ofstream::out);
  outfile<<"Maximum Power Consumption := "<<max_PowCons<<"W"<<endl;
  outfile<<"\n\nMultiband Power Consumption := "<<multib_PowCons<<"W"<<endl;
  outfile<<"\nMultiband Solution CPU-Time := "<<multibTime<<"s"<<endl;
  outfile<<"\n\nNominal Power Consumption := "<<nr_PowCons<<"W"<<endl;
  outfile<<"\nNominal Solution CPU-Time := "<<nrTime<<"s"<<endl;
  if(multib_feasFlag) {
    outfile<<"\n\nOptimal Future Power Consumption := "<<fut_PowCons<<"W"<<endl;
    outfile<<"\n\nMultiband Future Power Consumption := "<<mf_PowCons<<"W"<<endl;  
  }
  else
    outfile<<"\n\nMultiband Future Power Consumption := UNFEASIBLE"<<endl;
  //if(nr_feasFlag)
  if(multib_feasFlag && nr_feasFlag)
    outfile<<"\n\nNominal Future Power Consumption := "<<nrf_PowCons<<"W"<<endl;
  else
    outfile<<"\n\nNominal Future Power Consumption := UNFEASIBLE"<<endl;
  outfile<<"\n\nr_x[*][*] := "<<endl<<"\t\t";
  for (j = 0; j < numAPs; j++) 
      outfile<<cs_vect[j]<<" ";
  outfile<<endl<<"\t";
  for (i = 0; i < numUTs; i++) {
    outfile<<ut_vect[i]<<" ";
    for (j = 0; j < numAPs; j++) 
      outfile<<xSol[i][j]<<" ";
    outfile<<endl<<"\t";
  }
  outfile<<"\n\nr_y[*][*] := "<<endl;
  for (j = 0; j < numAPs; j++) {
    outfile<<cs_vect[j]<<" "<<ySol[j]<<endl;
  }
  outfile<<"\n\nnr_x[*][*] := "<<endl<<"\t\t";
  for (j = 0; j < numAPs; j++) 
      outfile<<cs_vect[j]<<" ";
  outfile<<endl<<"\t";
  for (i = 0; i < numUTs; i++) {
    outfile<<ut_vect[i]<<" ";
    for (j = 0; j < numAPs; j++) 
      outfile<<nr_xSol[i][j]<<" ";
    outfile<<endl<<"\t";
  }
  outfile<<"\n\nnr_y[*][*] := "<<endl;
  for (j = 0; j < numAPs; j++) {
    outfile<<cs_vect[j]<<" "<<nr_ySol[j]<<endl;
  }
  if(multib_feasFlag) {
    outfile<<"\n\nfut_x[*][*] := "<<endl<<"\t\t";
    for (j = 0; j < numAPs; j++) 
        outfile<<cs_vect[j]<<" ";
    outfile<<endl<<"\t";
    for (i = 0; i < numUTs; i++) {
      outfile<<ut_vect[i]<<" ";
      for (j = 0; j < numAPs; j++) 
        outfile<<fut_xSol[i][j]<<" ";
      outfile<<endl<<"\t";
    }
    outfile<<"\n\nf_y[*][*] := "<<endl;
    for (j = 0; j < numAPs; j++) {
      outfile<<cs_vect[j]<<" "<<fut_ySol[j]<<endl;
    }
  }
  outfile.close();


}// END printResults 




void readData(IloEnv env, string input_filename, IloInt& K, IloNumArray H, IloNum& rho, 
              vector<string>& cs_vect, vector<string>& ut_vect, IloNumArray w, 
              IloNumArray b, IloNumArray p_w, IloNumArray2D r_curr, 
              IloNumArray3D r_min, IloNumArray2D r_future)
{ 
  
  ifstream infile;
  string temp;  
  
  // Seek for b(j)
  infile.open (input_filename.c_str(), ifstream::in);
  string candidate_site;
  IloNum baseline_power;
  while(temp!="b"){
    infile >> temp;
    if(infile.eof()) {
      cerr<<"ERROR: Missed 'b' parameter in input file!"<<endl;
      break;
    }
  }
  infile >> temp; // read ":="
  
  do{
     infile >> candidate_site >> baseline_power;
     if(candidate_site!=";") {
       cs_vect.push_back(candidate_site);
       b.add(baseline_power);
     }
     else
       break;
  }while(true);
  infile.close();
  
  
  // Read p_w(j)
  infile.open (input_filename.c_str(), ifstream::in);
  IloNum wireless_power;
  while(temp!="p_w"){
    infile >> temp;
  }
  infile >> temp; // read ":="
  do{
     infile >> candidate_site >> wireless_power;
     if(candidate_site!=";") {
       p_w.add(wireless_power);
     }
     else
       break;
  }while(true);
  infile.close();
 

  // Read rho
  infile.open (input_filename.c_str(), ifstream::in);
  while(temp!="rho"){
    infile >> temp;
  }
  infile >> temp; // read ":="
  infile >> rho;
  infile.close();
  

  // Read w_i
  infile.open (input_filename.c_str(), ifstream::in);
  string user_terminal;
  IloNum dem;
  while(temp!="w"){
    infile >> temp;
  }
  infile >> temp; // read ":="
  do{
     infile >> user_terminal >> dem;
     if(user_terminal!=";")
       ut_vect.push_back(user_terminal);
       w.add(dem);
  }while(user_terminal!=";");
  infile.close();


  // Read r_curr(i,j)
  infile.open (input_filename.c_str(), ifstream::in);
  IloNum r_elem;
  IloInt UT_index = 0;
  while(temp!="r_curr:"){
    infile >> temp;
  }
  do{
    infile >> temp; // read ":="
    if(temp==":=")
    break;
  }while(true);
  do{
     infile >> user_terminal;
     if(user_terminal!=";"){
       r_curr.add(IloNumArray(env,cs_vect.size()));
       for(size_t i=0; i<cs_vect.size(); ++i){
         infile >> r_elem;
         r_curr[UT_index][i] = r_elem;
       }
       ++UT_index;
     }
  }while(user_terminal!=";");
  infile.close();
  

  
  // Read r_min-b(i,j)
  infile.open (input_filename.c_str(), ifstream::in);
  UT_index = 0;
  while(temp!="r_min-b"){
    infile >> temp;
  }
  do{
    infile >> temp; // read ":="
    if(temp==":=")
    	break;
  }while(true);
  
  do{
	 r_min.add(IloArray<IloNumArray>(env, cs_vect.size()));
     infile >> user_terminal;
     if(user_terminal!=";"){
    	 for(size_t i=0; i<cs_vect.size(); ++i){
    		 infile >> candidate_site;
    		 infile >> temp; // read "["
    		 IloNumArray auxiliary(env);
    		 do{
    			 infile>>temp;
    			 if(temp!="]") 
    				 auxiliary.add(stod(temp));
    		 }while(temp!="]");
    		 r_min[UT_index][i] = auxiliary;
    	 }
       ++UT_index;
     }
  }while(user_terminal!=";");
  
  infile.close();
 
  
  // Read r_future(i,j)
  infile.open (input_filename.c_str(), ifstream::in);
  UT_index = 0;
  while(temp!="r_future:"){
    infile >> temp;
  }
  do{
    infile >> temp; // read ":="
    if(temp==":=")
    break;
  }while(true);
  do{
     infile >> user_terminal;
     if(user_terminal!=";"){
       r_future.add(IloNumArray(env,cs_vect.size()));
       for(size_t i=0; i<cs_vect.size(); ++i){
         infile >> r_elem;
         r_future[UT_index][i] = r_elem;
       }
       ++UT_index;
     }
  }while(user_terminal!=";");
  infile.close();
 
  // Read K
  infile.open (input_filename.c_str(), ifstream::in);
  while(temp!="K"){
    infile >> temp;
  }
  infile >> temp; // read ":="
  infile >> K;
  infile.close();
  
  // Read H[b]
  infile.open (input_filename.c_str(), ifstream::in);
  IloNum cardinality;
  while(temp!="H[b]"){
    infile >> temp;
  }
  infile >> temp; // read ":="
  do{
     infile >> temp >> cardinality;
     if(temp!=";")
       H.add(cardinality);
  }while(temp!=";");
  infile.close();
  
}// END readData 

//      instance-maker.cc
//
//      Copyright 2016 Alessio Massidda <alessio.massidda@hotmail.it>
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.

/*
 * Generatore di istanze per MULTIBAND Green WLANs
 * Build command: g++ -std=c++0x -o inst-maker -I. ConfigFile.cpp instance-maker.cc
 */
#include <random>
#include <vector>
#include <iostream>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include "ConfigFile.h"

using namespace std;

///PI (pi greco) constant
#define PI 3.14159265
/// Speed of light
#define SL 2.99792458e+8
///Number of candidate positions
#define N_pos 3

///Maximum goodput per AP  [kbps]
const unsigned max_goodput = 18000;


///Not significant, for the moment being there is only 1 class
int n_AP_classes = 1;

/// NOT USED in this version of the algorithm
double *cWeight;
///Power value used to compute candidate sites
double *P0;

/// Maximum power transmission [W]
double Pmax;
/// Minimum power transmission [W]
double Pmin;
/// Number of power levels used. Currently it is ONLY 1.
int nPwrLevels;

///Power Step Mode
/**
* - if it is set to 0 it means that tabulated values are used (e.g. "PT_tab = value1,value2,...").
*\n USUALLY, when this mode is set, the number of values is the same of the number of power levels.

* - if it is set to 1, it is linear, between Pmin and Pmax.

* - if it is set to 2, it is geometric (linear in dB), between Pmin and Pmax.

* - if it is set to 3, it grows by steps of 3dB, starting from Pmin.

* @see PT_tab Pmin Pmax
*/
int psm = 0;

/// Transmission Power tabulation (meaningful if psm is equal to 0) [W]
/**
* @see psm
*/
string PT_tab;

///For each AP class, it stores a list of power levels.
/**
* For example, eta[i][j] is pointing to the (j+1)-th power level of the (i+1)-th AP class.
* It is initialized in the scope of "readInput()" function by converting "eta_tab" values (from string to double).
*\n Notice that currently there is JUST ONE AP CLASS and ONE POWER LEVEL considered.
* @see readInput
*/
double **eta;
/// The new mu[j] parameter, either a random value from a uniform distrib: mu[j] = mu_j_avg +/- mu_j_delta, or an indexed value
double *mu;
/// The power consumption for frame reception
double *P_demod;

///Transmission powers. <b> Meaningful only if psm=0. </b>
/** It is an array of length nPwrLevels, initialized in createFile() from the tabulation contained in PT_tab
* @see nPwrLevels PT_tab psm
*/
double *Pt;

/// Frequency band used (when the documentation has been written, it was 2400) [MHz]
double freq;
/// Antenna gain (when the documentation has been written, it was computed as pow(10,x/10), where 'x' is the value fetched from the ini file)
double antgain;
/// attenuation factor or p.l. exponent
double attenf;
/// distance used when computing path loss [m]
/** @see computeRefPL computePL_freespace */
double d0 = 1;
/// ref. Path Loss at d0: if it is set to -1 it is computed automatically(free space formula) [dB]
double PL0 = -1;
/// path loss model: 0 = indoor LOS (default), 1 = multi-wall NLOS
int pl_model = 0;
// path loss formula (log-distance, NLOS multi wall)
/// constant path loss multi-wall [dB]
double Lc = 0;

///Stands for "overall-multi-wall-loss-dB".
/**
* If present and non-negative will replace all terms about columns, walls and losses
* @see Nw Ncl wall_sep col_sep Lc Lw Lcl
*/
double LMW = -1;
/// wall loss [dB]
double Lw = 0;
/// column loss [dB]
double Lcl = 0;
/// number of crossed walls (-1 to compute automatically)
int Nw = -1;
/// number of crossed columns (-1 to compute automatically)
int Ncl = -1;
/// distance between walls [m]
double wall_sep = 0;
/// distance between columns [m]
double col_sep = 0;

/// Noise of the channel [dB]
double noise;
/// Parameter for rate computation [kb/s]
double theta;
/// Parameter for rate computation [kb/s]
double beta;
///  Parameter for rate computation [kb/s]
double delta;
/// Parameter used for Shannon formula
double zeta;
/// Parameter used for Shannon formula [kb/s]
double B;

/// Received power threshold [dB] (taken from Cisco Aironet 2600)
double pwrThr;
///SNIR threshold [dB] (taken from Cisco Aironet 2600)
double snirThr;
/// Maximum PHY throughput possible between every couple AP/UT [kbps]
double r_max;

/// Maximum goodput at IP/APP level. It is computed in the function "readInput()" as "r_max/oh_factor" [kbps] .
/**
* @see readInput
*/
double goodput_max;

///Optional scale factor for protocol overhead: used to convert goodput from MAC level to APP/IP level.
/** For example: from MAC to APP this value could be 3, because we assume that actual data is encapsulated by TCP, IP and MAC headers.*/
double oh_factor = 1;

///Encodes the formula used to compute the rate
/** \n 0 = linear approx. ( r = theta * SNR ) (default),
* 1 = logarithmic ( r = beta * SNR[dB] + delta )
* 2 = Shannon ( r = zeta * B * log(1+SNIR) )
*/
int formula = 0;

///Width of the place considered for the simulation [m].
/**
* ATTENTION: width and height are set EQUAL by the scripting system
* because in the file "gwlan-configs.dat" currently you have to provide just one dimension.
* \n So, THE PLACE IS SQUARED !
*/
double fieldSizeX = 1;

///Height of the place considered for the simulation [m].
/**
* ATTENTION: width and height are set EQUAL by the scripting system
* because in the file "gwlan-configs.dat" currently you have to provide just one dimension.
* \n So, THE PLACE IS SQUARED !
*/
double fieldSizeY = 1;
/// Minimum distance allowed between each connected AP and UT [m]. If the value is 0 it means that this constraint is neglected.
double minDist = 0;

/// APs and UTs deployment strategy: 0 means full random, 1 means big squares, 2 means 2-density zones (based on squares)
/**
* @see createFullRandom createWithRegularSquares createWithTwoDensityZones
*/
int deployMode = 0;

/// Ratio considered as a High-Density Region (hdr) in the area used for the simulation
double hdrRatio = 0.5;

///The fraction of users in the High Density Region (applies iff "deployMode" is equal to 2).
/**
* NOTE: if for instance frac_UT_hdr=0.7 , it means that 70% of users are in 30% of area.
* @see deployMode
* @see hdrRatio
*/
double frac_UT_hdr = 0.5;

/// Maximum traffic demand of users(IP/APP level) [kbps]
double maxDemand = 0;
/// Minimum traffic demand of users(IP/APP level) [kbps]
double minDemand = 0;
/// Set to 1 to autogenerate the demand interval (maxDemand and minDemand are neglected); if set to 0 maxDemand and minDemand are considered.
int autoGenDemand = 0;
///Saturation factor of the network. In other words it is the reduction of the average threshold value self-computed (vml)
double utilFactor = 1.0;

///Never used.
/** It is multiplied to the vml (autoGenDemand variable) to have a factor of the size of the auto-generated interval.
* It must be between 0 and 2 (0 = all equal to vml, 2 = variable between 0 and 2*vml)
*/
double expFactor = 1.0;

///Service Data Unit (just PAYLOAD) at MAC LEVEL measured in BITs
int frameLength;
///PAYLOAD at Application level , measured in BYTEs
int dataLength;
///Protocol overhead, from APP to LLC (e.g. HEADERS), measured in BYTEs
int ohLength;
///  Computed in "readInput()" after having fetched "frameLength" , "dataLength" and "ohLength" [bytes] .
double msdu = 1500;
///dl fraction
double dlFraction = 0.75;
///Maximum airtime for APs (as known as rho)
double max_airtime = 1;
/// Number of Access Points given in this configuration.
int nAPs;
/// Number of User Terminals given in this configuration.
int nUsers;
///Number of users belonging to each band (H[0]= number of static users, H[1]= number of low-mobility users...)
int * H;
/// Important parameter for the multiband model (number of links subject to fluctuation).
int K;
///Random variable. It will be initalized properly in the main function.
int seed = 0;

///Number of bands for user mobility in Multiband GWLAN
int Nbands;
///When a user moves, he has an upper bound in the number of meters: it is computed as 'maxPercMov' percentage of fieldSizeX
float maxPercMov;
///Mapping UT-band, valid for computing future rates
vector<int> mapping;

///Generator of random numbers, properly initialized in the main function
default_random_engine generator;

/// Maximum speed of the users who move [m/s]
double s_max = 1;
/// Delta T (it depends on the configuration taken into account)
//double delta_t = 1;
/// Variance of the current position (it depends on the configuration taken into account)
double varianceCurrent = 0.2;
/// Variance of the future position (it depends on the configuration taken into account)
double varianceFuture = 0.1;
/// Quantile set for the Rayleigh distribution (it depends on the configuration taken into account)
double quantileRayleigh = 0.025;
/// Sigma for the Rayleigh distribution (it depends on the configuration taken into account)
double sigmaRayleigh = 1/sqrt(2);



///Data structure to represent positions of objects (i.e. users and APs) in the area considered for the simulation.
struct coord {
    	///Position of the object in the x-axis
	double x;
	///Position of the object in the y-axis
	double y;
};

/// Generic user's coordinates.
/** @see coord
*/
struct coord *user;

/// Generic AP's coordinates. See definition of the type "coord".
/** @see coord
*/
struct coord *ap;

///Stores the distance between each UT and each AP
/** For example dist[i][j] contains the distance between user i and AP j */
double **dist;
/// path gain [dB]. It is NEVER SET TO SOME SPECIFIC VALUE (it has always the default value 0)
/** Please recall that <i>path_gain[dB] = - path_loss[dB]</i> */
double **alpha_dB;
///UT-AP transmission rate matrix [kbps]. It is initialized in computeMaxPhyRates()
double **r;

/// for each AP , the closest UT is stored.
int *closest;

///Currently it is USELESS since there is just one AP class.
/** It will be an array that will represent for each AP the class it has been randomly assigned to it
* (this will be done by "readInput()" ).
*/
int *apclass;

///This variable will host a function for computing the path loss, depending on pl_model (path loss model).
/**
* It is a pointer to a function, which returns a double and takes as only parameter a double. \n
* The pointer to the proper function is assigned in the main() and will be either
* computePL_logdistNLOS(double d) or computePL_freespace(double d).
* @see pl_model
*/
double (*computePathLoss)(double);

///Flag which is used for debugging purposes.
int debug = 0;

///Output file which will contain data about the generated instance.
fstream fu;
///File read by this program and generated by the calling script.
const char* iniName = "gwlan-input.ini";
///Name for the output file
/**  @see fu  */
const char* fuName = "gwlan-input.dat";

const char* fcName = "extra.dat";

///Reads some parameters from "gwlan-input.ini" and assignes values to the proper global variables.
/**
* There are also some inner variables in this function, please check them in the source code. The variables initialized have a link below.
* @see n_AP_classes P0 mu P_demod Pmax Pmin nPwrLevels psm PT_tab freq antgain attenf pl_model d0 PL0 Lc LMW Lw Lcl Nw Ncl wall_sep
* @see col_sep noise theta beta delta pwrThr snirThr r_max oh_factor formula fieldSizeX fieldSizeY minDist frac_UT_hdr hdrRatio
* @see maxDemand minDemand autoGenDemand utilFactor expFactor frameLength dataLength ohLength dlFraction max_airtime nAPs nUsers
* @see deployMode zeta B s_max delta_t varianceCurrent varianceFuture quantileRayleigh sigmaRayleigh expFactor msdu goodput_max
* @see apclass eta K H
*/
void readInput()
{
    	/* a = mu_j_avg, but that value in the ini file is COMMENTED! a IS NOT INITIALIZED */
	double a;
	/* b = mu_j_delta, but that value in the ini file is COMMENTED! b IS NOT INITIALIZED */
	double b;
	// Antenna gain: depending on the model of AP it ranges from 2 to 4/5 dBi
	double x;
    string s;
	char c[25];
    	// This, when present, is in place of eta. It is a collection of 4 eta(s)
	string *eta_tab;
	char *u, *t;

	try
	{
		ConfigFile config( iniName );

        /* It has been already initialized to 1 in line 47 */
		config.readInto(n_AP_classes, "n_AP_classes");

		// Array of power values for computing candidate sites (one array-index per AP class)
		P0 = new double [n_AP_classes] ();
		// Array of eta_tab (one array-index per AP class)
		eta_tab = new string [n_AP_classes] ();
		// Array of mu (one array-index per AP class)
		mu = new double [n_AP_classes] ();
		// Array of P_demod (one array-index per AP class)
		P_demod = new double [n_AP_classes] ();

		x = 0;

		/* For each AP_class, fill the arrays for "P_zero" , "eta_tab" , "mu" , "P_demod" */
		for(int i = 0; i<n_AP_classes; i++) {
			sprintf(c, "P_zero[%d]", i+1);
			config.readInto(P0[i], c);
			sprintf(c, "eta_tab[%d]", i+1);
			config.readInto(eta_tab[i], c);
			sprintf(c, "mu[%d]", i+1);
			config.readInto(mu[i], c);
			sprintf(c, "P_demod[%d]", i+1);
			config.readInto(P_demod[i], c);
		}
		
		config.readInto(Nbands, "n_of_bands");
		config.readInto(maxPercMov, "max_perc_mov");
		config.readInto(Pmax, "PT_max");
		config.readInto(Pmin, "PT_min");
		config.readInto(nPwrLevels, "n_Pwr-Level");
		config.readInto(psm, "pwr-step-mode");
		config.readInto(PT_tab, "PT_tab");
		config.readInto(a, "mu_j_avg");
		config.readInto(b, "mu_j_delta");
		
		config.readInto(freq, "central_frequency-MHz");
		config.readInto(x, "antenna_gain-dB");
		antgain = pow(10, x/10);
		config.readInto(attenf, "attenuation_factor");
		config.readInto(pl_model, "PL-model");
		config.readInto(d0, "ref-distance");
		config.readInto(PL0, "ref-path-loss-dB");
		config.readInto(Lc, "const-path-loss-mw-dB");
		config.readInto(LMW, "overall-multi-wall-loss-dB");
		config.readInto(Lw, "wall-loss-dB");
		config.readInto(Lcl, "column-loss-dB");
		config.readInto(Nw, "num-of-walls");
		config.readInto(Ncl, "num-of-columns");
		config.readInto(wall_sep, "wall-separation");
		config.readInto(col_sep, "column-separation");

		config.readInto(noise, "noise-dB");
		config.readInto(theta, "theta");
		config.readInto(beta, "beta");
		config.readInto(delta, "delta");

		config.readInto(pwrThr, "recv-pwr_threshold-dB");
		config.readInto(snirThr, "SNIR_threshold-dB");
		config.readInto(r_max, "max-rate");
		config.readInto(oh_factor, "overhead_factor");
		config.readInto(formula, "rate-formula");

		config.readInto(fieldSizeX, "field-size_X");
		config.readInto(fieldSizeY, "field-size_Y");
		config.readInto(minDist, "min_dist_AP-user");
		config.readInto(frac_UT_hdr, "frac_UT_hdr");
		hdrRatio = 1 - frac_UT_hdr;

		config.readInto(maxDemand, "max_demand");
		config.readInto(minDemand, "min_demand");

		config.readInto(autoGenDemand, "autogen_demand");
		config.readInto(utilFactor, "util_factor");
		config.readInto(expFactor, "exp_factor");

		// This variable must be greater or equal than 0
		if(expFactor<0)
			expFactor=0;
        	//This variable must be less or equal than 2
		if(expFactor>2)
			expFactor=2;

		config.readInto(frameLength, "frame_len");
		config.readInto(dataLength, "datagram_len");
		config.readInto(ohLength, "overhead_len");
		config.readInto(dlFraction, "dl_fraction");
		config.readInto(max_airtime, "max_airtime");

		config.readInto(nAPs, "n_AP");
		config.readInto(nUsers, "n_User");

		config.readInto(deployMode, "deploy_strat");
		config.readInto(zeta, "zeta");
		config.readInto(B, "B");

		config.readInto(s_max, "max_speed");
		//config.readInto(delta_t, "delta_t");
		config.readInto(varianceCurrent, "variance_present_pos");
		config.readInto(varianceFuture, "variance_future_pos");
		config.readInto(quantileRayleigh, "quantile");
		//config.readInto(fracMovUT, "frac_moving_UT");
		config.readInto(sigmaRayleigh, "sigma_rayleigh");
	}
	catch( ConfigFile::file_not_found& e ) {
		cout << "Errore: file di configurazione '" << e.filename << "' non trovato.\n";
		exit(1);
	}
	catch( ConfigFile::key_not_found& e ) {
		cout << "Errore: chiave '" << e.key << "' non trovata.\n";
		exit(1);
	}

	// set frame length at the MAC layer (i.e. the MSDU). These variables are declared in the global scope
	if ( frameLength > 0 )
		msdu = frameLength/8;
	else if ( dataLength > 0)
		msdu = dataLength + ohLength;

	// set max IP/app rate
	goodput_max = r_max/oh_factor;

	// immediately assign a class to each AP. In practice it is USELESS, because when this comment has been written there was just 1 class
	apclass = new int [nAPs] ();
	for(int i = 0; i<nAPs; i++) {
        //EVERY AP GETS THE SAME CLASS, IT IS ONLY ONE.
		apclass[i] = (int)(rand() % n_AP_classes);
	}

	// Transform the strings (eta_tab) into doubles (eta)
	eta = new double* [n_AP_classes];
	int i;
	for(int j = 0; j<n_AP_classes; j++) {
		eta[j] = new double [nPwrLevels] ();
		u = (char*)eta_tab[j].c_str();
		t = strtok(u, ",");
		i=0;
		while (t && i<nPwrLevels) {
			eta[j][i] = atof(t);
			i++;
			t = strtok(NULL, ",");
		}
		if(i<nPwrLevels) {
			for(int k=i ; k<nPwrLevels; k++)
				eta[j][k] = eta[j][i];
		}
	}// We end up with just one value (e.g. 75) because we have 1 AP class and 1 power level.


    // Assigning H and K: important parameters for the robust model discussed in the paper.
    K = floor(0.9*nUsers);
    //H will be initialized in the createFile() function
    H = new int[Nbands+1] ();

}


///For each UT and each AP assigns a completely random position within the area used for the simulation
void createFullRandom()
{
	int i,j;

	for(i=0; i < nUsers; i++) {
		user[i].x = (double)rand()/RAND_MAX * fieldSizeX;
		user[i].y = (double)rand()/RAND_MAX * fieldSizeY;
	}

	for(j=0; j < nAPs; j++) {
		ap[j].x = (double)rand()/RAND_MAX * fieldSizeX;
		ap[j].y = (double)rand()/RAND_MAX * fieldSizeY;
	}
}

///The area used for the simulation is mentally divided in subregions and APs and UTs are evenly distributed among such subregions.
/** How? First of all, 2 divisors are computed:
* divisore1 is the square root of the number of AP in the current configuration minus 0.5 ;
* divisore2 is the number of AP in the current configuration divided by divisore1. For example, if |J| is 4, then divisore1 and
* divisore2 will both be equal to 2.
* Now, the space is partitioned in divisore1*divisore2 subareas whose width is fieldSizeX/divisore1 and whose height is
* fieldSizeY/divisore2. For sure, every subarea will be covered by one AP. It could be the case that |J| is more than the number
* of generated sub-areas (divisore1*divisore2) , so the remaining |J| - (divisore1*divisore2) APs will be assigned to a random position
* in the full area.
* The same argument is then applied for users, but instead of placing 1 per subarea, they are placed by groups of
* |I|/(divisore1*divisore2) per subarea.
* @see fieldSizeX fieldSizeY nAPs nUsers
*/
void createWithRegularSquares()
{
	int i,j,k,m;
	double subFieldStepX, subFieldStepY;

	int divisore1 = ceil(sqrt(nAPs) - 0.5); //ceil rounds upward (per eccesso)
	int divisore2 = floor(nAPs / divisore1); //floor rounds downward (per difetto)
	int resto = nAPs - divisore1 * divisore2; //number of APs which will be installed completely randomly after that every sub-area has 1 AP.
	if(debug) cout<<"nAPs="<<nAPs<<", divisore1="<<divisore1<<", divisore2="<<divisore2<<", resto="<<resto<<"\n";

	subFieldStepX = fieldSizeX / divisore1;
	subFieldStepY = fieldSizeY / divisore2;

	if(debug) cout<<"subFieldStepX="<<subFieldStepX<<", subFieldStepY="<<subFieldStepY<<", fieldSizeX="<<fieldSizeX<<", fieldSizeY="<<fieldSizeY<<"\n";

	/* BEGINNING OF THE PART REGARDING APs */
	for(k=0; k < divisore1; k++){
		for(m=0; m < divisore2; m++){
			j = k*divisore2 + m;
			ap[j].x = k * subFieldStepX + (double)rand()/RAND_MAX * subFieldStepX;
			ap[j].y = m * subFieldStepY + (double)rand()/RAND_MAX * subFieldStepY;
		}
	}
    // the remaining APs are spred completely randomly
	for(j=k*m; j < nAPs; j++) {
		ap[j].x = (double)rand()/RAND_MAX * fieldSizeX;
		ap[j].y = (double)rand()/RAND_MAX * fieldSizeY;
	}
    /* END OF THE PART REGARDING APs */

    // place nUserGrp users per subarea
	int nUsrGrp = floor(nUsers / (divisore1 * divisore2));

    /* BEGINNING OF THE PART REGARDING UTs */
    //assign to each subarea nUserGrp users
	for(k=0; k < divisore1; k++){
		for(m=0; m < divisore2; m++){
			for(j=0; j < nUsrGrp; j++) {
				//i = k*divisore2*nUsrGrp + m*nUsrGrp + j;
				i = (k*divisore2 + m)*nUsrGrp + j;
				user[i].x = k * subFieldStepX + (double)rand()/RAND_MAX * subFieldStepX;
				user[i].y = m * subFieldStepY + (double)rand()/RAND_MAX * subFieldStepY;
			}
		}
	}
    //remaining users are placed fully randomly
	for(i=j*k*m; i < nUsers; i++) {
		user[i].x = (double)rand()/RAND_MAX * fieldSizeX;
		user[i].y = (double)rand()/RAND_MAX * fieldSizeY;
	}
}

///Distribute the users un-evenly among a periphereal low-density region and a central high-density region.
/**
* APs are placed following exactly the same argument in createWithRegularSquares().
* UTs are instead split in two groups: the first one is distributed according to the same argument in createWithRegularSquares();
* the second one is distributed only in the high density region (which is central with respect to the whole area).
* Parameters such as hdrRatio and frac_UT_hdr here are used. To see further details about the implementation, please refer
* to the comments on the source code.
* @see createWithRegularSquares() hdrRatio frac_UT_hdr
*/
 void createWithTwoDensityZones()
{
    /* The part regarding the APs is the same as in createWithRegularSquares() */
	int i,j,k,l,m;
	double subFieldStepX, subFieldStepY;

	int divisore1 = ceil(sqrt(nAPs) - 0.5);
	int divisore2 = floor(nAPs / divisore1);
	int resto = nAPs - divisore1 * divisore2;

	if(debug) cout<<"nAPs="<<nAPs<<", divisore1="<<divisore1<<", divisore2="<<divisore2<<", resto="<<resto<<"\n";

	subFieldStepX = fieldSizeX / divisore1;
	subFieldStepY = fieldSizeY / divisore2;

    //assign to each sub-area 1 AP
	for(k=0; k < divisore1; k++){
		for(m=0; m < divisore2; m++){
			j = k*divisore2 + m;
			ap[j].x = k * subFieldStepX + (double)rand()/RAND_MAX * subFieldStepX;
			ap[j].y = m * subFieldStepY + (double)rand()/RAND_MAX * subFieldStepY;
		}
	}

    // the remaining APs are spread completely randomly
	for(j=k*m; j < nAPs; j++) {
		ap[j].x = (double)rand()/RAND_MAX * fieldSizeX;
		ap[j].y = (double)rand()/RAND_MAX * fieldSizeY;
	}
    /* ********************************* END OF THE PART REGARDING APs ********************************** */


    /* ********************************* STARTING POINT FOR PLACING UTs ********************************** */
    //number of users in the high density region
	int nUsersHighDens = floor(nUsers * frac_UT_hdr + 0.5);
	//number of users in the low density region
	int nUsersLowDens = nUsers - nUsersHighDens;
	//width of the high density region
	double hdrSizeX = fieldSizeX * sqrt(hdrRatio);
	//heigh of the high density region
	double hdrSizeY = fieldSizeY * sqrt(hdrRatio);
	//total area of the high density region
	double hdrArea = hdrSizeX * hdrSizeY;
	//total area of the low density region
	double ldrArea = fieldSizeX * fieldSizeY - hdrArea;
	//density of UTs in high density region
	double densHDR = nUsersHighDens/hdrArea;
	//density of UTs in low density region
	double densLDR = nUsersLowDens/ldrArea;
    //UTs are placed in 2 passes: the first is the same of createWithRegularSquares()
	int nUsersFirstPass = floor(nUsersLowDens + densLDR * hdrArea + 0.5); //why densLDR * hdrArea ???


	if(debug) cout<<"nUsersHighDens="<<nUsersHighDens<<", nUsersLowDens="<<nUsersLowDens<<", hdrSizeX="<<hdrSizeX<<", hdrArea="<<hdrArea<<", densHDR="<<densHDR<<", densLDR="<<densLDR<<", nUsersFirstPass="<<nUsersFirstPass<<"\n";

    //users are placed in the subareas in groups. Given that 2 passes are performed,
    // this is the dimension of the groups for the 1st pass
	int nUsrGrp = floor(nUsersFirstPass / (divisore1 * divisore2));

    /* BEGINNING OF 1st PASS */
	for(k=0; k < divisore1; k++){
		for(m=0; m < divisore2; m++){
			for(j=0; j < nUsrGrp; j++) {
				i = (k*divisore2 + m)*nUsrGrp + j;
				user[i].x = k * subFieldStepX + (double)rand()/RAND_MAX * subFieldStepX;
				user[i].y = m * subFieldStepY + (double)rand()/RAND_MAX * subFieldStepY;
			}
		}
	}

	for(i=j*k*m; i < nUsersFirstPass; i++) {
		user[i].x = (double)rand()/RAND_MAX * fieldSizeX;
		user[i].y = (double)rand()/RAND_MAX * fieldSizeY;
	}
    /* END OF 1st PASS */

    //users are placed in 2 passes, this is the second and the last.
	int nUsersSecondPass = nUsers - nUsersFirstPass;
    //number of APs in the high density region
	int nAPsInHDR = floor(nAPs * hdrRatio + 0.5);

	/* The following paramenters are all RECOMPUTED because of the second pass to perform. */
	divisore1 = ceil(sqrt(nAPsInHDR) - 0.5);
	divisore2 = floor(nAPsInHDR / divisore1);
	resto = nAPsInHDR - divisore1 * divisore2;
    subFieldStepX = hdrSizeX / divisore1; //also the high density region is divided in subareas :)
	subFieldStepY = hdrSizeY / divisore2; //also the high density region is divided in subareas :)
    nUsrGrp = floor(nUsersSecondPass / (divisore1 * divisore2));

    //necessary to have as a starting x the border of the high density region
	double hdrShiftX = (fieldSizeX - hdrSizeX) / 2;
	//necessary to have as a starting y the border of the high density region
	double hdrShiftY = (fieldSizeY - hdrSizeY) / 2;

	if(debug) cout<<"nAPsInHDR="<<nAPsInHDR<<", divisore1="<<divisore1<<", divisore2="<<divisore2<<", subFieldStepY="<<subFieldStepY<<", hdrShiftX="<<hdrShiftX<<", resto="<<resto<<"\n";

    /* BEGINNING OF THE 2nd PASS */
    //UTs are distributed in groups in each subregion of the high density region
	for(k=0; k < divisore1; k++){
		for(m=0; m < divisore2; m++){
			for(j=0; j < nUsrGrp; j++) {
				i = (k*divisore2 + m)*nUsrGrp + j + nUsersFirstPass;
				user[i].x = k * subFieldStepX + (double)rand()/RAND_MAX * subFieldStepX + hdrShiftX;
				user[i].y = m * subFieldStepY + (double)rand()/RAND_MAX * subFieldStepY + hdrShiftY;
			}
		}
	}
    //the UTs have not yet been assign to any subregion of the high density region, are placed randomly in the HDR
	for(i=j*k*m; i < nUsersSecondPass; i++) {
		l = i + nUsersFirstPass;
		user[l].x = (double)rand()/RAND_MAX * hdrSizeX + hdrShiftX;
		user[l].y = (double)rand()/RAND_MAX * hdrSizeY + hdrShiftY;
	}
}

///UTs too close to the APs are moved.
/**
 * This function is called only if the variable minDist is greater than 0. If it is the case, the distance between the user i and
 * the AP j is adjusted accordingly.
 * @see minDist dist user
 */
void adjustDistance()
{
	double x_app, y_app; // how far to move from the AP
	int i,j;

	for(i=0; i < nUsers; i++) {
		for(j=0; j < nAPs; j++) {
			if ( dist[i][j] < minDist ) {
                //distance between user i and AP j in the x-axis
				x_app = abs(user[i].x - ap[j].x);
				//distance between user i and AP j in the y-axis
				y_app = abs(user[i].y - ap[j].y);
				//now the variable stores the distance on the x-axis for the user from the AP
				x_app *= minDist / dist[i][j];
				//now the variable stores the distance on the y-axis for the user from the AP
				y_app *= minDist / dist[i][j];

				/* Apply new coordinates to the user */
				user[i].x = x_app + ap[j].x;
				user[i].y = y_app + ap[j].y;

				/* Checks whether the new coordinates go beyond the area of the simulation.
				* If it is the case, the user is moved on the opposite side while keeping the new distance.
				*/
				if ( user[i].x > fieldSizeX)
					user[i].x = ap[j].x - x_app;
				if ( user[i].y > fieldSizeY)
					user[i].y = ap[j].y - y_app;
			}
		}
	}
}


///Generates the position of UTs and APs and compute the distance between each couple UT-AP.
/**
 * The strategy for deploying UTs and APs is given by deployMode , passed as actual parameter.
 * @see deployMode
 */
void computeInstance(int mode /*!< formal parameter for deployMode */)
{
	int i,j;

	// allocate memory for data structures
	user = new struct coord [nUsers] ();
	ap = new struct coord [nAPs] ();
	dist = new double* [nUsers];
	alpha_dB = new double* [nUsers];
	r = new double* [nUsers];
	//allocating the second dimension for the 2-dimension arrays dist , alpha_dB and r
	for(i=0; i < nUsers; i++) {
		dist[i] = new double [nAPs] ();
		alpha_dB[i] = new double [nAPs] ();
		r[i] = new double [nAPs] ();
	}

	closest =  new int [nAPs] ();

	// invoke the correct discipline of deployment for UTs and APs
	switch(mode)
	{
	case 2:
		createWithTwoDensityZones();
		break;
	case 1:
		createWithRegularSquares();
		break;
	case 0:
	default:
		createFullRandom();
	}

	// compute distance
	for(j=0; j < nAPs; j++) {
		closest[j]=0;
		for(i=0; i < nUsers; i++) {
            dist[i][j] = sqrt( pow(user[i].x - ap[j].x, 2) + pow(user[i].y - ap[j].y, 2) );
            //looking for the closest user to AP j
			if ( dist[i][j] < dist[closest[j]][j] )
				closest[j] = i;
		}
	}
}


///Returns the reference path loss at distance d using "propagation in open space" formula.
/**
*The formula used is \n p = -10 * log10 ( antgain * antgain * pow( SL / (4 * PI * freq * 1000000 * d0), 2) )
* @see d0 antgain SL PI freq */
double computeRefPL(double d /*!< formal parameter for d0 */ )
{
	double g=antgain;
	//g = pow(10, antgain/10);
	double p = -10 * log10 ( g * g * pow( SL / (4 * PI * freq * 1000000 * d0), 2) );
	if(debug) cout<<"PL0 di riferimento a d0 = "<<p<<" dB\n";
	return p;
}

///Returns the freespace path loss.
/**
 * Log Distance Model is used. More precisely the formula used is <h5>-10 * log10 ( propF1 / pow (d, attenf) )</h5> \n where
 *  the propagation factor propF1 is computed as <h5> antgain * antgain * pow( SL / (4 * PI * freq * 1000000), attenf ) </h5>
 * @see d0 antgain SL PI freq attenf
 */
double computePL_freespace(double d /*!< formal parameter for d0 */ )
{
	double g = antgain;
	static double propF1 = 0;
	if (propF1 == 0) {
		propF1 = g * g * pow( SL / (4 * PI * freq * 1000000), attenf );
		if(debug) cout<<"fattore di propagazione = "<<propF1<<"\n";
	}
	return -10 * log10 ( propF1 / pow (d, attenf) );
}


///Returns path loss value [dB]
/**
 *  It has been used the Log Distance Model for N-LOS with multiple walls. Short formulation is:
 *  <h5> PL = PL0 + Lc + 10 n log (d/d0) + Nw Lw + Ncl Lcl </h5>
 * @see PL0 Lc d0 Nw Lw Ncl Lcl
 */
double computePL_logdistNLOS(double d /*!< it is a distance */ )
{
    if (PL0 < 0)
		PL0 = computeRefPL(d0);

	if (! LMW < 0)
		return PL0 + Lc + 10 * attenf * log10 ( d / d0 ) + LMW;
	else {
		int nw = Nw;
		int ncl = Ncl;
		if (nw < 0 || ncl < 0) {
			nw = (int)(d / wall_sep);
			ncl = (int)(d / col_sep);
		}
		return PL0 + Lc + 10 * attenf * log10 ( d / d0 ) + nw * Lw + ncl * Lcl;
	}
}


///Returns path loss value [dB]. It is used in the computation of the application-level data rate.
 /**
 * It uses Rayleigh Fading with unitary power for N-LOS environment. Short formulation is:
 *  <h5> PL = |h|^2 / d^gamma </h5>
 * For further details about the formula above please check the source code.
 */
double computePL_rayleigh(double d /*!< it is a distance */)
{
  //Rayleigh realization generation
  double u = (double)rand()/RAND_MAX;
  double h = sqrt(-2*pow(sigmaRayleigh,2) * log(u));
  double PL_dB = 20 * log10 (h) - computePathLoss(d);

  return PL_dB;

}
///Returns the average path loss value [dB]. It is used in the computation of the application-level data rate.
 /**
 *  It uses Rayleigh Fading with unitary power for N-LOS environment. \n
 * For further details about this function please check the source code.
 */
double computePL_rayleighAVG(double d /*!< it is a distance */ )
{
  double avg =  20 * log10(sigmaRayleigh) + 10 * log10 (PI/2) - computePathLoss(d);
  return avg;
}

///Returns the min95 path loss value [dB]. It is used in the computation of the application-level data rate.
 /**
 *  It uses Rayleigh Fading with unitary power for N-LOS environment. \n
 * For further details about this function please check the source code.
 */
double computePL_rayleighMIN(double d)
{
  double min95 = 20 *log10(sigmaRayleigh * sqrt(-2 * log(1-quantileRayleigh))) - computePathLoss(d);
  return min95;
}

///Computes the maximum transmission rate between UTs and APs.
/**
 * It populates the matrix r. For a user i and an AP j, the rate computed here (in kbps) is given by:
 *  <h5> r = max (0, min (R*, Rmax)) </h5> \n
 *  where: <b>R* = theta · PR / noise</b> if <i>SNIR(i,j) > gammaS</i> and <i>PR(i,j) > gammaP</i> (otherwise 0) \n
 *  or: <b>R* = beta · (PR - noise) + delta</b> if <i>SNIR(i,j) > gammaS</i> and <i>PR(i,j) > gammaP</i> (otherwise 0) \n
 *  where <b> PR = 10 * log10(Pmax) + alpha_dB(i,j) </b>
 * \n \n <b> PLEASE NOTICE THAT alpha_dB for whichever couple <i>i</i> and <i>j</i> WILL ALWAYS BE EQUAL TO 0.</b>
 * @see r Pmax alpha_dB zeta B noise theta beta delta r_max
 */
void computeMaxPhyRates()
{
    double pr, ru;

	for(int i=0; i < nUsers; i++) {
		for(int j=0; j < nAPs; j++) {

            pr = 10 * log10(Pmax) + alpha_dB[i][j];
            if (formula == 2) {
                r[i][j] = zeta * B * log(1 + pow(10, (pr-noise)/10));
                continue;
            }
            if ( pr < pwrThr )
				continue;
			if ( pr - noise < snirThr )
				continue;
            if (formula == 1)
				ru = theta * pow( 10, (pr-noise)/10 );
			else
				ru = (beta * (pr-noise) + delta);
			if (ru < 0)
				continue;
			r[i][j] = (ru > r_max ? r_max : ru);
		}
	}
}



///Returns the maximum application-layer rate [kbps] between a specific couple AP-UT.
/**
 * First the SNIR is computed: <b> SNIR_dB = PL_dB + 10 * log10 (Pt[0]) - noise </b> \n
 * Then the physical data rate is computed; finally application data rate is returned. \n
 * For further details about how this function is implemented, please take a look at the source code.
 */
double computeAppRate(double PL_dB /*!< path loss between a couple UT-AP */)
{
	double pr, ru;


    double SNIR_dB = PL_dB + 10 * log10 (Pt[0]) - noise;

    //compute physical data rate
    double r_phy = 0;
	pr = PL_dB + 10 * log10(Pt[0]);
    if(pr > pwrThr) {
        ru = beta * SNIR_dB + delta;
        if(ru > 0) {
            r_phy = (ru > r_max ? r_max : ru);
        }
    }

    // compute application data rate (in kbps)
    double r_app = 0;
    if(r_phy > 0) {
        r_app = msdu*8 * 1000 /
        (157.5 + 4 * (156 + msdu*8 + 272) / (4 * r_phy / 1000));
    }

    return r_app;
}


///Returns a new position for a certain UT.
/**
* The new position is computed taking randomly a value from a gaussian curve with mean (equal to the user position passed as parameter)
* and variance provided as parameters of the function.
*/
struct coord bivariateNormal(struct coord mean_pos/*!< position of the user <i>i</i>  */,
                            double sigma/*!< variance */){

  //new position
  struct coord new_pos;

  normal_distribution<double> x_distribution (mean_pos.x,sigma);
  normal_distribution<double> y_distribution (mean_pos.y,sigma);

  new_pos.x = x_distribution(generator);
  new_pos.y = y_distribution(generator);

  return new_pos;
}



///Returns N_pos candidate positions for each band for a given user.
/**
* For the current coordinates of a user, Nbands*N_pos candidate positions are generated as a whole. 
*/
struct coord ** computeCandidatePositions(struct coord current_pos /*!< current position of the UT */ , 
										double * rho /*!< array of radiuses for multiband model */)
{
  struct coord ** candidate_pos;
  candidate_pos = new struct coord * [Nbands] ();
  struct coord reference_pos;
  
  for(int j=0; j<Nbands; j++){
	  double radius;
	  candidate_pos[j] = new struct coord [N_pos]();
	  
	  if (j==0) 
		  radius = (double) rand()/RAND_MAX * rho[0]; 
	  else 
		  radius = (double) rand()/RAND_MAX * rho[0] + rho[j-1];
	  
	  double direction = (double) rand()/RAND_MAX * (2*PI);
	  reference_pos.x = radius * cos(direction);
	  reference_pos.y = radius * sin(direction);
	  
	  for(int k=0;k<N_pos;k++){
		 bool radiusOk = false;
		 while(true){		  
		 	reference_pos = bivariateNormal(reference_pos,varianceCurrent);
		 	candidate_pos[j][k].x = reference_pos.x + current_pos.x;
		 	candidate_pos[j][k].y = reference_pos.y + current_pos.y;
		 	//check correctness of new position
		 	radiusOk = (candidate_pos[j][k].x>=0 && candidate_pos[j][k].x<=fieldSizeX && candidate_pos[j][k].y>=0 && candidate_pos[j][k].y<=fieldSizeY);
		 	//take another direction , since this one isn't good
		 	if(!radiusOk){
		 		direction = (double) rand()/RAND_MAX * (2*PI);
		 		reference_pos.x = radius * cos(direction);
		 		reference_pos.y = radius * sin(direction);
		 	}
		 	//the candidate position generated is compatible with the simulation area :)
		 	else 
		 		break;
		 }
	  }	  
  }
  return candidate_pos;
}



///Returns a UT-band mapping for the future.
/**
* Of course they are choosen randomly according to upper bounds H[b] for b=0,1,...,|B|
* @see H
*/
vector<int> assignBandstoUTs(){
	
	int * available_places = new int[Nbands+1]();
	for(int i=0;i<=Nbands;i++){
		available_places[i] = H[i];
	}
	
	vector<int> user_set;
	for(int i=0;i<nUsers;i++){
		bool bandOk = false;
		int selected_band=-1;	
		do{
			selected_band = abs(rand() + seed)%(Nbands+1);
			if(available_places[selected_band]>0) bandOk = true;
		}while(!bandOk);
		
		available_places[selected_band]-=1;
		user_set.push_back(selected_band);
	}
	return user_set;
}


///Returns the minimum application data-rate given a set of candidate positions (relative to the band specified as input) and a specific AP.
/**
* The distance between the i-th candidate position and the associated AP (passed as parameter) is computed and their APP/data-rate is achieved. 
* \n After all candidate positions are considered, the minimum rate is returned ( aka r_min-b[i][j] ).
* @see computeAppRate
*/
double computeDevminRate(struct coord ** candidate_pos/*!< set of candidate positions for a certain user */, 
						 int band /*!<  */,
                         int AP_id/*!< the AP used for that set of candidate positions */) {
	
  //a huge value which represent infinity, for computing the min rate
  double min_rate = 99999999;
  //checking the minimum data rate among candidate positions
  for(int i=0; i<N_pos; i++) { 
    //distance between the candidate position and the AP provided as parameter 
    double distance = sqrt(pow(ap[AP_id].x - candidate_pos[band][i].x,2) + pow(ap[AP_id].y - candidate_pos[band][i].y,2));
    //computing the application data-rate given the distance computed above
    double potential_rate = computeAppRate(computePL_rayleighMIN(distance));
    if (potential_rate < min_rate)
      min_rate=potential_rate;
  }
  return min_rate;
}


///Selects a random position for a user among his/her candidate positions, for a certain band.
struct coord computeFuturePosition(struct coord * candidate_pos){ 

  struct coord future_pos;
  int selected_pos = rand() % N_pos;
  future_pos = bivariateNormal(candidate_pos[selected_pos],varianceFuture);

  if(future_pos.x>=0 && future_pos.x<=fieldSizeX && future_pos.y>=0 && future_pos.y<=fieldSizeY)
	  return future_pos;
  
  else
	  return candidate_pos[selected_pos]; 
}



///Creates the file containing all information about the generated instance. It works ONLY FOR <b>PSM=0</b>
/**
* The first computation done in the function is the normalization factor: <i> <b>(P0 + eta * Pt + mu * 3.750) / (eta * Pt) </b> </i>.
* Then, some information is computed and then reported in the output file "gwlan-input.dat": rho(max airtime); H[b] for each band (even for b=0); K; 
* parameters "b" and "p_w" (base and airtime power consumption) for each candidate site, computed respectively as 
* <i> <b>P0*norm_factor </b></i> and <i><b> (eta * dlFraction + P_demod * (1-dlFraction)) * norm_factor[j]</b> </i> ;
* demand values of each UT; data rate values for UTs in each possible candidate site, according to the 
* <i> <b>computePL_rayleighAVG </b> </i> function.
* \n Then, for each triple <i>(UT,AP,band)</i>, we select the minimum data rate among the candidate positions for the band 'b', 
* <i>rate_min-B(i,j)</i>, by calling <i> <b>computeDevminRate() </b> </i>. 
* \n Then, for each band b (even 0) choose H[b] users, computing related new positions and future data rates (using the 
* <i> <b>computePL_rayleigh </b></i> function).
* 
* @see P0 eta Pt mu max_airtime H K dlFraction P_demod computePL_rayleighAVG computePL_rayleighMIN computeDevminRate computePL_rayleigh
*/
void createFile()
{
	int i,j,tmp;
	double pwrStep;
	double p,p0,p1;
	double * rho; //will contain upper radius for each band b (rho[0] would be the max radius of the lowest mob. band)
	double maxavgload;
	char *s, *t;
	
    //Initializing some variables depending on the chosen power step mode
    switch (psm) {

        //parsing transmission power tabulation and storing in Pt
        case 0:{
            //transform a C++ string in a C string(array of characters)
            s = (char*)PT_tab.c_str();
            //create tokens for substrings separated by comma (remember that tabulation are something like "x,y,z,...")
            t = strtok(s, ",");
            i=1;
            Pt = new double [nPwrLevels] ();
            if(debug) cout<<"power tabulation is: ";

            while (t && i<=nPwrLevels) {
                Pt[i-1] = atof(t);
                if(debug) cout<<Pt[i-1]<<" ";
                i++;
                t = strtok(NULL, ",");
            }
            if(i <= nPwrLevels) {
                cout << "Errore: numero di valori PT_tab inferiore al numero di livelli richiesto!\n";
                exit(2);
            }
            if(debug) cout<<endl;
            break;
		}
		default: {
			cout << "Power step mode not allowed (psm="<<psm<<")! Implementation not yet provided for values different from 0.\n";
			exit(99);
        }
    }

    // Computing normalization factor
    double* norm_factor;
    norm_factor = new double [nAPs] ();
    for(int j=0; j < nAPs; j++){
        norm_factor[j] = (P0[apclass[j]] + eta[apclass[j]][0] * Pt[0] + mu[apclass[j]] * 3.750) / (eta[apclass[j]][0] * Pt[0]);
    }
    
    //computing different radiuses for the multiband model (rho[1] , rho[2] , ... )
    rho = new double [Nbands] ();
	double radiusFrac = (maxPercMov*fieldSizeX)/Nbands;
    for(int i=0;i<Nbands;i++){
    	rho[i] = radiusFrac * (i+1);
    }

    //computing cardinalities of bands (H[0] contains 'x' users , H[1] contains 'y' users , . . . )
    H[0]=ceil(0.5*nUsers); 
    int remaining_users = nUsers - H[0];
    normal_distribution<float> frac_distr(50,35);//average value=50 , variance=35
    //number of users in each band is completely random
    while(remaining_users>0){
    	float chosen_frac = abs(frac_distr(generator));
    	if(chosen_frac>50) continue;
    	int chosen_band = ((int)chosen_frac + rand())%Nbands +1;
    	int supplement = ceil(chosen_frac * remaining_users/100);
    	H[chosen_band] += supplement;
    	remaining_users -= supplement;
    }
    
    
    
    //Opening the file to write (name is "gwlan-input.dat")
	fu.open(fuName, ios::out);
	if(fu.fail()){
		cout << "Errore nella creazione del file "<< fuName <<"!\n";
		exit(2);
	}
	fu<<"### Green WLAN Input Parameters ###\n\n";

	fu<<"param rho := "<<max_airtime<<" ;\n\n";
	fu<<"param K := "<<K<<" ;\n\n";
	
	
	fu<<"param H[b] :=";
	for(j=0;j<=Nbands;j++){
		fu<<"\n\t"<<"H"<<j<<"\t"<<H[j];
	}
	fu<<" ;\n\n";


	fu<<"param: CANDIDATE_SITES: b := ";
	for(j=0; j < nAPs; j++) {
		fu<<"\n\t"<<"CS"<<j<<"\t"<<P0[apclass[j-1]]*norm_factor[j];
	}
	fu<<" ;\n\n";


	fu<<"param p_w :=";
	for(j=0; j < nAPs; j++) {
		fu<<"\n\t"<<"CS"<<j<<"\t"<<
        (eta[apclass[j]][0] * Pt[0] * dlFraction + P_demod[apclass[j]] * (1-dlFraction)) * norm_factor[j];
  }
	fu<<" ;\n\n";


	fu<<"# (Application level) Demand values w[i] are now in kb/s"<<"\n";
	fu<<"param: USER_TERMINALS: w :=";
	tmp = maxDemand - minDemand;
	for(i=1; i <= nUsers; i++) {
		p = minDemand + (rand() % tmp);
		p = p * oh_factor;	// convert to PHY rate
		fu<<"\n\t"<<"UT"<<i<<"\t"<<p;
	}
	fu<<" ;\n\n";


    fu<<"# (Application-level) Data rate values r_curr[i][j] are now in kb/s"<<"\n";
	fu<<"param r_curr:";
	for(i=1; i <= nAPs; i++)
		fu<<"\t CS"<<i;
	fu<<"\t:=";
	for(i=0; i < nUsers; i++) {
		fu<<"\n\t"<<"UT"<<i+1;
		for(j=0; j < nAPs; j++) {
			fu<<"\t"<< computeAppRate(computePL_rayleighAVG(dist[i][j]));
		}
	}
	fu<<" ;\n\n";
	

	//Candidate positions are computed for ALL USERS and all bands! 
    struct coord *** CandidatePositions = new struct coord** [nUsers];
    for(i=0; i < nUsers; i++) 
    		CandidatePositions[i] = computeCandidatePositions(user[i], rho);
    


    //This is computed assuming all UT move!!!
    fu<<"# (Application-level) Data rate values r_min-b[i][j] are now in kb/s"<<"\n";
	fu<<"param r_min-b :=";
	for(i=0; i < nUsers; i++) {
		fu<<"\n\t"<<"UT"<<i+1;
		for(j=0; j < nAPs; j++) {
			//r(b=0,f1)[i][j]
			fu<<"\n\t\tCS"<<j+1<<" [ "<<computeAppRate(computePL_rayleighMIN(dist[i][j]))<<" ";
			for(int b=0;b<Nbands;b++){
				//r(b>0,f1)[i][j]
				fu<<computeDevminRate(CandidatePositions[i], b, j);
				fu<<" ";
			}
			fu<<"]";
		}
	}
	fu<<" ;\n\n";

	
	//randomly assign each user a band (0th is included) according to H[b]
    mapping = assignBandstoUTs(); 
	
    //Here depending on the band of each user, the future rates are computed
    fu<<"# (Application-level) Data rate values r_future[i][j] are now in kb/s"<<"\n";
	fu<<"param r_future:";
	for(i=1; i <= nAPs; i++)
		fu<<"\t CS"<<i;
	fu<<"\t:=";
	for(i=0; i < nUsers; i++) {
		fu<<"\n\t"<<"UT"<<i+1;
		//checking whether this user is static 
		if(mapping[i]==0){
			for(int j=0;j<nAPs;j++) fu<<"\t"<<computeAppRate(computePL_rayleigh(dist[i][j]));
		}
		//depending on the band, choose the proper candidate position
		else{
			struct coord fut_pos = computeFuturePosition(CandidatePositions[i][mapping[i]-1]);
			for(int j=0;j<nAPs;j++){
				double fut_dist = sqrt(pow(ap[j].x-fut_pos.x,2) + pow(ap[j].y-fut_pos.y,2));
				fu<<"\t"<< computeAppRate(computePL_rayleigh(fut_dist));
			}
		}
    }
	fu<<" ;\n\n";
	
	fu.close(); 
}



///Saves in "extra.dat" coordinates of UTs; coordinates of APs; distances between i-th UT and j-th AP; mappings UT-band for the future; maximum radius of the greatest band
/** @see computeMaxPhyRates */
void saveExtraInfo()
{
	fstream fc;
	const char* fcName = "extra.dat";
	int i,j;

	fc.open(fcName, ios::out);
	if(fc.fail()){
		cout << "Errore nella creazione del file "<< fcName <<"!\n";
		return;
	}
	fc << "Field size = "<<fieldSizeX<<" x "<<fieldSizeY<<"\n\n";

	fc << "*** Coordinates of the users ***\n\n";
	for(i=0; i < nUsers; i++) {
		//fc<<"user["<<i+1<<"]: x="<<user[i].x<<", y="<<user[i].y<<"\n";
		fc<<"user["<<i+1<<"]: "<<user[i].x<<" "<<user[i].y<<"\n"; //GN
	}
	fc << "\n\n";

	fc << "*** Coordinates of the APs ***\n\n";
	for(j=0; j < nAPs; j++) {
		//fc<<"ap["<<j+1<<"]: x="<<ap[j].x<<", y="<<ap[j].y<<"\n";
		fc<<"ap["<<j+1<<"]: "<<ap[j].x<<" "<<ap[j].y<<"\n"; //GN
	}
	fc << "\n\n";

	fc << "*** Distances user-AP ***\n\n";
	for(j=0; j < nAPs; j++) {
		for(i=0; i < nUsers; i++) {
			fc<<"ap["<<j+1<<"]-user["<<i+1<<"] = "<<dist[i][j]<<"\n";
		}
	}
	fc << "\n\n";

	fc <<"*** Mapping users-bands ***\n\n";
	for(int i=0;i<nUsers;i++){
		fc<<"user["<<i+1<<"] --> band "<<mapping[i]<<"\n";
	}
	fc<<"\n\n";
	
	fc<<"Maximum radius of the greatest band = "<<maxPercMov*fieldSizeX<<" m";
	
	fc.close();
}


///Generates the random seed, parses configuration data, checks if the configuration is consistent, generates a scenario and saves the instance.
/**
 * Parameters for the main are:
 * - -d = debug mode, which displays verbose information
 * - -s 'integer' = passing an integer for the seed generation (helps good randomization)
 * \n\n For reading and parsing configurations we use the <b> readInput </b> function. Depending on <b> pl_model </b>, <b>computePathLoss</b> is 
 * initialized with the proper PL computing function (either <b>computePL_logdistNLOS</b> or <b>computePL_freespace</b>). The instance is generated
 * by calling computeInstance() and finally <b> createFile() </b> and <b> saveExtraInfo() </b> are invoked.
 * 
 * @see readInput pl_model computePathLoss computePL_logdistNLOS computePL_freespace computeInstance createFile saveExtraInfo
 */
int main(int argc, char** argv)
{
	int i,j;


	int c=1, inst;
	while (c < argc) {
		if(strcmp(argv[c], "-d") == 0)
			debug = 1;
		else if (strcmp(argv[c], "-s") == 0)
			inst = atoi(argv[++c]);
		else {
			cout << "Unknown parameter: " << argv[c] << endl;
			return 1;
		}
		c++;
	}
	srand(time(NULL));
	seed = abs(inst * rand());
	default_random_engine random_number(seed);
	generator = random_number;
	if(debug) cout << "Seed = "<< seed <<endl;

	readInput();

  	// check that problem size is greater than zero
	if ( nUsers < 1 || nAPs < 1 || nPwrLevels < 1) {
		cout << "La dimensione di una o più delle varibili del problema è minore di uno! Esecuzione terminata!\n";
		exit(3);
	}

		// check traffic constraint: do not create an instance that will be unfeasible for sure
	if ( nUsers*minDemand > nAPs*goodput_max ) {
		cout << "La domanda è maggiore della capacità del sistema, lo scenario non è fattibile!\n";
		exit(3);
	}

	// check multi-wall parameters
	if (( Nw < 0 || Ncl < 0) && !(wall_sep > 0 && col_sep > 0)) {
		cout << "Errore: parametri MW-NLOS non consistenti! (Nw, Ncl, wall_sep, col_sep = "<<Nw<<Ncl<<wall_sep<<col_sep<<")\n";
		exit(3);
	}

	// set propagation model

	if ( pl_model = 1)
		computePathLoss = computePL_logdistNLOS;
	else
		computePathLoss = computePL_freespace;

	// start computations
	computeInstance(deployMode);

	// check closeness
	if (minDist > 0)
		adjustDistance();

	createFile();

	saveExtraInfo();

	//cout << "**********************Istanza creata con successo!***********************\n";
	return 0;
}

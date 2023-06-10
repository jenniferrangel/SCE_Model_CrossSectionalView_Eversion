//============================================================================
// Name        : Main.cpp
// Author      : Wenzhao Sun, Ali Nematbakhsh, Kevin Tsai
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

//#include "MeshGen.h"
#include "commonData.h"
#include "CellInitHelper.h"
#include <vector>
#include "SimulationDomainGPU.h"

using namespace std;

GlobalConfigVars globalConfigVars;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort =
		true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
				line);
		if (abort)
			exit(code);
	}
}
//test here

/////////////////////////////////////////////////////////////
/////////////// Reading configuration file //////////////////
/////////////////////////////////////////////////////////////
void initializeSlurmConfig(int argc, char* argv[]) {
	// read configuration.
	ConfigParser parser;
	std::string configFileNameDefault = "./resources/disc_M.cfg";
	globalConfigVars = parser.parseConfigFile(configFileNameDefault);
	std::string configFileNameBaseL = "./resources/disc_";
	std::string configFileNameBaseR = ".cfg";

	// Unknown number of input arguments.
	if (argc != 1 && argc != 3) {
		std::cout << "ERROR: Incorrect input argument count.\n"
				<< "Expect either no command line argument or three arguments"
				<< std::endl;
		exit(0);
	}
	// one input argument. It has to be "-slurm".
	else if (argc == 3) {
		if (strcmp(argv[1], "-slurm") != 0) {
			std::cout
					<< "ERROR: one argument received from commandline but it's not recognized.\n"
					<< "Currently, the argument value must be -slurm"
					<< std::endl;
			exit(0);
		} else {
			std::string configFileNameM(argv[2]);
			std::string configFileNameCombined = configFileNameBaseL
					+ configFileNameM + configFileNameBaseR;
			parser.updateConfigFile(globalConfigVars, configFileNameCombined);
		        int myDeviceID =
				globalConfigVars.getConfigValue("GPUDeviceNumber").toInt();
		        gpuErrchk(cudaSetDevice(myDeviceID));
		}
	}
	// no input argument. Take default.
	else if (argc == 1) {

		// set GPU device.
		int myDeviceID =
				globalConfigVars.getConfigValue("GPUDeviceNumber").toInt();
		gpuErrchk(cudaSetDevice(myDeviceID));
	}
}


////////////////////////////////////////////////////
////////// Calculate division threshold ////////////
////////////////////////////////////////////////////
// Used in the apical view model to determine the threshold value for division to occur
// Called in the current script but the information is not used in the current simulation
void updateDivThres(double& curDivThred, uint& i, double& curTime,  //Ali
		double& decayCoeff, double& divThreshold) {
        
        //cout<<"The value of initial time stage in updateDivThres is"<<curTime<<endl ;  
	double decay = exp(-curTime * decayCoeff);
	//curDivThred = 1.0 - (1.0 - divThreshold) * decay;
	curDivThred = divThreshold ;
}


/////////////////////////////////////////////////////
//////////// Start of the simulation ////////////////
/////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
	// initialize random seed.
	srand(time(NULL));

	// Slurm is computer-cluster management system.
	initializeSlurmConfig(argc, argv);

        cout<< "I am in main file after slurm "<<endl; 
	// initialize simulation control related parameters from config file.
	SimulationGlobalParameter mainPara;

        cout<< "I am in main file after simulation global parameter "<<endl; 
	mainPara.initFromConfig();

        cout<< "I am in main file before Cell IniHelper instance creation"<<endl; 
	// initialize simulation initialization helper.
	CellInitHelper initHelper;
        cout<< "I am in main file after Cell IniHelper instance creation"<<endl; 

	// initialize simulation domain.
	SimulationDomainGPU simuDomain;

        cout<< "I am in main file after simulationDomainGPU instance creation"<<endl;

	// comment: initInput_M() goes in cellInitHelper.cpp. There the information are read from input txt files and saved in CPU of the computer.
	// Then initialize_v2_M will go to the simuDomain.cu and initially create some vector spaces for GPU values and then assigned the CPU values to GPU vectors either in SceNodes or SceCells. The next two functions below are the main functions for initialization of the code from txt and transfer them to GPU.
	SimulationInitData_V2_M initData = initHelper.initInput_M(); // it will go inside cellinithelper.cpp  Ali 

        cout<< "I am in main file after initInput_M creation"<<endl; 
	simuDomain.initialize_v2_M(initData,mainPara.InitTimeStage);
	
        cout<< "I am in main file after initInput_v2_M creation"<<endl; 
	std::string polyStatFileNameBase = globalConfigVars.getConfigValue(
			"PolygonStatFileName").toString();
	std::string uniqueSymbol =
			globalConfigVars.getConfigValue("UniqueSymbol").toString();
	std::string polyStatFileName = polyStatFileNameBase + uniqueSymbol + ".txt";

	std::remove(polyStatFileName.c_str());

	std::string detailStatFileNameBase = globalConfigVars.getConfigValue(
			"DetailStatFileNameBase").toString() + uniqueSymbol;
	double divThreshold =
			globalConfigVars.getConfigValue("DivThreshold").toDouble();
	double decayCoeff =
			globalConfigVars.getConfigValue("ProlifDecayCoeff").toDouble();
	double curDivThred;

	int maxStepTraceBack =
			globalConfigVars.getConfigValue("MaxStepTraceBack").toInt();

	uint aniFrame = 0;
	// main simulation steps.
    std::string stressStrainFileNameBase="StressStrain" ;
    std::string stressStrainFileName=stressStrainFileNameBase +uniqueSymbol+ ".CSV" ; 
    SingleCellData singleCellData(stressStrainFileName);
	double current_Time;
	double total_Time = mainPara.totalTimeSteps*mainPara.dt;
	std::cout<<"mainPara.totalTimeStpes = "<<mainPara.totalTimeSteps<<std::endl;
	std::cout<<"mainPara.dt = "<<mainPara.dt<<std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// These variables are needed for several functions in the current code but is not used. It was from an /
	// earlier version of cell division algorithm. Do not delete without making changes in other files and //
	// functions that still take this as an input. //////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	double timeRatio;
	double timeRatio_Crit_actomyo = 2.0;
	// std::cout<<"Critical timeRatio for actomyosin strength reduction = "<<timeRatio_Crit_actomyo<<std::endl;
	bool reduced_actomyo_triggered = false;
	double timeRatio_Crit_ECM = 2.0;
	// std::cout<<"Critical timeRatio for ecm strength reduction = "<<timeRatio_Crit_ECM<<std::endl;
	double timeRatio_Crit_Division = 0.2;
	// std::cout<<"Critical timeRatio for growth and cell division = "<<timeRatio_Crit_Division<<std::endl;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	bool reduced_ecm_triggered = false;
	// Initialization of the true-false marker for the utilization of anisotropic ECM stiffness.
	
	bool Division_triggered = false;
	// Initialization of the true-false marker for the occurrence of cell division.
	
	double volume_Increase_Target_Ratio = 1.4;
	std::cout<<"Target ratio for cell volume increase = "<<volume_Increase_Target_Ratio<<std::endl;
	// Target cell volume increase during mitotic rounding.

	double volume_Increase_Scale = 1.0;
	std::cout<<"How fast is the volume increase happening = x"<<volume_Increase_Scale<<" rate of change"<<std::endl;
	// A scale that determines how fast target cell volume during mitotic rounding will be achieved.

	double postDivision_restorationRateScale = 1.0;//0.5;//0.1;
	std::cout<<"How fast is the volume restoration happening post division = x"<<postDivision_restorationRateScale<<" rate of change"<<std::endl;
	// A scale that determines how fast the divided cell will go back to its original equilibrium cell volume.
	
	bool volume_restoration_rate_restore = false;

	double thresholdToIntroduceNewCell = 0.25;//0.5;//1.0;//0.25;//-1.0;
	std::cout<<"The likelihood (probability) a new cell will be introduced in the same cross section after cell division = "<<thresholdToIntroduceNewCell<<std::endl;

	double mitoticThreshold = 0.95;//0.8973;
	std::cout<<"mitoticThreshold = "<<mitoticThreshold<<std::endl;
	// The counter value that needs to be reached for division to occur. The max value of the counter is 1.

	double distFromNucleus_max = 0.0;//4.0;
	double distFromNucleus_min = 0.3;//-14.0;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////// These values determine the portion of cells occupied by the contractile springs /////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	double distFromNucleus_normalMax1 = 0.3;//0.275;//0.20;//0.355;//-8.0;
	double distFromNucleus_normalMax2 = 0.2645;//0.275;//0.20;//0.355;//-8.0;
	double distFromNucleus_normalMax3 = 0.2615;//0.275;//0.20;//0.355;//-8.0;
	double distFromNucleus_normalMax_apical1 = 0.2918;//0.175;//0.20;//0.295;//7.5;
	double distFromNucleus_normalMax_apical2 = 0.2883;//0.175;//0.20;//0.295;//7.5;
	double distFromNucleus_normalMax_apical3 = 0.2429;//0.175;//0.20;//0.295;//7.5;
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	double percentage_before_timeRatio_Crit_Division_scaling = 4.0; //No longer in use, but don't delete yet since it is still passed into several functions. //Kevin
	double mitoRndActomyoStrengthScaling = 5.0; //Please consult SceNodes.cu to see what are the scaling applied to non-dividing cells and disc_NXX_X.cfg file to see what is the corresponding default spring constant. //Kevin
	std::cout<<"Cell division requires the contractile spring to increase strength by "<<percentage_before_timeRatio_Crit_Division_scaling<<" fold."<<std::endl;
	std::cout<<"Contractile spring minimum at "<<distFromNucleus_min<<" and maximum at "<<distFromNucleus_max<<" away from the cell center."<<std::endl;
	std::cout<<"But under non-growth (stationary) circumstances, the basal contractile spring is set at "<<distFromNucleus_normalMax1<<", "<<distFromNucleus_normalMax2<<", "<<distFromNucleus_normalMax3<<"*cellheight away from the cell center."<<std::endl;
	std::cout<<"But under non-growth (stationary) circumstances, the apical contractile spring is set at "<<distFromNucleus_normalMax_apical1<<", "<<distFromNucleus_normalMax_apical2<<", "<<distFromNucleus_normalMax_apical3<<"*cellheight away from the cell center."<<std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////// Weights applied to the coefficients of the contractile springs ///////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
		//Example: to set the contractility higher in the medial region of the pouch as compared to the lateral sides (1 and 3), use the following values:
			//*Basal: 1.5 | 2.0 |1.5
			//*Apical: 0.0625| 1.333 | 0.0625
	double contractActomyo_multip_perCell1 = 1.5;  
	double contractActomyo_multip_perCell2 = 2.0;  
	double contractActomyo_multip_perCell3 = 1.5;  
	double contractActomyo_multip_perCell_apical1 = 0.0625;
	double contractActomyo_multip_perCell_apical2 = 1.333;
	double contractActomyo_multip_perCell_apical3 = 0.0625;
	std::cout<<"Basal_left_actomyo_weight = "<<contractActomyo_multip_perCell1<<std::endl;
	std::cout<<"Basal_mid_actomyo_weight = "<<contractActomyo_multip_perCell2<<std::endl;
	std::cout<<"Basal_right_actomyo_weight = "<<contractActomyo_multip_perCell3<<std::endl;
	std::cout<<"Apical_left_actomyo_weight = "<<contractActomyo_multip_perCell_apical1<<std::endl;
	std::cout<<"Apical_mid_actomyo_weight = "<<contractActomyo_multip_perCell_apical2<<std::endl;
	std::cout<<"Apical_right_actomyo_weight = "<<contractActomyo_multip_perCell_apical3<<std::endl;
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////

	double growthProgressSpeed = 0.5*3.84e-7;//3.84e-7;//(0.02)*0.25*(0.001*0.002); 
	// 0.002 is the default simulation time step size. The value above indicate how much progress is gained per time step.
	std::cout<<"Growth speed for cell cycle per time step is = "<<growthProgressSpeed<<std::endl;

	int maxApicalBasalNodeNum = 18;//9999;
	std::cout<<"Max number of apical and basal nodes, respectively, for columnar cells = "<<maxApicalBasalNodeNum<<std::endl;
	// int minApicalBasalNodeNum = 21;
	// std::cout<<"Min number of apical and basal nodes, respectively, for columnar cells = "<<minApicalBasalNodeNum<<std::endl;

	double maxLengthToAddMemNodes = 0.195;//0.26;
	std::cout<<"Max length for each edge to qualify for new node additions = "<<maxLengthToAddMemNodes<<std::endl;
	// This is restricted to apical and basal side of the cell only, due to the constraint of the number of nodes placed laterally.

	std::vector<int> cycle_vec;
	cycle_vec.push_back(-1);
	int cycle;

	double fixed_dt = mainPara.dt;
	for (int i = 0; i < cycle_vec.size(); i++){	
		Division_triggered = false;
		cycle = cycle_vec[i];
		for (uint i = 0; i <= (uint) (mainPara.totalTimeSteps); i++) {
			
			// current_Time = i*mainPara.dt + mainPara.InitTimeStage;
			current_Time = i*fixed_dt + mainPara.InitTimeStage;
			timeRatio = current_Time/total_Time;
			
			// this if is just for output data// 
			if (i % mainPara.aniAuxVar == 0) {
				double curTime=i*mainPara.dt + mainPara.InitTimeStage;  //Ali - Abu
				updateDivThres(curDivThred, i, curTime, decayCoeff,divThreshold);

				CellsStatsData polyData = simuDomain.outputPolyCountData();  //Ali comment
				std::cout<<"survived outputPolyCountData()?"<<std::endl;
				singleCellData=simuDomain.OutputStressStrain() ;              
				singleCellData.printStressStrainToFile(stressStrainFileName,curTime) ;
							
				// prints brief polygon counting statistics to file
				polyData.printPolyCountToFile(polyStatFileName, curDivThred);
				// prints detailed individual cell statistics to file
				polyData.printDetailStatsToFile(detailStatFileNameBase, aniFrame);
				// prints the animation frames to file. They can be open by Paraview
				simuDomain.outputVtkColorByCell_polySide(mainPara.animationNameBase,
						aniFrame, mainPara.aniCri);
				
				simuDomain.outputResumeData(aniFrame) ; 
				aniFrame++;
			}
			simuDomain.runAllLogic_M(mainPara.dt,mainPara.Damp_Coef,mainPara.InitTimeStage,
									timeRatio, timeRatio_Crit_actomyo, timeRatio_Crit_ECM, timeRatio_Crit_Division,
										volume_Increase_Target_Ratio, volume_Increase_Scale, postDivision_restorationRateScale, cycle,
										distFromNucleus_max, distFromNucleus_min, distFromNucleus_normalMax1, distFromNucleus_normalMax2, distFromNucleus_normalMax3,
										distFromNucleus_normalMax_apical1, distFromNucleus_normalMax_apical2,distFromNucleus_normalMax_apical3,
										percentage_before_timeRatio_Crit_Division_scaling, growthProgressSpeed, maxApicalBasalNodeNum, maxLengthToAddMemNodes,mitoRndActomyoStrengthScaling, thresholdToIntroduceNewCell,
										contractActomyo_multip_perCell1, contractActomyo_multip_perCell2, contractActomyo_multip_perCell3,
										contractActomyo_multip_perCell_apical1, contractActomyo_multip_perCell_apical2, contractActomyo_multip_perCell_apical3, mitoticThreshold);  //Ali //Kevin
		}
	}

	return 0;
}

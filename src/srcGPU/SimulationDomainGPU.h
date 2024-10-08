#ifndef SimulationDomainGPU_H_
#define SimulationDomainGPU_H_

//#include "SceNodes.h"
//#include "SceECM.h"
#include "SceCells.h"
#include "commonData.h"
#include "NetworkInfo.h"
#include "Solver.h"
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std; 
/**
 * This class will control the process of stabilizing cell centers.
 */
class StabPara {
public:
	bool isProcessStab;

	int outputFrameCount;
	int totalIterCount;
	double bdrySpacingRatio;
	//double dt;
	std::string outputAniName;
};

/**
 * This class is responsible for domain-wise highest level logic, e.g. output animation.
 */

class SimulationDomainGPU {
	/**
	 * Variable that contains information for nodes.
	 * Handles node level interaction logic.
	 */

	Solver solver ; 
	SceNodes nodes;

	SceECM   eCM ;  // Ali 
	/**
	 * Variable that contains information for cells.
	 * Handles cell level logics like growth and division.
	 */
	SceCells cells;

	NetworkInfo netInfo;

	std::vector<std::vector<PreT1State> > preT1Vec;

	std::set<int> t1CellSet;
	std::vector<double> cellColorVec;

	/**
	 * memory related parameters.
	 */
	SceMemPara memPara;

	/**
	 * domain related parameters.
	 */
	SceDomainPara domainPara;

	/**
	 * chemical related parameters.
	 */
	SceChemPara chemPara;

	/**
	 * parameters used for stabilize the initial cell positions
	 */
	StabPara stabPara;

	/**
	 * reads memory related parameters.
	 */
	void readMemPara();

	/**
	 * reads domain related parameters.
	 */
	void readDomainPara();

	/**
	 * reads all parameters by calling all other reading methods.
	 */
	void readAllParameters();

	/**
	 * Initializes data vectors by given vectors.
	 * improved from the previous version.
	 */
	void initializeNodes_M(std::vector<SceNodeType> &nodeTypes,
			std::vector<bool> &nodeIsActive, std::vector<CVector> &initNodesVec,std::vector<CVector> &initNodeMultip_actomyo, std::vector<CVector> &initNodeMultip_integrin,
			std::vector<uint> &numOfInitActiveEpiNodeCounts,
			std::vector<uint> &numOfInitActiveInternalNodeCounts,
			std::vector<double> &initGrowProgVec, 
			std::vector<ECellType> & eCellTypeV1,
			std::vector<double> & mDppV,
			std::vector<MembraneType1> & mTypeV
			,double InitTimeStage);

	NetworkInfo buildNetInfo(CellsStatsData &polyData);
	std::set<int> findT1Transition();

	void outputVtkGivenCellColor(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri, std::vector<double>& cellColorVec, std::vector<double> & cellsPerimeter);  //AliE
	std::vector<double> processT1Color();
	std::vector<double> processPolySideColor(std::vector<double> & cellsPerimeter); //AliE

public:
	/**
	 * Default constructor.
	 * Reads all the configuration values from file.
	 */

	SimulationDomainGPU();

	/**
	 * Domain initialization.
	 * Assigns values to the data fields in simulation domain.
	 * @param initData initial data set for simulation domain
	 */
	void initialize_v2_M(SimulationInitData_V2_M &initData, double InitTimeStage);

	/**
	 * Run one step of simulation in the domain.
	 * Contains cell level logics and node level logics.
	 * @param dt timestep
	 */
	void runAllLogic(double dt);

	/**
	 * Run one step of simulation in the domain.
	 * Contains cell level logics and node level logics.
	 * @param dt timestep
	 */
//Ali 	void runAllLogic_M(double dt);
	// void runAllLogic_M(double & dt,double Damp_Coef,double InitTimeStage, 
	// 						double timeRatio, double timeRatio_Crit_actomyo, double timeRatio_Crit_ECM, double timeRatio_Crit_Growth,
	// 							double volume_Increase_Target_Ratio, double volume_Increase_Scale, double postDivision_restorationRateScale, int cycle,
	// 							double distFromNucleus_max, double distFromNucleus_min, double distFromNucleus_normalMax, double distFromNucleus_normalMax_apical, double percentage_before_timeRatio_Crit_Division_scaling,
	// 							double growthProgressSpeed, int maxApicalBasalNodeNum, int minApicalBasalNodeNum, double maxLengthToAddMemNodes);
	void runAllLogic_M(double & dt,double Damp_Coef,double InitTimeStage, 
							double timeRatio, double timeRatio_Crit_actomyo, double timeRatio_Crit_ECM, double timeRatio_Crit_Growth,
								double volume_Increase_Target_Ratio, double volume_Increase_Scale, double postDivision_restorationRateScale, int cycle,
								double distFromNucleus_max, double distFromNucleus_min, double distFromNucleus_normalMax1, double distFromNucleus_normalMax2, double distFromNucleus_normalMax3,
								double distFromNucleus_normalMax_apical1, double distFromNucleus_normalMax_apical2, double distFromNucleus_normalMax_apical3,
								double percentage_before_timeRatio_Crit_Division_scaling,
								double growthProgressSpeed, int maxApicalBasalNodeNum, double maxLengthToAddMemNodes, double mitoRndActomyoStrengthScaling, double thresholdToIntroduceNewCell,
								double contractActomyo_multip_perCell1, double contractActomyo_multip_perCell2, double contractActomyo_multip_perCell3,
								double contractActomyo_multip_perCell_apical1, double contractActomyo_multip_perCell_apical2, double contractActomyo_multip_perCell_apical3, double mitoticThreshold);
	bool isDividing_ForAni();

	/**
	 * Method that animates the domain to VTK format.
	 * @param scriptNameBase name of the vtk animation series.
	 * @param rank frame sequence in the vtk animation series.
	 * @param aniCri criteria for outputing animation.
	 */
	void outputVtkFilesWithCri(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri);

	/**
	 * Method that animates the domain to VTK format.
	 * @param scriptNameBase name of the vtk animation series.
	 * @param rank frame sequence in the vtk animation series.
	 * @param aniCri criteria for outputing animation.
	 */
	void outputVtkFilesWithCri_M(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri);

	void outputVtkColorByCell_T1(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri);

	void outputVtkColorByCell_polySide(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri);

	void outputResumeData(uint frame);

	/**
	 * Method that output the simulation domain as label matrix.
	 * @param resultNameBase name of the labelMatrix series.
	 * @param rank frame sequence in the labelMatrix series.
	 * @param pixelPara criteria for labelMatrix generation.
	 */
	vector<vector<int> > outputLabelMatrix(std::string resultNameBase, int rank,
			PixelizePara &pixelPara);

	/**
	 * method that prints out the growth progress vector of cells to a file.
	 */
	void outputGrowthProgressAuxFile(int step);

	/**
	 * Post processing for the label matrix.
	 */
	void analyzeLabelMatrix(vector<vector<int> > &labelMatrix, int step,
			std::string &imageFileNameBase, std::string &statFileName);

	void performAblation(AblationEvent &ablEvent);

	CellsStatsData outputPolyCountData();
	SingleCellData OutputStressStrain();

	void processT1Info(int maxStepTraceBack, CellsStatsData &polyData);
};



#endif

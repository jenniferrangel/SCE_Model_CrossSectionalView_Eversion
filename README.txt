DOI: https://zenodo.org/badge/latestdoi/314704785

CUDA code for developmental biology using Subcellular Element Method

Hardware requirement: 
Nvidia video card that supports SM 2.0+ and CUDA 4.0 

Software environment requirement: 
CMAKE ----- Build system.
CUDA  ----- Provide runtime support for parallel GPU computation.
CGAL  ----- Computational geometry library.
Thrust ---- Build-in library of cuda, similar to STL of C++
Paraview -- (Optional) Visualization software for animation purpose. 

*******************************************
Location of configuration files:
 ./resources
 
Additional parameter-configurating files:
  (1) disc_NXX_X.cfg, X is some number, contains:
        animation (vtk) files storage designation.
        animation files name desgination.
        initial contractile spring coefficient.
  (2) disc_M.cfg
        initial parameter for membrane stiffness, bending stiffness, volume exclusion, simulation time, time step size, etc.
  (3) ECM_input.cfg
        ECM stiffness, adhesion with membrane node, volume exclusion with membrane node.
  (4) discMain_M.cpp
        growth speed, anisotropic contractile spring coefficients, etc.
  (5) SceCells.cu
        cell cycle length
        
*******************************************
How to identify the flow of algorithm calling and simulation proceeding?

1. Look for : int main(int argc, char* argv[])
2. Find : simuDomain.runAllLogic_M
3. Go to SimulationDomainGPU.cu and SimulationDomainGPU.h in /src/SceGPU subfolder 
4. Find : SimulationDomainGPU::runAllLogic_M(double & dt...)
5. You will see that in this function, the followings key algorithms are executed:
   (a) nodes.sceForcesDisc_M
     (i) Go to SceNodes.cu and SceNodes.h in /src/SceGPU subfolder for details
   (b)	cells.runAllCellLogicsDisc_M
     (i) Go to SceCells.cu and SceCells.h in /src/SceGPU subfolder for details


For (a):
  This function deals with: (a) Identifying node-to-node correspondence for contractile springs and cell-cell adhesion, (b) Volume exclusion (via Morse-like interaction) between membrane nodes.
  In SceNodes::sceForcesDisc_M, the following key algorithms are executed:
    (I) prepareSceForceComputation_M
    (II) Setting up the coefficients of the apical and basal contractile springs.
    (III) applySceForcesDisc_M
      (i) Identify membrane nodes connected by cell-cell adhesion forces.
      (ii) Identify the membrane nodes connected by contractile springs.
      (iii) Calculate forces due to membrane-to-membrane volume exclusion (or repulsion) using the function : AddForceDisc_M
    (IV) processMembrAdh_M
      (i) Calculate adhesion between membrane nodes using : applyMembrAdh_M

  For (b):
     (I) Assignment for the number of contractile spring for each cell.
     (II) Speed counter for growth progress and etc are set.
     In cells.runAllCellLogicsDisc_M, the following key algorithms are executed:
         (III) Populating nucleus node for dividing cells.
         (IV) computeApicalLoc
         (V) computeBasalLoc
         (VI) computeIndividualCellHeight
         (VII) eCMCellInteraction
         (VIII) computeCenterPos_M2: computes centers of cells.
         (IX) computeInternalAvgPos_M(): computes average nucleus position of each cell.
         (X) applySceCellDisc_M: volume exclusion between membrane and internal (nucleus) nodes.
         (XI) applyMembContraction2: calculate contractile spring force.
         (XII) applyMemForce_M
         (XIII) applyVolumeConstraint
         (XIV) growAtRandom_M_Ver2: advance the growth progress counter with specific conditions related to neighboring cells.
         (XV) allComponentsMove: update position of each node based on previously calculated forces.
 

*******************************************
To run simulation on UCR HPCC cluster: 
   After uploading all the files under this repository, follow the command below.
   (1) cd [The folder where all your files are placed in]
   (2) module load singluarity
   (3) module load centos
   (4) centos.sh
   (5) module load extra; module load GCC; module load cuda;
   (6) (IF it is the first time compiling) module load CGAL; module load cmake;
   (7) (IF it is the first time compiling) cmake .
   (8) make
   (9) After compilation, enter: exit
   (10) sbatch -p gpu --gres=gpu:1 --time=X:00:00 EpiScale_run.sh; (X here is the number of hours you want to keep the simulation running)
   
********************************************
[Data processing]
ExtractApiBasNucInfo.m : used to extract the apical node, basal node, and nucleus center information.
CalcApicalCurvature.m : used to calculate the apical curvature using the apical node of each individual cell.
CalcBasalCurvature.m : used to calculate the basal curvature using the basal node of each individual cell.
Both CalcApicalCurvature and CalcBasalCurvature include cell height calculation.
The curvature calculated here is the Menger curvature.

1. The input for cell information has to be extracted manually from files with a "detailedStat_" prefix. 
   Every triplet (row-wise, i.e. row1 ~ row3) represents the apical node, basal node, and nucleus center 
   position, respectively.
2. For simulations involving growth, the order of cells must be entered to acquire the correct curvature calculation. 
   The cell order can be found in the printout created (see EpiScale_run.sh). For the detail on what information is 
   presented in the files with "detailedStat_" prefix, see detailedStat_example.txt in the dataOutput folder.

CUDA code for developmental biology using Subcellular Element Method

Hardware requirement: 
Nvidia video card that supports SM 2.0+ and CUDA 4.0 

Software environment requirement: 
CMAKE ----- Build system.
CUDA  ----- Provide runtime support for parallel GPU computation.
CGAL  ----- Computational geometry library.
Thrust ---- Build-in library of cuda, similar to STL of C++
Paraview -- (Optional) Visualization software for animation purpose. 

To compile:
 (1) In project root folder, type "cmake ." ("sudo cmake ." preferred)
 (2) type "make" 
Please note that CMake, CUDA, CGAL, Thrust, are all required for compilation.  
%%% During this process, make sure SET(BUILD_DOCUMENTATION OFF) is set to OFF in the CMakeLists.txt to avoid compile error.

To run unit test from project root folder:
 Option 1: Simple run: type "make test"
 Option 2: See more details about unit test: type "./bin/UnitTest"
 
To run performance test from project root folder:
 In project root folder, type "./bin/PerfTest"

To run simulation:
 In project root folder, type "./bin/run***Simulation"
 Currently, two simulations are available: Beak and Disc.


************************
To run simulation on slurm cluster (acms-gpu is powered by slurm) 
 (1) In project root folder, cd ./scripts
 (2) sbatch *.sh, for example, sbatch discN01G02.sh means take 
     the first configuration file and then submit it to gpu02 compute node 
     so. The actual GPU device number that is going to run the program is 
     controled by slurm, and specified GPUDevice in the config file is ignored 
     if running on cluster.

Location of configuration files:
 ./resources

************************************************************************************************

For compilation on UCR hpcc:
  module load cmake; module load CGAL
are needed for "cmake ."

Module loading:
  module load extra; module load GCC; module load cuda
are needed for "make" compilation.
Compilation must be done on an interactive GPU session otherwise the there is insufficient 
computing power to compile the code.

Please refer to UCR HPCC "managing jobs" section for job submission.

================================================================================================
How to manipulate specific parameters?
1. Membrane stiffness is found in disc_M.cfg file in /resources
2. Default contractile spring coefficient is found in disc_NXX_X.cfg in /resources
3. Contractile spring coefficient multiplier is found in SceNodes.cu by searching : "infoVecs.contractActomyo_multip[i] ="
4. Equilibrium area, speed of increase in the number of contractile spring, and population of nucleus particles are found in SceCells.cu. These be found by searching : "ssspeed =" 
5. discMain_M.cpp also contains many parameters used in the simulations. Too many to be listed out and explained here effectively.
      a. Parameters that controls the number of apical and basal actomyosin contractile springs (via ratios with respect to individual cell heights).
      b. Strenght of contractile spring for generating a narrow basal section when a cell undergoes mitotic rounding.
      c. Max number of apical and basal node per columnar cells.
      d. Max extending of each edge needed before a new node needs to be inserted.

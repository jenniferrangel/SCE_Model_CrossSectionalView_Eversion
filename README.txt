CUDA code for developmental biology using Subcellular Element Method

Hardware requirement: 
Nvidia video card that supports SM 2.0+ and CUDA 4.0 

Software environment requirement: 
CMAKE ----- Build system.
CUDA  ----- Provide runtime support for parallel GPU computation.
CGAL  ----- Computational geometry library.
Thrust ---- Build-in library of cuda, similar to STL of C++
Paraview -- (Optional) Visualization software for animation purpose. 

************************************************************************************************
************************************************************************************************
************************************************************************************************
COMPILATION PROCESS FOR UCR HPCC COMPUTING SERVICE
After uploading all the files contained in this repository onto your HPCC account, proceed with the following commands:
1. Log into your HPCC account using either a terminal command prompt or MOBAXTERM
2. cd [The folder where you place the files from this repository]
3. module load singularity;
4. module load centos;
5. centos.sh;
6. module load extra; module load GCC; module load cuda;
7. (IF YOU HAVEN'T COMPILE THE CODE YET IN THIS DIRECTORY) module load cmake; module load CGAL;
8. (IF YOU HAVEN'T COMPILE THE CODE YET IN THIS DIRECTORY) cmake .
9. make
(9.5 Or you can use : make -j X , where X is an integer between 2 to 24)
10. exit; (To quit the centos container we called using centos.sh)
11. sbatch -p gpu --gres=gpu:1 --time=1:00:00 EpiScale_run.sh; (To submit a job, --time indicates the amount of time to be allocated for this simulation)

For detailed job managing process, see the UCR HPCC manual: https://hpcc.ucr.edu/manuals/hpc_cluster/jobs/

************************************************************************************************
************************************************************************************************
************************************************************************************************

Location of configuration files:
 ./resources

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
6. If you want to run cell proliferation simulations. Make sure to go into SceCells.cu and search for "quiescence1", "quiescence2", and "quiescence3". These three parameters are controlling the cell cycles so they need to be adjusted accordingly.
      
      
      
      
      
      
      
      
************************************************************************************************
************************************************************************************************
************************************************************************************************
************************************************************************************************
************************************************************************************************

LEGACY VERSION OF THE COMPILATION PROCESS
(THIS PART IS OBSOLETE)
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

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

Known issue: When the apical contractility is sufficiently strong (basal ratio: 0.25, apical ratio: 1.25 for instance), the simulation can run into seg fault after a significant duration of simulation.

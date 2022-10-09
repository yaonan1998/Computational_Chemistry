# Running & Analyzing `M`olecular `D`ynamics Simulations
## 0. NOTE
For now, the following procedures are all LAMMPS-based. (https://www.lammps.org/)
## 1. Submit Jobs
* **Preliminary files**
1. For running MD by LAMMPS, the **topology (\*.xyz)** and **force field parameter (\*.lmp)** files of 
the simulated species are required.

2. As the **input script (in.lammps)** of LAMMPS is barely changed (but sometimes it will be tailored for certain purposes), 
this file is also provided in advance.  

&emsp;You can find these necessary files in `Computational_Chemistry/MD/run/pre_files/`.

&emsp;If required files for new species are missing, you have to add them by yourself. :)
* **Run**
1. Please state the system information and simulation settings in `lammps_run.init`:

&emsp;(NOTE: The following JSON file is ONLY for instructions, so it contains annotations starting with #.)
```text
{
  # The directory where you put the preliminary files.
  "PublicFilesDir": "/gpfs01/zhq_work/yn/PUBLIC_FILES/LAMMPS/",
  
  # A list of dictionary(ies) stating the compositions and corresponding molar rations of species in a system.
  "mol_ratio": [
    {
      "DME": 1.4,
      "LiFSI": 1
    },
	{
      "DME": 3.8,
      "LiFSI": 1
    }
  ],
  
  # The expand factor multiplies the molar ratio equals the number of species in a system.
  # The number of elements in `expand_factor` list should equal that in `mol_ratio` list.
  "expand_factor": [
    200,
    200
  ],
  
  # The number of parallel runs for one system with different initial configurations.
  # If only 1 run is required, you can index 1 or false.
  "parallel_run": 2,
  
  # The seed keyword in Packmol.
  # If set to false, default seed=1001 will be used.
  # For parallel runs, the setting here will be ignored, and random seed between 1-1000 will be used.
  "seed": false,
  
  # The simulation and annealing temperature (unit=K).
  "temperature": [
    298.15,
    450.15
  ],
  
  # The simulation time for each fix command in in.lammps file (unit=fs).
  "runtime":[
    2000000,    # npt
    500000,     # npt, heating
    500000,     # npt
    500000,     # npt, cooling
    4000000,    # npt
    2000000,    # nvt
    5000000     # nvt, sampling
  ],
  
  # The compute command for specific output infomation.
  # Options: false (N/A) / "dm" (dipole moment) / "press" (pressure)/ "msd" (msd) or combinations like "dm msd"
  "compute_type": false
}
```
2. Now you can submit MD jobs to the queue in supercomputer(e5/wx) by executing the following command:
```shell
python lammps_run.py --sc e5 --path `pwd` --intemp in.lammps --xyzlmp xyzlmp
```
&emsp;For detailed explanations to the above 4 arguments, please execute `python lammps_run.py -h`.
## 2. Analyze Results
### Solvation structure
#### Radial distribution function (RDF) & Coordination number (CN)
#### Cluster statistics
* Single-cation-centered structures
* Cation–Anion aggregates
### Mean square displacement (MSD)
### Visualization


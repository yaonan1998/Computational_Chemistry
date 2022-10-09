import numpy as np
import math
import scipy.constants as const
import os
import json
import argparse
import random
import logging
from pymatgen.io.lammps.data import LammpsData, CombinedData
from typing import Dict, List, Tuple, Union
import string
from multiprocessing import Pool


def load_json_file(path: str) -> Union[Dict, List]:
    '''
    Read a text file and parse as JSON.

    Args:
    - path: Path to the Json file.

    Returns:
    - data: A Python dict of list of parsed data.
    '''
    with open(path, 'r') as fin:
        return json.load(fin)


def write_json_file(path: str, data: Union[Dict, List]):
    '''
    Write Jsoniazable data to the given path.

    Args:
    - path: Path to the Json file.
    - data: A dict of data.
    '''
    with open(path, 'w') as fout:
        json.dump(data, fout, indent=2)


def density() -> Dict:
    '''
    Get a dict of solvents and corresponding densities.

    Returns:
    - (A dict of solvents and corresponding densities, unit: g cm^-3)
    '''
    return {
        'DME': 0.87,
        'DOL': 1.06,
        'DMP': 0.868,
        'EC': 1.321,
        'DMC': 1.063,
        'FEC': 1.454,
        'HFE':1.53,
        'EA': 0.902,
        'DIPS': 0.81
    }


def molecular_weight() -> Dict:
    '''
    Get a dict of solvents and corresponding molecular weights.

    Returns:
    - (A dict of solvents and corresponding molecular weights, unit: g mol^-1)
    '''
    return {
        'DME': 90.12,
        'DOL': 74.08,
        'DMP': 104.148,
        'EC': 88.06,
        'DMC': 90.07,
        'FEC': 106.052,
        'HFE':232.07,
        'EA': 88.105,
        'DIPS': 118.24
    }


def MDmodel_info(mol_ratio: Dict, expand_factor: int) -> Dict:
    '''
    Estimate the theoretical mass/density/volume of the electrolyte and the length of the MD cell.

    Args:
    - sol_ratio: A dict of solvents/salts and corresponding molar ratios.
    - expand_factor: The expand_factor multiplies the molar ratio equals the number of solvents/salts.

    Returns:
    - Dict,
    {
        'mass (unit: g)': theo_mass,
        'density (unit: g cm^-3)': theo_density,
        'volume (unit: cm^3)': theo_volume,
        'concentration (unit: M)': theo_conc,
        'total_molecule_number': total_molecule_number,
        'composition': composition,
        'box_size (unit: Angstrom)': box_size
    }
    '''
    mol_weight = molecular_weight()
    mass = 0
    elements = []
    weights = []
    dens = density()
    composition = {}
    for mol, ratio in mol_ratio.items():
        mol_num = int(math.ceil(expand_factor * ratio))
        composition.update({mol: mol_num})
        if mol not in ['LiPF6', 'LiFSI', 'LiTFSI', 'LiNO3', 'Li2S4', 'Li2S6', 'Li2S8']:
            mass += mol_weight[mol] * mol_num
            elements.append(dens[mol])
            weights.append(ratio)
        else:
            salt_num = mol_num
            # suppose salts are DME solvents for estimating the system information
            mol = 'DME'
            mass += mol_weight[mol] * mol_num
            elements.append(dens[mol])
            weights.append(ratio)
    theo_mass = mass / const.Avogadro
    theo_density = np.average(elements, weights=weights)
    theo_volume = theo_mass / theo_density
    theo_conc = salt_num / const.Avogadro / (theo_volume * 10 ** (-3))
    box_size = (theo_volume * 10 ** 24) ** (1/3)
    # print(f"{theo_density:.4f}", f"{theo_mass:.4e}")
    return {
        'mass (unit: g)': theo_mass,
        'density (unit: g cm^-3)': theo_density,
        'volume (unit: cm^3)': theo_volume,
        'concentration (unit: M)': theo_conc,
        'total_molecule_number': int(sum(composition.values())),
        'composition': composition,
        'box_size (unit: Angstrom)': box_size
    }


def def_model(composition) -> Tuple[str, List[str], List[str]]:
    '''
    Define the model (i.e., detailed compositions) for the MD job.

    Args:
    - composition: A dict of molecules and corresponding numbers.

    Returns:
    - model: Combined name of the MD model.
    - mol_names: List of molecule names.
    - mol_nums: List of number (str format) of molecules.
    '''
    model = ''
    for mol, num in composition.items():
        if model != '':
            model = model + '+' + str(num) + mol
        else:
            model = str(num) + mol
    mol_names = list(composition.keys())
    mol_nums = [str(i) for i in list(composition.values())]
    return model, mol_names, mol_nums


def mkfile_modelxyz(dir_run: str, model: str, mol_names: List[str], mol_nums: List[str],
                    box_size: float, seed: int = None) -> str:
    '''
    Make the .xyz file of a model containing multiple molecules for the MD job by Packmol program.

    Args:
    - dir_run: Path of the folder to run LAMMPS.
    - model: Model name.
    - mol_names: List of molecule names.
    - mol_nums: Numbers of each used molecule.
    - box_size: Length of the cubic simulation box.
    - seed: Seed keyword in packmol.

    Returns:
    - xyz_path: Path of the `.xyz` file containing packed molecules.
    '''
    pkmin = os.path.join(dir_run, 'packmol.inp')
    pkmout = os.path.join(dir_run, 'packmol.oup')
    modelxyz = os.path.join(dir_run, model + '.xyz')
    if not seed:
        seed = 1001
    with open(pkmin, 'w') as fout:
        fout.write('seed %s\n' % str(seed))
        fout.write('tolerance 2.0\n')
        fout.write('filetype xyz\n')
        fout.write('output ' + modelxyz + '\n')
        for i in range(len(mol_names)):
            fout.write('structure ' +
                       os.path.join(init['PublicFilesDir'], args.xyzlmp, f'{mol_names[i]}.xyz') + '\n')
            fout.write('  number ' + mol_nums[i] + '\n')
            fout.write('  inside box 0.0 0.0 0.0 {l} {l} {l}\n'.format(l=box_size))
            fout.write('end structure\n')
    sh_command = 'packmol <%s >%s && echo packmol done' % (pkmin, pkmout)
    os.system(sh_command)
    with open(pkmout, 'r') as f1:
        info = f1.read()
    assert 'Success!' in info, 'Packmol reported failure!'
    return os.path.join(dir_run, modelxyz)


def mkfile_data(dir_run: str, model: str, modelxyz_path: str,
                mol_names: List[str], mol_nums: List[str]) -> str:
    '''
    Make a .data file for the MD job.

    Args:
    - dir_run: Path of the folder to run LAMMPS.
    - model: Model name.
    - modelxyz_path: Path of the `.xyz` file containing packed molecules.
    - mol_names: List of molecule names.
    - mol_nums: Numbers of each used molecule.

    Returns:
    - data_path: Path of the LAMMPS data file.
    '''
    os.chdir(dir_run)
    mols = []
    list_of_numbers = []
    for path, num in zip(
            [os.path.join(init['PublicFilesDir'], args.xyzlmp, i+'.lmp') for i in mol_names], mol_nums):
        mols.append(LammpsData.from_file(path))
        list_of_numbers.append(int(num))
    coordinates = CombinedData.parse_xyz(modelxyz_path)
    combined = CombinedData.from_lammpsdata(mols, mol_names, list_of_numbers, coordinates)
    data_path = os.path.join(dir_run, model + '.data')
    combined.write_file(data_path)
    return data_path


def lammps_fix_commands(compute_type: str):
    '''
    Get related `compute` & `fix` commands of the given type.

    Args:
    - compute_type: The compute keywords of the MD job, na:None/dm:dipole moment/press:pressure/msd:msd/com:com.

    Returns:
    - commands: Lammps commands.
    '''
    compute_command_info = {
        "na": [
            "\n"
        ],
        "dm": [
            "compute dpall1 all chunk/atom molecule\n",
            "compute dpall2 all dipole/chunk dpall1\n",
            "fix dpall all ave/time 1000 1 1000 c_dpall2[*] file $$ModelName_dipole.out mode vector\n"
        ],
        "press": [
            "variable pxy equal pxy\n",
            "variable pxz equal pxz\n",
            "variable pyz equal pyz\n",
            "fix pressure all ave/time 1 1 1 v_pxy v_pxz v_pyz file $$ModelName_pressure.out\n"
        ],
        "msd": [
            "group cation type $$CationType\n",
            "compute cation1 cation chunk/atom molecule\n",
            "compute cation2 cation msd/chunk cation1\n",
            "group anion type $$AnionType\n",
            "compute anion1 anion chunk/atom molecule\n",
            "compute anion2 anion msd/chunk anion1\n",
            "fix msdcation cation ave/time 1 1 1000 c_cation2[*] file Li_msd.out mode vector\n",
            "fix msdanion anion ave/time 1 1 1000 c_anion2[*] file $$AnionName_msd.out mode vector\n"
        ],
        "com": [
            "TOBEDETERMINED"
        ]
    }
    compute_command_list = []
    for c in compute_type.split():
        if c in compute_command_info:
            ccjoin = ''.join(compute_command_info[c])
            compute_command_list.append(ccjoin)
        else:
            logging.warning('Compute type %s NOT defined in compute_commands file!' % c)
            compute_command_list.append('\n')
    return '\n'.join(compute_command_list)


class LammpsTemplate(string.Template):
    """A string class for supporting $$-PascalStyleToken."""
    delimiter = '$$'
    idpattern = '[a-zA-Z][a-zA-Z0-9]*'


def modfile_inlammps(dir_run: str, in_template: str, model: str, compute_type: str) -> None:
    '''
    Modify `in.lammps` files for the MD job.

    Args:
    - dir_run: Path of the folder to run LAMMPS.
    - in_template: Name of the in.lammps template file.
    - model: Model name.
    - compute_type: The compute keywords of the MD job, 0:None/1:dipole/2:pressure/3:msd.
    '''
    compute_commands = lammps_fix_commands(compute_type)
    temperature1, temperature2, runtime = init['temperature'][0], init['temperature'][1], init['runtime']
    to_replace = {
        'ModelName': model,
        'Temperature1': str(temperature1),
        'Temperature2': str(temperature2)
    }
    to_replace['ComputeType'] = LammpsTemplate(compute_commands).substitute(**to_replace)
    for idx, item in enumerate(runtime):
        to_replace['Runtime' + str(idx)] = str(item)
    intemp = os.path.join(init['PublicFilesDir'], 'intemplate', in_template)
    with open(intemp, 'r') as fin:
        raw_content = fin.read()
    filled_content = LammpsTemplate(raw_content).substitute(**to_replace)
    in_ = os.path.join(dir_run, 'in.lammps')
    with open(in_, 'w') as fout:
        fout.write(filled_content)


def mkfile_subscript(dir_run: str) -> None:
    '''
    Make the submission script file.

    Args:
    - dir_run: Path of the folder to run LAMMPS.
    '''
    dir_sub = os.path.join(dir_run, 'lammps.sh')
    if args.sc == 'wx':
        with open(dir_sub, 'w+') as f2:
            f2.write('bsub -q q_x86_share  -N 1 -np 24 -o lammps.log lmp_mpi -in in.lammps\n')
        os.chdir(dir_run)
        os.system('chmod +x lammps.sh')
        # os.system('./lammps.sh')
    elif args.sc == 'e5':
        with open(dir_sub, 'w+') as f2:
            f2.write('#!/bin/bash\n'
                     '#SBATCH -J lmp_yn\n'
                     '#SBATCH -p cnall\n'
                     '#SBATCH -N 1\n'
                     '#SBATCH -o lammps.log\n'
                     '#SBATCH -e stderr.%j\n'
                     '#SBATCH --no-requeue\n'
                     '#SBATCH --ntasks-per-node=28\n'
                     'module load compiles/intel/2019/u4/config\n'
                     'mpiexec.hydra /home/zhq/lammps-29Sep2021/src/lmp_mpi -in in.lammps\n')
        os.chdir(dir_run)
        # os.system('chmod +x lammps.sh && sbatch lammps.sh')


def submit_job(index: int) -> None:
    '''
    Submit MD jobs to the supercomputer.

    Args:
    - index: Index of MD jobs to be submitted.
    '''
    mol_ratios = init['mol_ratio']
    expand_factors = init['expand_factor']
    para_run = init['parallel_run']
    seed = init['seed']
    compute_type = init['compute_type']
    if not compute_type:
        compute_type = 'na'
    model_info = MDmodel_info(mol_ratios[index], expand_factors[index])
    for i, j in model_info.items():
        print(f'{i:30}{j}')
    model, mol_names, mol_nums = def_model(model_info['composition'])
    model_info['model'] = model
    dir_run = os.path.join(args.path, '+'.join(mol_names) + '_' + '+'.join(mol_nums))
    if para_run and para_run != 1:
        dirs_run = [dir_run + f'_{para_index + 1}' for para_index in range(para_run)]
    else:
        dirs_run = [dir_run]
    for dir_run in dirs_run:
        if not os.path.exists(dir_run):
            os.mkdir(dir_run)
        if len(dirs_run) > 1:
            seed = random.randint(1, 1000)  # Ignore the 'seed' setting in .init file to ensure different seeds
        modelxyz_path = mkfile_modelxyz(dir_run, model, mol_names, mol_nums,
                                        float(np.around(model_info['box_size (unit: Angstrom)'], 2)),
                                        seed)
        mkfile_data(dir_run, model, modelxyz_path, mol_names, mol_nums)
        modfile_inlammps(dir_run, args.intemp, model, compute_type)
        mkfile_subscript(dir_run)
        write_json_file(os.path.join(dir_run, 'model'), model_info)


def main():
    global init
    init = load_json_file('lammps_run.init')
    pool = Pool(len(init['mol_ratio']))
    pool.map(submit_job, range(len(init['mol_ratio'])))
    pool.close()
    pool.join()
    # 1.4: 4.02 M, 3.8: 2.01 M, 8.2: 1.05 M


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sc', help='The supercomputer, wx or e5', type=str, default=None)
    parser.add_argument('--path', help='The path of the folder to run LAMMPS', type=str, default=None)
    parser.add_argument('--intemp', help='The name of the in.lammps template file', type=str, default=None)
    parser.add_argument('--xyzlmp', help='The name of the folder containing .xyz and .lmp files', type=str, default=None)
    args = parser.parse_args()
    main()


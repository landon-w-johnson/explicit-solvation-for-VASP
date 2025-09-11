import numpy as np
import re
import sys
import math
import random



atomic_names = ['H','He','Li','Be','B','C','N','O','F','Ne']
atomic_masses = np.array([1.008,4.003,6.941,9.012,10.811,12.011,14.007,15.999,18.998,20.180], dtype=np.double) # from https://sciencenotes.org/periodic-table-with-atomic-mass/



min_header_skip = 7
opt_header_skip = 0

poscarfile = open('POSCAR', 'r')
poscar = poscarfile.readlines()
poscarfile.close()

scaling_factor = re.split(' +', poscar[1].lstrip().rstrip())
sc_fac = np.ones(3, dtype=np.double)
if len(scaling_factor) == 1:
    if float(scaling_factor[0]) < 0:
        print(f'detected negative scaling factor of {scaling_factor[0]}')
        print('the calculated scaling factor will be printed after the lattice vectors')
        print()
    else:
        for i in range(3):
            sc_fac[i] = np.double(scaling_factor[0])
        print('scaling factor:')
        print(sc_fac)
        print()
elif len(scaling_factor) == 3:
    for i in range(3):
        sc_fac[i] = np.double(scaling_factor[i])
    print('scaling factor:')
    print(sc_fac)
    print()
else:
    print('ERROR! Scaling factor did not have either 1 or 3 elements!')
    print(f'Found {len(scaling_factor)} elements in the following line!')
    print(poscar[1])
    sys.exit()

lat_vec = np.zeros((3,3), dtype=np.double)
for i in range(3):
    vec_list = re.split(' +', poscar[i+2].lstrip().rstrip())
    lat_vec[i,0] = vec_list[0]
    lat_vec[i,1] = vec_list[1]
    lat_vec[i,2] = vec_list[2]
print('lattice vectors:')
print(lat_vec)
print()
lat_lens = np.zeros(3, dtype=np.double)
for i in range(3):
    lat_lens[i] = np.sqrt(lat_vec[i,0]**2 + lat_vec[i,1]**2 + lat_vec[i,2]**2)
if float(scaling_factor[0]) < 0:
    raw_vol = lat_lens[0]*lat_lens[1]*lat_lens[2]
    actual_sc_fac = -float(scaling_factor[0])/raw_vol
    for i in range(3):
        sc_fac[i] = actual_sc_fac
    print('scaling factor:')
    print(actual_sc_fac)
    print()

sc_lat_vec = np.array(lat_vec)
sc_lat_vec_len = np.zeros(3, dtype=np.double)
for i in range(3):
    sc_lat_vec[i,0] *= sc_fac[i]
    sc_lat_vec[i,1] *= sc_fac[i]
    sc_lat_vec[i,2] *= sc_fac[i]
    sc_lat_vec_len[i] = math.sqrt(sc_lat_vec[i,0]**2+sc_lat_vec[i,1]**2+sc_lat_vec[i,2]**2)
print('scaled lattice vectors:')
print(sc_lat_vec)
print()
print('scaled lattice vector lengths:')
print(sc_lat_vec_len)
print()

from_poscar_mat = np.zeros((3,3), dtype=np.double)
from_poscar_mat[:,0] = sc_lat_vec[0,:]
from_poscar_mat[:,1] = sc_lat_vec[1,:]
from_poscar_mat[:,2] = sc_lat_vec[2,:]
print('from_poscar_mat:')
print(from_poscar_mat)
print()
det_from_poscar_mat = from_poscar_mat[0,0]*(from_poscar_mat[1,1]*from_poscar_mat[2,2]-from_poscar_mat[1,2]*from_poscar_mat[2,1])\
    + from_poscar_mat[0,1]*(from_poscar_mat[1,2]*from_poscar_mat[2,0]-from_poscar_mat[1,0]*from_poscar_mat[2,2])\
    + from_poscar_mat[0,2]*(from_poscar_mat[1,0]*from_poscar_mat[2,1]-from_poscar_mat[1,1]*from_poscar_mat[2,0])
print('determinant of from_poscar_mat:')
print(det_from_poscar_mat)
print()
to_poscar_mat = np.zeros((3,3), dtype=np.double)
to_poscar_mat[0,:] = np.cross(from_poscar_mat[:,1],from_poscar_mat[:,2])
to_poscar_mat[1,:] = np.cross(from_poscar_mat[:,2],from_poscar_mat[:,0])
to_poscar_mat[2,:] = np.cross(from_poscar_mat[:,0],from_poscar_mat[:,1])
print('to_poscar_mat before dividing by determinant:')
print(to_poscar_mat)
print()
to_poscar_mat = np.divide(to_poscar_mat,det_from_poscar_mat)
print('to_poscar_mat after dividing by determinant:')
print(to_poscar_mat)
print()

current_line = 5
species_names = []
ions_per_species = re.split(' +', poscar[current_line].lstrip().rstrip())
species_line_ind = current_line
num_ions_line_ind = current_line
if re.search('[0-9]', ions_per_species[0][0]) == None:
    species_names = ions_per_species
    print('species names:')
    print(species_names)
    print()
    current_line += 1
    ions_per_species = re.split(' +', poscar[current_line].lstrip().rstrip())
    num_ions_line_ind = current_line
    for i in range(len(ions_per_species)):
        ions_per_species[i] = int(ions_per_species[i])
else:
    for i in range(len(ions_per_species)):
        ions_per_species[i] = int(ions_per_species[i])
print('ions per species:')
print(ions_per_species)
print()
num_ions = 0
for i in range(len(ions_per_species)):
    num_ions += ions_per_species[i]
print('total number of ions:')
print(num_ions)
print()

current_line += 1
sel_dyn = False
coord_scheme = poscar[current_line].lstrip().rstrip()
if coord_scheme[0].lower() == 's':
    sel_dyn = True
    current_line += 1
    coord_scheme = poscar[current_line].lstrip().rstrip()
if coord_scheme[0].lower() == 'c' or coord_scheme[0].lower() == 'k':
    coord_scheme = 'Cartesian'
else:
    coord_scheme = 'Direct'
print('selective dynamics?')
print(sel_dyn)
print()
print('coordinate scheme:')
print(coord_scheme)
print()
header_skip = current_line + 1

ion_positions = np.zeros((num_ions,3), dtype=np.double)
i = 0
for line_ind in range(header_skip,header_skip+num_ions):
    pos = re.split(' +', poscar[line_ind].lstrip().rstrip())
    ion_positions[i,0] = pos[0]
    ion_positions[i,1] = pos[1]
    ion_positions[i,2] = pos[2]
    i += 1
print('ion positions:')
print(ion_positions)
print()

car_ion_pos = np.zeros((num_ions,3), dtype=np.double)
if coord_scheme == 'Direct':
    for i in range(num_ions):
        car_ion_pos[i,:] = np.matmul(from_poscar_mat,ion_positions[i,:])
else:
    for i in range(num_ions):
        car_ion_pos[i,:] = np.multiply(ion_positions[i,:],sc_fac)
print('ion positions in standard Cartesian reference frame:')
print(car_ion_pos)
print()

tot_vol = np.double(lat_lens[0]*sc_fac[0]*lat_lens[1]*sc_fac[1]*lat_lens[2]*sc_fac[2])
print('tot_vol:')
print(tot_vol)
print()





solvent_file = open('SOLVENT', 'r')
solvent = solvent_file.readlines()
solvent_file.close()

sol_scaling_factor = re.split(' +', solvent[1].lstrip().rstrip())
sol_sc_fac = np.ones(3, dtype=np.double)
if len(sol_scaling_factor) == 1:
    if float(sol_scaling_factor[0]) < 0:
        print(f'detected negative scaling factor in solvent of {scaling_factor[0]}')
        print('the calculated solvent scaling factor will be printed after the lattice vectors')
        print()
    else:
        for i in range(3):
            sol_sc_fac[i] = np.double(sol_scaling_factor[0])
        print('solvent scaling factor:')
        print(sol_sc_fac)
        print()
elif len(sol_scaling_factor) == 3:
    for i in range(3):
        sol_sc_fac[i] = np.double(sol_scaling_factor[i])
    print('solvent scaling factor:')
    print(sol_sc_fac)
    print()
else:
    print('ERROR! Solvent scaling factor did not have either 1 or 3 elements!')
    print(f'Found {len(scaling_factor)} elements in the following line!')
    print(solvent[1])
    sys.exit()

sol_lat_vec = np.zeros((3,3), dtype=np.double)
for i in range(3):
    sol_vec_list = re.split(' +', solvent[i+2].lstrip().rstrip())
    sol_lat_vec[i,0] = sol_vec_list[0]
    sol_lat_vec[i,1] = sol_vec_list[1]
    sol_lat_vec[i,2] = sol_vec_list[2]
print('solvent lattice vectors:')
print(sol_lat_vec)
print()

sol_lat_lens = np.zeros(3, dtype=np.double)
for i in range(3):
    sol_lat_lens[i] = np.sqrt(sol_lat_vec[i,0]**2 + sol_lat_vec[i,1]**2 + sol_lat_vec[i,2]**2)
if float(sol_scaling_factor[0]) < 0:
    sol_raw_vol = sol_lat_lens[0]*sol_lat_lens[1]*sol_lat_lens[2]
    sol_actual_sc_fac = -float(sol_scaling_factor[0])/sol_raw_vol
    for i in range(3):
        sol_sc_fac[i] = sol_actual_sc_fac
    print('solvent scaling factor:')
    print(sol_actual_sc_fac)
    print()

sol_sc_lat_vec = np.array(sol_lat_vec)
sol_sc_lat_vec_len = np.zeros(3, dtype=np.double)
for i in range(3):
    sol_sc_lat_vec[i,0] *= sol_sc_fac[i]
    sol_sc_lat_vec[i,1] *= sol_sc_fac[i]
    sol_sc_lat_vec[i,2] *= sol_sc_fac[i]
    sol_sc_lat_vec_len[i] = math.sqrt(sol_sc_lat_vec[i,0]**2+sol_sc_lat_vec[i,1]**2+sol_sc_lat_vec[i,2]**2)
print('solvent scaled lattice vectors:')
print(sol_sc_lat_vec)
print()
print('solvent scaled lattice vector lengths:')
print(sol_sc_lat_vec_len)
print()

from_solvent_mat = np.zeros((3,3), dtype=np.double)
from_solvent_mat[:,0] = sol_sc_lat_vec[0,:]
from_solvent_mat[:,1] = sol_sc_lat_vec[1,:]
from_solvent_mat[:,2] = sol_sc_lat_vec[2,:]
print('from_solvent_mat:')
print(from_solvent_mat)
print()

sol_current_line = 5
sol_species_names = []
sol_ions_per_species = re.split(' +', solvent[sol_current_line].lstrip().rstrip())
if re.search('[0-9]', sol_ions_per_species[0][0]) == None:
    sol_current_line += 1
    sol_species_names = sol_ions_per_species
    print('solvent species names:')
    print(sol_species_names)
    print()
    sol_ions_per_species = re.split(' +', solvent[sol_current_line].lstrip().rstrip())
    for i in range(len(sol_ions_per_species)):
        sol_ions_per_species[i] = int(sol_ions_per_species[i])
else:
    print('ERROR! Solvent file does not include species names!')
    print('')
    sys.exit()
print('solvent ions per species:')
print(sol_ions_per_species)
print()
sol_num_ions = 0
for i in range(len(sol_ions_per_species)):
    sol_num_ions += sol_ions_per_species[i]
print('total number of solvent ions:')
print(sol_num_ions)
print()

sol_current_line += 1
sol_sel_dyn = False
sol_coord_scheme = solvent[sol_current_line].lstrip().rstrip()
if sol_coord_scheme[0].lower() == 's':
    sol_sel_dyn = True
    sol_current_line += 1
    sol_coord_scheme = solvent[sol_current_line].lstrip().rstrip()
if sol_coord_scheme[0].lower() == 'c' or sol_coord_scheme[0].lower() == 'k':
    sol_coord_scheme = 'Cartesian'
else:
    sol_coord_scheme = 'Direct'
print('selective dynamics for solvent?')
print(sol_sel_dyn)
print()
print('solvent coordinate scheme:')
print(sol_coord_scheme)
print()
sol_header_skip = current_line + 1

sol_ion_positions = np.zeros((sol_num_ions,3), dtype=np.double)
i = 0
for sol_line_ind in range(sol_header_skip,sol_header_skip+sol_num_ions):
    sol_pos = re.split(' +', solvent[sol_line_ind].lstrip().rstrip())
    sol_ion_positions[i,0] = sol_pos[0]
    sol_ion_positions[i,1] = sol_pos[1]
    sol_ion_positions[i,2] = sol_pos[2]
    i += 1
print('solvent ion positions:')
print(sol_ion_positions)
print()

sol_sc_ion_pos = np.array(sol_ion_positions, dtype=np.double)
sol_center_pos = np.zeros(3, dtype=np.double)
if sol_coord_scheme == 'Direct':
    for i in range(sol_ion_positions.shape[0]):
        sol_sc_ion_pos[i,:] = np.matmul(from_solvent_mat,sol_ion_positions[i,:])
        sol_center_pos = np.add(sol_center_pos,sol_sc_ion_pos[i,:])
else: # Cartesian
    print('WARNING! The Cartesian scaling code block has not been tested or debugged!')
    for i in range(sol_ion_positions.shape[0]):
        sol_sc_ion_pos[i,:] = np.multiply(sol_ion_positions[i,:],sol_sc_fac)
        sol_center_pos = np.add(sol_center_pos,sol_sc_ion_pos[i,:])
sol_center_pos[0] /= sol_num_ions
sol_center_pos[1] /= sol_num_ions
sol_center_pos[2] /= sol_num_ions
print('scaled solvent ion positions:')
print(sol_sc_ion_pos)
print()
print('solvent center position:')
print(sol_center_pos)
print()
sol_rel_pos = np.zeros(sol_ion_positions.shape, dtype=np.double)
for i in range(sol_sc_ion_pos.shape[0]):
    sol_rel_pos[i,:] = np.subtract(sol_sc_ion_pos[i,:],sol_center_pos)
print('solvent relative positions:')
print(sol_rel_pos)
print()

sol_mass = np.double(0)
for i in range(len(sol_species_names)):
    sol_mass += sol_ions_per_species[i]*atomic_masses[atomic_names.index(sol_species_names[i])]
print('solvent mass:')
print(f'{sol_mass} atomic mass units')
print()

print('Enter your desired solvent density in g/mL:')
density = np.double(input())
print()

density *= 0.6022136652 # conversion from g/mL to amu/Angstrom^3
volume_per_sol = sol_mass/density
print('volume per solvent molecule:')
print(volume_per_sol)
print()

too_big = True
partitions = 1
part_len = 1
mini_vol = np.double(tot_vol)
while too_big:
    partitions += 1
    mini_vol = tot_vol/(partitions**3)
    if np.less_equal(mini_vol, volume_per_sol):
        too_big = False
        part_len = np.double(1/partitions)
print('partitions:')
print(partitions)
print()
print('partition length:')
print(part_len)
print()
print('mini volume:')
print(mini_vol)
print()

occ_rat = np.double(mini_vol/volume_per_sol)
print('occupation ratio:')
print(occ_rat)
print()

########## REMOVE? ##########
mini_cells = np.zeros((3,partitions), dtype=np.double)
for i in range(partitions):
    scaling = i/partitions
    mini_cells[0,i] = scaling*lat_lens[0]
    mini_cells[1,i] = scaling*lat_lens[1]
    mini_cells[2,i] = scaling*lat_lens[2]
print('mini_cells:')
print(mini_cells)
print()
#############################

use_cell = np.ones((partitions,partitions,partitions), dtype=bool)
if coord_scheme == 'Direct':
    for ion in range(ion_positions.shape[0]):
        x_cell = partitions-1
        y_cell = partitions-1
        z_cell = partitions-1
        for x in range(1,partitions):
            if ion_positions[ion,0] < x/partitions:
                x_cell = x - 1
                break
        for y in range(1,partitions):
            if ion_positions[ion,1] < y/partitions:
                y_cell = y - 1
                break
        for z in range(1,partitions):
            if ion_positions[ion,2] < z/partitions:
                z_cell = z - 1
                break
        use_cell[x_cell,y_cell,z_cell] = False
else: # Cartesian
    print('Cartesian coordinate scheme not yet implemented')
    print('Please convert your POSCAR to direct coordinate scheme')
    print('Exiting program')
    sys.exit()
print('cells used:')
print(use_cell)
print()

num_usable_mini_cells = 0
avail_cells = []
for x in range(partitions):
    for y in range(partitions):
        for z in range(partitions):
            if use_cell[x,y,z]:
                num_usable_mini_cells += 1
                avail_cells.append((x,y,z))
print('number of usable mini cells before shuffle and pop:')
print(num_usable_mini_cells)
print()
print('available cells before shuffle and pop:')
print(avail_cells)
print()

num_unocc = math.floor((1-occ_rat)*num_usable_mini_cells)
print('number of unoccupied mini cells')
print(num_unocc)
print()

random.shuffle(avail_cells)
for i in range(num_unocc):
    avail_cells.pop()
num_used_cells = len(avail_cells)
print('number of used cells after shuffle and pop:')
print(num_used_cells)
print()
print('used cells after shuffle and pop:')
print(avail_cells)
print()

sol_ion_identifier = []
for i,name in enumerate(sol_species_names):
    for j in range(sol_ions_per_species[i]):
        sol_ion_identifier.append(name)
sep_sol_pos = np.zeros((num_used_cells*sol_num_ions, 3), dtype=np.double)
rot_mat = np.zeros((3,3), dtype=np.double)
ind = 0
for cell_ind,cell in enumerate(avail_cells):
    x = random.random()
    y = random.random()
    z = random.random()
    trans_array = np.array([(cell[0]+x)*part_len, (cell[1]+y)*part_len, (cell[2]+z)*part_len], dtype=np.double)
    print('cell index:')
    print(cell_ind)
    print()
    print('cell:')
    print(cell)
    print()
    print('translation array:')
    print(trans_array)
    print()
    alpha = random.random()*2*math.pi
    beta = random.random()*2*math.pi
    gamma = random.random()*2*math.pi
    rot_mat[0,0] = math.cos(alpha)*math.cos(beta)
    rot_mat[0,1] = math.cos(alpha)*math.sin(beta)*math.sin(gamma)\
        - math.sin(alpha)*math.cos(gamma)
    rot_mat[0,2] = math.cos(alpha)*math.sin(beta)*math.cos(gamma)\
        + math.sin(alpha)*math.sin(gamma)
    rot_mat[1,0] = math.sin(alpha)*math.cos(beta)
    rot_mat[1,1] = math.sin(alpha)*math.sin(beta)*math.sin(gamma)\
        + math.cos(alpha)*math.cos(gamma)
    rot_mat[1,2] = math.sin(alpha)*math.sin(beta)*math.cos(gamma)\
        - math.cos(alpha)*math.sin(gamma)
    rot_mat[2,0] = -math.sin(beta)
    rot_mat[2,1] = math.cos(beta)*math.sin(gamma)
    rot_mat[2,2] = math.cos(beta)*math.cos(gamma)
    print('rotation matrix:')
    print(rot_mat)
    print()
    for sol_ion in range(sol_num_ions):
        rot_pos = np.matmul(rot_mat,sol_rel_pos[sol_ion,:])
        print(f'rotated ion position for ion index {sol_ion}:')
        print(rot_pos)
        print()
        new_basis_pos = np.matmul(to_poscar_mat, rot_pos)
        trans_pos = np.add(new_basis_pos, trans_array)
        print(f'transalted ion positions for ion index {sol_ion}:')
        print(trans_pos)
        print()
        sep_sol_pos[ind, :] = trans_pos
        ind += 1
print('species separated solvent ion positions:')
print(sep_sol_pos)
print()

new_sol_pos = np.zeros((num_used_cells*sol_num_ions,3), dtype=np.double)
ind = 0
ion_shift = 0
for ions in sol_ions_per_species:
    for mol in range(num_used_cells):
        for ion in range(ions):
            sep_ind = mol*sol_num_ions + ion_shift + ion
            new_sol_pos[ind,:] = sep_sol_pos[sep_ind,:]
            ind += 1
    ion_shift += ions
print('new solvent ion positions after separating by species:')
print(new_sol_pos)
print()




new_poscar = open('POSCAR_WITH_SOLVENT', 'w')
for i in range(species_line_ind):
    new_poscar.write(poscar[i])
if species_names:
    species_names_output = re.sub('\n', '', poscar[species_line_ind])
    for species in sol_species_names:
        species_names_output += f' {species} '
    species_names_output += '\n'
    print('species_names_output:')
    print(species_names_output)
    print()
    new_poscar.write(species_names_output)
num_ions_output = re.sub('\n', '', poscar[num_ions_line_ind])
for ions in sol_ions_per_species:
    num_ions_output += f' {num_used_cells*ions} '
num_ions_output += '\n'
print('num_ions_output:')
print(num_ions_output)
print()
new_poscar.write(num_ions_output)
for line_ind in range(num_ions_line_ind+1,len(poscar)):
    if re.search('[a-zA-Z0-9.]+', poscar[line_ind]) == None:
        break
    new_poscar.write(poscar[line_ind])
if sel_dyn:
    for sol_ion in range(len(new_sol_pos)):
        x = new_sol_pos[sol_ion,0]
        y = new_sol_pos[sol_ion,1]
        z = new_sol_pos[sol_ion,2]
        new_poscar.write(f'  {x:.16f}  {y:.16f}  {z:.17f}   T   T   T\n')
else:
    for sol_ion in range(len(new_sol_pos)):
        x = new_sol_pos[sol_ion,0]
        y = new_sol_pos[sol_ion,1]
        z = new_sol_pos[sol_ion,2]
        new_poscar.write(f'  {x:.16f}  {y:.16f}  {z:.17f}\n')
new_poscar.close()





print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('!REMEMBER TO UPDATE YOUR POTCAR!')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print()

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
elif len(scaling_factor) == 3:
    for i in range(3):
        sc_fac[i] = np.double(scaling_factor[i])
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
lat_lens = np.zeros(3, dtype=np.double)
for i in range(3):
    lat_lens[i] = np.sqrt(lat_vec[i,0]**2 + lat_vec[i,1]**2 + lat_vec[i,2]**2)
if float(scaling_factor[0]) < 0:
    raw_vol = lat_lens[0]*lat_lens[1]*lat_lens[2]
    actual_sc_fac = -float(scaling_factor[0])/raw_vol
    for i in range(3):
        sc_fac[i] = actual_sc_fac

sc_lat_vec = np.array(lat_vec)
sc_lat_vec_len = np.zeros(3, dtype=np.double)
for i in range(3):
    sc_lat_vec[i,0] *= sc_fac[i]
    sc_lat_vec[i,1] *= sc_fac[i]
    sc_lat_vec[i,2] *= sc_fac[i]
    sc_lat_vec_len[i] = math.sqrt(sc_lat_vec[i,0]**2+sc_lat_vec[i,1]**2+sc_lat_vec[i,2]**2)

from_poscar_mat = np.zeros((3,3), dtype=np.double)
from_poscar_mat[:,0] = sc_lat_vec[0,:]
from_poscar_mat[:,1] = sc_lat_vec[1,:]
from_poscar_mat[:,2] = sc_lat_vec[2,:]
det_from_poscar_mat = from_poscar_mat[0,0]*(from_poscar_mat[1,1]*from_poscar_mat[2,2]-from_poscar_mat[1,2]*from_poscar_mat[2,1])\
    + from_poscar_mat[0,1]*(from_poscar_mat[1,2]*from_poscar_mat[2,0]-from_poscar_mat[1,0]*from_poscar_mat[2,2])\
    + from_poscar_mat[0,2]*(from_poscar_mat[1,0]*from_poscar_mat[2,1]-from_poscar_mat[1,1]*from_poscar_mat[2,0])

to_poscar_mat = np.zeros((3,3), dtype=np.double)
to_poscar_mat[0,:] = np.cross(from_poscar_mat[:,1],from_poscar_mat[:,2])
to_poscar_mat[1,:] = np.cross(from_poscar_mat[:,2],from_poscar_mat[:,0])
to_poscar_mat[2,:] = np.cross(from_poscar_mat[:,0],from_poscar_mat[:,1])
to_poscar_mat = np.divide(to_poscar_mat,det_from_poscar_mat)

current_line = 5
species_names = []
ions_per_species = re.split(' +', poscar[current_line].lstrip().rstrip())
species_line_ind = current_line
num_ions_line_ind = current_line
if re.search('[0-9]', ions_per_species[0][0]) == None:
    species_names = ions_per_species
    current_line += 1
    ions_per_species = re.split(' +', poscar[current_line].lstrip().rstrip())
    num_ions_line_ind = current_line
    for i in range(len(ions_per_species)):
        ions_per_species[i] = int(ions_per_species[i])
else:
    for i in range(len(ions_per_species)):
        ions_per_species[i] = int(ions_per_species[i])
num_ions = 0
for i in range(len(ions_per_species)):
    num_ions += ions_per_species[i]

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

header_skip = current_line + 1
ion_positions = np.zeros((num_ions,3), dtype=np.double)
i = 0
for line_ind in range(header_skip,header_skip+num_ions):
    pos = re.split(' +', poscar[line_ind].lstrip().rstrip())
    ion_positions[i,0] = pos[0]
    ion_positions[i,1] = pos[1]
    ion_positions[i,2] = pos[2]
    i += 1

car_ion_pos = np.zeros((num_ions,3), dtype=np.double)
if coord_scheme == 'Direct':
    for i in range(num_ions):
        car_ion_pos[i,:] = np.matmul(from_poscar_mat,ion_positions[i,:])
else:
    for i in range(num_ions):
        car_ion_pos[i,:] = np.multiply(ion_positions[i,:],sc_fac)

tot_vol = np.double(lat_lens[0]*sc_fac[0]*lat_lens[1]*sc_fac[1]*lat_lens[2]*sc_fac[2])





solvent_file = open('SOLVENT', 'r')
solvent = solvent_file.readlines()
solvent_file.close()

sol_scaling_factor = re.split(' +', solvent[1].lstrip().rstrip())
sol_sc_fac = np.ones(3, dtype=np.double)
if len(sol_scaling_factor) == 1:
    if float(sol_scaling_factor[0]) < 0:
        print('found negative scaling factor')
    else:
        for i in range(3):
            sol_sc_fac[i] = np.double(sol_scaling_factor[0])
elif len(sol_scaling_factor) == 3:
    for i in range(3):
        sol_sc_fac[i] = np.double(sol_scaling_factor[i])
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

sol_lat_lens = np.zeros(3, dtype=np.double)
for i in range(3):
    sol_lat_lens[i] = np.sqrt(sol_lat_vec[i,0]**2 + sol_lat_vec[i,1]**2 + sol_lat_vec[i,2]**2)
if float(sol_scaling_factor[0]) < 0:
    sol_raw_vol = sol_lat_lens[0]*sol_lat_lens[1]*sol_lat_lens[2]
    sol_actual_sc_fac = -float(sol_scaling_factor[0])/sol_raw_vol
    for i in range(3):
        sol_sc_fac[i] = sol_actual_sc_fac

sol_sc_lat_vec = np.array(sol_lat_vec)
sol_sc_lat_vec_len = np.zeros(3, dtype=np.double)
for i in range(3):
    sol_sc_lat_vec[i,0] *= sol_sc_fac[i]
    sol_sc_lat_vec[i,1] *= sol_sc_fac[i]
    sol_sc_lat_vec[i,2] *= sol_sc_fac[i]
    sol_sc_lat_vec_len[i] = math.sqrt(sol_sc_lat_vec[i,0]**2+sol_sc_lat_vec[i,1]**2+sol_sc_lat_vec[i,2]**2)

from_solvent_mat = np.zeros((3,3), dtype=np.double)
from_solvent_mat[:,0] = sol_sc_lat_vec[0,:]
from_solvent_mat[:,1] = sol_sc_lat_vec[1,:]
from_solvent_mat[:,2] = sol_sc_lat_vec[2,:]

sol_current_line = 5
sol_species_names = []
sol_ions_per_species = re.split(' +', solvent[sol_current_line].lstrip().rstrip())
if re.search('[0-9]', sol_ions_per_species[0][0]) == None:
    sol_current_line += 1
    sol_species_names = sol_ions_per_species
    sol_ions_per_species = re.split(' +', solvent[sol_current_line].lstrip().rstrip())
    for i in range(len(sol_ions_per_species)):
        sol_ions_per_species[i] = int(sol_ions_per_species[i])
else:
    print('ERROR! Solvent file does not include species names!')
    print('')
    sys.exit()
sol_num_ions = 0
for i in range(len(sol_ions_per_species)):
    sol_num_ions += sol_ions_per_species[i]

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

sol_header_skip = sol_current_line + 1
sol_ion_positions = np.zeros((sol_num_ions,3), dtype=np.double)
i = 0
for sol_line_ind in range(sol_header_skip,sol_header_skip+sol_num_ions):
    sol_pos = re.split(' +', solvent[sol_line_ind].lstrip().rstrip())
    sol_ion_positions[i,0] = sol_pos[0]
    sol_ion_positions[i,1] = sol_pos[1]
    sol_ion_positions[i,2] = sol_pos[2]
    i += 1

sol_sc_ion_pos = np.array(sol_ion_positions, dtype=np.double)
sol_center_pos = np.zeros(3, dtype=np.double)
if sol_coord_scheme == 'Direct':
    for i in range(sol_ion_positions.shape[0]):
        sol_sc_ion_pos[i,:] = np.matmul(from_solvent_mat,sol_ion_positions[i,:])
        sol_center_pos = np.add(sol_center_pos,sol_sc_ion_pos[i,:])
else: # Cartesian
    for i in range(sol_ion_positions.shape[0]):
        sol_sc_ion_pos[i,:] = np.multiply(sol_ion_positions[i,:],sol_sc_fac)
        sol_center_pos = np.add(sol_center_pos,sol_sc_ion_pos[i,:])
sol_center_pos[0] /= sol_num_ions
sol_center_pos[1] /= sol_num_ions
sol_center_pos[2] /= sol_num_ions

sol_rel_pos = np.zeros(sol_ion_positions.shape, dtype=np.double)
for i in range(sol_sc_ion_pos.shape[0]):
    sol_rel_pos[i,:] = np.subtract(sol_sc_ion_pos[i,:],sol_center_pos)

sol_mass = np.double(0)
for i in range(len(sol_species_names)):
    sol_mass += sol_ions_per_species[i]*atomic_masses[atomic_names.index(sol_species_names[i])]

print('Enter your desired solvent density in g/mL:')
density = np.double(input())
print()

density *= 0.6022136652 # conversion from g/mL to amu/Angstrom^3
volume_per_sol = sol_mass/density

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

occ_rat = np.double(mini_vol/volume_per_sol)

atoms_in_mini_cell = []
for x in range(partitions):
    atoms_in_mini_cell.append([])
    for y in range(partitions):
        atoms_in_mini_cell[x].append([])
        for z in range(partitions):
            atoms_in_mini_cell[x][y].append([])

use_cell = np.ones((partitions,partitions,partitions), dtype=bool)
if coord_scheme == 'Direct':
    for ion in range(num_ions):
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
        atoms_in_mini_cell[x_cell][y_cell][z_cell].append(np.array(car_ion_pos[ion,:], dtype=np.double))
else: # Cartesian
    for ion in range(num_ions):
        direct_pos = np.matmul(to_poscar_mat,car_ion_pos[ion,:])
        x_cell = partitions-1
        y_cell = partitions-1
        z_cell = partitions-1
        for x in range(1,partitions):
            if direct_pos[0] < x/partitions:
                x_cell = x - 1
                break
        for y in range(1,partitions):
            if direct_pos[1] < y/partitions:
                y_cell = y - 1
                break
        for z in range(1,partitions):
            if direct_pos[2] < z/partitions:
                z_cell = z - 1
                break
        use_cell[x_cell,y_cell,z_cell] = False
        atoms_in_mini_cell[x_cell][y_cell][z_cell].append(np.array(car_ion_pos[ion,:], dtype=np.double))

num_usable_mini_cells = 0
avail_cells = []
for x in range(partitions):
    for y in range(partitions):
        for z in range(partitions):
            if use_cell[x,y,z]:
                num_usable_mini_cells += 1
                avail_cells.append((x,y,z))

num_unocc = math.floor((1-occ_rat)*num_usable_mini_cells)

random.shuffle(avail_cells)
for i in range(num_unocc):
    avail_cells.pop()
num_used_cells = len(avail_cells)

sol_ion_identifier = []
for i,name in enumerate(sol_species_names):
    for j in range(sol_ions_per_species[i]):
        sol_ion_identifier.append(name)
sep_sol_pos = np.zeros((num_used_cells*sol_num_ions, 3), dtype=np.double)
rot_mat = np.zeros((3,3), dtype=np.double)
ind = sol_num_ions
if coord_scheme == 'Direct':
    for cell_ind,cell in enumerate(avail_cells):
        x = random.random()
        y = random.random()
        z = random.random()
        trans_array = np.matmul(from_poscar_mat, np.array([(cell[0]+x)*part_len, (cell[1]+y)*part_len, (cell[2]+z)*part_len], dtype=np.double))
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
        too_close = True
        while too_close:
            too_close = False
            shift_vec = np.zeros(3, dtype=np.double)
            tmp_ion_storage = np.zeros((sol_num_ions,3), dtype=np.double)
            ind -= sol_num_ions
            for sol_ion in range(sol_num_ions):
                rot_pos = np.matmul(rot_mat,sol_rel_pos[sol_ion,:])
                trans_pos = np.add(rot_pos, trans_array)
                for x_cell in [(cell[0]-1)%partitions, cell[0], (cell[0]+1)%partitions]:
                    for y_cell in [(cell[1]-1)%partitions, cell[1], (cell[1]+1)%partitions]:
                        for z_cell in [(cell[2]-1)%partitions, cell[2], (cell[2]+1)%partitions]:
                            for already_existing_atom_pos in atoms_in_mini_cell[x_cell][y_cell][z_cell]:
                                pbc_pos = already_existing_atom_pos
                                if x_cell==0 and cell[0]==partitions-1:
                                    pbc_pos = np.add(pbc_pos, sc_lat_vec[0,:])
                                elif x_cell==partitions-1 and cell[0]==0:
                                    pbc_pos = np.subtract(pbc_pos, sc_lat_vec[0,:])
                                if y_cell==0 and cell[1]==partitions-1:
                                    pbc_pos = np.add(pbc_pos, sc_lat_vec[1,:])
                                elif y_cell==partitions-1 and cell[1]==0:
                                    pbc_pos = np.subtract(pbc_pos, sc_lat_vec[1,:])
                                if z_cell==0 and cell[2]==partitions-1:
                                    pbc_pos = np.add(pbc_pos, sc_lat_vec[2,:])
                                elif z_cell==partitions-1 and cell[2]==0:
                                    pbc_pos = np.subtract(pbc_pos, sc_lat_vec[2,:])
                                diff_vec = np.subtract(trans_pos, pbc_pos)
                                diff_len = math.sqrt(diff_vec[0]**2 + diff_vec[1]**2 + diff_vec[2]**2)
                                if diff_len < 2:
                                    too_close = True
                                    shift_vec = np.add(shift_vec, diff_vec)
                tmp_ion_storage[sol_ion,:] = trans_pos
                ind += 1
            trans_array = np.add(trans_array, shift_vec)
        for i in range(sol_num_ions):
            sep_sol_pos[ind-sol_num_ions+i,:] = np.matmul(to_poscar_mat, tmp_ion_storage[i,:])
            atoms_in_mini_cell[cell[0]][cell[1]][cell[2]].append(np.array(tmp_ion_storage[i,:]))
        ind += sol_num_ions
else: # Cartesian
    for cell_ind,cell in enumerate(avail_cells):
        x = random.random()
        y = random.random()
        z = random.random()
        trans_array = np.zeros(3, dtype=np.double)
        #trans_array[0] = (cell[0]+x)*part_len*(sc_lat_vec[0,0]+sc_lat_vec[1,0]+sc_lat_vec[2,0])
        #trans_array[1] = (cell[1]+y)*part_len*(sc_lat_vec[0,1]+sc_lat_vec[1,1]+sc_lat_vec[2,1])
        #trans_array[2] = (cell[2]+z)*part_len*(sc_lat_vec[0,2]+sc_lat_vec[1,2]+sc_lat_vec[2,2])
        trans_array[0] = (cell[0]+x)*part_len*sc_lat_vec[0,0]\
            + (cell[1]+y)*part_len*sc_lat_vec[1,0]\
            + (cell[2]+z)*part_len*sc_lat_vec[2,0]
        trans_array[1] = (cell[0]+x)*part_len*sc_lat_vec[0,1]\
            + (cell[1]+y)*part_len*sc_lat_vec[1,1]\
            + (cell[2]+z)*part_len*sc_lat_vec[2,1]
        trans_array[2] = (cell[0]+x)*part_len*sc_lat_vec[0,2]\
            + (cell[1]+y)*part_len*sc_lat_vec[1,2]\
            + (cell[2]+z)*part_len*sc_lat_vec[2,2]
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
        too_close = True
        while too_close:
            too_close = False
            shift_vec = np.zeros(3, dtype=np.double)
            tmp_ion_storage = np.zeros((sol_num_ions,3), dtype=np.double)
            ind -= sol_num_ions
            for sol_ion in range(sol_num_ions):
                rot_pos = np.matmul(rot_mat,sol_rel_pos[sol_ion,:])
                trans_pos = np.add(rot_pos, trans_array)
                for x_cell in [(cell[0]-1)%partitions, cell[0], (cell[0]+1)%partitions]:
                    for y_cell in [(cell[1]-1)%partitions, cell[1], (cell[1]+1)%partitions]:
                        for z_cell in [(cell[2]-1)%partitions, cell[2], (cell[2]+1)%partitions]:
                            for already_existing_atom_pos in atoms_in_mini_cell[x_cell][y_cell][z_cell]:
                                pbc_pos = already_existing_atom_pos
                                if x_cell==0 and cell[0]==partitions-1:
                                    pbc_pos = np.add(pbc_pos, sc_lat_vec[0,:])
                                elif x_cell==partitions-1 and cell[0]==0:
                                    pbc_pos = np.subtract(pbc_pos, sc_lat_vec[0,:])
                                if y_cell==0 and cell[1]==partitions-1:
                                    pbc_pos = np.add(pbc_pos, sc_lat_vec[1,:])
                                elif y_cell==partitions-1 and cell[1]==0:
                                    pbc_pos = np.subtract(pbc_pos, sc_lat_vec[1,:])
                                if z_cell==0 and cell[2]==partitions-1:
                                    pbc_pos = np.add(pbc_pos, sc_lat_vec[2,:])
                                elif z_cell==partitions-1 and cell[2]==0:
                                    pbc_pos = np.subtract(pbc_pos, sc_lat_vec[2,:])
                                diff_vec = np.subtract(trans_pos, pbc_pos)
                                diff_len = math.sqrt(diff_vec[0]**2 + diff_vec[1]**2 + diff_vec[2]**2)
                                if diff_len < 2:
                                    too_close = True
                                    shift_vec = np.add(shift_vec, diff_vec)
                tmp_ion_storage[sol_ion,:] = trans_pos
                ind += 1
                trans_array = np.add(trans_array, shift_vec)
        for i in range(sol_num_ions):
            sep_sol_pos[ind-sol_num_ions+i,:] = np.divide(tmp_ion_storage[i,:], sc_fac)
            atoms_in_mini_cell[cell[0]][cell[1]][cell[2]].append(np.array(tmp_ion_storage[i,:]))
        ind += sol_num_ions

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




new_poscar = open('POSCAR_WITH_SOLVENT', 'w')
for i in range(species_line_ind):
    new_poscar.write(poscar[i])
if species_names:
    species_names_output = re.sub('\n', '', poscar[species_line_ind])
    for species in sol_species_names:
        species_names_output += f' {species} '
    species_names_output += '\n'
    new_poscar.write(species_names_output)
num_ions_output = re.sub('\n', '', poscar[num_ions_line_ind])
for ions in sol_ions_per_species:
    num_ions_output += f' {num_used_cells*ions} '
num_ions_output += '\n'
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

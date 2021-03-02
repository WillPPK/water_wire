import sys

import MDAnalysis as mda

import MDAnalysis.analysis.hbonds

import numpy as np



u = mda.Universe('3_sort.pdb','md_skip_end.xtc')



def analysis(current, output, u, order):

        # decompose the first hydrogen bond.

    sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle = current[0]

        # decompose the last hydrogen bond.

    atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle = current[-1]

        # expand the atom index to the resname, resid, atom names

    sele1 = u.atoms[sele1_index]

    sele2 = u.atoms[sele2_index]

    (s1_resname, s1_resid, s1_name) = (sele1.resname, sele1.resid, sele1.name)

    (s2_resname, s2_resid, s2_name) = (sele2.resname, sele2.resid, sele2.name)

        # if the residue name is ASP and the atom name is OD2 or OD1,

        # the atom name is changed to OD

    if s1_resname == 'LYS' and (s1_name == 'HZ1' or s1_name == 'HZ2' or s1_name == 'HZ3'):

        s1_name = 'HZ'

    if s2_resname == 'GLU' and (s2_name == 'OE1' or s2_name == 'OE2'):

        s2_name = 'OE'

    key = (s1_resname, s1_resid, s1_name, s2_resname, s2_resid, s2_name)

    # Order = 0 is hydrogen bond which has a current length of 1

    if len(current) <= order + 1:

        output[key] = 1



f = open('water2.txt', 'w')

w = MDAnalysis.analysis.hbonds.WaterBridgeAnalysis(u, 'resid 405 and segid A',

                                                   'resid 984 and segid A',

                                                   order=6)

w.run(verbose=True)

for i in range(7):

    f.write('This is order: ' + str(i))

    try:

        result = w.count_by_type(analysis_func=analysis, order=i)

        print(result)

        f.write(str(result) + '\n')

    except:

        continue

w = MDAnalysis.analysis.hbonds.WaterBridgeAnalysis(u, 'resid 405 and segid B', 'resid 984 and segid B', order=6)

w.run(verbose=True)

for j in range(7):

    f.write('This is order: '+str(j))

    try:

        result = w.count_by_type(analysis_func=analysis, order=j)

        print(result)

        f.write(str(result) + '\n')

    except:

        continue

w = MDAnalysis.analysis.hbonds.WaterBridgeAnalysis(u, 'resid 405 and segid C', 'resid 984 and segid C', order=6)

w.run(verbose=True)

for k in range(7):

    f.write('This is order: '+str(k))

    try:

        result = w.count_by_type(analysis_func=analysis, order=k)

        print(result)

        f.write(str(result) + '\n')

    except:

        continue

f.close()


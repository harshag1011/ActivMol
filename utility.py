import math

import numpy as np
from Bio.PDB import PDBList, PDBParser, NeighborSearch, Selection

# collecting metal atoms using Element Symbol method
Metal_Symbols = ["SC", "TI", " V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
                 " Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD",
                 "HF", "TA", " W", "RE", "OS", "IR", "PT", "AU", "HG",
                 "LI", "NA", " K", "RB", "CS", "BE", "MG", "CA", "SR", "BA",
                 "AL", "GA", "IN", "TL", "SN", "PB"]
radius_dic = {'C': 0.76, 'O': 0.66, 'N': 0.71, 'P': 1.06, 'S': 1.05, 'FE': 1.17,
              'ZN': 1.22, 'MG': 1.45, 'CA': 1.94, 'NA': 1.86, 'MN': 1.39, ' K': 2.27,
              'NI': 1.24, 'CU': 1.28, 'CO': 1.26, 'CD': 1.44, 'HG': 1.49, 'PT': 1.28,
              'MO': 1.39, 'AL': 1.18, 'BE': 0.9, 'BA': 2.15, 'RU': 1.34, 'V': 1.53, 'SR': 2.19,
              'CS': 2.24}


def parse_pdb(file_name_id):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(file_name_id, pdir='.', file_format='pdb')
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = parser.get_structure(file_name_id, "pdb{}.ent".format(file_name_id))
    return structure


def atoms_name(file_name_id):
    ATOM_name_dic = {}
    HETATM_name_dic = {}
    fp = open('pdb/pdb{}/pdb{}.ent'.format(file_name_id,file_name_id))
    for line in fp:
        if line.startswith('ATOM'):
            ATOM_name_dic[line[12:16].strip()] = line[76:78].strip()
        elif line.startswith('HETATM'):
            HETATM_name_dic[line[12:16].strip()] = line[76:78].strip()
    fp.close()
    return ATOM_name_dic, HETATM_name_dic


def find_metal_atoms(structure, file_name_id):
    ls = []
    # iterate over lines in pdb and getting the data from line starting from 'HETATM'
    fp = open('pdb/pdb{}/pdb{}.ent'.format(file_name_id,file_name_id))
    for line in fp:
        if line.startswith('HETATM'):
            if line[76:78] in Metal_Symbols:
                atom_name = line[12:16].strip()
                chain_id = line[21:22]
                res_id = int(line[22:26].strip())
                ls.append((atom_name, chain_id, res_id))
    fp.close()
    metal_atoms = []
    atom_all = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    id = atom.get_full_id()
                    # print(id)
                    for x in ls:
                        # print((x[0], x[1], x[2]))
                        # print((id[4][0], id[2], int(id[3][1])))
                        if x[0] == id[4][0] and x[1] == id[2] and x[2] == id[3][1]:
                            metal_atoms.append(atom)
                    atom_all.append(atom)
    return metal_atoms, atom_all


def find_pos_sites(metal_atoms, atom_all):
    # creating neighborSearch object using all atoms in a structure
    ns = NeighborSearch(atom_all)
    # searching neighbour residues in 3A radius around metal center coordinates
    pos_sites = []  # contains list(metal_atom,neighbour_residues)
    radius = 3.0
    for atom in metal_atoms:
        center = atom.get_coord()
        ns_res = ns.search(center, radius, level="R")
        pos_sites.append([atom, ns_res])
    return pos_sites


def find_active_sites(pos_sites):
    active_sites = []
    temp = []
    for itr in range(len(pos_sites)):
        temp.append(1)  # 1 represent not appended any pos_sites inside any list
    # logic for polynuclear active site
    for i in range(0, len(pos_sites)):
        ls = []
        if temp[i] == 1:
            ls.append(pos_sites[i])
            temp[i] = 0  # 0 represent appended inside list ls
        for j in range(i + 1, len(pos_sites)):
            chain_id_i = pos_sites[i][0].get_full_id()[2]
            chain_id_j = pos_sites[j][0].get_full_id()[2]
            # checking 2 list each containing tuple(atom object,neighbour residues)
            if chain_id_i == chain_id_j:
                set_i = set()
                set_j = set()
                for x in pos_sites[i][1]:
                    # print(x.get_id())
                    set_i.add((chain_id_i, x.get_id()[1]))
                for y in pos_sites[j][1]:
                    set_j.add((chain_id_j, y.get_id()[1]))
                # print(set_i)
                # print(set_j)
                if len(set_i.intersection(set_j)) != 0:  # for poly-nuclear active site
                    if temp[j] == 1:
                        ls.append(pos_sites[j])
                        temp[j] = 0
        if len(ls) != 0:
            active_sites.append(ls)
    return active_sites


def print_active_sites(active_sites,file_id):
    with open('pdb/pdb{}/{}.txt'.format(file_id,file_id),'w') as f:
        i = 1
        for site in active_sites:
            #print("*" * 20,file=f)
            print("site{}".format(i),file=f)
            # print(site)
            # print("------------------------------------------")
            for x in site:
                print("Metal atom full id {} have neighbours:-".format(x[0].get_full_id()),file=f)
                for y in x[1]:
                    print(y, end="  ",file=f)
                print('\n', file=f)
            i += 1


def get_atoms_of_residue(res):
    atom_list = Selection.unfold_entities(res, 'A')

    water_atom = None
    #print("Atom list for {} residue is\n{}".format(res, atom_list))

    if len(atom_list) == 1:  # O of h20 or metal atom - het res
        res_id = res.get_id()
        if res_id[0][0] == "W":
            water_atom = atom_list

    return atom_list, water_atom
    # return atom_list


def normalize_box(atom_list):
    min_x = min(atom.get_coord()[0] for atom in atom_list)
    max_x = max(atom.get_coord()[0] for atom in atom_list)
    max_y = max(atom.get_coord()[1] for atom in atom_list)
    min_y = min(atom.get_coord()[1] for atom in atom_list)
    min_z = min(atom.get_coord()[2] for atom in atom_list)
    max_z = max(atom.get_coord()[2] for atom in atom_list)
    # Compute the center and size of the bounding box
    center = [(min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2]
    size = [max_x - min_x, max_y - min_y, max_z - min_z]
    return center, size


sf = 10


def modify_variable(new_value):
    global sf
    sf = new_value


def normalize_metal_atoms(metal_atoms, center, size, HETATM_name_dic):
    # scaling factor is 10 here
    # sf = 10
    metal_coords = [(round(sf * ((atom.get_coord()[0] - center[0]) / size[0]) - 1, 2),
                     round(sf * ((atom.get_coord()[1] - center[1]) / size[1]) - 1, 2),
                     round(sf * ((atom.get_coord()[2] - center[2]) / size[2]) - 1, 2),
                     HETATM_name_dic[atom.get_full_id()[4][0]]) for atom in metal_atoms]
    # normal_coords = [(atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2]) for atom in metal_atoms]
    # print(normal_coords)
    return metal_coords

def normalize_atom(atom,center,size):
    atom_coord = (round(sf * ((atom.get_coord()[0] - center[0]) / size[0]) - 1, 2),
                     round(sf * ((atom.get_coord()[1] - center[1]) / size[1]) - 1, 2),
                     round(sf * ((atom.get_coord()[2] - center[2]) / size[2]) - 1, 2))
    return atom_coord

def normalize_residue(bonded_coords, center, size):
    # scaling factor is 3.5 here    this is variable
    # sf = 10
    b_coords = [((round(10 * ((x[0][0] - center[0]) / size[0]) - 1, 2),
                  round(sf * ((x[0][1] - center[1]) / size[1]) - 1, 2),
                  round(sf * ((x[0][2] - center[2]) / size[2]) - 1, 2)),
                 (round(sf * ((x[1][0] - center[0]) / size[0]) - 1, 2),
                  round(sf * ((x[1][1] - center[1]) / size[1]) - 1, 2),
                  round(sf * ((x[1][2] - center[2]) / size[2]) - 1, 2)),
                 x[2]) for x in bonded_coords]  # [(x1,y1,z1),(x2,y2,z2),(ele1,ele2)]
    return b_coords


def check_het(atom):
    res_id = atom.get_full_id()[3]
    res_type = res_id[0]
    if res_type[0] == "H":
        return True
    else:
        return False


magic = 0.2


# def modify_magic(new_value):
#     global magic
#     magic = new_value


def find_bonded_atoms(atom_list, HETATM_name_dic, ATOM_name_dic):
    bonded_atoms = []
    for i in range(0, len(atom_list)):
        for j in range(i + 1, len(atom_list)):
            atom1 = atom_list[i]
            atom2 = atom_list[j]
            if check_het(atom1):
                atom1_name = HETATM_name_dic[atom1.get_full_id()[4][0]]
            else:
                atom1_name = ATOM_name_dic[atom1.get_full_id()[4][0]]
            if check_het(atom2):
                atom2_name = HETATM_name_dic[atom2.get_full_id()[4][0]]
            else:
                atom2_name = ATOM_name_dic[atom2.get_full_id()[4][0]]

            x1, y1, z1 = atom1.get_coord()
            x2, y2, z2 = atom2.get_coord()
            d = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
            # print(d, "   ", atom1, "   ", atom2)
            r1 = radius_dic[atom1_name]
            r2 = radius_dic[atom2_name]
            radial_d = r1 + r2 + magic  # first it was 0.5
            if (d <= radial_d):
                #print(d, "   ", atom1, "   ", atom2)
                bonded_atoms.append((atom1, atom2))
                # bonded_atoms_names.append((name_dic[atom1.get_full_id()[4][0]],name_dic[atom2.get_full_id()[4][0]]))
    # return bonded_atoms,bonded_atoms_names
    return bonded_atoms


# donor_range = 3.0
#
#
# def modify_donor_range(new_value):
#     global donor_range
#     donor_range = new_value


def find_interaction_atoms(metal_atom, atom_list,HETATM_name_dic,ATOM_name_dic):
    ns = NeighborSearch(atom_list)
    # searching neighbour residues in 3A radius around metal center coordinates
    #donors = []  # contains list(metal_atom,neighbour_residues)
    donor_range = 3.0  # variable
    #print(metal_atom)
    center = metal_atom.get_coord()
    ns_atoms = ns.search(center, donor_range, level="A")

    donor_atoms = {'N', 'O', 'S', 'P'}
    # print(ns_atoms)
    # remove other atoms apart from donor_atoms
    temp = []
    for i in ns_atoms:
        if check_het(i):
            name = HETATM_name_dic[i.get_full_id()[4][0]]
        else:
            name = ATOM_name_dic[i.get_full_id()[4][0]]
        # print(name)
        if name in donor_atoms:
            temp.append(i)
    # print("temp is ", temp)
    return temp


def get_bonded_coords(bonded_atoms, HETATM_name_dic, ATOM_name_dic):
    b_coords_list = []
    for b_a in bonded_atoms:
        atom1 = b_a[0]
        atom2 = b_a[1]
        if check_het(atom1):
            atom1_name = HETATM_name_dic[atom1.get_full_id()[4][0]]
        else:
            atom1_name = ATOM_name_dic[atom1.get_full_id()[4][0]]
        if check_het(atom2):
            atom2_name = HETATM_name_dic[atom2.get_full_id()[4][0]]
        else:
            atom2_name = ATOM_name_dic[atom2.get_full_id()[4][0]]

        b_coords_list.append(((b_a[0].get_coord()[0],
                               b_a[0].get_coord()[1],
                               b_a[0].get_coord()[2]),
                              (b_a[1].get_coord()[0],
                               b_a[1].get_coord()[1],
                               b_a[1].get_coord()[2]),
                              (atom1_name, atom2_name)))
    return b_coords_list


def flat(lis):
    flatList = []
    # Iterate with outer list
    for element in lis:
        if type(element) is list:
            # Check if type is list than iterate through the sublist
            for item in element:
                flatList.append(item)
        else:
            flatList.append(element)
    return flatList

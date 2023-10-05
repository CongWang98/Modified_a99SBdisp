import MDAnalysis as mda

def save_pairs(pairs, file_name, param=None, atom_1_name='OB', atom_2_name='HB'):
    pairs = sorted(list(pairs))
    with open(file_name, 'w') as f:
        if not param:
            f.write(f'{atom_1_name}\t{atom_2_name}\n')
        else:
            assert len(pairs) == len(param)
            f.write(f'{atom_1_name}\t{atom_2_name}\tepsilon\tQQ\tsigma\n')
        for i, (at1, at2) in enumerate(pairs):
            if not param:
                f.write(f'{at1}\t{at2}\n')
            else:
                f.write(f'{at1}\t{at2}\t{param[i][0]}\t{param[i][1]}\t{param[i][2]}\n')
    return

def find_1_2_pairs(mda_universe, atom_type1, atom_type2):
    pairs_1_2 = set()
    for bond in mda_universe.bonds:
        at1, at2 = bond[0], bond[1]

        if atom_type1 == atom_type2 == at1.type == at2.type:
            pairs_1_2.add((min(at1.index, at2.index), max(at1.index, at2.index)))

        elif at1.type == atom_type1 and at2.type == atom_type2:
            pairs_1_2.add((at1.index, at2.index))
        elif at1.type == atom_type2 and at2.type == atom_type1:
            pairs_1_2.add((at2.index, at1.index))
    return pairs_1_2


def find_1_3_pairs(mda_universe, atom_type1, atom_type2):
    pairs_1_3 = set()
    for angle in mda_universe.angles:
        at1, at3 = angle[0], angle[2]
        if atom_type1 == atom_type2 == at1.type == at3.type:
            pairs_1_3.add((min(at1.index, at3.index), max(at1.index, at3.index)))
        elif at1.type == atom_type1 and at3.type ==  atom_type2:
            pairs_1_3.add((at1.index, at3.index))
        elif at1.type ==  atom_type2 and at3.type == atom_type1:
            pairs_1_3.add((at3.index, at1.index))
    return pairs_1_3


def find_1_4_pairs(mda_universe, atom_type1, atom_type2):
    pairs_1_4 = set()
    for dihedral in mda_universe.dihedrals:
        at1, at4 = dihedral[0], dihedral[3]
        if atom_type1 == atom_type2 == at1.type == at4.type:
            pairs_1_4.add((min(at1.index, at4.index), max(at1.index, at4.index)))
        elif at1.type == atom_type1 and at4.type == atom_type2:
            pairs_1_4.add((at1.index, at4.index))
        elif at1.type == atom_type2 and at4.type == atom_type1:
            pairs_1_4.add((at4.index, at1.index))
    return pairs_1_4


def find_all_pairs(mda_universe, atom_type1, atom_type2):
    if atom_type1 == atom_type2:
        at_lis = [at.index for at in mda_universe.atoms if at.type == atom_type1]
        return set([(at_lis[i], at_lis[j]) for i in range(len(at_lis)) for j in range(i + 1, len(at_lis))])
    at1_lis = []
    at2_lis = []
    for at in mda_universe.atoms:
        if at.type == atom_type1:
            at1_lis.append(at.index)
        elif at.type == atom_type2:
            at2_lis.append(at.index)
    return set([(at1, at2) for at1 in at1_lis for at2 in at2_lis])



def get_nb_param(pairs, combined_sigma, combined_epsilon, charge_list, pair_type, fudgeLJ=0.5, fudgeQQ=0.8333):
    if pair_type == 'LR':
        return [(combined_epsilon, charge_list[at1]*charge_list[at2], combined_sigma) for at1, at2 in pairs]
    elif pair_type == 'SR':
        return [(0, 0, 1) for _ in pairs]
    elif pair_type == '1-4':
        return [(fudgeLJ*combined_epsilon, fudgeQQ*charge_list[at1]*charge_list[at2], combined_sigma) for at1, at2 in pairs]
    else:
        raise TypeError('Pair type can only be "LR", "SR", or "1-4".')


if __name__ == "__main__":
    #[ nonbond_params ]
    #; i    j    funct    sigma   epsilon
    #  OB   HB     1      0.150   1.2552
    top_file = 'disp_dry.top'
    gmx_path='/home/gridsan/congwang/software/gromacs/share/gromacs/top'
    u = mda.Universe(top_file, topology_format='ITP',infer_system=True, include_dir=gmx_path)

    
    atp1, atp2 = 'OB', 'HB'
    combined_sigma = 0.150
    combined_epsilon = 1.2552

    pair_1_2 = find_1_2_pairs(u, atp1, atp2)
    pair_1_3 = find_1_3_pairs(u, atp1, atp2)
    pair_1_4 = find_1_4_pairs(u, atp1, atp2)
    pair_all = find_all_pairs(u, atp1, atp2)

    pair_short = sorted(list(pair_1_2 | pair_1_3))
    pair_long = sorted(list(pair_all - pair_1_2 - pair_1_3 - pair_1_4))
    pair_1_4 = sorted(list(pair_1_4))


    save_pairs(pair_short, 'pairs_short_mtd3.txt')
    save_pairs(pair_long, 'pairs_long_mtd3.txt')
    save_pairs(pair_1_4, 'pairs_1_4_mtd3.txt')

    param_short = get_nb_param(pair_short, combined_sigma, combined_epsilon, u.atoms.charges, 'SR')
    param_long = get_nb_param(pair_long, combined_sigma, combined_epsilon, u.atoms.charges, 'LR')
    param_1_4 = get_nb_param(pair_1_4, combined_sigma, combined_epsilon, u.atoms.charges, '1-4')

    save_pairs(pair_short, 'pairs_short_param.txt', param_short, atp1, atp2)
    save_pairs(pair_long, 'pairs_long_param.txt', param_long, atp1, atp2)
    save_pairs(pair_1_4, 'pairs_1_4_param.txt', param_1_4, atp1, atp2)






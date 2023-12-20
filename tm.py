from Bio.PDB import PDBParser

#slownik z przybliżonymi promieniami van der Waalsa
vdw_radii = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.85,
    "P": 1.90,
}


def calculate_clash_score(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    clash_score = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_name = atom.get_id()[0]
                    vdw_radius = vdw_radii.get(atom_name, 1.5)  # Domyślny promień 1.5 A

                    #promień clash to promień + 0.4A (r)
                    clash_radius = vdw_radius + 0.4

                    for other_model in structure:
                        for other_chain in other_model:
                            for other_residue in other_chain:
                                for other_atom in other_residue:
                                    if atom != other_atom:
                                        distance = atom - other_atom
                                        if distance < clash_radius:
                                            clash_score += 1

    return clash_score


pdb_filename = "R1107_reference.pdb"
score = calculate_clash_score(pdb_filename)
print("Clash score dla pliku: ", pdb_filename, score)

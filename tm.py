from Bio.PDB import *

def calculate_clash_score(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("my_structure", pdb_file)

    clash_count = 0
    total_atoms = 0

    vdw_radii = {
        "H": 1.20,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "S": 1.85,
        "P": 1.90,
    }

    for model in structure:
        for chain in model:
            ns = NeighborSearch(list(chain.get_atoms()))
            for residue in chain:
                atoms_to_check = [atom for atom in residue.get_atoms() if atom.element in vdw_radii]

                for atom in atoms_to_check:
                    total_atoms += 1
                    atom_coord = atom.get_coord()

                    close_atoms = ns.search(atom_coord, vdw_radii[atom.element])
                    for close_atom in close_atoms:
                        if close_atom.get_parent().get_full_id() != residue.get_full_id() and close_atom.element in vdw_radii:
                            vdw_sum = vdw_radii.get(atom.element) + vdw_radii.get(close_atom.element)
                            atom_distance = atom - close_atom

                            if atom_distance < vdw_sum - 0.4:
                                clash_count += 1

    return 1000 * (clash_count / total_atoms)

# Przykładowe użycie
pdb_file_path = "C:/Users/karol/PycharmProjects/RSMD/R1107_reference.pdb"  # Zmień na właściwą ścieżkę
#pdb_file_path = "C:/Users/karol/PycharmProjects/RSMD/4.pdb"  # Zmień na właściwą ścieżkę
#pdb_file_path = "C:/Users/karol/PycharmProjects/RSMD/2.pdb"  # Zmień na właściwą ścieżkę
#pdb_file_path = "C:/Users/karol/PycharmProjects/RSMD/test.txt"  # Zmień na właściwą ścieżkę



score = calculate_clash_score(pdb_file_path)
print(f"Clash score: {score}")

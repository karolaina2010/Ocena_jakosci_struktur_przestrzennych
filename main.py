from Bio.PDB import *
import numpy as np
import io

#liczenie rmsd
def calculate_rmsd(atom_list_1, atom_list_2):
    if len(atom_list_1) != len(atom_list_2):
        raise ValueError("Listy atomow maja rozna dlugosc")

    squared_distances = [(np.linalg.norm(atom_1.get_coord() - atom_2.get_coord())) ** 2
                         for atom_1, atom_2 in zip(atom_list_1, atom_list_2)]

    rmsd = np.sqrt(np.mean(squared_distances))
    return rmsd

#superpozucja
def process_model(model_data, model_number):
    model_structure = parser.get_structure("model", io.StringIO(''.join(model_data)))
    super_imposer = Superimposer()
    ref_atoms = []
    model_atoms = []

    for ref_chain in reference_structure.get_chains():
        model_chain_id = ref_chain.id
        try:
            model_chain = model_structure[0][model_chain_id]
        except KeyError:
            continue

        for ref_residue, model_residue in zip(ref_chain, model_chain):
            for ref_atom, model_atom in zip(ref_residue.get_atoms(), model_residue.get_atoms()):
                ref_atoms.append(ref_atom)
                model_atoms.append(model_atom)

    if len(ref_atoms) > 0 and len(model_atoms) > 0:
        if len(ref_atoms) != len(model_atoms):
            min_length = min(len(ref_atoms), len(model_atoms))
            ref_atoms = ref_atoms[:min_length]
            model_atoms = model_atoms[:min_length]

        super_imposer.set_atoms(ref_atoms, model_atoms)
        super_imposer.apply(model_structure[0].get_atoms())

        calculated_rmsd = calculate_rmsd(ref_atoms, model_atoms)
        print(f"RMSD dla modelu {model_number}: {calculated_rmsd}")
    else:
        print(f"Brak atomów do superpozycji dla modelu {model_number}")

#odczytywanie modeli
parser = PDBParser(QUIET=True)
reference_structure = parser.get_structure("reference", "R1107_reference.pdb")

models_file = "R1107TS081.pdb"

#glowna petla
with open(models_file, 'r') as file:
    model_started = False
    current_model = []
    model_number = 1

    for line in file:
        if line.startswith("MODEL"):
            if current_model:
                process_model(current_model, model_number)
                model_number += 1
            current_model = []

            model_started = True
        elif line.startswith("ENDMDL"):
            model_started = False
        elif model_started:
            current_model.append(line)

#przetworzenie ostatniego modelu
if current_model:
    process_model(current_model, model_number)

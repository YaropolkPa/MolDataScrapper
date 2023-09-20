import csv
import os
import re
import requests
from bs4 import BeautifulSoup

dataset_file_name = 'datasets/test.csv'
list_allowed_atoms = ["H", "C", "N", 'O', "S", "F", "P", "Cl", 'Br', "I"]
max_heavy_atoms = 50
max_Fsp3 = 1.0
max_iter = 20000000


checkpoint = 0
if not os.path.isfile(dataset_file_name):
    with open(dataset_file_name, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ZINC', "SMILES", "FORMULA", "mol_weigth", "num_heavy_atoms", "Fsp3"])
else:
    with open(dataset_file_name, 'r') as file:
        csvreader = list(csv.reader(file))
        checkpoint = int(csvreader[-1][0])


class MoleculeInfo:
    def __init__(self, ZINC, SMILES, FORMULA, MWT, LOGP, H_A, FSP3):
        self.zinc = str(ZINC)
        self.smiles = str(SMILES)
        self.formula = str(FORMULA)
        self.molecular_weight = float(MWT)
        self.logP = float(LOGP)
        self.num_heavy_atoms = int(H_A)
        self.Fsp3 = float(FSP3)

        self.atoms_list = split_chemical_formula(self.formula)

    def is_allowed(self, allowed_atoms, max_h_a, max_fsp3):
        if max_h_a >= self.num_heavy_atoms and \
                max_fsp3 >= self.Fsp3 and set(self.atoms_list).issubset(set(allowed_atoms)):
            return True
        return False

    def __str__(self):
        return f'''
        ZINC: {self.zinc}
        SMILES: {self.smiles}
        MWT: {self.molecular_weight}
        HEAVY_ATOMS: {self.num_heavy_atoms}
        FSP3: {self.Fsp3}     
        '''


def split_chemical_formula(input_formula):
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, input_formula)

    atoms = []
    for match in matches:
        atom, count = match
        atoms.append(atom)

    return atoms


for idx in range(checkpoint+1, max_iter):
    url = f'https://zinc.docking.org/substances/ZINC{"{:0>12}".format(idx)}/'
    page = requests.get(url)

    soup = BeautifulSoup(page.text, features="html.parser")

    if soup.title.text == "Page Not Found":
        continue

    # Mwt = 3, logP = 4, Formula = 6, Heavy Atoms = 8, Fsp3 = 10
    all_table_args = soup.find_all('td')
    lst_all_table_args = [item.text for item in all_table_args]

    # Molecule parameters
    mol_weight = lst_all_table_args[3]
    logP = lst_all_table_args[4]
    formula = lst_all_table_args[6]
    heavy_atoms = lst_all_table_args[8]
    fsp3 = lst_all_table_args[10]
    smiles = soup.input['value']

    mol_obj = MoleculeInfo(
        ZINC="{:0>12}".format(idx),
        SMILES=smiles,
        FORMULA=formula,
        MWT=mol_weight,
        LOGP=logP,
        H_A=heavy_atoms,
        FSP3=fsp3
    )
    print(mol_obj)

    if mol_obj.is_allowed(list_allowed_atoms, max_heavy_atoms, max_Fsp3):
        print("Molecule added! \n")
        with open(dataset_file_name, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([mol_obj.zinc,
                             mol_obj.smiles,
                             mol_obj.formula,
                             mol_obj.molecular_weight,
                             mol_obj.num_heavy_atoms,
                             mol_obj.Fsp3])

print("Done")

import csv
import os
import re
import requests
from bs4 import BeautifulSoup
from multiprocessing import Pool
import time


dataset_file_name = 'datasets/zinc_drugs_2k_v2.csv'
source = "https://zinc.docking.org/substances/subsets/in-vivo/"
list_allowed_atoms = ["H", "C", "N", 'O', "S", "F", "P", "Cl", 'Br', "I"]
max_heavy_atoms = 50
max_Fsp3 = 1.0
# max_iter = 20000000


checkpoint = 0
if not os.path.isfile(dataset_file_name):
    with open(dataset_file_name, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ZINC', "SMILES", "FORMULA", "mol_weigth", "num_heavy_atoms", "logP", "Fsp3"])


# else:
#     with open(dataset_file_name, 'r') as file:
#         csvreader = list(csv.reader(file))
#         checkpoint = int(csvreader[-1][0])


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
        FORMULA: {self.formula}
        MWT: {self.molecular_weight}
        HEAVY_ATOMS: {self.num_heavy_atoms}
        logP: {self.logP}
        FSP3: {self.Fsp3}     
        '''


def parse(mol_url):
    full_url = "https://zinc.docking.org" + mol_url
    dat = {'q': 'goog'}
    mol_info = []

    # Mwt = 3, logP = 4, Formula = 6, Heavy Atoms = 8, Fsp3 = 10
    try:
        page = requests.get(full_url, params=dat, headers={'User-Agent': 'Mozilla/5.0'}, timeout=10)
        time.sleep(2)
        if page.status_code == 200:
            soup = BeautifulSoup(page.text, features="html.parser")
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
                ZINC=mol_url.split('/')[-2],
                SMILES=smiles,
                FORMULA=formula,
                MWT=mol_weight,
                LOGP=logP,
                H_A=heavy_atoms,
                FSP3=fsp3
            )
            if mol_obj.is_allowed(list_allowed_atoms, max_heavy_atoms, max_Fsp3):
                print(mol_obj)
                print("Molecule added! \n")
                mol_info = [mol_obj.zinc,
                            mol_obj.smiles,
                            mol_obj.formula,
                            mol_obj.molecular_weight,
                            mol_obj.num_heavy_atoms,
                            mol_obj.logP,
                            mol_obj.Fsp3]
    except Exception as e:
        print(e)
    finally:
        if len(mol_info) != 0:
            return mol_info
        else:
            return None


def url_list_from_page(page_url):
    page = requests.get(page_url, headers={'User-Agent': 'Mozilla/5.0'})
    soup = BeautifulSoup(page.text, features="html.parser")
    urls_list = [item.find('a').get('href') for item in soup.find_all("h4")]
    return urls_list


def split_chemical_formula(input_formula):
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, input_formula)

    atoms = []
    for match in matches:
        atom, count = match
        atoms.append(atom)

    return atoms


if __name__ == "__main__":
    for page in range(19, 1147):
        url = f'{source}?page={page}'
        molecules_urls = url_list_from_page(url)
        p = Pool(20)
        records = p.map(parse, molecules_urls)
        p.terminate()
        p.join()
        with open(dataset_file_name, 'a', newline='') as f:
            for item in records:
                writer = csv.writer(f)
                if item is None:
                    continue
                else:
                    writer.writerow(list(item))

    print("Done")

# if __name__ == "__main__":
#     for page in range(1, 30):
#         url = f'{source}?page={page}'
#         molecules_urls = url_list_from_page(url)
#         for mol_url in molecules_urls:
#             mol_info_list = parse(mol_url)
#             with open(dataset_file_name, 'a', newline='') as f:
#                     writer = csv.writer(f)
#                     writer.writerow(mol_info_list)

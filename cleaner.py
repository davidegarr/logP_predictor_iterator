import csv
from rdkit import Chem
from rdkit.Chem import Descriptors

def main():
    with open("raw_data.csv", "r") as raw_data:
        csvreader = csv.DictReader(raw_data, delimiter = ";")
        
        #opens target file and creates the writer
        with open("clean_data.csv", "w", newline='') as clean_data:
            headers = ["ChEMBL_ID", "SMILES", "AlogP", "MolWt", "NumHAcceptors", "NumHDonors","NumHeteroatoms", "NumRotatableBonds", "NumAromaticCarbocycles", "NumAromaticHeterocycles", "NumSaturatedHeterocycles", "NumSaturatedCarbocycles", "FractionCSP3", "TPSA"]
            csvwriter = csv.DictWriter(clean_data, fieldnames=headers)
            
            csvwriter.writeheader()
            #copies only the indicated values in the new row if all of them are present in source
            n = 0

            for row in csvreader:
                if row["AlogP"] == "None" or row["Smiles"] == "":
                    continue

                #creates molecule object
                molecule = Chem.MolFromSmiles(str(row["Smiles"]))
                #checks that the molecule is created correctly and that it is organic (has at least one carbon atom)
                if molecule is None:
                    continue
                if not has_carbon(molecule):
                    continue

                features = compute_descriptors(molecule)

                #writes into the new file only relevant fields.
                csvwriter.writerow({
                    "ChEMBL_ID": str(row["ChEMBL ID"]), 
                    "SMILES": str(row["Smiles"]),
                    "AlogP": float(row["AlogP"]),
                    "MolWt": features["MolWt"],
                    "NumHAcceptors" : features["NumHAcceptors"],
                    "NumHDonors": features["NumHDonors"],
                    "NumHeteroatoms": features["NumHeteroatoms"],
                    "NumRotatableBonds": features["NumRotatableBonds"],
                    "NumAromaticCarbocycles": features["NumAromaticCarbocycles"],
                    "NumAromaticHeterocycles": features["NumAromaticHeterocycles"],
                    "NumSaturatedHeterocycles": features["NumSaturatedHeterocycles"],
                    "NumSaturatedCarbocycles": features["NumSaturatedCarbocycles"],
                    "FractionCSP3": features["FractionCSP3"],
                    "TPSA": features["TPSA"]
                    

                })

def has_carbon(molecule):
    #checks if at least there is a carbon atom in the molecule
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 6:
            return True
    return False

def compute_descriptors(molecule):
    
    features = {}

    features["MolWt"] = float(Descriptors.MolWt(molecule))
    features["NumHAcceptors"] = float(Descriptors.NumHAcceptors(molecule))
    features["NumHDonors"] = float(Descriptors.NumHDonors(molecule))
    features["NumHeteroatoms"] = float(Descriptors.NumHeteroatoms(molecule))
    features["NumRotatableBonds"] = float(Descriptors.NumRotatableBonds(molecule))
    features["NumAromaticCarbocycles"] = float(Descriptors.NumAromaticCarbocycles(molecule))
    features["NumAromaticHeterocycles"] = float(Descriptors.NumAromaticHeterocycles(molecule))
    features["NumSaturatedHeterocycles"] = float(Descriptors.NumSaturatedHeterocycles(molecule))
    features["NumSaturatedCarbocycles"] = float(Descriptors.NumSaturatedCarbocycles(molecule))
    features["FractionCSP3"] = float(Descriptors.FractionCSP3(molecule))
    features["TPSA"] = float(Descriptors.TPSA(molecule))

    return features

main()

from sklearn.preprocessing import StandardScaler
from tensorflow.keras.models import load_model
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import numpy as np
import pandas as pd
from joblib import load
from collections import deque

def main():
    user_smiles = input("Enter SMILES: ")
    initial_molecule = Chem.MolFromSmiles(user_smiles)

    if initial_molecule is None:
        print("Invalid SMILES string.")

    initial_logP = perform_prediction(initial_molecule)
    
    print(f"Predicted AlogP: {initial_logP[0][0]:.3f}")

    direction = input("Would you like to increase or decrease the logP?: ")
    
    top_3_molecules = tree_search(initial_molecule, initial_logP, direction)
    print("According to your direction, the top 3 molecules are:", top_3_molecules)


def tree_search(initial_molecule, initial_logP, direction, depth=5):
    visited = set()
    queue = deque([(initial_molecule, initial_logP, 0)])
    best_molecules = []

    while queue:
        current_molecule, current_logP, current_depth = queue.popleft()
        
        # Check if this molecule is better than previously seen molecules
        if is_better(current_logP, initial_logP, direction):
            best_molecules.append((current_molecule, current_logP))

        # If we've reached our desired depth, stop expanding this branch
        if current_depth >= depth:
            continue
        
        # Expand the current molecule by applying all transformations
        transformed_molecules = apply_all_transformations(current_molecule)
        for smiles in transformed_molecules:
            if smiles not in visited and Chem.MolFromSmiles(smiles) != None:
                visited.add(smiles)
                logP = perform_prediction(Chem.MolFromSmiles(smiles))
                queue.append((smiles, logP, current_depth + 1))

    # Return the top 3 molecules based on their logP values, according to the desired direction
    top_3_molecules = sorted(best_molecules, key=lambda x: x[1], reverse=(direction == 'increase'))[:3]
    return top_3_molecules


def is_better(new_logP, best_logP, direction):
    if direction == "increase" and new_logP > best_logP:
        return True
    elif direction == "decrease" and new_logP < best_logP:
        return True
    return False


def perform_prediction(molecule):
    #get descriptors from string
    features = compute_descriptors(molecule)

    scaler = load("scaler.joblib")
    loaded_model = load_model("./predictor_v1")

    feature_df = pd.DataFrame([features])
    feature_scaled = scaler.transform(feature_df)

    # Make the prediction
    return loaded_model.predict(feature_scaled)


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


def oxidize_primary_alcohol(molecule):
    rxn = AllChem.ReactionFromSmarts('[CH2:1][O:2]>>[CH:1]=[O:2]')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)
    

def reduce_aldehyde(molecule):
    rxn = AllChem.ReactionFromSmarts('[CH:1]=[O:2]>>[CH2:1][O:2]')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)
    

def oxidize_secondary_alcohol(molecule):
    rxn = AllChem.ReactionFromSmarts('[CH:1][O:2]>>[C:1]=[O:2]')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)
    

def reduce_ketone(molecule):
    rxn = AllChem.ReactionFromSmarts('[CH0:1]=[O:2]>>[CH1:1][O:2]')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)


def acetylate_primary_amine(molecule):
    rxn = AllChem.ReactionFromSmarts('[NH2:1]>>[NH:1]-C(=O)C')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)


def acylate_secondary_amine(molecule):
    rxn = AllChem.ReactionFromSmarts('[NH:1]>>[N:1]-C(=O)C')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)


def alcohol_benzylation(molecule):
    rxn = AllChem.ReactionFromSmarts('[C:1][O:2]>>[C:1]-[O:2]-CC1=CC=CC=C1')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)


def benzylate_primary_amine(molecule):
    rxn = AllChem.ReactionFromSmarts('[NH2:1]>>[NH:1]-CC1=CC=CC=C1')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)


def benzylate_secondary_amine(molecule):
    rxn = AllChem.ReactionFromSmarts('[NH:1]>>[N:1]-CC1=CC=CC=C1')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)


def ester_to_ether(molecule):
    rxn = AllChem.ReactionFromSmarts('[C:1][C:2](=[O:3])[O:4][C:5]>>[C:1][C:2][O:4][C:5]')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)
    

def ether_to_ester(molecule):
    rxn = AllChem.ReactionFromSmarts('[C:1][C:2][O:4][C:5]>>[C:1][C:2](=[O:3])[O:4][C:5]') 
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)
    

def amide_to_amine(molecule):
    rxn = AllChem.ReactionFromSmarts('[*:1][C:2](=[O:3])[N:4]>>[*:1][C:2][N:4]')
    ps = rxn.RunReactants((molecule,))
    for x in ps:
        product = x[0]
        return Chem.MolToSmiles(product)


def apply_all_transformations(molecule):
    if type(molecule) == str:
        molecule = Chem.MolFromSmiles(molecule)

    products = []

    reaction_functions = [
        ("Oxidize Primary Alcohol", oxidize_primary_alcohol),
        ("Reduce Aldehyde", reduce_aldehyde),
        ("Oxidize Secondary Alcohol", oxidize_secondary_alcohol),
        ("Reduce Ketone", reduce_ketone),
        ("Acetylate Primary Amine", acetylate_primary_amine),
        ("Acylate Secondary Amine", acylate_secondary_amine),
        ("Alcohol Benzylation", alcohol_benzylation),
        ("Benzylate Primary Amine", benzylate_primary_amine),
        ("Benzylate Secondary Amine", benzylate_secondary_amine),
        ("Ester to Ether", ester_to_ether),
        ("Ether to Ester", ether_to_ester),
        ("Amide to Amine", amide_to_amine)
    ]
    
    for name, func in reaction_functions:
        product = func(molecule)
        if product is not None:  # skip Nones
            products.append(product)
    
    return products


if __name__ == "__main__":
    main()
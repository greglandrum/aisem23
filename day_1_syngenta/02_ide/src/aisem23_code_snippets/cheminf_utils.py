from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw


def get_tanimoto_similarity(mol1, mol2, useChiral=False):

    if mol1 is None or mol2 is None:
        return 0.0

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024,
                                                useChirality=useChiral)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024,
                                                useChirality=useChiral)

    tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)

    return round(tanimoto, 2)


def get_sim_from_smiles(smiles1, smiles2, useChiral=False):

    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        return 0.0

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024,
                                                useChirality=useChiral)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024,
                                                useChirality=useChiral)

    tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)

    return round(tanimoto, 2)


def draw_molecule(rdk_mol):
    if rdk_mol is None:
        return None

    return Draw.MolToImage(rdk_mol)


def draw_molecule_from_smiles(smiles):
    rdk_mol = Chem.MolFromSmiles(smiles)
    if rdk_mol is None:
        return None

    return Draw.MolToImage(rdk_mol)

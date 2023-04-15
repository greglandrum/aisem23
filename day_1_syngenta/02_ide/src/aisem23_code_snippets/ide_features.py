from cheminf_utils import draw_molecule_from_smiles, get_sim_from_smiles
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import sys


# rename to get_tanimoto_similarity (using 'Rename Symbol')
# rename test to fp1 & reference to fp2 (using 'Rename Symbol')
# add chirality as an agrument (edit using multicursor)
def get_similarity(smiles1, smiles2):

    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        return 0.0

    test = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024,
                                                 useChirality=False)
    reference = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024,
                                                      useChirality=False)

    tanimoto = DataStructs.TanimotoSimilarity(test, reference)

    return round(tanimoto, 2)


if __name__ == "__main__":

    print(sys.version)

    isomers = ["F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2",
               "F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2",
               "F[C@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2",
               "F[C@@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2"]

    self_similarity = get_similarity(isomers[0], isomers[0])

    print(f"A molecule has {self_similarity:.0%} similarity to itself")

    # data structure to hold similarity results
    sim_no_stereo = [[0]*len(isomers) for i in range(len(isomers))]

    for i in range(len(isomers)):
        for j in range(len(isomers)):
            sim_no_stereo[i][j] = get_similarity(isomers[i], isomers[j])

    print("Similarity matrix for al the 4 stereoisomers:")
    for row in sim_no_stereo:
        print(row)

    # display the same matrix but with useChirality=True

    


    # display all 4 isomers
    columns = 2
    rows = 2

    fig = plt.figure(figsize=(8, 8))
    for i in range(0, len(isomers)):
        img = draw_molecule_from_smiles(isomers[i])
        fig.add_subplot(rows, columns, i + 1)
        plt.axis('off')
        plt.imshow(img)
    plt.show()

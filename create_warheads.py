from rdkit import Chem
from rdkit.Chem import BRICS
from typing import List, Optional

def create_pairwise_warheads(smiles: str) -> Optional[List[str]]:
    """
    Create pairwise warheads from a given SMILES string using BRICS decomposition.

    Args:
        smiles (str): The input SMILES string.

    Returns:
        Optional[List[str]]: A list of pairwise warheads, or None if an error occurs.

    Raises:
        ValueError: If the input SMILES is invalid or BRICS decomposition fails.
    """
    try:
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        # Perform BRICS decomposition
        brics_fragments = list(BRICS.BRICSDecompose(mol, singlePass=True))
        
        # Check if we have at least two fragments
        if len(brics_fragments) < 2:
            return []

        # Remove the last element (original molecule) and create pairwise warheads
        brics_fragments = brics_fragments[:-1]
        warheads = []

        for i, frag1 in enumerate(brics_fragments):
            for frag2 in brics_fragments[i+1:]:
                if frag1 != frag2:
                    warhead = f"{frag1}|{frag2}"
                    warheads.append(warhead)

        return warheads

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None

def write_warheads_to_file(warheads: List[str], filename: str = "warheads.smi"):
    """
    Write warheads to a file, one per line.

    Args:
        warheads (List[str]): List of warhead strings.
        filename (str): Name of the output file.
    """
    with open(filename, "w") as f:
        for warhead in warheads:
            f.write(f"{warhead}\n")
    print(f"Warheads written to {filename}")

# Example usage:
smiles = "Cc1cc(C(F)(F)F)nn1CC(=O)N1CCC(c2nc(C3=NOC(c4c(F)cccc4F)C3)cs2)CC1"
result = create_pairwise_warheads(smiles)

if result is not None:
    print(f"Number of warheads: {len(result)}")
    print("First 5 warheads:")
    for warhead in result[:5]:
        print(warhead)
    
    # Write warheads to file
    write_warheads_to_file(result)
else:
    print("Failed to create warheads.")
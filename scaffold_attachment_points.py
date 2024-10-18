from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
import py3Dmol
from IPython.display import display
import ipywidgets as widgets

def extract_scaffold_with_controlled_attachments(smiles, attachment_indices=None):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    
    if attachment_indices is None:
        # If no specific indices are provided, use the original method
        exit_vectors = []
        for atom in scaffold.GetAtoms():
            if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() > 0:
                exit_vectors.append(atom.GetIdx())
        attachment_indices = exit_vectors[:2]
    
    scaffold = Chem.RWMol(scaffold)
    
    for i, idx in enumerate(attachment_indices, start=1):
        if idx < scaffold.GetNumAtoms():
            scaffold.ReplaceAtom(idx, Chem.Atom(0))
            scaffold.GetAtomWithIdx(idx).SetIsotope(i)
        else:
            print(f"Warning: Atom index {idx} is out of range. Skipping this attachment point.")

    # Use kekuleSmiles=True to avoid aromatic notation
    scaffold_smiles = Chem.MolToSmiles(scaffold, kekuleSmiles=True, canonical=True, isomericSmiles=False)

    scaffold_smiles = scaffold_smiles.replace('*:1', '[*:1]').replace('*:2', '[*:2]')

    return scaffold_smiles

# Example usage
molecule_smiles = "Fc1cccc(F)c1C1CC(c2csc(N3CC=C(c4ccc(C(F)(F)F)cn4)CC3)n2)=NO1"

# Use the function without specifying attachment points (original behavior)
result1 = extract_scaffold_with_controlled_attachments(molecule_smiles)
print("Result with automatic attachment points:", result1)

# Use the function with specified attachment points
result2 = extract_scaffold_with_controlled_attachments(molecule_smiles, attachment_indices=[10, 20])
print("Result with specified attachment points:", result2)



def smiles_to_3d_interactive_with_selection(smiles, width=400, height=400):
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    
    # Add hydrogens and generate 3D coordinates
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    # Convert to PDB format
    pdb = Chem.MolToPDBBlock(mol)
    
    # Create py3Dmol view
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb, 'pdb')
    
    # Style the molecule
    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.3}})
    
    # Create a list to store selected atoms
    selected_atoms = []
    
    # Create output widget for displaying selected atoms
    output = widgets.Output()
    
    # Function to update selected atoms display
    def update_selected_atoms():
        with output:
            output.clear_output()
            print(f"Selected atoms: {', '.join(map(str, selected_atoms))}")
    
    # Add click events for atoms
    view.setClickable({
        'clickable': 'true',
        'callback': '''function(atom,viewer) {
            var elem = atom.elem;
            var index = atom.index;
            var selected = !atom.selected;
            atom.selected = selected;
            if (selected) {
                viewer.setStyle({serial: atom.serial}, {stick: {color: 'red'}, sphere: {color: 'red'}});
            } else {
                viewer.setStyle({serial: atom.serial}, {stick: {color: 'gray'}, sphere: {color: 'gray'}});
            }
            viewer.render();
            Jupyter.notebook.kernel.execute('selected_atoms = ' + JSON.stringify(viewer.getModel().selectedAtoms()) + '\nupdate_selected_atoms()');
        }'''
    })
    
    view.zoomTo()
    
    # Display the view and the output widget
    display(view, output)
    
    return view, output

# Example usage
smiles = "CC1=C(C(=O)NC2=C(C=CC(=C2)C#N)F)C(=C(N1)C3=CC=CC=C3)O"
view, output = smiles_to_3d_interactive_with_selection(smiles)
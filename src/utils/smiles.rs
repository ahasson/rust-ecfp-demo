use chemcore::{daylight::{ read_smiles, SmilesInputError }, molecule::DefaultMolecule};
use std::collections::HashMap;

pub fn smiles_to_mol(smiles: &str, map:Option<&mut HashMap<usize, u16>> ) -> Result<DefaultMolecule, SmilesInputError> {
    let molecule = read_smiles(smiles, map)?;
    
    Ok(molecule)
}
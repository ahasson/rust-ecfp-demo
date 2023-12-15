use chemcore::daylight::read_smiles;
use rust::fingerprints::ecfp;
use rust::utils::bitvect;
use std::collections::HashMap;
use std::io::{self, BufRead};

fn main() {
    let stdin = io::stdin();

    for line in stdin.lock().lines() {
        let smiles = line.expect("Failed to read line");
        let mut map = HashMap::new();
        let molecule = read_smiles(&smiles, Some(&mut map));
        if molecule.is_ok() {
            let fingerprints = ecfp::mol_to_features(molecule.unwrap(), 2);
            let mut _bitvect = bitvect::hash_to_bitvect(fingerprints, 64);
            // println!("{:?}", _bitvect);
        } else {
            eprintln!("Error: {:?} - {:?}", molecule.err(), smiles);
        }
    }
}

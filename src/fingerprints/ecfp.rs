use chemcore::molecule::{DefaultMolecule, Molecule};
use gamma::graph::Graph;
use std::collections::HashSet;
use std::fmt::Display;

use crate::utils::atom::invariant_identifier;
use crate::utils::hash::{hash_64, HashBaseType};

#[derive(Debug, Clone, Hash)]
struct Pair {
    atom: HashBaseType,
    bond: u8,
}

fn sort_pairs(pairs: &mut Vec<Pair>) {
    // Sort the pairs by atom type, then bond type, when SORT_BY_ATOM
    pairs.sort_by(|a, b| {
        if a.bond == b.bond {
            a.atom.cmp(&b.atom)
        } else {
            a.bond.cmp(&b.bond)
        }
    });
}

fn pairs_to_hash(pairs: &Vec<Pair>) -> HashBaseType {
    // Sort the pairs, and then hash the resulting array
    let mut sorted_pairs = pairs.to_vec();
    sort_pairs(&mut sorted_pairs);

    hash_64(&sorted_pairs)
}
#[derive(Debug, Clone, Eq, PartialEq)]
struct BondSet {
    bonds: HashSet<u8>, // The bonds in the set
    depth: u8,          // Depth the set was added at
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ECFPFeature {
    fingerprint: HashBaseType,
    bondset: BondSet,
}

impl Display for ECFPFeature {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "F: {:?}, Bondset: {:?}:@{:?}", self.fingerprint, self.bondset.bonds, self.bondset.depth)
    }
}
/*

   1. If the features were generated from different numbers of iterations,
    then the feature with the larger number of iterations is rejected.
   2. If the features were generated from the same number of iterations,
    then the larger hashed identifier value (interpreted as an integer) is rejected.
   (Rogers, D. R.; Hahn, M. W. (2006).
    "Extended-connectivity fingerprints".
    Journal of Chemical Information and Modeling. 46 (6): 1420–1425. doi:10.1021/ci060265f. PMID 17181570.)
*/
fn all_extend_bondset(features: Vec<ECFPFeature>, pairs: &Vec<Pair>, bondset: &mut BondSet){
    for feature in pairs {
        if let Some(found_feature) = features.iter().find(|f| f.fingerprint == feature.atom) {
            for bond in found_feature.bondset.bonds.iter() {
                bondset.bonds.insert(*bond);
            }
        }
    }
}
fn add_feature(
    features: &mut Vec<ECFPFeature>,
    current_features: &mut Vec<ECFPFeature>,
    pairs: &Vec<Pair>,
    mut bondset: BondSet,
) -> HashBaseType {
    // Extend the bondset of the source feature, by adding the bondsets of any neighbouring features.
    // This is used to identify unique substructures.
    
    let mut combined_features = features.clone();
    combined_features.extend(current_features.clone());
    
    all_extend_bondset(combined_features, pairs, &mut bondset);
    let bonds = bondset.bonds.clone();
    let depth = bondset.depth;
    

    let hash = pairs_to_hash(pairs);
        
    if current_features.iter().any(|f| f.fingerprint == hash) {
        // If the feature is already present, then return the hash,
        // Don't add it again
        return hash;
    }

    if bondset.bonds.is_empty() {
        // If the bondset is empty, then this is the first features
        current_features.push(ECFPFeature {
            fingerprint: hash,
            bondset,
        });
        
        return hash;
    }

    let mut to_remove = false;
    let mut j = 0;
    // Check if the bondset is unique, etc.
    for (i,feature) in features.iter().enumerate() {
        if feature.bondset.bonds.is_empty() {
            continue;
        }

        if feature.bondset.bonds == bonds {
            // If any bondsets overlap, structurally the same.
            //println!("Found a duplicate bondset: {:?} {:?}", feature.bondset, bondset);
            if feature.bondset.depth == depth {
                if feature.fingerprint > hash {
                    to_remove = true;
                    j = i;
                    break;
                } else {
                    if feature.fingerprint == hash {
                        panic!(
                            "The same fingerprint was added twice! {:?} {:?}",
                            feature.fingerprint, hash
                        )
                    }
                    // The feature is already present, but with a lower hash, retain the original, don't add this fp
                    return hash;
                }
            } else if feature.bondset.depth > depth {
                // The feature is already present, but at a higher depth, so remove it
                to_remove = true;
                j = i;
                break;
            } else {
                // feature.bondset.depth < bondset.depth
                // The feature is already present, but at a lower depth, retain the original, don't add this fp
                return hash;
            }
        }
    }

    if to_remove {
        features.swap_remove(j); // 
        
    }
    // Otherwise, the bondset is unique, so add it to the list of features
    current_features.push(ECFPFeature {
        fingerprint: hash,
        bondset,
    });

    hash
}

struct BondId(usize);

impl BondId {
    fn new(sid: usize, tid: usize) -> Self {
        // Add the two ids together as strings, independent of order
        let mut bond_id = String::new();
        if sid < tid {
            bond_id.push_str(&sid.to_string());
            bond_id.push_str(&tid.to_string());
        } else {
            bond_id.push_str(&tid.to_string());
            bond_id.push_str(&sid.to_string());
        }
        BondId(bond_id.parse().unwrap())
    }
}

pub fn mol_to_features(molecule: DefaultMolecule, radius: u8) -> HashSet<HashBaseType> {
    let mut features = Vec::new();
    // Initialise the atoms identifiers.
    let mut molecular_identifiers: Vec<HashBaseType> = invariant_identifier(&molecule);
    let mut next_identifiers: Vec<HashBaseType> = molecular_identifiers.clone();
    let mut current_iteration: u8 = 0;

    //: let sid:usize; // Source ID
    //: let tid:usize; // Target ID
    while current_iteration <= radius {
        let mut current_features: Vec<ECFPFeature> = Vec::new();
        for sid in molecule.ids() {
            let mut pairs = Vec::new();

            if current_iteration == 0 {
                // Add root pair
                let atom_identifier = molecular_identifiers[sid];
                pairs.push(Pair {
                    atom: atom_identifier,
                    bond: 0,
                });
                molecular_identifiers[sid] = add_feature(
                    &mut features,
                    &mut current_features,
                    &pairs,
                    BondSet {
                        bonds: HashSet::<u8>::new(),
                        depth: 0,
                    },
                );
                
            } else {
                // Else

                let neighbours = molecule.neighbors(sid).unwrap();
                let mut current_bonds = HashSet::<u8>::new();
                for tid in neighbours {
                    let bond = molecule.get_bond(sid, tid).unwrap();
                    let bond_type = bond.order() as u8;
                    let atom_identifier = molecular_identifiers[tid];

                    pairs.push(Pair {
                        atom: atom_identifier,
                        bond: bond_type,
                    });

                    current_bonds.insert(BondId::new(sid, tid).0 as u8);
                }
                
                
                // Add feature, and update the identifier
                next_identifiers[sid] = add_feature(
                    &mut features,
                    &mut current_features,
                    &pairs,
                    BondSet {
                        bonds: current_bonds,
                        depth: current_iteration,
                    },
                );
            }
        }
        
        features.append(&mut current_features);
        molecular_identifiers = next_identifiers.clone();
        
        current_iteration += 1;

        
    }

    features.iter().map(|feature| feature.fingerprint).collect()
}

// A number of functions that are used to calculate information
// about atoms. 

use chemcore::molecule::Molecule;

use crate::utils::hash::{hash_64, HashBaseType};

pub fn invariant_identifier(molecule: &dyn Molecule) -> Vec<HashBaseType> {
    // This function is used to set the identifier of the atom
    // This is based on the Daylight atomic invariants.
    // "...are six properties of an atom in a molecule that do not depend on initial atom numbering." (Rogers, 2010)
    // : 1. Immediate Heavy Atom Neighbours, 2. Valence - No. Hydrogens, 
    // : 3. Atomic Number, 4. Atomic Mass, 
    // : 5. Atomic Charge, 6. Number of Hydrogens (Both Implicit and Explicit)

    // NOTE: "We include one additional property: whether the atom is contained in at least one ring." (Rogers, 2010)
    // : 7. Ring Membership

    let mut identifiers = Vec::new();

    for i in molecule.ids() {
        let atom = molecule.atom(i).unwrap().clone();
        let element = atom.element.unwrap();

        let n_hydrogens = atom.hydrogens;
        let valence = element.valence_electrons();
        let atomic_number = element.atomic_number();

        let n_neigbours = molecule.neighbors(i).unwrap().count();

        let charge = atom.electrons;

        let hash = hash_64(&vec![
            n_neigbours as i64,
            valence as i64 - n_hydrogens as i64,
            atomic_number as i64,
            charge as i64,
            n_hydrogens as i64,
        ]);
        
        identifiers.push(hash as HashBaseType);
    }
    identifiers
}
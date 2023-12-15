use std::collections::HashSet;

use bitvec::prelude::*;

use super::hash::HashBaseType;

pub fn hash_to_bitvect(hash: HashSet<HashBaseType>, size: usize) -> BitVec {
    let mut bitvect = BitVec::repeat(false, size);
    for hash in hash {
        // Fold the hash to the size of the bitvect
        let hash = hash % size as HashBaseType;
        
        bitvect.set(hash as usize, true);
    }
    bitvect
}

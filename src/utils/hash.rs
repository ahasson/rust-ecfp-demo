use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

pub type HashBaseType = u64;

pub fn hash_64<T>(obj: T) -> HashBaseType
where
    T: Hash,
{
    let mut hasher = DefaultHasher::new();
    obj.hash(&mut hasher);
    hasher.finish() as HashBaseType
}
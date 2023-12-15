use super::Parity;

#[derive(Debug, PartialEq)]
pub struct Bond {
    pub electrons: u8,
    pub parity: Option<Parity>,
    pub tid: usize,
    pub idx: usize,
}

impl Bond {
    pub fn new(electrons: u8, parity: Option<Parity>, tid: usize, idx: usize) -> Self {
        Self {
            electrons,
            parity,
            tid,
            idx,
        }
    }

    pub fn order(&self) -> f32 {
        self.electrons as f32 / 2 as f32
    }
}

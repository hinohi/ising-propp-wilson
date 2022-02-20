use rand::Rng;
use rand_pcg::Mcg128Xsl64;

fn gibbs_sampler(beta: f64, e0: f64, e1: f64) -> f64 {
    let p0 = (-e0 * beta).exp();
    let p1 = (-e1 * beta).exp();
    p1 / (p0 + p1)
}

pub struct ProppWilsonIsing {
    n2: usize,
    local_prob2: f64,
    local_prob4: f64,
    neighbors: Vec<[usize; 4]>,
    /// Manage Rng
    pointer: usize,
    /// Generated Random number
    /// lower 29-bits are index, higher 3-bits are cell-flip.
    data: Vec<u32>,
    cell_up_start: Vec<u8>,
    cell_down_start: Vec<u8>,
}

const CELL_INDEX_BIT: u32 = 29;

impl ProppWilsonIsing {
    pub fn new(n: usize, beta: f64) -> ProppWilsonIsing {
        assert!(n * n < 2usize.pow(CELL_INDEX_BIT));
        let mut neighbors = Vec::with_capacity(n * n);
        for x in 0..n {
            for y in 0..n {
                neighbors.push([
                    (y + n - 1) % n * n + x,
                    y * n + (x + n - 1) % n,
                    y * n + (x + 1) % n,
                    (y + 1) % n * n + x,
                ]);
            }
        }
        ProppWilsonIsing {
            n2: n * n,
            local_prob2: gibbs_sampler(beta, -2.0, 2.0),
            local_prob4: gibbs_sampler(beta, -4.0, 4.0),
            neighbors,
            pointer: 0,
            data: Vec::new(),
            cell_up_start: vec![1; n * n],
            cell_down_start: vec![0; n * n],
        }
    }

    pub fn rerun<R: Rng>(&mut self, n: usize, rng: &mut R) {
        assert_eq!(n % self.n2, 0);
        for _ in 0..n {
            let i = rng.gen_range(0..self.n2) as u32;
            let c = match rng.gen::<f64>() {
                p if p < self.local_prob4 => 5,
                p if p < self.local_prob2 => 4,
                p if p < 0.5 => 3,
                p if p < 1.0 - self.local_prob2 => 2,
                p if p < 1.0 - self.local_prob4 => 1,
                _ => 0,
            };
            self.data.push(i + (c << CELL_INDEX_BIT));
        }
        self.pointer = self.data.len();
        self.cell_up_start = vec![1; self.n2];
        self.cell_down_start = vec![0; self.n2];
    }

    pub fn update_one(&mut self) {
        debug_assert!(self.pointer > 0);
        self.pointer -= 1;
        let i = self.data[self.pointer];
        let index = (i & ((1 << CELL_INDEX_BIT) - 1)) as usize;
        let threshold = (i >> CELL_INDEX_BIT) as u8;
        match threshold {
            0 => {
                self.cell_up_start[index] = 1;
                self.cell_down_start[index] = 1;
                return;
            }
            5 => {
                self.cell_up_start[index] = 0;
                self.cell_down_start[index] = 0;
                return;
            }
            _ => (),
        }
        let mut m_up_start = 0;
        let mut m_down_start = 0;
        for j in self.neighbors[index] {
            m_up_start += self.cell_up_start[j];
            m_down_start += self.cell_down_start[j];
        }
        self.cell_up_start[index] = if threshold <= m_up_start { 1 } else { 0 };
        self.cell_down_start[index] = if threshold <= m_down_start { 1 } else { 0 };
    }

    pub fn update_n(&mut self) -> bool {
        assert!(self.pointer >= self.n2);
        for _ in 0..self.n2 {
            self.update_one();
        }
        self.pointer > 0
    }

    pub fn distance(&self) -> u32 {
        let mut s = 0;
        for (&c0, &c1) in self.cell_up_start.iter().zip(self.cell_down_start.iter()) {
            if c0 != c1 {
                s += 1;
            }
        }
        s
    }
}

fn main() {
    let mut rng = Mcg128Xsl64::new(1);
    let mut ising = ProppWilsonIsing::new(20, 1.0 / 2.1);
    let mut n = ising.n2;
    'OUT: loop {
        ising.rerun(n, &mut rng);
        let mut i = 0;
        while ising.update_n() {
            i += 1;
            let dist = ising.distance();
            if dist == 0 {
                println!("{}", i);
                break 'OUT;
            }
        }
        n *= 2;
    }
}

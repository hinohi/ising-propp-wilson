use rand::Rng;

const CELL_INDEX_BIT: u32 = 29;

fn gibbs_sampler(beta: f64, e0: f64, e1: f64) -> f64 {
    let p0 = (-e0 * beta).exp();
    let p1 = (-e1 * beta).exp();
    p1 / (p0 + p1)
}

pub struct Ising {
    n2: usize,
    local_prob2: f64,
    local_prob4: f64,
    neighbors: Vec<[usize; 4]>,
    /// Manage Rng
    pointer: usize,
    /// Generated Random number
    /// lower 29-bits are index, higher 3-bits are cell-flip.
    generated_random_number: Vec<u32>,
    cell_up_start: Vec<u8>,
    cell_down_start: Vec<u8>,
}

impl Ising {
    pub fn new(n: usize, beta: f64) -> Ising {
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
        Ising {
            n2: n * n,
            local_prob2: gibbs_sampler(beta, -2.0, 2.0),
            local_prob4: gibbs_sampler(beta, -4.0, 4.0),
            neighbors,
            pointer: 0,
            generated_random_number: Vec::new(),
            cell_up_start: vec![1; n * n],
            cell_down_start: vec![0; n * n],
        }
    }

    fn setup_next<R: Rng>(&mut self, rng: &mut R) {
        let n = if self.generated_random_number.is_empty() {
            1
        } else {
            // n *= 2;
            self.generated_random_number.len()
        };
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
            self.generated_random_number.push(i + (c << CELL_INDEX_BIT));
        }
        self.pointer = n;
        self.cell_up_start = vec![1; self.n2];
        self.cell_down_start = vec![0; self.n2];
    }

    fn update_one(&mut self) {
        self.pointer -= 1;
        let i = self.generated_random_number[self.pointer];
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

    pub fn rerun<R: Rng>(&mut self, rng: &mut R) -> bool {
        self.setup_next(rng);
        while self.pointer > 0 {
            self.update_one();
        }
        self.cell_up_start == self.cell_down_start
    }

    pub fn magnetization(&self) -> i64 {
        let mut m = 0;
        for &c in self.cell_up_start.iter() {
            m += if c == 0 { -1 } else { 1 };
        }
        m
    }

    pub fn energy(&self) -> i64 {
        let mut e = 0;
        for (&c, nei) in self.cell_up_start.iter().zip(self.neighbors.iter()) {
            let mut local_m = 0;
            for &j in nei {
                local_m += self.cell_up_start[j];
            }
            e += (c as i64 * 2 - 1) * (local_m as i64 - 2);
        }
        e
    }
}

pub fn run<R: Rng>(rng: &mut R, n: usize, temperature: f64, limit: u8) -> Option<(u8, i64, i64)> {
    let mut ising = Ising::new(n, 1.0 / temperature);
    let mut loop_count = 0;
    let done = loop {
        if loop_count == limit {
            break false;
        }
        if ising.rerun(rng) {
            break true;
        }
        loop_count += 1;
    };
    if done {
        Some((loop_count, ising.magnetization(), ising.energy()))
    } else {
        None
    }
}

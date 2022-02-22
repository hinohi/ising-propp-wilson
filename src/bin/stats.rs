use std::{fs::read_to_string, path::PathBuf};

use clap::Parser;

#[derive(Debug, Parser)]
struct Args {
    data: PathBuf,
    #[clap(short)]
    n: Option<u32>,
}

#[derive(Debug, Default, Copy, Clone)]
struct S {
    n: f64,
    s1: f64,
    s2: f64,
    s4: f64,
}

impl S {
    fn add(&mut self, s: f64) {
        self.n += 1.0;
        self.s1 += s;
        let s2 = s * s;
        self.s2 += s2;
        self.s4 += s2 * s2;
    }

    fn mean(&self) -> f64 {
        self.s1 / self.n
    }

    fn sig2(&self) -> f64 {
        self.s2 / (self.n - 1.0) - self.s1 * self.s1 / self.n / (self.n - 1.0)
    }

    fn binder(&self) -> f64 {
        1.0 - self.s4 / self.s2 / self.s2 / 3.0 * self.n
    }
}

fn main() {
    let args: Args = Args::parse();
    let n = if let Some(n) = args.n {
        n
    } else {
        let name = args.data.file_name().unwrap().to_string_lossy();
        let name = name.splitn(2, '.').next().unwrap();
        name.parse().unwrap()
    };
    let n = (n * n) as f64;
    let mut pre_t = None;
    let mut s_tau = S::default();
    let mut s_m = S::default();
    let mut s_e = S::default();
    for line in read_to_string(args.data).unwrap().lines() {
        let mut words = line.split(' ');
        let (t, dt, tau, m, e) = (
            words.next().unwrap().parse::<f64>().unwrap(),
            words.next().unwrap().parse::<f64>().unwrap(),
            words.next().unwrap().parse::<f64>().unwrap(),
            words.next().unwrap().parse::<f64>().unwrap(),
            words.next().unwrap().parse::<f64>().unwrap(),
        );
        if pre_t.is_none() {
            pre_t = Some(t);
        }
        if pre_t != Some(t) {
            println!(
                "{} {} {} {} {} {} {} {} {} {}",
                t,
                dt,
                s_m.n,
                s_m.mean(),
                s_m.sig2().sqrt(),
                s_m.binder(),
                s_e.mean(),
                s_e.sig2().sqrt(),
                s_tau.mean(),
                s_tau.sig2().sqrt(),
            );
            pre_t = Some(t);
            s_tau = S::default();
            s_m = S::default();
            s_e = S::default();
        }
        s_m.add(m.abs() / n);
        s_e.add(e / n);
        s_tau.add(tau);
    }
}

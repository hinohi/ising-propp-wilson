use clap::Parser;
use rand_pcg::Mcg128Xsl64;

use ising_propp_wilson::run;

#[derive(Debug, Parser)]
struct Args {
    /// The size of one side of a two-dimensional square lattice
    #[clap(short)]
    n: usize,
    /// Propp-Wilson algorithm iteration limit
    #[clap(long, default_value_t = 30)]
    limit: u8,
    #[clap(long, short, default_value_t = 100)]
    samples: usize,
    #[clap(long, default_value_t = 1)]
    seed: u128,
}

fn main() {
    let tc = 2.0 / (1.0 + std::f64::consts::SQRT_2).ln();
    let sample_targets = {
        let mut t = Vec::new();
        // 9,8,・・・,1
        for i in (1..=9).rev() {
            t.push(i as f64);
        }
        // 0.9,0.8,・・・,0.1
        for i in (1..=9).rev() {
            t.push(i as f64 / 10.0);
        }
        // 0.099,0.098,・・・,0.01,0,-0.01,-0.02,・・・,-0.99,-1,・・・
        for i in (-2000..=99).rev() {
            t.push(i as f64 / 1000.0);
        }
        t
    };
    let args: Args = Args::parse();
    let mut rng = Mcg128Xsl64::new(args.seed * 2);
    for dt in sample_targets {
        let t = tc + dt;
        let mut ok = 0;
        for _ in 0..args.samples {
            if let Some((tau, m, e)) = run(&mut rng, args.n, t, args.limit) {
                ok += 1;
                println!("{} {} {} {} {}", t, dt, tau, m, e);
            } else {
                eprintln!("NG: t={} dt={}", t, dt);
            }
        }
        eprintln!("stats: t={} dt={} ok={}", t, dt, ok);
        if ok * 2 < args.samples {
            break;
        }
    }
}

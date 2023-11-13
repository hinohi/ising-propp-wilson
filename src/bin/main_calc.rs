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
    let args: Args = Args::parse();
    let mut rng = Mcg128Xsl64::new(args.seed * 2);
    for dt in (-100..=100).rev() {
        let dt = dt as f64 / args.n as f64 / 4.0;
        let t = tc + dt;
        if t <= 0.0 {
            break;
        }
        let mut ok = 0;
        for _ in 0..args.samples {
            if let Some((c, ising)) = run(&mut rng, args.n, t, args.limit) {
                ok += 1;
                println!(
                    "{} {} {} {} {}",
                    t,
                    dt,
                    c,
                    ising.magnetization(),
                    ising.energy()
                );
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

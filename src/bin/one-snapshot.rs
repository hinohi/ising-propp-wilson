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
    #[clap(long, default_value_t = 1)]
    seed: u128,
    #[clap(long, default_value_t = 0.0)]
    dt: f64,
}

fn main() {
    let tc = 2.0 / (1.0 + std::f64::consts::SQRT_2).ln();
    let args: Args = Args::parse();
    let mut rng = Mcg128Xsl64::new(args.seed * 2);
    let t = tc + args.dt;
    eprintln!("t={}", t);
    if let Some((c, ising)) = run(&mut rng, args.n, t, args.limit) {
        eprintln!(
            "loop_count={} M={} E={}",
            c,
            ising.magnetization(),
            ising.energy(),
        );
        println!("{}", ising.spin_snapshot());
    }
}

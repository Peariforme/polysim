use bigsmiles::parse;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use polysim_core::builder::{linear::LinearBuilder, BuildStrategy};

fn bench_polyethylene(c: &mut Criterion) {
    let mut group = c.benchmark_group("homopolymer/polyethylene");

    for n in [10usize, 100, 1_000, 10_000] {
        let bs = parse("{[]CC[]}").unwrap();
        let builder = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n));
        // Throughput = nombre d'unités de répétition générées par seconde
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &builder, |b, builder| {
            b.iter(|| builder.homopolymer().unwrap());
        });
    }
    group.finish();
}

fn bench_polystyrene(c: &mut Criterion) {
    // Polystyrène : ring renumbering activé (ring 1 par unité)
    let mut group = c.benchmark_group("homopolymer/polystyrene");

    for n in [10usize, 100, 1_000] {
        let bs = parse("{[]CC(c1ccccc1)[]}").unwrap();
        let builder = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n));
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &builder, |b, builder| {
            b.iter(|| builder.homopolymer().unwrap());
        });
    }
    group.finish();
}

criterion_group!(benches, bench_polyethylene, bench_polystyrene);
criterion_main!(benches);

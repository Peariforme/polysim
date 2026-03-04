use bigsmiles::parse;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::mechanical::density,
};

fn build_chain(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).unwrap();
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .unwrap()
}

fn bench_density(c: &mut Criterion) {
    // PS n=100 : chaine avec phenyle, anneau aromatique --> plus couteux
    let mut group = c.benchmark_group("density/ps");

    for n in [10usize, 100, 1_000] {
        let chain = build_chain("{[]CC(c1ccccc1)[]}", n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &chain, |b, chain| {
            b.iter(|| density(chain).unwrap());
        });
    }
    group.finish();
}

criterion_group!(benches, bench_density);
criterion_main!(benches);

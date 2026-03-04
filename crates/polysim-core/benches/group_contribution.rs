use bigsmiles::parse;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::group_contribution::GroupDatabase,
};

fn build_chain(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).unwrap();
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .unwrap()
}

fn bench_group_decomposition_pe(c: &mut Criterion) {
    let mut group = c.benchmark_group("group_contribution/pe");

    for n in [10usize, 100, 1_000] {
        let chain = build_chain("{[]CC[]}", n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &chain, |b, chain| {
            b.iter(|| GroupDatabase::decompose(chain).unwrap());
        });
    }
    group.finish();
}

fn bench_group_decomposition_ps(c: &mut Criterion) {
    // PS : ring renumbering + detection des phényles → plus coûteux que PE
    let mut group = c.benchmark_group("group_contribution/ps");

    for n in [10usize, 100, 1_000] {
        let chain = build_chain("{[]CC(c1ccccc1)[]}", n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &chain, |b, chain| {
            b.iter(|| GroupDatabase::decompose(chain).unwrap());
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_group_decomposition_pe,
    bench_group_decomposition_ps
);
criterion_main!(benches);

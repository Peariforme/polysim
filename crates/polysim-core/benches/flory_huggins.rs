use bigsmiles::parse;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::solubility::flory_huggins_chi,
};

fn build_chain(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).unwrap();
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .unwrap()
}

fn bench_chi_pe_pvc(c: &mut Criterion) {
    let mut group = c.benchmark_group("flory_huggins/pe_pvc");

    for n in [10usize, 100, 1_000] {
        let pe = build_chain("{[]CC[]}", n);
        let pvc = build_chain("{[]C(Cl)C[]}", n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(n),
            &(&pe, &pvc),
            |b, (pe, pvc)| {
                b.iter(|| flory_huggins_chi(pe, pvc, 298.0).unwrap());
            },
        );
    }
    group.finish();
}

fn bench_chi_pe_ps(c: &mut Criterion) {
    // PS : décomposition phényle plus coûteuse → mesure l'impact sur chi
    let mut group = c.benchmark_group("flory_huggins/pe_ps");

    for n in [10usize, 100, 1_000] {
        let pe = build_chain("{[]CC[]}", n);
        let ps = build_chain("{[]CC(c1ccccc1)[]}", n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(n),
            &(&pe, &ps),
            |b, (pe, ps)| {
                b.iter(|| flory_huggins_chi(pe, ps, 298.0).unwrap());
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_chi_pe_pvc, bench_chi_pe_ps);
criterion_main!(benches);

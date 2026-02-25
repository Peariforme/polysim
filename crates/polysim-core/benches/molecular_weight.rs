use bigsmiles::parse;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::molecular_weight::{average_mass, monoisotopic_mass},
};

fn bench_average_mass(c: &mut Criterion) {
    let mut group = c.benchmark_group("molecular_weight/average_mass");

    for n in [10usize, 100, 1_000] {
        // Pré-construire la chaîne — on benche uniquement le calcul de masse
        let bs = parse("{[]CC[]}").unwrap();
        let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
            .homopolymer()
            .unwrap();
        // Throughput = atomes lourds dans la chaîne (2n carbons)
        group.throughput(Throughput::Elements(2 * n as u64));
        group.bench_with_input(BenchmarkId::new("polyethylene", n), &chain, |b, chain| {
            b.iter(|| average_mass(chain));
        });
    }

    // Polystyrène : atomes aromatiques, plus complexe à parser
    for n in [10usize, 100] {
        let bs = parse("{[]CC(c1ccccc1)[]}").unwrap();
        let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
            .homopolymer()
            .unwrap();
        // 8 atomes lourds par unité (2C aliphatique + 6C aromatique)
        group.throughput(Throughput::Elements(8 * n as u64));
        group.bench_with_input(BenchmarkId::new("polystyrene", n), &chain, |b, chain| {
            b.iter(|| average_mass(chain));
        });
    }

    group.finish();
}

fn bench_monoisotopic_mass(c: &mut Criterion) {
    let mut group = c.benchmark_group("molecular_weight/monoisotopic_mass");

    for n in [10usize, 100, 1_000] {
        let bs = parse("{[]CC[]}").unwrap();
        let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
            .homopolymer()
            .unwrap();
        group.throughput(Throughput::Elements(2 * n as u64));
        group.bench_with_input(BenchmarkId::new("polyethylene", n), &chain, |b, chain| {
            b.iter(|| monoisotopic_mass(chain));
        });
    }

    group.finish();
}

fn bench_by_target_mn(c: &mut Criterion) {
    let mut group = c.benchmark_group("molecular_weight/by_target_mn");

    // Bench complet : résolution de n + construction de la chaîne + calcul MW
    for target in [282.554f64, 2825.54, 28255.4] {
        let bs = parse("{[]CC[]}").unwrap();
        group.bench_with_input(
            BenchmarkId::new("polyethylene", target as usize),
            &(bs, target),
            |b, (bs, target)| {
                b.iter(|| {
                    let bs2 = parse("{[]CC[]}").unwrap();
                    let _ = bs2;
                    LinearBuilder::new(bs.clone(), BuildStrategy::ByTargetMn(*target))
                        .homopolymer()
                        .unwrap()
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_average_mass,
    bench_monoisotopic_mass,
    bench_by_target_mn
);
criterion_main!(benches);

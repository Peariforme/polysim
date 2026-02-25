.PHONY: install-hooks test fmt lint bench

## Install git hooks from .githooks/ into the local repository.
## Run this once after cloning: make install-hooks
install-hooks:
	git config core.hooksPath .githooks
	chmod +x .githooks/pre-commit .githooks/pre-push
	@echo "Git hooks installed (pre-commit: fmt+clippy, pre-push: test)."

## Run the full test suite.
test:
	cargo test --all-features

## Check and apply code formatting.
fmt:
	cargo fmt --all

## Run clippy with all features, treating warnings as errors.
lint:
	cargo clippy --all-targets --all-features -- -D warnings

## Compile benchmarks without running them.
bench:
	cargo bench --no-run --all-features

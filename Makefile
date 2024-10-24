BENCHMARKS_DIR := benchmarks
BENCHMARK_SCRIPTS := $(shell find $(BENCHMARKS_DIR) -name "*.jl")

build:
	julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.resolve()'

doc:
	make -C docs/

test:
	julia --project=. -e 'using Pkg; Pkg.test()'

clean:
	rm -rf docs/build/

preview:
	julia -e 'using LiveServer; serve(dir="docs/build")'

benchmark:
	@echo "Running all benchmark scripts in $(BENCHMARKS_DIR)..."
	@for benchmark_script in $(BENCHMARK_SCRIPTS); do \
		echo "==============================="; \
		echo "$$benchmark_script" ; \
		echo "-------------------------------"; \
	  julia --project=. $$benchmark_script; \
		echo "-------------------------------"; \
		echo; \
	done

.PHONY: build doc test clean preview benchmark

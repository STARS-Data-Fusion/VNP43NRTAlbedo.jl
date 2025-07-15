# Makefile for VNP43NRTAlbedo.jl

JULIA ?= julia

.PHONY: test

test:
	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

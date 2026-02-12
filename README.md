# PositiveChartsToricVarieties


PositiveChartsToricVarieties is a Julia package to compute positive charts of toric varieties. A detailed documentation is available at the Zenodo page [Zenodo page PositiveChartsToricVarieties](10.5281/zenodo.18613405).

## Installation

```julia
using Pkg
Pkg.add("PositiveChartsToricVarieties")
```
PositiveChartsToricVarieties is supported on Julia 1.12 and later.

## Example

The following is a minimum working example.

```julia
using PositiveChartsToricVarieties
using Oscar 

P = cube(3)

nefCube, F = nef_cone_modulo_lineality(P)
M, F = unimod_matrix(P)
f, M, F = unimod_nef_polynomials(P)
J = Y_variety(P)
u = parametrization_Y(P)
test_sat_conjecture(P)

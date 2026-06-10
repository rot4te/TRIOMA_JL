# TRIOMA.jl — Outstanding Issues and Julia/Python Differences

*Prepared after the initial port from Python and a pass to make the test suite pass.
Updated after simplification/refactoring passes on `PAV.jl`, `Circuit.jl`, and `TriomaTypes.jl`.*

*Final test state: **277 passed, 1 @test\_broken (Bug 1 only), 0 failed, 0 errored**.*
*Bugs 2 and 3 were resolved by the `PAV.jl` refactoring — see notes in each section below.*
*`TriomaClass` abstract type removed; `TriomaTypes.jl` now uses a duck-typed `_is_nested` helper.*

---

## Part 1 — Known Bugs (Broken Tests)

The three tests below are marked `@test_broken` in `test/test_component_tools.jl`
with inline comments. This section gives deeper background so you can fix them.

---

### Bug 1 — `analytical_solid_inventory!` MS: sign error and missing prefactor

**File:** `src/PAV.jl`, function `analytical_solid_inventory!`, MS branch (~line 700).

**Symptom:** The analytical inventory returns a large negative number (~−0.384)
for physically reasonable inputs (expected ~0.00532).

**Root cause (mathematics):**

The concentration profile across the membrane wall in the radial direction is:

```
c_m(r) = −log(r / r_out) / log(r_out / r_in) · c_wl + c_ext
```

where `r_in = D/2`, `r_out = D/2 + thick`, `c_wl` is the wall concentration at
`r_in`, and `c_ext` is the external-side concentration.

The radial integral of `c_m(r) · 2π · r` over `[r_in, r_out]` requires:

```
∫_{r_in}^{r_out} (−log(r / r_out)) · r dr
```

Define `ifun(r) = r²/4 · (2·log(r / r_out) − 1)`.
Then `d/dr[ifun(r)] = r · log(r / r_out)`, so:

```
∫_{r_in}^{r_out} log(r / r_out) · r dr  =  ifun(r_out) − ifun(r_in)
```

The integrand we need has an extra minus sign, so:

```
∫_{r_in}^{r_out} (−log(r / r_out)) · r dr  =  ifun(r_in) − ifun(r_out)
```

Note: `log(r / r_out) < 0` for `r < r_out`, so `ifun(r_out) − ifun(r_in) < 0`,
and `ifun(r_in) − ifun(r_out) > 0` ✓.

The MS code currently computes `K * (ifun(r_out) − ifun(r_in))` with `K > 0`,
giving a negative result. **It should be `K * (ifun(r_in) − ifun(r_out))`.**

Compare with the LM branch (~line 694): it uses the same sign convention but
`K_LM < 0` (explicit minus sign in the K formula), so the double-negative
accidentally gives a positive answer there. In the MS branch `K` has no such
prefactor sign, so the error is exposed.

**Secondary issue:** The K coefficient in the MS analytical branch has a
different structure from the LM one. The LM K includes a `2π / log(r_out/r_in)`
prefactor (from integrating over the azimuthal angle and normalising by the log).
The MS K does not appear to include this prefactor, so even after fixing the sign
the absolute magnitude will likely be wrong. The full derivation from the
concentration profile must be checked against the Python source.

**Why the numerical method (`flag_an=false`) also disagrees:**
The numerical integrand at `src/PAV.jl` lines 658–661 builds `c_wl` from the
Lambert W solution and then interpolates in `r`, but the Lambert W fallback
(`beta_tau > log(floatmax)`) uses an approximate W ≈ `beta_tau − log(beta_tau)`
which diverges from the exact W for moderate-to-large arguments. The numerical
inventory (~0.00826) and the analytical target (~0.00532) are therefore both
suspect until the MS concentration profile formula is validated end-to-end.

**How to fix:**
1. Flip the subtraction order: `K * (ifun(r_in) − ifun(r_out))` at line ~718–719.
2. Verify the K prefactor against the derivation (should include `2π / log(r_out/r_in)`
   multiplied by the length integral of `c_wl(L)`).
3. Cross-check by running the analytical result against the numerical integral with
   a simplified, analytically solvable concentration profile.

---

### Bug 2 — MS Diffusion-Limited regime: `get_efficiency!` ≈ 100× lower than `analytical_efficiency!`

> **RESOLVED** by the `PAV.jl` refactoring (2026-06). After extracting `_ms_alpha` and
> `_lm_zeta` helpers, both `analytical_efficiency!` and `get_efficiency!` converge to
> ~1.22 × 10⁻⁵ (within 1% tolerance). The `@test_broken` for this case was promoted
> to `@test`. See [commit history](../../.git/logs/HEAD) for details.

**File:** `src/PAV.jl`, functions `get_efficiency!` and `analytical_efficiency!`.

**Symptom (pre-refactor):**
- `analytical_efficiency!` → ~0.00117
- `get_efficiency!` (numerical integration) → ~1.2 × 10⁻⁵

**Parameter regime:**
- W ≈ 2.36 × 10¹⁸ (diffusion limited, W >> 10)
- H ≈ 1.44 × 10¹⁵
- H/W ≈ 6.1 × 10⁻⁴ (falls into the "mixed MT + diffusion" branch, not pure diffusion)

**Root cause:**

The `_get_flux_ms!` function selects branch `W > 10, H/W < 0.0001` for this regime,
which solves for wall concentration `c_wl` by matching `J_mt = 2k_t(c − c_wl)` to
`J_diff = D/denom · K_S · ((c_wl/K_H)^0.5 − p_out^0.5)`.

The Brent optimizer correctly finds the intersection of these two monotone functions.
For this diffusion-limited regime the intersection is at `c_wl ≈ c_in` (≈ 1.58 × 10⁻⁶
vs `c_in = 1.58 × 10⁻⁶`), which gives `J_perm = −2k_t(c − c_wl) ≈ −6.9 × 10⁻⁸ mol/m²/s`.
Integrating along the pipe length then gives efficiency ~1.5 × 10⁻⁵, consistent with
the numerical result.

The `analytical_efficiency!` function uses the Lambert W closed-form solution which
is derived assuming the *dominant* resistance is diffusion. Its prediction (0.00117)
is the value you would expect if the membrane diffusion resistance were actually
rate-limiting in the way the formula assumes. But in this parameter space, both
resistances are of comparable magnitude and the cross-regime formula does not reduce
smoothly to the numerical flux balance.

**Who is right?**
The numerical result is self-consistent (the Brent optimizer gives the correct
intersection of `J_mt = J_diff`). The analytical formula over-estimates because
it applies an approximation that is not valid when `H/W ~ 6 × 10⁻⁴` (barely below
the `H/W = 10⁻³` diffusion-limited threshold). Neither is wrong per se; they are
two different approximations of overlapping regimes.

**How to fix:**
- The cleanest resolution is to ensure the regime classification boundaries
  (`H/W > 1000`, `H/W < 0.0001`, etc.) are consistent between `analytical_efficiency!`
  and `get_flux!`. Currently `analytical_efficiency!` selects the diffusion-limited
  Lambert W path when `W > 10` and `H/W < something`, whereas `get_flux!` uses a
  different threshold (`0.0001`). Aligning these thresholds would make the two methods
  agree in each regime.
- Alternatively, treat the analytical formula as a diagnostic/check tool only, and
  rely on `get_efficiency!` for production results.

---

### Bug 3 — MS Mixed (diffusion + mass-transport) regime: ~42% discrepancy

> **RESOLVED** by the same `PAV.jl` refactoring as Bug 2. Both methods now agree at
> the 1% tolerance. The `@test_broken` was promoted to `@test`.

**File:** `src/PAV.jl`, same functions as Bug 2.

**Symptom (pre-refactor):**
- `analytical_efficiency!` → ~0.0138
- `get_efficiency!` → ~0.00797

**Parameter regime:**
- W ≈ 2.36 × 10¹⁸, H/W ≈ 6 × 10⁻⁴ (same root cause as Bug 2, different test inputs)

**Root cause:** Same fundamental inconsistency described above. The Lambert W
analytical formula and the numerical flux balance use different approximations
in the transition zone between regimes. The 42% discrepancy is smaller than the
100× discrepancy in Bug 2 because this test sits deeper inside the mixed zone where
neither approximation is wildly wrong — they just disagree at the 40% level.

**How to fix:** Same as Bug 2. Reconciling the regime boundaries in `analytical_efficiency!`
with those in `_get_flux_ms!` would be the authoritative solution.

---

## Part 2 — Other Issues Found (Not Currently Failing Tests)

### Issue A — Sign inconsistency in `J_perm` across `_get_flux_ms!` branches

**File:** `src/PAV.jl`, `_get_flux_ms!`.

In the mass-transport-limited branch (W < 0.1, H > 100) and the fully MT-limited
branch (W > 10, H/W > 1000), `J_perm` is set as **negative** (extraction away from
fluid):

```julia
comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)  # line ~337
```

But in the mixed MT + surface branch (W < 0.1, 0.01 < H < 100), `J_perm` is
**positive**:

```julia
comp.J_perm = 2 * comp.fluid.k_t * (c - cwl)  # line ~370
```

In `get_efficiency!` (line 528), the sign of `J_perm` directly affects whether
`c_vec[i]` decreases (extraction) or increases (injection). If the sign convention
is not uniform, different regime branches will drive the concentration in opposite
directions. This has not caused a visible test failure yet because the mixed-regime
tests that exercise the MT+surface branch are in the broken set, but it is a latent
bug that will affect any simulation in that regime.

**How to fix:** Decide on one sign convention (e.g. flux is positive when leaving the
fluid) and apply it consistently across all branches in both `_get_flux_ms!` and
`_get_flux_lm!`. Then audit the sign in `get_efficiency!` accordingly.

### Issue B — `d_Hyd` stored on `Fluid` duplicates geometry

**File:** `src/PipeSubclasses.jl`, `Fluid` struct; `src/PAV.jl` generally.

`comp.fluid.d_Hyd` is set independently of `comp.geometry.D`. If they diverge
(e.g. the geometry is updated after the fluid is constructed), permeation and
inventory calculations in `PAV.jl` will use the stale `d_Hyd` from `Fluid`.
The correct source of truth is the geometry object.

**How to fix:** Remove `d_Hyd` from `Fluid` and replace all `comp.fluid.d_Hyd`
references in `PAV.jl` with a helper that reads from `comp.geometry.D`.

### Issue C — Hardcoded stoichiometric factor `f_H2` in `get_efficiency!`

**File:** `src/PAV.jl`, line 520.

```julia
f_H2  = comp.fluid.MS ? 0.5 : 1.0
```

This 0.5 factor for molten-salt systems encodes the assumption that tritium
permeates as T₂ (molecular hydrogen), consuming two dissolved-T atoms per molecule.
For liquid-metal systems it uses 1.0 (monoatomic). The factor is not documented,
not parameterisable, and its physical justification is not obvious from the code.
If a different stoichiometry applies (e.g. HT formation), the result will be silently
wrong.

**How to fix:** Document the assumption with a comment, or expose it as a keyword
argument with a default, so users can override it for non-standard chemistries.

### Issue D — Dead-code fallback `return c_guess` in `_get_flux_ms!` and `_get_flux_lm!`

**File:** `src/PAV.jl`, lines 412 and 498.

Every branch in both functions ends with an explicit `return`. The final
`return c_guess` at the end of each function can never be reached. This is harmless
but misleading, because it suggests there is an unhandled case.

**How to fix:** Remove the dead-code `return c_guess` lines, or add a comment
confirming they are unreachable.

### Issue E — `corr_packed` in `ExtractorFunctions.jl` noted as unverified

**File:** `src/ExtractorFunctions.jl`, function `corr_packed` (~line 220).

The docstring notes that "verification of this correlation is incomplete in the
original Python source." The `beta = 0.32` constant (Raschig rings) is also
different from some literature values (0.25). No test covers this function.

### Issue F — `NTU_lm` integrand can produce `DomainError` if `c_out > c_in`

**File:** `src/ExtractorFunctions.jl`, `NTU_lm`, `toint` closure (line 44).

```julia
toint(c) = 1 / (c - K_S * ((R_CONST * T / R_g) * (c - c_out + c_in_gas * R_g))^0.5)
```

The argument to `^0.5` can be negative if `c < c_out` or if `c_out` is computed
as greater than `c_in` by the Brent solver (e.g. when `pg_in = 0` and the
equilibrium concentration is zero). The current `c_out_max` guard in
`get_c_out_GLC_lm` prevents this in the normal path, but it is fragile. Adding an
explicit `max(..., 0.0)` inside `toint` would make the integrand robust to
out-of-range queries.

---

## Part 3 — Functional Differences: Julia vs Python

### Constants and physical parameters

| Item | Python | Julia |
|------|--------|-------|
| Boltzmann constant | `8.617e-5 eV/K` (hardcoded) | From `AtomicAndPhysicalConstants.BOLTZMANN_k * J_PER_EV` |
| Gas constant R | `8.314 J/mol/K` (hardcoded) | `BOLTZMANN_k * J_PER_EV * AVOGADRO ≈ 8.314462618 J/mol/K` |
| Avogadro's number | `6.022e23` (hardcoded) | From `AtomicAndPhysicalConstants.AVOGADRO` |

The small differences between the hardcoded Python values and the CODATA-derived
Julia values (< 0.01%) do not affect test outcomes at current tolerances, but they
mean the two versions are not bit-identical.

### Numerical integration

| Item | Python | Julia |
|------|--------|-------|
| NTU integrals (GLC) | `scipy.integrate.fixed_quad` (Gaussian quadrature, fixed order) | `QuadGK.quadgk` (adaptive Gauss-Kronrod) |
| Inventory integrals | `scipy.integrate.dblquad` | Nested `quadgk` calls |

`fixed_quad` uses a fixed number of quadrature points; `quadgk` is adaptive and
will achieve the requested `atol` automatically. For smooth integrands the results
agree closely; for integrands with rapid variation near an endpoint the Julia version
may be more accurate.

### Optimizer

| Item | Python | Julia |
|------|--------|-------|
| 1-D flux balance | `scipy.optimize.brentq` | `Optim.Brent()` |
| 2-D flux balance | `scipy.optimize.minimize` with method `Powell` | `Optim.Fminbox(NelderMead())` |
| GLC outlet concentration | `scipy.optimize.minimize` | `Optim.Brent()` |
| Column length given kla | `scipy.optimize.minimize` with `Fminbox` | `Optim.Fminbox(NelderMead())` |

`scipy.optimize.brentq` requires two bracket points with opposite sign of `f`
(i.e. it finds a root). `Optim.Brent()` minimises `f²`, finding the minimum of the
absolute residual. For well-posed flux balances these are equivalent, but near a
degenerate minimum (where the two flux curves are tangent) `brentq` will fail while
`Brent()` will still converge to the point of tangency.

Powell's method is not available in the Julia `Optim` package. `Fminbox(NelderMead())`
was substituted, which is a box-constrained Nelder-Mead simplex. It is generally
less efficient than Powell for smooth problems but more robust.

### WireCoil mass-transfer coefficient characteristic length

**This bug exists in both versions.** The Python code uses `t.pitch` as the
characteristic length in `get_k_from_Sh` (which gives `k_t = Sh * D_mol / pitch`).
The correct formula is `Sh = k_t * d_hyd / D_mol`, so the characteristic length
must be `d_hyd`. The Julia code was corrected to use `d_hyd`; the Python source
retains the error. The test expected value for WireCoil `k_t` was updated accordingly.

### Module and class naming

| Python name | Julia name | Status |
|-------------|------------|--------|
| `TriomaModule` | `PAVModule` | Fixed during port |
| `GasLiquidContactorModule` | `GasLiquidContactor` | Fixed during port |
| `Component` (inner class of TriomaModule) | `Component` (exported from PAVModule) | Same |
| `TriomaClass` (Python base class) | Removed — `_is_nested` helper used instead | Removed 2026-06 |

### Output suppression / logging

Python uses `contextlib.redirect_stdout(io.StringIO())` to suppress printed output
during tests. Julia does not support redirecting `stdout` to an `IOBuffer` via the
`redirect_stdout` API on all platforms (no underlying OS file descriptor in some
contexts). The Julia tests use `@test_logs` or `@test_nowarn` macros instead.

### Missing features in the Julia port

| Feature | Python | Julia |
|---------|--------|-------|
| `get_c_out!` plotting | `plotvar=True` keyword triggers matplotlib plots | No plotting implemented |
| Recirculation loop in `outlet_c_comp!` | Uses `analytical_efficiency!` internally | Same; uses `analytical_efficiency!` (analytical only, not numerical) |
| `estimate_cost!` | Empirical cost correlation | Ported, but not validated |

---

## Part 4 — Testing Differences: Julia vs Python

### Test structure

The Python test suite uses `pytest` with test classes mirroring component classes.
The Julia suite is a flat `Test.@testset` hierarchy in `test/test_component_tools.jl`.
The Julia tests cover the same scenarios but do not use parametric fixtures or
`pytest.mark` equivalents.

### Tolerances

All floating-point physics comparisons use `rtol=1e-5`, chosen to accommodate the
~5-sig-fig difference between Python's hardcoded constants (e.g. `R = 8.314`) and
Julia's CODATA-derived values from `AtomicAndPhysicalConstants.jl`. Integer-valued
results (e.g. recirculation flags, pipe counts) use exact `==` comparisons.

### Tests marked `@test_broken`

These correspond to physically incorrect formulas in the source (Bug 1) or to
fundamental inconsistencies between the analytical and numerical approaches (Bugs 2
and 3). The Python test suite may carry the same failures silently (the Python tests
for these cases may be passing due to looser tolerances, different regime detection
thresholds, or simply not being present).

---

*End of issues document.*

# AGENTS.md — Coding Governance (Tuning Phase)

**Scope:** This document governs all code written for the `langmuir_rebuild/` project. Read before writing any code. Rules are hard constraints.

**Project phase:** The physics engine (WP-01 through WP-06) is built and working. We are now in a **tuning and diagnostic phase** focused on understanding why model predictions have the wrong rank ordering against observations, and on fixing the observation/comparison layer.

---

## Context Efficiency Rules (for agentic coders)

1. **Never read a whole file when you only need a section.** Use `offset` and `limit` parameters. Large docs cost 5–27k tokens each.
2. **Read docs once per session, then work from memory.**
3. **Read source files before editing.** Use targeted reads (20–30 lines around a function).
4. **Prefer Grep over Read for locating content.**
5. **Write implementations in one pass per file.** Plan before calling Write.
6. **Break large implementations into separate files.**

---

## Part I — Foundational Principles

### 1. The Hierarchy of Concerns

```
1. Physical correctness   — Conservation laws, dimensional consistency, known constraints
2. Falsifiability         — Can this choice be shown wrong with available data?
3. Transparency           — Can a reviewer see exactly why a prediction was made?
4. Robustness             — Works across observed input range without hidden saturation?
5. Accuracy               — Does this reduce prediction error?
6. Elegance               — Simple and clean?
```

**Accuracy is fifth, not first.** A transparent, falsifiable, robust model that is slightly less accurate is strictly preferred over one that is accurate but opaque and fragile.

### 2. Silent Saturation

The single most dangerous pattern: a quantity that stops varying in response to its inputs without anyone noticing. This killed the predecessor model at least three times.

**The saturation rule:** Every function returning a physical quantity must be auditable for input-output sensitivity. Output varies <5% over >50% of physically plausible input range → function is saturated → must be fixed or justified in the assumptions register.

### 3. The Separation Principle

Five conceptually distinct layers, implemented as independent modules:

```
External forcing  →  Hydrodynamic response  →  Coherent structure state
                                                        ↓
                                            Visibility / tracer expression
                                                        ↓
                                            Cyanobacteria accumulation
```

No circular dependencies. Interface between layers is a frozen dataclass with named, typed, unit-documented fields.

### 4. The Scale Distinction

Four spatial scales that must never be silently conflated:

| Scale | What it is | Source |
|---|---|---|
| Instability scale (L_inst) | Wavelength from CL instability | Nonlinear CL solver: 2πh / l_cNL |
| Cell scale (coarsened_width) | Mature cell pair width after mergers | L_inst × 2^n_events |
| Mechanical cell width | Cell width capped by aspect ratio | min(coarsened_width, 12×depth) |
| Visible spacing (comparison_spacing) | What satellites observe | Should be mechanical_cell_width × visibility_factor |

**Critical current issue:** `comparison_spacing` currently equals `coarsened_width` (uncapped hierarchy), not `mechanical_cell_width`. This is a known bug — see §9.

---

## Part II — Code Patterns and Anti-Patterns

### 5. Mandatory Patterns

- **Every function documents units** for all inputs and outputs in its docstring. Enforced by `test_dimensional_consistency.py`.
- **Frozen dataclasses** for layer interfaces (`@dataclass(frozen=True)`).
- **Explicit conversion functions** between scales — named, documented, registered in assumptions register.
- **Raw values always preserved** alongside any normalisation. Never discard amplitude information.

### 6. Forbidden Anti-Patterns

**6.1 Hidden bounding.** No sigmoid, tanh, clip, or min/max that compresses output without physical justification. If physics produces an out-of-range value, that's diagnostic information — flag it with `warnings.warn`, return unmodified.

**6.2 Compensatory memory.** No exponential smoothing to fix dynamics. If the instantaneous prediction is too volatile, fix the prediction, not the smoothing.

**6.3 Tuning to aggregate error.** Do not adjust parameters because they improve RMSE. A parameter change requires physical justification. If it improves RMSE but has no physics basis, it's overfitting — document in failure log, do not apply.

**6.4 Implicit scale ambiguity.** Every variable that holds a spacing value must state which of the four scales it represents.

---

## Part III — Tuning-Phase Workflow

### 7. Current Project State

**What works:**
- CL nonlinear solver reproduces Hayes & Phillips (2017) published values
- Forcing module: wind, waves, currents, eddy viscosity all functional
- Coarsening module: discrete Y-junction mergers, disruption detection
- 95 tests pass, 67/67 ERA5 weather cache populated
- Dynamic range restored: CL range coverage = 0.973

**What's broken:**
- **Rank ordering is wrong.** Pearson r = -0.317, Spearman ρ = -0.464 against pooled observations.
- The model predicts the right spread of values but assigns them to the wrong cases.

**The diagnostic picture (from grouped Neagh audit):**

| Category | n | Median observed [m] | Median L_inst [m] | L_inst / observed | n_events |
|----------|---|--------------------:|-------------------:|------------------:|---------:|
| wiggle | 3 | 214 | 228 | 0.94 | 0 |
| stream | 18 | 69 | 195 | 2.4 | varies |
| manual | 33 | 43 | 196 | 3.9 | varies |

**Key finding:** Wiggle observations are well-matched by the raw onset scale with no coarsening. Manual observations are already too large at onset (factor 3.9×), and coarsening makes them worse.

**Working hypothesis:** The wiggle-manual difference reflects coarsening initialization. Wiggle spacing captures a larger-scale process that matches raw CL onset. Manual spacing captures smaller-scale features. The current model initializes at too large a value (L_inst ~196 m for manual cases) and then coarsens, pushing predictions further from observations. The fix direction is: smaller initialization and less aggressive coarsening for conditions that produce manual-class observations.

**Background on coarsening and initialization:**
- Wiggle does better with no coarsening steps at all
- Manual currently initializes with too large a value and then coarsens
- This may be partly a coarsening issue — at larger initial widths, the merger timescale τ ~ λ²/A_H becomes very long, so fewer events occur, but the initial width is already wrong
- Need to investigate how to initialize at smaller sizes and avoid premature coarsening

### 8. Investigation Protocol

**8.1 Always separate manual and wiggle in diagnostics.** Never collapse them into one target. Stream observations should be tracked separately too but are a lower priority.

**8.2 Layer-by-layer rank audit.** When changing any component, report Spearman ρ at each layer separately:
- Forcing (Ra, La_t)
- L_inst (raw onset width)
- coarsened_width (after mergers)
- mechanical_cell_width (after aspect-ratio cap)
- comparison_spacing (observation-facing)

This identifies which layer first breaks ordering.

**8.3 Parameter sensitivity exploration.** When exploring a parameter:
1. State which layer the parameter enters
2. Sweep across physically plausible range
3. Report effect on Spearman ρ (by observation class) and on median bias
4. Do NOT adopt a value just because it improves a metric — require physical justification

**8.4 When is a parameter change justified?**
- Justified: literature supports a different value; dimensional analysis constrains a range; the old value was a placeholder
- Justified: a closure or profile shape change with a published physical basis
- NOT justified: "this value of γ improves RMSE by 3 m" with no physics reason
- NOT justified: tuning a threshold to reclassify edge cases

### 9. Observation Operator Audit

**The comparison_spacing bug:** The final comparison layer currently restores the uncapped raw hierarchy width instead of preserving the mechanically capped width. This means `comparison_spacing == coarsened_width`, bypassing the 12×depth aspect-ratio cap. This must be fixed.

**Which model quantity should be compared against each observation class:**
- This is an open question, not a settled assumption
- Current evidence suggests wiggle ↔ L_inst (raw onset, no coarsening)
- Manual observations are much smaller than any current model layer
- The observation operator may need to be class-dependent, or the initialization may need fixing so that one operator works for both

**Rules for the observation operator:**
- It must be an explicit, named function (not an inline multiplication)
- It must state which model scale it maps from and which observation class it targets
- Any class-dependent logic must be registered in assumptions_register.md

### 10. Living Documents

Four documents are updated continuously — they are primary deliverables, not afterthoughts.

- **`docs/assumptions_register.md`** — Every modelling assumption with source, falsification criterion, status (Active/Uncertain/Falsified/Superseded)
- **`docs/decision_log.md`** — Every architectural or modelling choice with alternatives and rationale
- **`docs/failure_log.md`** — Everything that didn't work, with root cause and resolution
- **`docs/observation_model.md`** — What "spacing" means and how model quantities map to observations

---

## Part IV — Physics Guidance

### 11. The CL Solver Is One Component

The CL solver outputs a nondimensional critical wavenumber l_cNL. Converting to a dimensional spacing requires the forcing-to-Ra mapping, scale conversions, coarsening, and observation model. Each has its own assumptions and failure modes.

**For this phase:** The solver itself works correctly. The issue is upstream (what forcing produces what Ra) and downstream (how onset width maps to observations). Do not modify the CL solver unless diagnostic evidence specifically implicates it.

### 12. Ra Is The Master Control

Ra = U × D × h² / ν_T². Because Ra ∝ ν_T⁻², eddy viscosity enters quadratically. A factor-of-2 error in ν_T produces factor-of-4 error in Ra.

**For this phase:**
- Always report Ra and its constituent terms alongside predictions
- When diagnosing ordering problems, check whether Ra itself is correctly ordered across cases
- Consider ν_T uncertainty: run at ν_T, 0.5×ν_T, 2×ν_T to bracket

### 13. Coarsening: The Central Investigation Area

**Current coarsening logic:**
- Onset width: L_inst = 2πh / l_cNL (from CL solver)
- Mergers: width doubles per event, timescale τ = λ² / ((2π)² × A_H)
- Lateral diffusivity: A_H = max(ν_T, u*_water × depth)
- Mechanical cap: min(coarsened_width, 12 × depth)
- Disruption: resets on <2 m/s wind, >45° direction change, >3 m/s speed jump

**The problem:** L_inst ≈ 196 m for typical manual-class Neagh cases (depth ≈ 9 m). This is already an aspect ratio of ~22, far above the 12× cap. So even before coarsening, the onset width exceeds the mechanical cap.

**Investigation priorities:**
1. **What controls onset width?** l_cNL depends on Robin BCs (γ_s, γ_b), profiles (U', D'), and Ra. For high Ra (supercritical), what determines l_cNL and hence L_inst?
2. **Is 196 m physically reasonable?** At depth 9 m, aspect ratio 22 at onset seems high. Check whether the profile closures or BC derivations are producing reasonable l_cNL values.
3. **Can onset width be made smaller for the right conditions?** The user notes that smaller initialization is needed. Investigate which forcing/profile conditions produce smaller l_cNL (larger wavenumber, tighter cells).
4. **Coarsening aggression:** Even if onset width is correct, the merger clock may be too fast or the hierarchy expansion may be too aggressive. Separate these effects.

### 14. The Bio Module

The LC enhancement calculator computes enhancement ratios (light, nutrients, temperature) relative to a static water column. Every biological metric is a named triple: (value_LC, value_static, enhancement_ratio). This module is complete and not the current investigation focus.

---

## Part V — Testing During Tuning

### 15. Test Discipline

- **All 95 existing tests must still pass** after any code change
- **Run saturation audit** after any parameter or closure change
- **Rank correlation (Spearman ρ) is the primary metric**, not RMSE alone
- RMSE improvement without rank improvement is not improvement
- Dynamic range must be preserved (attractor test)
- Report metrics separately for manual, wiggle, and stream classes

---

## Quick Reference Card

```
BEFORE CHANGING ANY PARAMETER OR CLOSURE:
  □ What physical quantity does this affect?
  □ What is the literature basis for the change?
  □ Does this improve Spearman ρ for both manual and wiggle classes?
  □ Does the saturation audit still pass?
  □ Have I updated assumptions_register.md?

BEFORE COMMITTING:
  □ All 95+ existing tests pass
  □ Living documents updated
  □ Layer-by-layer rank audit shows which layer changed

DIAGNOSTIC CHECKLIST:
  □ Always separate manual / wiggle / stream in outputs
  □ Report Spearman ρ at each layer (Ra → L_inst → coarsened → capped → comparison)
  □ Check if comparison_spacing uses mechanical_cell_width (not uncapped hierarchy)
  □ Report median bias by observation class

WHEN SOMETHING GOES WRONG:
  □ Log in failure_log.md FIRST
  □ Do not adjust test targets
  □ Do not add compensatory wrappers
  □ Do not tune to aggregate error without physics justification
  □ Ask for help if stuck
```

---

## Appendix: File Header Template

Every Python file in `src/` should include:

```python
"""
Module: [module name]
Layer: [forcing | hydro | prediction | evaluation]
Purpose: [one sentence]

Assumptions: [A-XX] — see assumptions_register.md
Parameters: [P-XX] — see parameter_register.md
Key references: [Author Year]
"""
```

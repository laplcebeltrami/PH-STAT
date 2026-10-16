# PH-STAT workflow animation

The animation is a visual walkthrough of methods in [`SCRIPT.mlx`](https://github.com/laplcebeltrami/PH-STAT/blob/835d5eecc760181ecef06883afbbb42c24bcfb0e/SCRIPT.mlx), using **computed synthetic examples**.

**Format:** 1600 × 1000 pixels, 48.7 seconds, 20 frames per second, continuous loop. **GIF size:** 3,016,488 bytes (3.016 MB), below the 25 MB limit.

## Source and reproducibility

- Toolbox source: [PH-STAT commit `835d5ee`](https://github.com/laplcebeltrami/PH-STAT/tree/835d5eecc760181ecef06883afbbb42c24bcfb0e), inspected October 16, 2026.
- Computation: MATLAB R2026b with Statistics and Machine Learning Toolbox.
- Fixed seeds: `20261016` for the new demonstrations; `1` and `2` for the Live Script's two-group simulation.
- Data export: [`export_workflow_data.m`](../../tools/workflow/export_workflow_data.m).
- Saved numeric outputs: [`data.json`](../../tools/workflow/data.json).
- Renderer: [`render_workflow.py`](../../tools/workflow/render_workflow.py), using NumPy and Pillow. Text and graphics are rendered at twice the final resolution and downsampled. A shared palette preserves colors across frames.

The selected examples were executed and checked. This does not constitute an end-to-end test of all sections of the Live Script.

## What each sequence shows

| Sequence | Computation and interpretation |
| --- | --- |
| **Morse filtration** | The Live Script's signal on `[0,1]`, Gaussian noise with standard deviation 0.2, cosine degree 100, and smoothing bandwidth 0.0001. The displayed finite pairs are returned by `PH_morse1D`. |
| **Rips complexes** | A new 16-point noisy ring, with radius increasing from 0.05 to 1.9. `PH_rips` constructs edges and filled triangles; `PH_boundary_betti` computes components and loops. The loop forms and then fills. |
| **Graph filtration** | A new 12-node weighted network, with weights derived from spatial proximity plus small symmetric jitter. `PH_betti` uses the strict rule `weight > threshold`. Increasing the threshold removes edges: beta-0 is nondecreasing and beta-1 is nonincreasing. No triangles are filled in this graph sequence. |
| **Inference** | The Live Script's 50-node, 4-versus-6 simulation. `WS_pdist2` supplies combined squared Wasserstein distances, and `WS_ratio` supplies the ratio of between-group to within-group distance sums. All 210 label allocations are enumerated for this small example: the displayed upper-tail p-value is exact conditional on these simulated networks. A separate 10,000-step `WS_transpositions` run is also saved for reference. The animation does not present the enumeration as a transposition run. |
| **Clustering and embedding** | `simulate_circle_difference` generates four groups of eight point clouds, with 24 points per circle and noise standard deviation 0.18. Pairwise Euclidean distances define the weighted networks. `WS_cluster` operates on full birth/death features; `WS_embed` supplies the displayed mean birth and death coordinates. Colors indicate known generating groups, not fitted cluster labels. Two-dimensional proximity is not an exact representation of the full Wasserstein distance. |
| **Hodge decomposition** | The Live Script's five-node, six-edge example with one filled triangle. `Hodge_decompose` returns gradient, curl, and harmonic flows. Edge indexing is checked against the adjacency matrix before plotting. The same arrow-width scale is used across the four panels. Moving markers indicate direction; they are not particle simulations. |
| **Edge diffusion** | The same filled-triangle example uses `L1 = B1' * B1 + B2 * B2'`. The graph term comes from `Hodge_laplacian`; the face term is added explicitly. Spectral diffusion `exp(-t*L1)*X` follows the formula in the Live Script. The time parameter is a mathematical diffusion scale. Gradient and curl components decay toward the harmonic projection. |

The numerical checks include monotonic graph Betti curves, symmetry and additivity of the squared-distance matrices, and Hodge reconstruction, divergence, curl, and orthogonality residuals. The largest Hodge residual was below `1.2e-15` for this example. The inference example produced `R = 50.279597` and `p = 1/210 = 0.0047619`; these are synthetic demonstration values, not evidence about a clinical group.

## Regenerate

To render from the saved outputs, run the following from the repository root in an environment with Python 3, NumPy, and Pillow:

```sh
python3 tools/workflow/render_workflow.py
```

To regenerate the numeric data, place a checkout of the pinned PH-STAT commit in `reference/github-PH-STAT`, then run from the repository root in MATLAB:

```matlab
addpath(fullfile(pwd,'tools','workflow'))
export_workflow_data(fullfile(pwd,'reference','github-PH-STAT'), ...
    fullfile(pwd,'tools','workflow','data.json'))
```

The renderer uses Arial fonts on macOS, or DejaVu Sans on Linux. Font substitution can change text appearance. It writes the GIF and still overview to `assets/workflow/`, along with intermediate previews and validation metadata in `build/workflow/`. The export enforces a maximum GIF size of **25,000,000 bytes**.

# Closed-Form Moment-Curvature Model

## Overview
This repository provides a GitHub remonstrator for the **Closed-Form Moment-Curvature Model**, which simulates the complete flexural response of fiber-reinforced concrete (FRC). This analytical model extends previous work on moment-curvature relationships by incorporating multilinear constitutive models applicable to hybrid FRC, textile-reinforced concrete, and structural steel sections.

## Background
The Closed-Form Moment-Curvature Model was first introduced by **Soranakom and Mobasher**, who derived equations based on a bilinear constitutive stress-strain law for both tension and compression. Later advancements by **Yao et al.** and **Pleesudjai et al.** generalized the approach for various reinforcement configurations, allowing analytical predictions of rotations, deflections, and crack widths in hybrid systems.

This model uses:
- **Bilinear and tri-linear representations** for reinforcing bars in singly and doubly reinforced configurations.
- **Hybrid FRC modeling**, extended to Ultra-High-Performance Concrete (UHPC), incorporating:
  - A quadrilinear tensile model.
  - A stiffness parameter in the post-yield compression zone.
  - A stiffness parameter for transverse reinforcement to capture hardening effects.

## Features
- **Parametric Model:**
  - Simulates beam cross-sections under flexure.
  - Defines strain and stress distributions at ultimate limit states.
- **Closed-Form Equations:**
  - Provides material models as piecewise linear equations.
  - Uses dimensionless parameters for scalable computation.
- **Hybrid System Compatibility:**
  - Supports different reinforcement materials (FRP, steel, strengthening solutions).
  - Captures strain-hardening and softening behavior.
- **Hierarchical Stage Approach:**
  - Identifies transition points between material zones.
  - Determines moment-curvature relationships iteratively.

## Implementation
The model is based on a **flexural hinge approach** using **Euler-Bernoulli beam theory**, ensuring precise calculation of stress and strain states at any given strain level. The neutral axis depth is determined by solving an internal force equilibrium equation in closed-form.

The strain incrementally increases, and materials transition between their defined stages, updating stress-strain relationships accordingly. The hierarchical derivation results in an envelope of moment-curvature response curves, where each stage’s boundary defines transition points for the next phase.

## Key Parameters
- **Tension Model:** Defined by stiffness parameters and strain parameters for each post-crack zone.
- **Compression Model:** Includes post-yield stiffness and strain limit parameters.
- **Reinforcement Model:** Accounts for post-yield hardening effects in transverse reinforcement.
- **Neutral Axis Depth:** Computed based on force equilibrium conditions.

## Results & Validation
- The model generates **closed-form solutions** for all possible stage combinations.
- It is applicable to **any combination of material characteristics**, including hybrid systems.
- The model establishes an **envelope for the moment-curvature response**, allowing direct determination of the most critical force path.
- The normalized moment-curvature relationship is solved analytically, ensuring accurate prediction of flexural performance.

## Example Outputs
- **Neutral axis depth (k)** vs. **Normalized Moment (M/Mcr)**.
- **Moment-curvature response envelopes** for different reinforcement and material combinations.
- **Color-graded tension zones**, distinguishing transitions between different stages.

## References
- Soranakom and Mobasher (2008)

## Future Work
- Implementation of **numerical solvers** to complement closed-form solutions.
- Expansion to **3D beam models** for more complex structural elements.
- Integration with **finite element analysis tools** for broader applications.

---
This GitHub repository aims to facilitate further research and practical application of the Closed-Form Moment-Curvature Model in advanced FRC and UHPC design. Contributions and collaborations are welcome!

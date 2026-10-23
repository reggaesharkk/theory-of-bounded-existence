import Lake
open Lake DSL

package «boundaryVerification» where
  leanOptions := #[
    ⟨`pp.unicode.fun, true⟩,
    ⟨`autoImplicit, false⟩
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "10156080d0e7997b71f3ef86bd6f0211b42ab90d"

@[default_target]
lean_lib «BoundaryVerification» where
  roots := #[`BoundaryVerification]

@[default_target]
lean_lib «ExistenceVerification» where
  roots := #[`ExistenceVerification]

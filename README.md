# sml-nbody

Barnes-Hut quadtree gravitational N-body simulation in pure Standard ML — a
mass/center-of-mass quadtree, an opening-angle (`theta`) force approximation, a
direct O(n²) sum for cross-checking, and a symplectic velocity-Verlet
(leapfrog) integrator over all bodies. Built on
[`sml-glm`](https://github.com/sjqtentacles/sml-glm) for the `Vec2` 2D vector
type. No FFI, no external dependencies, and **deterministic**, byte-identically
under both [MLton](http://mlton.org/) and [Poly/ML](https://www.polyml.org/).

## Status

- 31 assertions, green on MLton and Poly/ML (byte-identical pass count).
- Basis-library + vendored `sml-glm` only; deterministic across compilers.
- Vendors `sml-glm` (Layout B), so the repo builds standalone.

## Install

With [`smlpkg`](https://github.com/diku-dk/smlpkg):

```
smlpkg add github.com/sjqtentacles/sml-nbody
smlpkg sync
```

Include the MLB from your own (it pulls in the vendored `sml-glm`):

```
local
  $(SML_LIB)/basis/basis.mlb
  lib/github.com/sjqtentacles/sml-nbody/... (via smlpkg)
in
  ...
end
```

This brings `structure NBody` (and the vendored `Glm`) into scope.

## Quick start

```sml
structure V = NBody.V

(* two equal masses in a mutual circular orbit (separation 2, G = m = 1) *)
val bodies =
  [ { pos = V.v (~1.0, 0.0), vel = V.v (0.0, ~0.5), mass = 1.0 }
  , { pos = V.v ( 1.0, 0.0), vel = V.v (0.0,  0.5), mass = 1.0 } ]

(* one symplectic leapfrog step with the exact direct sum ... *)
val b1 = NBody.stepDirect 1.0 0.001 bodies
(* ... or with the Barnes-Hut approximation (opening angle theta) *)
val b2 = NBody.stepBH { g = 1.0, theta = 0.5 } 0.001 bodies
val bN = NBody.stepNBH { g = 1.0, theta = 0.5 } 0.001 1000 bodies

(* gravitational acceleration at a point: a = G m / r^2 *)
val t = NBody.build bodies
val a = NBody.bhAccel { g = 1.0, theta = 0.5 } t (V.v (0.0, 3.0))

(* conserved quantities *)
val p = NBody.totalMomentum bodies
val e = NBody.totalEnergy 1.0 bodies
```

## API (`signature NBODY`)

```sml
structure V : sig          (* a thin view over Glm.Vec2 *)
  type t
  val v : real * real -> t   val x : t -> real   val y : t -> real
  val zero : t
  val add : t * t -> t       val sub : t * t -> t
  val scale : real * t -> t  val dot : t * t -> real
  val length : t -> real     val lengthSq : t -> real   val dist : t * t -> real
end

type body = { pos : V.t, vel : V.t, mass : real }

(* quadtree *)
type tree
val build    : body list -> tree
val treeMass : tree -> real
val treeCom  : tree -> V.t

(* gravitational acceleration at a point *)
val directAccel : real -> body list -> V.t -> V.t
val bhAccel     : { g : real, theta : real } -> tree -> V.t -> V.t

(* symplectic (velocity-Verlet / leapfrog) stepping over all bodies *)
val stepDirect  : real -> real -> body list -> body list
val stepBH      : { g : real, theta : real } -> real -> body list -> body list
val stepNDirect : real -> real -> int -> body list -> body list
val stepNBH     : { g : real, theta : real } -> real -> int
                  -> body list -> body list

(* diagnostics *)
val totalMomentum   : body list -> V.t
val centerOfMass    : body list -> V.t
val kineticEnergy   : body list -> real
val potentialEnergy : real -> body list -> real
val totalEnergy     : real -> body list -> real
```

### Conventions

- Gravity is Newtonian with a caller-supplied constant `g`: the acceleration at
  a point is `Σⱼ g·mⱼ·(xⱼ − point) / |xⱼ − point|³`. A mass at zero distance
  contributes nothing, so a body never attracts itself and coincident masses are
  safe.
- `build` computes the bounding square automatically and partitions the bodies
  into the four quadrants recursively; each node stores its total mass, center
  of mass, and side length.
- `bhAccel {g, theta}` treats a node as a single point mass at its center of
  mass when its width `s` and the distance `d` to the query point satisfy
  `s/d < theta`, otherwise it recurses. **`theta = 0` recurses fully and equals
  `directAccel` to rounding**; larger `theta` trades accuracy for speed.
- `stepDirect` / `stepBH` advance every body one velocity-Verlet (leapfrog)
  step. The method is **symplectic** (energy stays bounded) and, because the
  pairwise forces are equal and opposite, **conserves total momentum** to
  rounding. `stepBH` rebuilds the tree each step.
- `potentialEnergy g` sums `−g·mᵢ·mⱼ / rᵢⱼ` over distinct pairs;
  `totalEnergy = kineticEnergy + potentialEnergy`.
- Everything is pure and deterministic: bodies are given explicitly (no RNG),
  so the same inputs produce the same reals on every run, machine, and compiler.

## Build & test

```
make test        # MLton
make test-poly   # Poly/ML
make all-tests   # both
make example     # build + run examples/demo.sml (writes assets/nbody.txt)
make clean
```

Both compilers run the same strict-TDD suite, seeded with closed-form vectors:
the quadtree mass / center-of-mass aggregation, the two-body law `a = G m / r²`,
the exact `theta = 0` Barnes-Hut / direct-sum equality and small-`theta`
agreement, and a two-body circular orbit whose separation stays constant, whose
total momentum is conserved to rounding, and whose energy (`E = −0.25` by
construction) stays bounded under leapfrog.

## Example

`make example` integrates a heavy central mass plus a ring of eight unit masses
on near-circular orbits with Barnes-Hut leapfrog, printing an
energy/momentum-vs-step table and an ASCII scatter of the final positions
(output is byte-identical under MLton and Poly/ML, and committed to
[`assets/nbody.txt`](assets/nbody.txt)):

```
=== sml-nbody demo ==========================================

Central mass 1000.0 at origin, plus a ring of 8 unit masses
at radius 5.0 on circular orbits (v = 14.1421).
Barnes-Hut theta = 0.50, leapfrog dt = 0.0005.

   step       t          KE            PE           E       |P|
      0   0.000     800.0000    -1604.4878   -804.4878  0.0000
    250   0.125     800.2784    -1604.7662   -804.4878  0.0000
    500   0.250     801.0810    -1605.5688   -804.4878  0.0001
    750   0.375     802.3092    -1606.7969   -804.4878  0.0001
   1000   0.500     803.8131    -1608.3008   -804.4878  0.0001
   1250   0.625     805.4061    -1609.8938   -804.4878  0.0001
   1500   0.750     806.8889    -1611.3767   -804.4878  0.0002
   1750   0.875     808.0741    -1612.5619   -804.4878  0.0002
   2000   1.000     808.8098    -1613.2976   -804.4878  0.0002

Final positions after 2000 steps  (O = central, * = ring body):
  +-----------------------------------------+
  |                                         |
  |                                         |
  |                                         |
  |                                         |
  |                        *                |
  |              *                          |
  |                                         |
  |                               *         |
  |        *                                |
  |                                         |
  |                    O                    |
  |                                         |
  |                                *        |
  |         *                               |
  |                                         |
  |                          *              |
  |                *                        |
  |                                         |
  |                                         |
  |                                         |
  |                                         |
  +-----------------------------------------+

============================================================
```

The total energy `E` holds at `−804.4878` across the run (symplectic leapfrog)
and the total momentum magnitude `|P|` stays at the round-off floor — the ring
is momentum-balanced by construction.

### Poly/ML note

CI builds Poly/ML 5.9.1 from source rather than using the Ubuntu package
(Poly/ML 5.7.1), whose X86 code generator crashes (`asGenReg raised while
compiling`) on heavy real-arithmetic code. See `.github/workflows/ci.yml`.

## License

MIT — see [LICENSE](LICENSE).

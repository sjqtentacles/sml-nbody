(* nbody.sig

   Barnes-Hut quadtree gravitational N-body simulation in pure Standard ML,
   built on the `Glm.Vec2` 2D vector type from sml-glm.

   The pieces:

     - A `tree`: a quadtree over the bodies' bounding square that aggregates the
       total mass and center of mass of each region.  `build` partitions the
       bodies into the four quadrants recursively.

     - Force approximation: `bhAccel {g, theta}` walks the tree, treating a
       region as a single point mass at its center of mass when its width `s`
       and the distance `d` to the query point satisfy `s/d < theta` (the
       opening angle), and otherwise recursing.  `theta = 0` recurses fully and
       so reproduces the direct O(n^2) sum (`directAccel`) to rounding.

     - A symplectic velocity-Verlet (leapfrog) integrator over all bodies,
       available both with the Barnes-Hut approximation (`stepBH`) and the exact
       direct sum (`stepDirect`).  Pairwise forces are equal and opposite, so
       total momentum is conserved.

   Acceleration is computed at a *point*; a mass at zero distance (a body's own
   contribution, or coincident masses) contributes nothing, so self-interaction
   is excluded automatically.  Gravity uses a caller-supplied constant `g`.

   Everything is pure and deterministic (no FFI, randomness, or wall-clock):
   bodies are given explicitly, so results are byte-identical across MLton and
   Poly/ML.  Test comparisons go through an explicit epsilon. *)

signature NBODY =
sig
  structure V : sig
    type t
    val v        : real * real -> t
    val x        : t -> real
    val y        : t -> real
    val zero     : t
    val add      : t * t -> t
    val sub      : t * t -> t
    val scale    : real * t -> t
    val dot      : t * t -> real
    val length   : t -> real
    val lengthSq : t -> real
    val dist     : t * t -> real
  end

  type body = { pos : V.t, vel : V.t, mass : real }

  (* --- quadtree --- *)
  type tree

  (* Build a quadtree from a list of bodies (the bounding square is computed
     automatically).  An empty list yields an empty tree. *)
  val build : body list -> tree

  (* Aggregate mass and mass-weighted center of mass of a (sub)tree. *)
  val treeMass : tree -> real
  val treeCom  : tree -> V.t

  (* --- gravitational acceleration at a point --- *)

  (* Direct O(n^2) sum: acceleration at `point` due to every body, using
     Newtonian gravity a = sum_j g m_j (x_j - point) / |x_j - point|^3.  A body
     at `point` (distance 0) contributes nothing. *)
  val directAccel : real -> body list -> V.t -> V.t

  (* Barnes-Hut approximation of the same acceleration with opening angle
     `theta`.  `theta = 0` recurses fully (== directAccel to rounding). *)
  val bhAccel : { g : real, theta : real } -> tree -> V.t -> V.t

  (* --- symplectic (velocity-Verlet / leapfrog) stepping over all bodies --- *)

  (* One leapfrog step of size dt using the exact direct sum. *)
  val stepDirect : real -> real -> body list -> body list
  (* One leapfrog step of size dt using Barnes-Hut (tree rebuilt internally). *)
  val stepBH     : { g : real, theta : real } -> real -> body list -> body list

  val stepNDirect : real -> real -> int -> body list -> body list
  val stepNBH     : { g : real, theta : real } -> real -> int
                    -> body list -> body list

  (* --- diagnostics --- *)
  val totalMomentum : body list -> V.t
  val centerOfMass  : body list -> V.t
  val kineticEnergy : body list -> real
  (* Sum over distinct pairs of -g m_i m_j / r_ij. *)
  val potentialEnergy : real -> body list -> real
  val totalEnergy     : real -> body list -> real
end

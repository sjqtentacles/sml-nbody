(* test_force.sml -- gravitational acceleration: direct sum and Barnes-Hut.

   The two-body law a = G m / r^2 is pinned exactly.  Barnes-Hut with theta = 0
   recurses fully and so reproduces the direct sum to rounding; with a small
   positive theta it approximates the direct sum within a tight tolerance when
   the query point is well separated from a cluster. *)

structure ForceTests =
struct
  open Support
  structure N = NBody

  fun run () =
    let
      val () = Harness.section "force: two-body a = G m / r^2"

      (* a mass m = 5 at the origin; acceleration on a point at (2,0), g = 1.
         magnitude = g m / r^2 = 5/4 = 1.25, directed toward the mass (-x). *)
      val a = N.directAccel 1.0 [ body (0.0, 0.0, 0.0, 0.0, 5.0) ] (V.v (2.0, 0.0))
      val () = checkVecTol 1E~12 "accel = (-1.25, 0)" (V.v (~1.25, 0.0), a)
      val () = checkApproxTol 1E~12 "|accel| = G m / r^2 = 1.25"
                 (1.25, V.length a)

      (* gravitational constant scales the force linearly *)
      val a2 = N.directAccel 2.0 [ body (0.0, 0.0, 0.0, 0.0, 5.0) ] (V.v (2.0, 0.0))
      val () = checkApproxTol 1E~12 "doubling g doubles accel"
                 (2.5, V.length a2)

      (* a body does not attract itself: query at a body's own position skips it *)
      val aSelf = N.directAccel 1.0 [ body (1.0, 1.0, 0.0, 0.0, 9.0) ] (V.v (1.0, 1.0))
      val () = checkVecTol eps "no self-interaction" (V.zero, aSelf)

      val () = Harness.section "Barnes-Hut theta=0 reproduces the direct sum"

      val bs =
        [ body (0.0, 0.0, 0.0, 0.0, 1.0)
        , body (1.0, 0.0, 0.0, 0.0, 2.0)
        , body (0.0, 1.0, 0.0, 0.0, 1.0)
        , body (1.0, 1.0, 0.0, 0.0, 3.0)
        , body (0.5, 2.0, 0.0, 0.0, 1.0) ]
      val t = N.build bs
      val g = 1.0
      val () =
        List.app
          (fn b =>
             let
               val d  = N.directAccel g bs (#pos b)
               val bh = N.bhAccel { g = g, theta = 0.0 } t (#pos b)
             in checkVecTol 1E~9 "BH(theta=0) == direct at a body" (d, bh) end)
          bs

      val () = Harness.section "Barnes-Hut small theta approximates direct"

      (* a compact cluster, queried from far away: the opening criterion lumps
         the cluster into its center of mass, with negligible error *)
      val cluster =
        [ body (0.0, 0.0, 0.0, 0.0, 1.0)
        , body (0.1, 0.0, 0.0, 0.0, 1.0)
        , body (0.0, 0.1, 0.0, 0.0, 1.0)
        , body (0.1, 0.1, 0.0, 0.0, 1.0) ]
      val tc = N.build cluster
      val q = V.v (10.0, 10.0)
      val dFar  = N.directAccel g cluster q
      val bhFar = N.bhAccel { g = g, theta = 0.5 } tc q
      val () = checkVecTol 1E~6 "BH(theta=0.5) ~ direct for a far cluster"
                 (dFar, bhFar)
      (* and a larger theta is still within a loose tolerance *)
      val bhFar2 = N.bhAccel { g = g, theta = 1.0 } tc q
      val () = checkVecTol 1E~4 "BH(theta=1.0) ~ direct for a far cluster"
                 (dFar, bhFar2)
    in
      ()
    end
end

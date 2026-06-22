(* test_orbit.sml -- leapfrog two-body dynamics.

   Two equal masses on opposite sides of the origin in a mutual circular orbit:
   with separation 2 (orbit radius 1), G = m = 1, the gravitational
   acceleration on each is G m / r_sep^2 = 1/4, which equals the centripetal
   v^2 / R for v = 1/2.  Under the symplectic velocity-Verlet integrator the
   orbit stays near-circular (separation ~ constant over a period), total
   momentum is conserved exactly, and total energy stays bounded. *)

structure OrbitTests =
struct
  open Support
  structure N = NBody

  fun nth (xs, i) = List.nth (xs, i)
  fun sep bs = V.dist (#pos (nth (bs, 0)), #pos (nth (bs, 1)))

  fun run () =
    let
      val g = 1.0
      val v = 0.5
      val bs0 =
        [ body (~1.0, 0.0, 0.0, ~v, 1.0)
        , body ( 1.0, 0.0, 0.0,  v, 1.0) ]

      val period = 4.0 * Math.pi          (* 2*pi*R/v = 2*pi*1/0.5 *)
      val steps = 4000
      val dt = period / real steps

      val () = Harness.section "orbit: separation stays ~ constant (near-circular)"

      val q1 = N.stepNDirect g dt (steps div 4) bs0
      val q2 = N.stepNDirect g dt (steps div 4) q1
      val q3 = N.stepNDirect g dt (steps div 4) q2
      val q4 = N.stepNDirect g dt (steps div 4) q3
      val () = checkApproxTol 2E~2 "separation at T/4 ~ 2" (2.0, sep q1)
      val () = checkApproxTol 2E~2 "separation at T/2 ~ 2" (2.0, sep q2)
      val () = checkApproxTol 2E~2 "separation at 3T/4 ~ 2" (2.0, sep q3)
      val () = checkApproxTol 2E~2 "separation at T ~ 2" (2.0, sep q4)

      val () = Harness.section "orbit: momentum and energy"

      (* total momentum starts at zero and is conserved to rounding *)
      val () = checkVecTol eps "initial total momentum = 0"
                 (V.zero, N.totalMomentum bs0)
      val () = checkVecTol 1E~9 "total momentum conserved over a period"
                 (V.zero, N.totalMomentum q4)

      (* KE = 2 * 1/2 m v^2 = 0.25 ; PE = -G m m / r_sep = -0.5 ; E = -0.25 *)
      val e0 = N.totalEnergy g bs0
      val () = checkApproxTol 1E~12 "initial energy = -0.25" (~0.25, e0)
      val () = checkApproxTol 1E~3 "energy conserved over a period"
                 (e0, N.totalEnergy g q4)

      val () = Harness.section "orbit: Barnes-Hut agrees with direct stepping"

      (* for two well-separated bodies theta=0 BH and the direct sum coincide *)
      val sD = N.stepDirect g dt bs0
      val sB = N.stepBH { g = g, theta = 0.0 } dt bs0
      val () = checkVecTol 1E~9 "stepBH(theta=0) == stepDirect, body 0 pos"
                 (#pos (nth (sD, 0)), #pos (nth (sB, 0)))
      val () = checkVecTol 1E~9 "stepBH(theta=0) == stepDirect, body 1 vel"
                 (#vel (nth (sD, 1)), #vel (nth (sB, 1)))

      (* a small theta over a quarter period stays close to the direct result *)
      val bhQ = N.stepNBH { g = g, theta = 0.3 } dt (steps div 4) bs0
      val () = checkApproxTol 2E~2 "BH(theta=0.3) keeps near-circular separation"
                 (2.0, sep bhQ)
    in
      ()
    end
end

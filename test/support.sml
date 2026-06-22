(* support.sml -- shared helpers for the sml-nbody tests.

   Every quantity is floating point, so the suite compares against analytic
   values through an explicit epsilon (`approx`) rather than string or
   structural equality: `Real.toString` differs between MLton and Poly/ML, and
   numerical methods only match the closed form up to truncation/rounding.  A
   tight `eps` (1e-12) pins exact algebraic identities; orbit / approximation
   checks pass their own, looser tolerances inline. *)

structure Support =
struct
  structure V = NBody.V

  val eps = 1E~12

  fun approx (a, b) = Real.abs (a - b) <= eps
  fun approxTol tol (a, b) = Real.abs (a - b) <= tol

  fun checkApprox name (expected, actual) =
    Harness.check name (approx (expected, actual))

  fun checkApproxTol tol name (expected, actual) =
    Harness.check name (approxTol tol (expected, actual))

  fun vApproxTol tol (a, b) =
    Real.abs (V.x a - V.x b) <= tol andalso Real.abs (V.y a - V.y b) <= tol

  fun checkVecTol tol name (expected, actual) =
    Harness.check name (vApproxTol tol (expected, actual))

  fun body (px, py, vx, vy, m) : NBody.body =
    { pos = V.v (px, py), vel = V.v (vx, vy), mass = m }
end

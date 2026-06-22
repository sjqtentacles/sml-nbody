(* test_quadtree.sml -- the quadtree aggregates mass and center of mass.

   A node's total mass is the sum of its bodies' masses and its center of mass
   is their mass-weighted average; these are pinned against hand calculation. *)

structure QuadtreeTests =
struct
  open Support
  structure N = NBody

  fun run () =
    let
      val () = Harness.section "quadtree: empty and singleton"

      val () = checkApprox "empty tree has zero mass" (0.0, N.treeMass (N.build []))
      val () = checkVecTol eps "empty tree com is zero"
                 (V.zero, N.treeCom (N.build []))

      val one = N.build [ body (3.0, ~4.0, 0.0, 0.0, 2.5) ]
      val () = checkApprox "singleton mass = body mass" (2.5, N.treeMass one)
      val () = checkVecTol eps "singleton com = body pos"
                 (V.v (3.0, ~4.0), N.treeCom one)

      val () = Harness.section "quadtree: mass and center of mass aggregate"

      (* four bodies, one per quadrant *)
      val bs =
        [ body (1.0,  1.0, 0.0, 0.0, 1.0)    (* NE *)
        , body (~1.0, 1.0, 0.0, 0.0, 1.0)    (* NW *)
        , body (~1.0, ~1.0, 0.0, 0.0, 2.0)   (* SW *)
        , body (1.0,  ~1.0, 0.0, 0.0, 4.0) ] (* SE *)
      val t = N.build bs
      val () = checkApprox "total mass = 8" (8.0, N.treeMass t)
      (* com_x = (1 - 1 - 2 + 4)/8 = 0.25 ; com_y = (1 + 1 - 2 - 4)/8 = -0.5 *)
      val () = checkVecTol 1E~12 "center of mass = (0.25, -0.5)"
                 (V.v (0.25, ~0.5), N.treeCom t)

      (* the standalone centerOfMass helper agrees with the tree *)
      val () = checkVecTol 1E~12 "centerOfMass helper agrees with tree"
                 (N.treeCom t, N.centerOfMass bs)

      val () = Harness.section "quadtree: deep subdivision still aggregates"

      (* tightly clustered (forces several levels of subdivision) plus a far one *)
      val cl =
        [ body (0.01, 0.01, 0.0, 0.0, 1.0)
        , body (0.02, 0.01, 0.0, 0.0, 1.0)
        , body (0.01, 0.02, 0.0, 0.0, 1.0)
        , body (5.0,  5.0,  0.0, 0.0, 1.0) ]
      val tc = N.build cl
      val () = checkApprox "clustered total mass = 4" (4.0, N.treeMass tc)
      val () = checkVecTol 1E~12 "clustered com matches helper"
                 (N.centerOfMass cl, N.treeCom tc)
    in
      ()
    end
end

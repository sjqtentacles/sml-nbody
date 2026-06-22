(* nbody.sml

   Implementation of NBODY on top of Glm.Vec2.  Pure and basis-only beyond the
   vendored sml-glm. *)

structure NBody :> NBODY =
struct
  structure V = struct
    type t = Glm.Vec2.t
    val v        = Glm.Vec2.v
    val x        = Glm.Vec2.x
    val y        = Glm.Vec2.y
    val zero     = Glm.Vec2.zero
    val add      = Glm.Vec2.add
    val sub      = Glm.Vec2.sub
    val scale    = Glm.Vec2.scale
    val dot      = Glm.Vec2.dot
    val length   = Glm.Vec2.length
    val lengthSq = Glm.Vec2.lengthSq
    val dist     = Glm.Vec2.dist
  end

  type body = { pos : V.t, vel : V.t, mass : real }

  (* A node carries its aggregate mass, center of mass, and side length
     (`size`, used for the opening-angle test). *)
  datatype tree =
      Empty
    | Leaf of { com : V.t, mass : real }
    | Node of { com : V.t, mass : real, size : real,
                nw : tree, ne : tree, sw : tree, se : tree }

  fun treeMass Empty = 0.0
    | treeMass (Leaf { mass, ... }) = mass
    | treeMass (Node { mass, ... }) = mass

  fun treeCom Empty = V.zero
    | treeCom (Leaf { com, ... }) = com
    | treeCom (Node { com, ... }) = com

  (* mass-weighted combination of two (com, mass) aggregates *)
  fun combine ((c1, m1), (c2, m2)) =
    let val m = m1 + m2
    in if Real.== (m, 0.0) then (V.zero, 0.0)
       else (V.scale (1.0 / m, V.add (V.scale (m1, c1), V.scale (m2, c2))), m)
    end

  (* Bounding square (center, half-side) covering all bodies, with a little pad
     so boundary points fall strictly inside. *)
  fun boundsOf [] = (V.zero, 1.0)
    | boundsOf (b :: bs) =
        let
          val p0 = #pos b
          fun go ([], minx, maxx, miny, maxy) = (minx, maxx, miny, maxy)
            | go (b' :: r, minx, maxx, miny, maxy) =
                let val p = #pos b'
                in go (r, Real.min (minx, V.x p), Real.max (maxx, V.x p),
                          Real.min (miny, V.y p), Real.max (maxy, V.y p)) end
          val (minx, maxx, miny, maxy) =
            go (bs, V.x p0, V.x p0, V.y p0, V.y p0)
          val cx = 0.5 * (minx + maxx)
          val cy = 0.5 * (miny + maxy)
          val span = Real.max (maxx - minx, maxy - miny)
          val half = 0.5 * span + 1.0E~9
          val half = if half <= 0.0 then 1.0 else half
        in (V.v (cx, cy), half) end

  val minHalf = 1.0E~9

  (* Recursively partition `bodies` within the square (center, half). *)
  fun buildIn (center, half) bodies =
    case bodies of
        [] => Empty
      | [b] => Leaf { com = #pos b, mass = #mass b }
      | _ =>
          if half < minHalf then
            (* coincident / sub-resolution cluster: aggregate into one leaf *)
            let
              val (com, mass) =
                List.foldl
                  (fn (b, acc) => combine (acc, (#pos b, #mass b)))
                  (V.zero, 0.0) bodies
            in Leaf { com = com, mass = mass } end
          else
            let
              val cx = V.x center and cy = V.y center
              val h2 = 0.5 * half
              fun quad b =
                let val p = #pos b
                in if V.x p >= cx
                   then (if V.y p >= cy then 1 else 3)   (* NE : SE *)
                   else (if V.y p >= cy then 0 else 2)   (* NW : SW *)
                end
              fun pick q = List.filter (fn b => quad b = q) bodies
              val nw = buildIn (V.v (cx - h2, cy + h2), h2) (pick 0)
              val ne = buildIn (V.v (cx + h2, cy + h2), h2) (pick 1)
              val sw = buildIn (V.v (cx - h2, cy - h2), h2) (pick 2)
              val se = buildIn (V.v (cx + h2, cy - h2), h2) (pick 3)
              val agg =
                List.foldl combine (V.zero, 0.0)
                  [ (treeCom nw, treeMass nw), (treeCom ne, treeMass ne)
                  , (treeCom sw, treeMass sw), (treeCom se, treeMass se) ]
              val (com, mass) = agg
            in Node { com = com, mass = mass, size = 2.0 * half
                    , nw = nw, ne = ne, sw = sw, se = se } end

  fun build bodies = buildIn (boundsOf bodies) bodies

  (* Acceleration at `point` from a point mass `m` at `com`. *)
  fun pull g (com, m) point =
    let
      val r = V.sub (com, point)
      val r2 = V.lengthSq r
    in
      if Real.== (r2, 0.0) then V.zero
      else V.scale (g * m / (r2 * Math.sqrt r2), r)
    end

  fun directAccel g bodies point =
    List.foldl
      (fn (b, acc) => V.add (acc, pull g (#pos b, #mass b) point))
      V.zero bodies

  fun bhAccel { g, theta } tree point =
    let
      fun acc Empty = V.zero
        | acc (Leaf { com, mass }) = pull g (com, mass) point
        | acc (Node { com, mass, size, nw, ne, sw, se }) =
            let val d = V.dist (point, com)
            in
              if d > 0.0 andalso size / d < theta
              then pull g (com, mass) point
              else V.add (V.add (acc nw, acc ne), V.add (acc sw, acc se))
            end
    in acc tree end

  (* ---- velocity Verlet (leapfrog) over all bodies ---- *)

  (* mkField builds, from the current bodies, the acceleration field a(point). *)
  fun vverlet mkField dt bodies =
    let
      val f0 = mkField bodies
      val a0 = List.map (fn b => f0 (#pos b)) bodies
      val bodies1 =
        ListPair.map
          (fn (b, a) =>
             { pos = V.add (#pos b,
                       V.add (V.scale (dt, #vel b), V.scale (0.5 * dt * dt, a)))
             , vel = #vel b, mass = #mass b })
          (bodies, a0)
      val f1 = mkField bodies1
      val a1 = List.map (fn b => f1 (#pos b)) bodies1
      fun finish (b1 :: bs1, a0i :: a0s, a1i :: a1s) =
            { pos = #pos b1
            , vel = V.add (#vel b1, V.scale (0.5 * dt, V.add (a0i, a1i)))
            , mass = #mass b1 } :: finish (bs1, a0s, a1s)
        | finish _ = []
    in finish (bodies1, a0, a1) end

  fun stepDirect g dt bodies =
    vverlet (fn bs => directAccel g bs) dt bodies

  fun stepBH { g, theta } dt bodies =
    vverlet (fn bs => bhAccel { g = g, theta = theta } (build bs)) dt bodies

  fun stepNDirect g dt n bodies =
    if n <= 0 then bodies else stepNDirect g dt (n - 1) (stepDirect g dt bodies)

  fun stepNBH params dt n bodies =
    if n <= 0 then bodies
    else stepNBH params dt (n - 1) (stepBH params dt bodies)

  (* ---- diagnostics ---- *)

  fun totalMomentum (bodies : body list) =
    List.foldl (fn (b, acc) => V.add (acc, V.scale (#mass b, #vel b)))
      V.zero bodies

  fun centerOfMass (bodies : body list) =
    let
      val (com, _) =
        List.foldl (fn (b, acc) => combine (acc, (#pos b, #mass b)))
          (V.zero, 0.0) bodies
    in com end

  fun kineticEnergy (bodies : body list) =
    List.foldl (fn (b, acc) => acc + 0.5 * #mass b * V.lengthSq (#vel b))
      0.0 bodies

  fun potentialEnergy g (bodies : body list) =
    let
      fun pairs ([], acc) = acc
        | pairs ((b : body) :: rest, acc) =
            let
              val s =
                List.foldl
                  (fn (b', a) =>
                     a - g * #mass b * #mass b' / V.dist (#pos b, #pos b'))
                  0.0 rest
            in pairs (rest, acc + s) end
    in pairs (bodies, 0.0) end

  fun totalEnergy g bodies = kineticEnergy bodies + potentialEnergy g bodies
end

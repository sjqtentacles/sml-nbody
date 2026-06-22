(* demo.sml

   A deterministic Barnes-Hut N-body run: a heavy central mass plus a ring of
   eight equal bodies started on near-circular orbits (positions and velocities
   from a fixed closed formula -- no RNG).  The demo integrates the system with
   the symplectic leapfrog stepper using the Barnes-Hut force approximation,
   prints an energy / momentum-vs-step table, and renders the final body
   positions as an ASCII scatter.  Output is byte-identical under MLton and
   Poly/ML.

   The full report is written to assets/nbody.txt and echoed to stdout.

   Build and run with `make example`. *)

structure N = NBody
structure V = NBody.V

fun fmt k x = Real.fmt (StringCvt.FIX (SOME k)) x
fun fmtD k x =
  let val s = fmt k x
  in if String.isPrefix "~" s then "-" ^ String.extract (s, 1, NONE) else s end

val buf : string list ref = ref []
fun line s = buf := (s ^ "\n") :: !buf

(* ---- the system: central mass + a ring of 8 bodies ---- *)
val g = 1.0
val centralM = 1000.0
val k = 8
val ringR = 5.0
val vcirc = Math.sqrt (g * centralM / ringR)

fun mk i =
  let
    val ang = 2.0 * Math.pi * real i / real k
    val px = ringR * Math.cos ang
    val py = ringR * Math.sin ang
    val vx = ~vcirc * Math.sin ang     (* tangential, counter-clockwise *)
    val vy =  vcirc * Math.cos ang
  in { pos = V.v (px, py), vel = V.v (vx, vy), mass = 1.0 } end

val bodies0 : N.body list =
  { pos = V.zero, vel = V.zero, mass = centralM } :: List.tabulate (k, mk)

val theta = 0.5
val dt = 0.0005

val () = line "=== sml-nbody demo =========================================="
val () = line ""
val () = line ("Central mass " ^ fmt 1 centralM ^ " at origin, plus a ring of "
               ^ Int.toString k ^ " unit masses")
val () = line ("at radius " ^ fmt 1 ringR ^ " on circular orbits (v = "
               ^ fmt 4 vcirc ^ ").")
val () = line ("Barnes-Hut theta = " ^ fmt 2 theta ^ ", leapfrog dt = "
               ^ fmt 4 dt ^ ".")
val () = line ""

(* ---- energy / momentum table ---- *)
val () = line "   step       t          KE            PE           E       |P|"
fun report (stepNo, bs) =
  let
    val t = real stepNo * dt
    val ke = N.kineticEnergy bs
    val pe = N.potentialEnergy g bs
    val e = ke + pe
    val p = V.length (N.totalMomentum bs)
  in
    line ("  " ^ StringCvt.padLeft #" " 5 (Int.toString stepNo)
          ^ "  " ^ StringCvt.padLeft #" " 6 (fmtD 3 t)
          ^ "  " ^ StringCvt.padLeft #" " 11 (fmtD 4 ke)
          ^ "  " ^ StringCvt.padLeft #" " 12 (fmtD 4 pe)
          ^ "  " ^ StringCvt.padLeft #" " 10 (fmtD 4 e)
          ^ "  " ^ StringCvt.padLeft #" " 6 (fmtD 4 p))
  end

val blocks = 8
val perBlock = 250
fun runTable (stepNo, bs, j) =
  if j > blocks then bs
  else
    ( report (stepNo, bs)
    ; if j = blocks then bs
      else runTable (stepNo + perBlock,
                     N.stepNBH { g = g, theta = theta } dt perBlock bs, j + 1) )
val finalBodies = runTable (0, bodies0, 0)
val () = line ""

(* ---- ASCII scatter of the final positions ---- *)
val Wc = 41
val Hc = 21
val halfWin = 8.0
val grid = Array.array (Wc * Hc, #" ")
fun plot (ch, p) =
  let
    val col = Real.round ((V.x p + halfWin) / (2.0 * halfWin) * real (Wc - 1))
    val row = Real.round ((halfWin - V.y p) / (2.0 * halfWin) * real (Hc - 1))
  in
    if col >= 0 andalso col < Wc andalso row >= 0 andalso row < Hc
    then Array.update (grid, row * Wc + col, ch) else ()
  end

val () =
  let
    fun go (_, []) = ()
      | go (i, (b : N.body) :: rest) =
          (plot (if i = 0 then #"O" else #"*", #pos b); go (i + 1, rest))
  in go (0, finalBodies) end

val () = line ("Final positions after " ^ Int.toString (blocks * perBlock)
               ^ " steps  (O = central, * = ring body):")
val () = line ("  +" ^ CharVector.tabulate (Wc, fn _ => #"-") ^ "+")
val () =
  let
    fun rows r =
      if r >= Hc then ()
      else
        let
          val cs = CharVector.tabulate (Wc, fn c => Array.sub (grid, r * Wc + c))
        in line ("  |" ^ cs ^ "|"); rows (r + 1) end
  in rows 0 end
val () = line ("  +" ^ CharVector.tabulate (Wc, fn _ => #"-") ^ "+")
val () = line ""
val () = line "============================================================"

val report = String.concat (List.rev (!buf))
val () =
  let val os = TextIO.openOut "assets/nbody.txt"
  in TextIO.output (os, report); TextIO.closeOut os end
val () = print report

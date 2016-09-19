(ns fv.coeff)
(use 'clojure.core.matrix)

(def τ 6.28318530718)
(defn fvHelix2P4
  "
 -- calculate track momentum 4-vector and 4x4 error matrix
 -- from helix parameters parameters h and Ch
 -- use w2pt for calculating pt from curvature w
 -- calculate energy assuming mass m
  "
  [h Ch m w2pt]
  (let [
        [w tl psi0 _ _]  (vec h)
        sph   (Math/sin psi0)
        cph   (Math/cos psi0)
        pt    (/ w2pt (Math/abs w))
        px    (* pt cph)
        py    (* pt sph)
        pz    (* pt tl)
        E     (Math/sqrt (+ (* m m) (* px px) (* py py) (* pz pz)))
        ps    (/ w2pt w)
        dpdk  (/ (* ps ps) w2pt)
        xy    (* 2.0 ps dpdk cph sph (mget Ch 0 2))
        sxx   (+ (* dpdk dpdk cph cph (mget Ch 0 0))
                 (* ps ps sph sph (mget Ch 2 2))
                 xy)
        sxy   (+ (* cph sph (- (* dpdk dpdk (mget Ch 0 0))
                               (* ps ps (mget Ch 2 2))))
                 (* ps dpdk (- (* sph sph) (* cph cph)) (mget Ch 0 2)))
        syy   (+ (* dpdk dpdk sph sph (mget Ch 0 0))
                 (* ps ps cph cph (mget Ch 2 2))
                 (- xy))
        sxz   (- (* dpdk dpdk cph tl (mget Ch 0 0))
                 (* ps dpdk (- (* cph (mget Ch 0 1))
                               (* sph tl (mget Ch 0 2))))
                 (* ps ps sph (mget Ch 1 2)))
        syz   (+ (* dpdk dpdk sph tl (mget Ch 0 0))
                 (- (* ps dpdk (+ (* sph (mget Ch 0 1))
                                  (* cph tl (mget Ch 0 2)))))
                 (* ps ps cph (mget Ch 1 2)))
        szz   (+  (* dpdk dpdk tl tl (mget Ch 0 0))
                 (* ps ps (mget Ch 1 1))
                 (- (* 2.0 ps dpdk tl (mget Ch 0 1))))
        sxe   (/ (+ (* px sxx)
                    (* py sxy)
                    (* pz sxz))
                 E)
        sye   (/ (+ (* px sxy)
                    (* py syy)
                    (* pz syz))
                 E)
        sze   (/ (+ (* px sxz)
                    (* py syz)
                    (* pz szz))
                 E)
        see   (/ (+ (* px px sxx)
                    (* py py syy)
                    (* pz pz szz)
                    (* 2.0 (+ (* px py sxy)
                              (* px pz sxz)
                              (* py pz syz))))
                 (* E E))

        p     (array [px py pz E])
        Cp    (matrix [[sxx sxy sxz sxe] [sxy syy syz sye] [sxz syz szz sze]
                       [sxe sye sze see]])
        ] [p Cp]))

(defn fvQ2P3
  "
 -- calculate track momentum 3-vector its error matrix
 -- from helix 3-vector q = [w tl psi0] and covariance matrix Q
 -- use w2pt for calculating pt from curvature w
  "
  [q Q w2pt]
  (let [
        [w tl psi0]  (vec q)
        sph   (Math/sin psi0)
        cph   (Math/cos psi0)
        pt    (/ w2pt (Math/abs w))
        px    (* pt cph)
        py    (* pt sph)
        pz    (* pt tl)
        ps    (/ w2pt w)
        dpdk  (/ (* ps ps) w2pt)
        xy    (* 2.0 ps dpdk cph sph (mget Q 0 2))
        sxx   (+ (* dpdk dpdk cph cph (mget Q 0 0))
                 (* ps ps sph sph (mget Q 2 2))
                 xy)
        sxy   (+ (* cph sph (- (* dpdk dpdk (mget Q 0 0))
                               (* ps ps (mget Q 2 2))))
                 (* ps dpdk (- (* sph sph) (* cph cph)) (mget Q 0 2)))
        syy   (+ (* dpdk dpdk sph sph (mget Q 0 0))
                 (* ps ps cph cph (mget Q 2 2))
                 (- xy))
        sxz   (- (* dpdk dpdk cph tl (mget Q 0 0))
                 (* ps dpdk (- (* cph (mget Q 0 1))
                               (* sph tl (mget Q 0 2))))
                 (* ps ps sph (mget Q 1 2)))
        syz    (+ (* dpdk dpdk sph tl (mget Q 0 0))
                  (- (* ps dpdk (+ (* sph (mget Q 0 1))
                                   (* cph tl (mget Q 0 2)))))
                  (* ps ps cph (mget Q 1 2)))
        szz    (+  (* dpdk dpdk tl tl (mget Q 0 0))
                  (* ps ps (mget Q 1 1))
                  (- (* 2.0 ps dpdk tl (mget Q 0 1))))
        p      (array [px py pz])
        Cp     (matrix [[sxx sxy sxz] [sxy syy syz] [sxz syz szz]])
        ] [p Cp]))

(defn fvq
  "
   -- calculate 3-vector q at a vetex position v (currently not used)
   -- corresponding to the 5-vector h for a track helix
   -- the more general case would calculate q taking into account
   -- a displaced vertex position
  "
  [h v] (subvector h 0 3))


(defn fvQ
  "
   -- calculate 3x3 cov matrix Q for a q vector,
   -- corresponding to the 5x5 cov matrix H
   -- for a helix parameter vector h
  "
  [H] (submatrix H 0 3 0 3))

(defn fvh
  "
  -- calculate and return helix parameters h = [w tl psi0 d0 z0]
  -- from v = [vx vy vz] and q = [w tl psi]
  "
  [v q]
  (let [
        [xx yy z]    (vec v)
        r            (sqrt (+ (* xx xx) (* yy yy)))
        phi    (Math/atan2 yy xx)

        [w tl psi]    (vec q)
        oow  (/ 1.0 w)
        xi   (- psi phi)
        cxi  (Math/cos xi)
        sxi  (Math/sin xi)
        gamma  (Math/atan (/ (* r cxi) (- oow (* r sxi))))
        h    (if (zero? w)
               (array [w
                       tl
                       (- psi gamma)
                       (- oow (/ (- oow (* r sxi)) (Math/cos gamma)))
                       (-  z  (* (/ gamma w) tl))
                       ])
               (array [w, tl, psi, (* r sxi), z])) ;;?? check this out
        ] h))

(defn fvABh0
  "
--
-- calculate coefficients of taylor expansion for estimated
-- track parameters p around point (v0,q0):
--
--   p    =    h(v,q) + eps
--        ≈    h(v0,q0) + A.(v-v0) + B.(q-q0) + eps
--        =    h0 + A.v + B.q + eps
-- where
--   A  =  dp/dv(v0,q0),
--   B  =  dp/dq(v0,q0), and
--   h0 =  h(v0,q0) - A.v0 - B.q0 are the
-- derivatives of track parameters w/r to vertex position
-- and track 3-momentum at space point v0 and estimated momentum q0
--
-- called with v0, q0 and returns A, B, h0
  "
  [v0 q0]
  (let [
    [xx yy z]   (vec v0)
    r           (Math/sqrt (+ (* xx xx) (* yy yy)))
    phi         (Math/atan2 yy xx)
    [w tl psi]  (vec q0)
    ;; some more derived quantities
    xi          (mod (- psi phi) τ)
    cxi         (Math/cos xi)
    sxi         (Math/sin xi)
    oow         (/ 1.0 w)
    rw          (* r w)

    gamma       (Math/atan (/ (* r cxi) (- oow (* r sxi))))
    sg          (Math/sin gamma)
    cg          (Math/cos gamma)

    ;; calculate transformed quantities
    psi0        (- psi gamma)
    d0          (- oow (/ (- oow (* r sxi)) cg))
    z0          (- z (* (/ gamma w) tl))

    ;; calc Jacobian
    [drdx drdy rdxidx rdxidy]   (if (not= 0 r)
                                   [(/ xx r) (/ yy r) (/ yy r) (- (/ xx r))]
                                   [0.0 0.0 0.0 0.0])
    dgdvar0     (/ 1.0 (+ 1.0 (* rw rw) (* -2.0 rw sxi)))
    dgdx        (* dgdvar0 (+ (* w cxi drdx) (* w (- rw sxi) rdxidx)))
    dgdy        (* dgdvar0 (+ (* w cxi drdy) (* w (- rw sxi) rdxidy)))
    dgdw        (* dgdvar0 r cxi)
    dgdpsi      (* dgdvar0 rw (- rw sxi))

    ;; fill matrix:
    ;; -- d w / d r, d phi, d z
    [ a11 a12 a13 ]     [ 0 0 0 ]
    ;; -- d tl / d x, d y, d z
    [ a21 a22 a23 ]     [ 0 0 0 ]
    ;; -- d psi0 / d x, d y, d z
    [ a31 a32 a33 ]     [ (- dgdx) (- dgdy) 0 ]
    ;; -- d d0 / d x, d y, d z
    [ a41 a42 a43 ]     [ (- (+ (/ (* cxi rdxidx) cg) (/ (* sxi drdx) cg))
                             (/ (* (- oow (* r sxi)) sg dgdx) cg cg))
                          (- (+ (/ (* cxi rdxidy) cg) (/ (* sxi drdy) cg))
                             (/ (* (- oow (* r sxi)) sg dgdy) cg cg))
                          0
                          ]
    ;; -- d z0 / d x, d y, d z
    [ a51 a52 a53 ]     [ (- (* (/ tl w) dgdx)) (- (* (/ tl w) dgdy)) 1.0 ]
    ;; B
    ;; -- d w / d w, d tl, d psi
    [ b11 b12 b13 ]    [ 1.0 0 0 ]
    ;; -- d tl / d w, d tl, d psi
    [ b21 b22 b23 ]    [ 0 1.0 0 ]
    ;; -- d psi0 / d w, d tl, d psi
    [ b31 b32 b33 ]    [ (- dgdw) 0 (- 1.0 dgdpsi) ]
    ;; -- d d0 / d w, d tl, d psi
    [ b41 b42 b43 ]    [ (- (+ (* oow oow (- 1.0 (/ 1.0 cg)))
                               (/ (* (- oow (* r sxi)) sg dgdw) cg cg)))
                        0
                        (- (* r (/ cxi cg))
                           (/ (* (- oow (* r sxi)) sg dgdpsi) cg cg)) ]
    ;; -- d z0 / d w, d tl, d psi
    [ b51 b52 b53 ]    [ (* (/ tl w) (- (/ gamma w) dgdw))
                         (- (/ gamma w))
                         (- (* (/ tl w) dgdpsi))
                       ]

    [ v01 v02 v03 ]       (vec v0)
    [ q01 q02 q03 ]       (vec q0)
    [ h01 h02 h03 h04 h05 ]    [ 0 0
                 (- psi0 (* a31 v01) (* a32 v02) (* b31 q01) (* b33 q03))
                 (- d0 (* a41 v01) (* a42 v02) (* b41 q01) (* b43 q03))
                 (- z0 (* a51 v01) (* a52 v02) (* a53 v03) (* b51 q01)
                       (* b52 q02) (* b53 q03)) ]
    A      (matrix [[a11 a12 a13] [a21 a22 a23] [a31 a32 a33]
                    [a41 a42 a43] [a51 a52 a53]])
    B      (matrix [[b11 b12 b13] [b21 b22 b23] [b31 b32 b33]
                    [b41 b42 b43] [b51 b52 b53]])
    h0     [h01 h02 h03 h04 h05]

    R      (vector A B h0)

    ] R))


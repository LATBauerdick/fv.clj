(ns fv.fv
  (:require [fv.coeff :refer :all]
            [clojure.core.matrix :refer :all]))

(set-current-implementation :vectorz)

;; a record if Helices contains a vertex, a 1x3 matrix of coord x,y,z
;;  and a list of tracks hl, 1x5 helices of helix params w, tl, psi0, d0, z0
;; together with the corresponding list of covariance matrices Hl
(defrecord Helices [v0 V0 hl Hl])
;; a Track has a helix h and its cov matrix H
(defrecord Track [h H])

;; A Prong record contains info about a n-prong tuple of 
;; a vertex v, a 1x3 matrix of coords x,y,z with its cov matrix V
;; a list of momentum vectors qQl coming from this vertex,
;; and a list of chi2 values chi2l for those momenta fitting the vertex
(defrecord Prong [v V qQl chi2l t#l])
;; A Momentum record has a 3-vector q given in [w, tl, psi]
;; together with the corresponding 3x3 covariance matrix Q
(defrecord Momentum [q Q])

(def fvLog false)
(def never false)

(defn fvPMerr "pretty-print vector and error" [s p P]
  (let [e    (-> P diagonal sqrt )]
    (print s "-->")
    (doseq [[x y] (map list p e)] (print (format "%9.3g ±%9.3g" x y)))
    (println )))
;;  (print s " ") (pm p)
;;  (print "±% ") (pm (-> P diagonal sqrt (div p) (scale 100.0) abs )))

(defn fvInverse ;;;;;LATB probably needs SVD inversion
  [M]
  (if-let
    [M-1 (inverse M)] ;;returns nil if can't invert, probably should use fvSVDinv
    (do (if (> (emax (div (sub (mmul M M-1)
                               (identity-matrix (row-count M)))
                          M)) 0.002)
          (if fvLog
            '(pm (scale (div (sub (mmul M M-1)
                                 (identity-matrix (row-count M))) M) 100.0))
            (print "🚩")))
        M-1)
    (do (println "fvInverse cannot invert"
                 (row-count M) "X" (row-count M) "matrix")
        (identity-matrix (row-count M)))))
(defn- fvAB [A B] (mmul A B))
(defn- fvsABAT [A B] (mmul A B (transpose A)))
(defn- fvsATBA [A B] (mmul (transpose A) B A))
(defn- fvAPB [A B] (add A B))
(defn- fvAMB [A B] (sub A B))
(defn- fvATB [A B] (mmul (transpose A) B))
(defn- fvNegA [A] (negate A))
(defn- fvATBC [A B C] (mmul (transpose A) B C))
(defn- fvABCT [A B C] (mmul A B (transpose C)))
(defn- fvABTCT [A B C] (mmul A (transpose B) (transpose C)))

(def 𝜒2cut 0.5)
(defn- goodEnough? [𝜒20 𝜒2] 
  (-> 𝜒2 (- 𝜒20) abs (< 𝜒2cut)))

(defn fvFilterer [v0 U0 h H] ;; we start with the "fvFilterer" implementation
  (let [
        [A B h0]     (fvABh0 v0 (fvq h v0))
        G            (fvInverse H)
        ]
    (loop [v0 v0, U0 U0, A A, B B, h0 h0, 𝜒20 1e10, iter 0]
      (let [
            ;; -- W = (B^T.G.B)^(-1)
            W            (fvInverse (fvsATBA B G))
            ;; -- GB = G - G.B.W.B^T.G^T
;;            GB           (fvAMB G (fvsABAT G (fvsABAT B W)))
            GB           (fvAMB G (fvsABAT (fvAB G B) W))
            ;; -- C = U^-1 = (U0 + A^T.GB.A)^-1
            U            (fvAPB U0 (fvsATBA A GB))
            C            (fvInverse U)   ;; check for editing of singular values?
            ;; -- m = h - h0
            m      (fvAMB h h0)
            ;; -- v = C. (U0.v0 + A^T.GB.(h-h0) )
            v      (fvAB C (fvAPB (fvAB U0 v0) (fvATB A (fvAB GB m))))
            ;; -- dm = h - h0 - A.v
            dm     (fvAMB m (fvAB A v))
            ;; -- q = W.BT.(G.dm)
            q      (fvAB  W  (fvATB  B  (fvAB  G  dm)))
            ;; -- D
            D      (fvAPB W (fvsATBA W (fvsATBA B (fvsATBA G (fvsABAT A C)))))
            ;; -- E
            E     (fvNegA  (fvATBC  W  (fvATBC  B  G  A)  C))
            ;; -- chi2 - (dm - B.q)T.G. (dm - B.q) + (v - v0)T.U0. (v-v0)
            𝜒2    (+ (scalar  (fvsATBA  (fvAMB  dm  (fvAB  B  q))  G))
                     (scalar  (fvsATBA  (fvAMB  v  v0)  U0)))
            ]
        (when-not fvLog (print ".")) ;; progress bar
        (when fvLog
          (print "Filter iteration " iter "yields chi2 " (format "%9.3g" 𝜒2))
          (fvPMerr " at x,y,z," v C))
        (if-not (goodEnough? 𝜒20 𝜒2)
          (do
            (let [ [A1 B1 h01]     (fvABh0 v (fvq h v))] ;; recalc derivs at v
              (recur v0 U0 A1 B1 h01 𝜒2 (inc iter))))
          [v U q (fvQ H) 𝜒2])))))

(defn fvSmoother [v U h H]
  (let [
          [A B h0]     (fvABh0 v (fvq h v))
;; -- m = h - h0
          m    (fvAMB h h0)
;; -- dm = h - h0 - A.v
          dm   (fvAMB m (fvAB A v))
          G    (fvInverse H)
;; -- W = (B^T.G.B)^(-1)
          W    (fvInverse (fvsATBA B G))
;; -- q = W.BT.(G.dm)
          q    (fvAB  W  (fvATB  B  (fvAB  G  dm)))
;; -- chi2 = (h - h(v,q))^T.Gh(v,q).(h - h(v,q))
;; -- where we calculate Gh(v,q) from the fit-result covariance
;; -- matrices C,D,E for the smoothed results v,q
;;
          C    (fvInverse U)     ;; check for editing of singular values?
          D    (fvAPB W (fvsATBA W (fvsATBA B (fvsATBA G (fvsABAT A C)))))
          E    (fvNegA  (fvATBC  W  (fvATBC  B  G  A)  C))
;;           Gh    (fvInverse (fvCh v q C D E))
;;  -- Ch = A.C.A^T + B.E.A^T + A.E^T.B^T + B.D.B^T
          Ch   (fvAPB (fvAPB (fvAPB (fvsABAT A C) (fvABCT B E A))
                             (fvABTCT A E B))
                       (fvsABAT B D))
;;          Gh   (fvInverse Ch)
          Gh   G  ;; using simpler method for the moment
          𝜒2   (scalar (fvsATBA (fvAMB h (fvh v q)) Gh))
          ]
          [ q D 𝜒2]))



(defn fvFilter
  "
  -- run kalman filter function ƒ recursively over each helix in list hl
  "
  [ƒ v0 U0 hl Hl]
  (loop
    [v0 v0, U0 U0, hl hl, Hl Hl, ql [], Ql [], ih 0]
    (if (empty? hl)
      (do (when-not fvLog (println))
          [v0 U0 ql Ql])   ;; return list of q vectors and final v
      (let [h (first hl) H (first Hl)
            [v U q Q 𝜒2] (ƒ v0 U0 h H)
            ]
        (when-not fvLog (print (str "h" ih)))
        (when never
          (println "Filter for track " ih " -------------------------"
                   ", final chi2 " (format "%9.3g" 𝜒2))
          (fvPMerr "v0" v0 (fvInverse U0))
          (fvPMerr "v " v (fvInverse U))
          (fvPMerr "h " h H)
          (fvPMerr "q " q Q)
          (fvPMerr "dq" (sub (fvq h v0) q) Q))
        (recur v U (next hl) (next Hl) (conj ql q) (conj Ql Q) (inc ih))))))

(defn fvSmooth [ƒ v U hl Hl]
  (loop
    [hl hl, Hl Hl, qQl [], ql [], Ql [], 𝜒2l [], ih 0]
    (if (empty? hl)
      [qQl ql Ql 𝜒2l]   ;; return list of q vectors and 𝜒2
      (let [h0 (first hl) H0 (first Hl)
            qQ        (ƒ v U h0 H0)
            [q Q 𝜒2]  qQ
            ]
        (when fvLog
          (println "Smoother for track" ih
                   ", 𝜒2 " (format "%9.3g" 𝜒2))
          '(fvPMerr "h" h0 H0)
          '(fvPMerr "q" q Q)
          )
        (recur (next hl) (next Hl) (conj qQl (->Momentum q Q))
               (conj ql q) (conj Ql Q) (conj 𝜒2l 𝜒2) (inc ih))))))

(defn fvFit
  "
 -- calculate vertex x and its covariance matrix Cx (x = {x,y,z}),
 -- the list of momentum vectors ql and their covariance matrices
 -- Cql (q = {w, tl, psi}) and the total chi-square chi2t
 -- for nt tracks with numbers in tList
 -- from initial value for vertex x0, Cx0 and
 -- from list of track parameters hl and list of
 -- inverse of covariance matrices of track parameters Ghl
 -- 3-momenta q are returned in {w, tl, psi}-system
 --   (we do not know the transformation of w->pt at this point!)
 --
 -- input:
 --     x0({1..3})     : initial guess of vertex position x0{x,y,z}
 --     Cx0({1..3},{1..3}) : covariance matrix of x0
 --     hl({1..5}, tList(i)) : track parameters
 --                      h={w,tl,psi0,d0,z0} of track i
 --     Ghl({1..5},{1..5}, tList(i)) : inverse of covariance matrix
 --                      (cov(h(i)))^(-1) of h of track i
 -- output:
 --     x({1..3})              <- v({x,y,z}), coords of fitted vertx
 --     Cx({1..3})             <- cov{v} covariance matrix of v
 --     ql({1..3}, tList(i))   <- q({w,tl,psi})
 --     Cql({1..3},{1..3},tList(i)) <- covariance matrix of q
 --     chi2l(tList(i))        <- chi2 for this track belonging to v

 -- Notation for helix paramter 5-vector h
 -- w:  omega, 1/R curvature >0 if helix turns anti-clockwise
 -- tl: tan(lambda), tangens of dip angle
 -- psi: azimuth angle (in r-phi plane)
 -- psi0: psi @ point of closest approach to coordinate origin
 -- d0:   distance of helix at point of closest approach to origin in r-phi plane
 --       sign convention of angular momentum Lz
 -- z0:   z coordinate of point of closest approach of helix to origin
  "
  [hel]
  (let [
        U0           (fvInverse (:V0 hel))
        [v U _ _]    (fvFilter fvFilterer (:v0 hel) U0 (:hl hel) (:Hl hel))
;;        v  (array [ 0.954     0.994      3.55])
;;        U (fvInverse (matrix [[0.469E-01 0.504E-01 0.826E-01]
;;                              [0.504E-01 0.542E-01 0.888E-01]
;;                              [0.826E-01 0.888E-01 0.150]]))
        V            (fvInverse U)
        [qQl ql Ql 𝜒2l]  (fvSmooth fvSmoother v U (:hl hel) (:Hl hel))
      ]
      (->Prong v V qQl 𝜒2l (range 0 (count 𝜒2l)))))

(defn fvRetlif [v V q h H]
  (let [
        U  (fvInverse V)
        G  (fvInverse H)
        [A B h0]     (fvABh0 v q)
 ;; -- W = (B^T.G.B)^(-1)
        W     (fvInverse (fvsATBA B G))
 ;; -- GB = G - G.B.W.B^T.G^T
        GB    (fvAMB G (fvsABAT (fvAB G B) W))
 ;; -- Gvp = U - A^T.GB.A
        Gvp    (fvAMB U (fvsATBA A GB))
 ;; -- Cp = Gvp^(-1)
        Cp     (fvInverse Gvp)
 ;; -- m = p0 - h0
        m   (fvAMB h h0)
 ;; -- vp = Cp.(U.v - A^T.GB.(p0-h0))
        vp  (fvAB Cp (fvAMB (fvAB U v) (fvATBC A GB m)))
 ;; -- Fruehwirth CERN 90-06 says:
 ;;-- distance of track from the new vertex is expressed by the chi-square
 ;; -- of the smoothed residuals
 ;; -- chi2 = (p0-h0-A.v-B.q)^T.G.(p0-h0-A.v-B.q) + (v-vp)^T.Gvp.(v-vp)
        chi2 (+ (scalar (fvsATBA (fvAMB (fvAMB m (fvAB  A v)) (fvAB B q)) G))
                (scalar (fvsATBA (fvAMB v vp) Gvp)))

        ]chi2)  )

(defn fvRemove [prong hel]
  ;; currently just calculates chi2s and prints them
  (let [
        chi2l (map #(fvRetlif (:v prong) (:V prong) (:q %1) %2 %3)
                   (:qQl prong) (:hl hel) (:Hl hel))
        ]
    (print "fvRemove: list of chi2s")
    (doall (map #(print (format "%9.3g" %)) chi2l))
    (println)
    prong))

(defn mass [pP]
  (let [
        sqr           (fn [x] (* x x))
        [px py pz e]  (vec (pP 0))
        P             (pP 1)
        m (sqrt (- (sqr e) (sqr px) (sqr py) (sqr pz)))
        sigm   (+ (* px px (mget P 0 0))
                  (* py py (mget P 1 1))
                  (* pz pz (mget P 2 2))
                  (* e e (mget P 3 3))
                  (* 2.0 (+ (* px (+ (* py (mget P 0 1))
                                     (* pz (mget P 0 2))
                                     (- (* e (mget P 0 3)))))
                            (* py (- (* pz (mget P 1 2))
                                     (* e (mget P 1 3))))
                            (- (* pz (mget P 2 3) e)))))
        ] [ m (/ (sqrt sigm) m) ]))

(defn- addp [aA pP]
  (let [
        a  (aA 0)
        p  (pP 0)
        px (+ (mget a 0) (mget p 0))
        py (+ (mget a 1) (mget p 1))
        pz (+ (mget a 2) (mget p 2))
        e  (+ (mget a 3) (mget p 3))
        A  (fvAPB (aA 1) (pP 1))  ;; sum up the covariance matrices
        ]
     [ [px py pz e] A]))
(defn invMass [prong]
  (let [
        pl (map fvQ2P4 (:qQl prong))
        printtt '(println pl)
        ptot  (reduce addp pl)
        ]
    (mass ptot) ))

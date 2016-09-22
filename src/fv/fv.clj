(ns fv.fv
  (:require [fv.coeff :refer :all]
            [clojure.core.matrix :refer :all]))

(set-current-implementation :vectorz)

(def fvLog false)

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

(defn fvFilterer [v0 U0 h H] ;; we start with the "fvFilterer" implementation
  (let [
          [A B h0]     (fvABh0 v0 (fvq h v0))
          G            (fvInverse H)
;; -- W = (B^T.G.B)^(-1)
          W            (fvInverse (fvsATBA B G))
;; -- GB = G - G.B.W.B^T.G^T
          GB           (fvAMB G (fvsABAT G (fvsABAT B W)));; or (G.B)W (G.B)T (fvABAT (fvAB G B) W)
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
           printxx      (do (println "------------------------------------fvFilterer-----------------------")
                            (print "v0-->") (pm v0)
                            (print "U0-->") (pm U0)
                            (pm A) (pm B) (pm h0) (pm W) (pm GB) (pm dm)
                            (pm D) (pm E)
                            (print "v-->") (pm v) (print "C-->") (pm C)
                            (print "chi2" 𝜒2 )
                            )
          ]
          [ v U q (fvQ (fvInverse G)) 𝜒2]))

(defn fvSVDfit [v0 U0 h G A B h0] (
      let [v v0
           V (fvInverse U0)
           q (fvq h v)
           Q (fvQ (fvInverse G))]
           [ v V U0 q Q ]))

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
;;;;;;;;          printtt (fvPMerr (sub q (fvq h v)))

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
          Gh   (fvInverse Ch)
;;          Gh   G  ;; using simpler method for the moment
          𝜒2   (scalar (fvsATBA (fvAMB h (fvh v q)) Gh))

          ]
          [ q (fvQ Ch) 𝜒2]))



(def 𝜒2cut 0.1)
(defn- goodEnough? [𝜒20 𝜒2] (-> 𝜒2 (- 𝜒20) (/ 𝜒2)  (#(* % %)) (< 𝜒2cut)))
(defn fvFilter
  "
  -- run kalman filter function ƒ recursively over each helix in list hl
  "
  [ƒ v0 U0 hl Hl]
  (loop
    [v0 v0, U0 U0, hl hl, Hl Hl, ql [], Ql [], 𝜒20 0, ih 0, iter 0]
    (if (empty? hl)
      (do (println) [v0 U0 ql Ql])   ;; return list of q vectors and final v
      (let [h (first hl) H (first Hl)
            [v U q Q 𝜒2] (ƒ v0 U0 h H)
            ]
        (when (zero? iter)
          (print (str "h" ih))) (print ".") ;; progress bar
        (when (and fvLog (zero? iter))
          (println "Filter for track " ih " -------------------------"))
        (when fvLog (println "Filter iteration " iter "yields chi2 "
                             (format "%9.3g" 𝜒2)))
        (if (goodEnough? 𝜒20 𝜒2) ;; if diff in 𝜒2 is large, recur with same h, H
          (do
            (when fvLog
              ;; print result of ƒ
              (println "Filter v0 and h result in v and q for track " ih 
                       " chi2 " (format "%9.3g" 𝜒2))
              (fvPMerr "v0" v0 (fvInverse U0))
              (fvPMerr "v " v (fvInverse U))
              (fvPMerr "h " h H)
              (fvPMerr "q " q Q)
              (fvPMerr "dq" (sub (fvq h v0) q) Q))
              (recur v U (next hl) (next Hl) (conj ql q) (conj Ql Q) 𝜒2 (inc ih) 0))
            (recur v U hl Hl ql Ql 𝜒2 ih (inc iter)))))))

(defn fvSmooth [ƒ v U hl Hl]
  (loop
    [hl hl, Hl Hl, ql [], Ql [], 𝜒2l []]
    (if (empty? hl)
      (do (println) [ql Ql 𝜒2l])   ;; return list of q vectors and 𝜒2
      (let [h0 (first hl) H0 (first Hl)
            [q Q 𝜒2] (ƒ v U h0 H0)
            ]
        (when fvLog
          (do ;; print result of ƒ
            (println)
            (fvPMerr "q" q Q)
            (println "𝜒2 " 𝜒2)))
        (recur (next hl) (next Hl)
               (conj ql q) (conj Ql Q) (conj 𝜒2l 𝜒2))))))


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
 --     chi2t                  <- total chi2 for the vertex fit

 -- Notation for helix paramter 5-vector h
 -- w:  omega, 1/R curvature >0 if helix turns anti-clockwise
 -- tl: tan(lambda), tangens of dip angle
 -- psi: azimuth angle (in r-phi plane)
 -- psi0: psi @ point of closest approach to coordinate origin
 -- d0:   distance of helix at point of closest approach to origin in r-phi plane
 --       sign convention of angular momentum Lz
 -- z0:   z coordinate of point of closest approach of helix to origin
  "
  [v0 V0 hl Hl]
  (let [
        U0           (fvInverse V0)
        [v U _ _]    (fvFilter fvFilterer v0 U0 hl Hl)
        [ql Ql 𝜒2l]  (fvSmooth fvSmoother v U hl Hl)
        V            (fvInverse U)
        𝜒2  (reduce + 𝜒2l) ;;LATB not sure about that...
      ]
      [v V ql Ql 𝜒2l 𝜒2]
   ))



;;(defrecord fvrec [v Cv 𝜒2v qs Cqs 𝜒2qs])
;;(->fvrec  tCx tnt (vec th) (vec tCh) tw2pt); add to rec

(ns fv.fvt)
(use 'clojure.core.matrix)
(set-current-implementation :vectorz)

(require 'fv.fv)
(refer 'fv.fv)


(def thisFile "dat/tr05129e001412.dat")
(defrecord fvtData [tx tCx tnt th tCh tw2pt])

(def mπ  0.1395675e0)


(defn- indexed [coll]  (map-indexed vector coll))

(defn fvtread
  "
  -- read from file with name df and return
  -- 
  "
  [df]
  (let [
        ;; receives list of numbers from file with name df
        ;;   tx(3), tCx(3,3), tw2pt, tnt, then tnt tracks th(5) and tCh (5,5)
        inp    (read-string (str "[" (slurp df) "]"))

        ;; aux functions
        get3   (take 9 (drop 3 inp))
 ;;       indexed (fn [coll]  (map-indexed vector coll))
        get5   (fn [i] (let [ sz 5 off (+ 14 (* i sz) (* i sz sz)) ]
                    (for [[idx elt] (indexed inp) :when (>= idx off) :when (< idx (+ off sz))]
                      elt)))
        get25  (fn [i] (let [ sz 25 off (+ 14 5 (* i 30)) ]
                    (for [[idx elt] (indexed inp) :when (>= idx off) :when (< idx (+ off sz))]
                      elt)))
        ;; fill data members
        tx     ;; "initial position" of vertex
               (vec (take 3 inp))
        tCx    ;; error matrix of vertex
               (matrix (vec (map vector (take 3 get3) (drop 3 get3) (drop 6 get3))))
        tw2pt  ;; factor to go from track curvature to pt, given magnetic field
               (first (drop 12 inp))
        tnt    ;; number of tracks associated to vertex
               (first (drop 13 inp))
        th     ;; list of track parameter vectors
               (for [i (range 0 tnt)] (vec (get5 i)))
        tCh    ;; list of track parameter error matrices
               (for [i (range 0 tnt)] (let [itCh (get25 i)] (
                  matrix (vec (map vector (take 5 itCh)
                                          (drop 5 itCh)
                                          (drop 10 itCh)
                                          (drop 15 itCh)
                                          (drop 20 itCh))))))

       ] (->fvtData tx tCx tnt (vec th) (vec tCh) tw2pt)))

(def fvtd (fvtread thisFile))

(pm (get (:tCx fvtd) 0))


(defn doFit []
  (let [
         v0       (:tx fvtd)
         Cv0      (emap #(* 10000.0 %) (:tCx fvtd))   ;;
         h        (:th fvtd); assuming trList is indeed [1 2 3 4 5 6], true for this run
         Ch       (:tCh fvtd)
         w2pt     (:tw2pt fvtd)
       ]
    (do
      (println "--------- input to vertext fit -------------------------")
      (println "initial vertex position v0 ")
      (pm v0)
      (print "±") (pm (sqrt (diagonal Cv0)))
      (print "1st track helix") (pm (first h))
      (print "±%") (pm (scale (div (sqrt (diagonal (first Ch))) (first h)) 100.0))
      (let [[p Cp] (fvHelix2P4 mπ w2pt (first h) (first Ch))]
        (print "1st track momentum ") (pm p)
        (print "±%") (pm (scale (div (sqrt (diagonal Cp)) p) 100.0))))
    (let [[x Cx ql Cql chi2l chi2t] (fvFit v0 Cv0 h Ch)]
      (println "--------- doFit result --------------------------------")
      (println "vertex fit converged, 𝜒2:" chi2t)
      (print "x: ") (pm x)
      (print "±") (pm (sqrt (diagonal Cx)))
      (print "Cx: ") (pm Cx)
      (println "--------- list of fitted q vectors---------------------")
      (doseq [[[iq q] Q chi2q] (map vector (indexed ql) Cql chi2l) ]
        (println "track#" iq ", 𝜒2: " chi2q ": " )
;;        (let [ [p Cp] (fvQ2P4 mπ w2pt q Q)]
;;         (pm p)
;;          (pm Q)
;;          (println ))
      )
      ;;         (println "ql: ") (doseq [ [it tt] (indexed ql) ] (print "track par. " it ": " ) (pm tt))
      )))
(doFit)


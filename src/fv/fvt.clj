(ns fv.fvt 
(:require [fv.coeff :refer :all]))

(use 'clojure.core.matrix)

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
         hl       (:th fvtd); assuming trList is indeed [1 2 3 4 5 6], true for this run
         Chl      (:tCh fvtd)
         w2pt     (:tw2pt fvtd)
       ]
    (do
      (when false (doFit))
      (println "--------- input to vertext fit -------------------------")
      (println "initial vertex position v0 ")
      (print "[x y z]=") (pm v0)
      (print "±") (pm (sqrt (diagonal Cv0)))
      (let [[v Cv ql Cql chi2l chi2t] (fvFit v0 Cv0 hl Chl)]
        (println "--------- doFit result --------------------------------")
        (println "vertex fit converged, 𝜒2:" chi2t)
        (print "v: ") (pm v)
        (print "±") (pm (sqrt (diagonal Cv)))
        (print "Cv: ") (pm Cv)
        (println "--------- list of fitted q vectors---------------------")
        (doseq [[[ih h] H q Q chi2q] (map vector (indexed hl) Chl ql Cql chi2l) ]
          (println "chi2 = " chi2q "prob" )
          (doseq [[x y] (map list h (sqrt (diagonal H)))]
            (print (format "%9.3g ±%9.3g" x y))) (println)
          (doseq [[x y] (map list q (sqrt (diagonal Q)))]
            (print (format "%9.3g ±%9.3g" x y))) (println)
          (println "track#" ih  ", 𝜒2: " chi2q ": ")
          (let [ [p P] (fvHelix2P4 h H mπ w2pt)]
            (fvPMerr "Helix [px py pz E]=" p P))
          (let [ [p P] (fvQ2P3 q Q w2pt)]
            (fvPMerr "Fit q [px py pz]  =" p P)))))))


;;(defrecord fvrec [v Cv 𝜒2v qs Cqs 𝜒2qs])
;;(->fvrec  tCx tnt (vec th) (vec tCh) tw2pt); add to rec

(ns fv.fvt
(:require [fv.coeff :refer :all]
          [clojure.core.matrix :refer :all]
          [fv.fv :refer :all]))

(def dataDirName "./dat")
(def data-dir (file-seq (clojure.java.io/file dataDirName)))
(defn only-files
  "Filter a sequence of files/directories by isFile property of java.io.File"
  [file-s]
  (filter #(.isFile %) file-s))
(defn names
  "Return the .getName property of a sequence of files"
  [file-s]
  (map #(str dataDirName "/" (.getName %)) file-s))
(def theseFiles (-> data-dir only-files names))

(def thisFile "dat/tr05129e001412.dat")
(def otherFile "dat/tr05158e004656.dat")
(def thirdFile "dat/tr05166e001984.dat")
(def tList [1 6 2 3 4 5])

(defrecord fvtData [tx tCx tnt th tCh tw2pt])


(defn- indexed [coll]  (map-indexed vector coll))

(defn- fvtread
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

(defn getHelices [fn]
  (let [
        fvtd  (fvtread fn)
        v0 (:tx fvtd)
        Cv0 (let [Cv00 (emap #(* 10000.0 %) (:tCx fvtd))]
              (array [[(mget Cv00 0 0) 0 0]
                      [0 (mget Cv00 1 1) 0]
                      [0 0 (mget Cv00 2 2)]]))
        w2pt (:tw2pt fvtd) ;; we don't do anything with this, but just use the default
        hl (let [tl tList]
             (vec (map first (sort-by second (map vector (:th fvtd) tl)))))
        Chl (let [tl tList]
              (vec (map first (sort-by second (map vector (:tCh fvtd) tl)))))
        ] (->Helices v0 Cv0 hl Chl)))

(def theseHelices (getHelices thisFile))
(def otherHelices (getHelices otherFile))

(defn printHeader [helices]
  (println "--------- input to vertext fit -------------------------")
  (println "vertex initial position")
  (fvPMerr "v0 [x y z]=" (:v0 helices) (:V0 helices)))

(defn printResults [prong hl Chl]
  (let [
        v (:v prong)
        Cv (:V prong)
        ql (:ql prong)
        Cql (:Ql prong)
        chi2l (:chi2l prong)
        ]
    (println "--------- doFit result --------------------------------")
    (println "vertex fit converged, 𝜒2:" (format  "%9.3g" (reduce + chi2l)))
    (fvPMerr "v [x y z]=" v Cv)
    (println "--------- list of fitted q vectors---------------------")
    (doseq [[[ih h] H q Q chi2q] (map vector (indexed hl) Chl ql Cql chi2l) ]
      (println "chi2=" (format  "%9.3g" chi2q)  "prob" )
      (fvPMerr "Helix params=" h H)
      (fvPMerr "q-vec params=" q Q)
      (println "track#" ih  ", 𝜒2: " (format  "%9.3g" chi2q) ": ")
      (let [ [p P] (fvHelix2P4 h H mπ)]
        (fvPMerr "Helix [px py pz E]=" p P))
      (let [ [p P] (fvQ2P3 q Q)]
        (fvPMerr "Fit q [px py pz]  =" p P)))))

(defn doFitTest [hel]
  (printHeader hel)
  (let [
        hl  (:hl hel)
        Hl  (:Hl hel)
        pr  (fvFit hel)
        p5  (->Prong (:v pr) (:V pr) (drop-last (:ql pr)) (drop-last (:Ql pr)) (drop-last (:chi2l pr)))
        ]
    (fvRemove pr hel)
    (println "Inv Mass " (invMass p5))
    (printResults pr hl Hl )))

(defn doFitTests [] (do
                  (doFitTest theseHelices); so I can call it from fireplace
                  (doFitTest otherHelices)
                  (println "----------------------------------------")
                  (println theseFiles)
                  (println "----------------------------------------")
                  (doall (take 10 (map #(doFitTest (getHelices %)) theseFiles)))))

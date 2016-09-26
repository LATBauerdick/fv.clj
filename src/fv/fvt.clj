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

(defn hSlurp [fn]
  (let [
        fvtd  (fvtread fn)
        v0    (:tx fvtd)
        Cv0   (let [Cv00 (emap #(* 10000.0 %) (:tCx fvtd))]
                (array [[(mget Cv00 0 0) 0 0]
                        [0 (mget Cv00 1 1) 0]
                        [0 0 (mget Cv00 2 2)]]))
        w2pt  (:tw2pt fvtd) ;; we don't do anything with this, but just use the default
        hl    (:th fvtd)
        Chl   (:tCh fvtd)
        ] (->Helices v0 Cv0 hl Chl)))
(defn hInto [ hel0 hel1]
  (let [
        hl (into (:hl hel0) (:hl hel1))
        Hl (into (:Hl hel0) (:Hl hel1))
        ;; here we should find the average vertex position.... 
        ;; instead we're using the one from hel0
        ] (->Helices (:v0 hel0) (:V0 hel0) hl Hl)))


(def theseHelices (hSlurp thisFile))
(def otherHelices (hSlurp otherFile))

(defn hFilter
  "
  -- return the Helices filtered with just the tracks that are in rng
  "
  [hel rng]
  (let [
        h#l (range 0 (count (:hl hel)))
        hl (map last (filter #(some #{(first %)}  rng)  (map vector  h#l  (:hl hel))))
        Hl (map last (filter #(some #{(first %)}  rng)  (map vector  h#l  (:Hl hel))))
        ] (->Helices (:v0 hel) (:V0 hel) hl Hl)))

(defn pFilter
  "
  -- return the Prong filtered with just the tracks that are in rng
  "
  [pr rng]
  (let [
        t#l (:t#l pr)
        qQl (map last (filter #(some #{(first %)}  rng)  (map vector  (:t#l pr)  (:qQl pr))))
        chi2l (map last (filter #(some #{(first %)}  rng)  (map vector  (:t#l pr)  (:chi2l pr))))
        ] (->Prong (:v pr) (:V pr) qQl chi2l rng)))

(defn printHeader [helices]
  (println "--------- input to vertext fit -------------------------")
  (println "vertex initial position")
  (fvPMerr "v0 [x y z]=" (:v0 helices) (:V0 helices)))

(defn printResults [prong hl Chl]
  (let [
        v (:v prong)
        Cv (:V prong)
        qQl (:qQl prong)
        chi2l (:chi2l prong)
        t#l (:t#l prong)
        ]
    (println "--------- doFit result --------------------------------")
    (println "vertex fit converged, 𝜒2:" (format  "%9.3g" (reduce + chi2l)))
    (fvPMerr "v [x y z]=" v Cv)
    (println "--------- list of fitted q vectors---------------------")
    (doseq [[h H qQ chi2 t#] (map vector hl Chl qQl chi2l t#l) ]
      (println "chi2=" (format  "%9.3g" chi2)  "prob" )
      (fvPMerr "Helix params=" h H)
      (fvPMerr "q-vec params=" (:q qQ) (:Q qQ))
      (println "track#" t#  ", 𝜒2: " (format  "%9.3g" chi2) ": ")
      (let [ [p P] (fvH2P4 h H :m mπ)]
        (fvPMerr "Helix [px py pz E]=" p P))
      (let [ [p P] (fvQ2P3 qQ)]
        (fvPMerr "Fit q [px py pz]  =" p P)))))

(defn doFitTest
  "
  -- do a series of test fits on helices hel, with optional 5-prong list l5
  -- fit and print inv mass for different configurations:
  -- mass of all helices, mass fitted, same with 5 prong, refitted 5 prong
  "
  [hel & {:keys [l5] :or {l5 (range 0 5)}}]

  (printHeader hel)
  (let [ pr  (fvFit hel) ]
    (fvRemove pr hel) ;; at the moment, just list contributing chi2 for each fitted helix
    (let [ pl (map #(fvH2P4 %1 %2) (:hl hel) (:Hl hel)) ]
      (->> pl
           invMass
           (apply format "%9.5g ±%9.3g GeV")
           (println "Inv Mass 6 helix")))
    (->> pr
         :qQl
         (map fvQ2P4)
         invMass
         (apply format "%9.5g ±%9.3g GeV")
         (println "Inv Mass 6 fit  "))
    (let [
          hel5  (hFilter hel l5)
          pl    (map #(fvH2P4 %1 %2) (:hl hel5) (:Hl hel5))
          ]
      (->> pl
           invMass
           (apply format "%9.5g ±%9.3g GeV")
           (println "Inv Mass 5 helix")))
    (->> (pFilter pr l5)
         :qQl
         (map fvQ2P4)
         invMass
         (apply format "%9.5g ±%9.3g GeV")
         (println "Inv Mass 5 fit  "))
    (->> (hFilter hel l5)
         fvFit
         :qQl
         (map fvQ2P4)
         invMass
         (apply format "%9.5g ±%9.3g GeV")
         (println "Inv Mass 5 refit"))
    (let [ hl (:hl hel),  Hl (:Hl hel)]
      (printResults pr hl Hl ))))

(def allHelices (reduce #(hInto %1 (hSlurp %2)) (hSlurp (first theseFiles))  (next theseFiles)))

(defn doFitTests [] (do
                  (doFitTest theseHelices :l5 [0 2 3 4 5]); so I can call it from fireplace
                  (doFitTest otherHelices :l5 [0 1 2 4 5])
                  (doFitTest (hInto theseHelices otherHelices) :l5 [0 2 3 4 5])
                  (doFitTest (hInto (hInto theseHelices otherHelices) (hSlurp thirdFile)) :l5 [0 2 3 4 5])
                  (doFitTest  allHelices :l5 [ 1 2 3 4 5 ])

                  (println "----------------------------------------")
                  (println theseFiles)
                  (println "----------------------------------------")
                  (doall (take 10 (map #(doFitTest (hSlurp %)) theseFiles)))))

(ns fv.fvt)
(use 'clojure.core.matrix)
(set-current-implementation :vectorz)

(def thisFile "dat/tr05129e001412.dat")
(defrecord fvtData [tx tCx tnt th tCh])


(defn fvtread [df]
  (let [
        ;; receives list of numbers from file with name df
        ;;   tx(3), tCx(3,3), tw2pt, tnt, then tnt tracks th(5) and tCh (5,5)
        inp    (read-string (str "[" (slurp df) "]"))

        ;; aux functions
        itx    (take 3 inp)
        itCx   (take 9 (drop 3 inp))
        itw2pt (first (drop 12 inp))
        indexed (fn [coll]  (map-indexed vector coll))
        get5   (fn [i] (let [ sz 5 off (+ 14 (* i sz) (* i sz sz)) ]
                    (for [[idx elt] (indexed inp) :when (>= idx off) :when (< idx (+ off sz))]
                      elt)))
        get25  (fn [i] (let [ sz 25 off (+ 14 5 (* i 30)) ]
                    (for [[idx elt] (indexed inp) :when (>= idx off) :when (< idx (+ off sz))]
                      elt)))
        ;; fill data members
        tx     ;; "initial position" of vertex
               (vec itx)
        tCx    ;; error matrix of vertex
               (matrix (vec (map vector (take 3 itCx) (drop 3 itCx) (drop 6 itCx))))
        tnt    ;; number of tracks associated to vertex
               (first (drop 13 inp))
        th     ;; list of track parameter vectors
               (for [i (range 0 (dec tnt))] (vec (get5 i)))
        tCh    ;; list of track parameter error matrices
               (for [i (range 0 (dec tnt))] (let [itCh (get25 i)]  (matrix (vec (map vector (take 5 itCh) (drop 5 itCh) (drop 10 itCh) (drop 15 itCh) (drop 20 itCh))))))

       ] (->fvtData tx tCx tnt th tCh)   ))

(def fvtd (fvtread thisFile))



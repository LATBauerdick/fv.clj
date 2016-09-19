(ns fv.try)

(defrecord my_rec [hl Hl ql Ql])
(def r (->my_rec [1 2 3 4 5] [10 20 30 40 50] [] []))

(def filt (fn [v V h H]
              (println "filt")
              (let [vf (+ v h)
                    Vf (+ V H)
                    q (* vf 1000)
                    Q (* Vf 1000)]
              [vf Vf q Q])))

(defn filterer [f v0 V0 hl Hl]
  (loop [v0 v0 V0 V0 hl hl Hl Hl ql [] Ql []]
    (if (empty? hl)
      [v0 V0 ql Ql]
      (do (println "filterer looping")
          (let [[v V q Q] (f v0 V0 (first hl) (first Hl))]
            (recur v V (next hl) (next Hl) (conj ql q) (conj Ql Q)))))))
(filterer filt 1000 10000 [1 2 3 4 5] [10 20 30 40 50])

(def filt1 (fn [[ v V] [h H]]
            (let [vf (+ v h) Vf (+ V H)]
                 (conj (:ql r) (* vf 1000))
                 (conj (:Ql r) (* Vf 1000))
                 (print h H (:ql r) " ")
                 [vf Vf])))
(reduce filt1 [100 1000] (map vector (:hl r) (:Hl r)))

(defn strict-map2 [f coll]
  (loop [coll coll acc []]
    (if (empty? coll)
      acc
      (recur (next coll)
             (conj acc (f (first coll)))))))
(strict-map2 - (range 5))


(defn my_index [coll]
  (cond
    (map? coll) (seq coll)
    (set? coll) (map vector coll coll)
    :else (map vector (iterate inc 0) coll)))

(defn my_pos [pred coll] ( for [[i v] (my_index coll) :when (pred v)] i))
(my_index {:a 1 :b 2 :c 3 :d 4})
(my_pos #{3 4} {:a 1 :b 2 :c 3 :d 4})
(my_pos #{3 4} [1 2 3 4 5])


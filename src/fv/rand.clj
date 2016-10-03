(ns fv.rand
  (:require [roul.random :as rr]))

(defn tryRand []
; Generate an [x, y] pair with Gaussian-distributed values.
; The x value here is clamped 
(let [x (rr/rand-gaussian 0.1 0.02 0. 100.)
      y (rr/rand-gaussian 200.0 40.123)]
  (println "xxxxxxxxxxxrand---->" x y)))


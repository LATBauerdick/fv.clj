(ns fv.matrix)
(use 'clojure.core.matrix.operators)
(set-currusent-implementation :vectorz)

(array (range 5))  ;; => [0 1 2 3 4]
(def M ( array ( for [i (range 3)] (for [j (range 3)] (+ i (* j 1000))))))
;; =>"[[0.0,1000.0,2000.0],\n[1.0,1001.0,2001.0],\n[2.0,1002.0,2002.  0]]"]
(mget M 1 2)
(slice M 0)
(apply + (slices M))
(map inc (slice M 1))
(ereduce + M)
(eseq M )
(zero-matrix 5 5)
(identity-matrix 5)
(permutation-matrix [3 1 0 2])

(transpose M)
(mmul [[9 2 7] [6 4 8]] [[2 8] [3 4] [5 9]])

(def π 3.141592653589793)
(def τ (* 2.0 π))
(defn rot [turns]
  (let [a (* τ turns)]
    [[  (cos a)  (sin a)]
     [(-(sin a)) (cos a)]]))
(mmul (rot 1/8) [3 4])

(defn proportions
  "Normalize a vector to a sum of 1.0"
  ([v] (/ v (esum v))))
(esum (proportions M))
(pm (proportions M))
; (mmul (inverse ( M )) M)

(def a (mutable [1 2]))
(add! a 1)

(+ [[1 2]
    [3 4]]
   (* (identity-matrix 2) 3.0))
; => [[4.0 2.0]
;     [3.0 7.0]]

(shape [[2 3 4] [5 6 7]]) ; => [2 3]

(mmul
  (array [[2 2] [3 3]])
  (array [[4 4] [5 5]])) ; => [[18 18] [27 27]]

(def v (array [1 2 3]))
v

; this is the left-hand-side of a linear system, 
; a 2x2 matrix.
(def A [[-0.66718   3.05353] [0.60487   2.18721]])


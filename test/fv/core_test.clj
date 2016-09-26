(ns fv.core-test
  (:require [clojure.test :refer :all]
            [fv.core :refer :all]
            [fv.fv :refer :all]
            [fv.coeff :refer :all]
            [fv.fvt :refer :all]))




(deftest chi2-test
  (testing "chi2 returned from fit on first data file."
    (is (> 0.0001 (Math/abs (- 334.34 (reduce + (:chi2l(fvFit(hSlurp thisFile))))))))))

(deftest invMass-test
  (testing "invariant mass of first 5 momenta returned from fit on first data file."
    (is (> 0.0001 (Math/abs (- 1.499 ((invMass (map fvQ2P4 (:qQl (fvFit(hFilter (hSlurp thisFile) [0 2 3 4 5]))))) 0)))))))

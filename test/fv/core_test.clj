(ns fv.core-test
  (:require [clojure.test :refer :all]
            [fv.core :refer :all]
            [fv.fv :refer :all]
            [fv.fvt :refer :all]))




(deftest chi2-test
  (testing "chi2 returned from fit on first data file."
    (is (> 0.0001 (Math/abs (- 334.34 (reduce + (:chi2l(fvFit(getHelices thisFile))))))))))

(ns fv.core
  (:gen-class))

(require 'fv.fvt)
(refer 'fv.fvt)

(defn -main
  [& args]
  (doFitTests))



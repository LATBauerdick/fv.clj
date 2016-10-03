(defproject fv "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [net.mikera/vectorz-clj "0.45.0"]
                 [net.mikera/core.matrix "0.55.0"]
                 [roul "0.2.0"]
                ]
  :main ^:skip-aot fv.core
  :target-path "target/%s"
  :profiles {:uberjar {:aot :all}})

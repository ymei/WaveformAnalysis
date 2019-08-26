;;;; -*- Mode: Lisp; -*-

(cl:in-package #:cl-user)
(defpackage #:wavana-system
  (:use #:cl #:asdf))
(in-package #:wavana-system)

(asdf:defsystem wavana
  :name "wavana"
  :version "0.1"
  :author "Yuan Mei"
  :licence "BSD"
  :components
  ((:file "wavana" :depends-on ("gnuclot" "mrdary"))
;  (:file "plot-pmt-pattern" :depends-on ("clppu"))
   (:file "gnuclot")
   (:file "mrdary"))
  :depends-on (cffi))

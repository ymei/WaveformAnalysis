;;;; The functions in this file reads gnuplot's datafile style data
;;;; into maxima's lists and matrices.
(in-package :maxima)

;; Read datafile into a list of matrices.  A contiguous data block
;; forms a matrix.  Matrices separated by blank lines in datafile are
;; stored into separate list elements.
(defun $read_matrix_list (stream-or-filename &rest args)
  (let ((sep-ch-flag (and (cdr args) (cadr args))))
    (if (streamp stream-or-filename)
        (read-nested-list-from-stream stream-or-filename sep-ch-flag)
        (let ((file-name (require-string stream-or-filename)))
          (with-open-file (in file-name :if-does-not-exist nil)
            (if (not (null in))
                (read-nested-list-from-stream in sep-ch-flag)
                (merror "read_matrix_list: no such file `~a'" file-name)))))))

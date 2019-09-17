;;;; The functions in this file reads/writes gnuplot's datafile style
;;;; data into maxima's lists and matrices.
(if (find-package 'maxima) ; A poor way of finding if it is in maxima.
    (in-package :maxima)
    (in-package :cl-user))

;;; Read datafile into a list of matrices.  A contiguous data block
;;; forms a matrix.  Matrices separated by blank lines in datafile are
;;; stored into separate list elements.
#|
(defun $read_matrix_list (stream-or-filename &rest args)
  (let ((sep-ch-flag (and (cdr args) (cadr args))))
    (if (streamp stream-or-filename)
        (read-nested-list-from-stream stream-or-filename sep-ch-flag)
        (let ((file-name (require-string stream-or-filename)))
          (with-open-file (in file-name :if-does-not-exist nil)
            (if (not (null in))
                (read-nested-list-from-stream in sep-ch-flag)
                (merror "read_matrix_list: no such file `~a'" file-name)))))))
|#

(declaim (optimize (speed 3) (debug 0) (safety 0)))
(defvar *commentschars* "#!%") ; the first element will be used in output.
(defvar *blankchars* '(#\Space #\Tab #\Linefeed #\Page #\Return))

(defun first-nonblank-is-commentschar-p (s)
  (let ((pos (position nil s :test (lambda (a b) (declare (ignore a))
                                     (not (member b *blankchars*))))))
    (if pos
        (if (find (char s pos) *commentschars*) t nil))))

(defun is-line-blank-p (s)
  (if (> (length s) 0)
      (if (find nil s :test (lambda (a b) (declare (ignore a))
                              (not (member b *blankchars*))))
          nil t)
      t))

(defun parse-string-to-list-of-numbers (s)
  (with-input-from-string (ss s)
    (let ((*read-default-float-format* 'double-float))
      (loop for x = (read ss nil nil) while x collect (coerce x 'double-float)))))

;;; Read gnuplot data file into a list of matrices.  A matrix is a
;;; block of consecutive lines of data.  A blank line indicates the
;;; end of a matrix and the next matrix is appended to the list, which
;;; is ultimately returned.
(defun read-gpdata-file (fname)
  (unless (probe-file fname)
    (error "File \"~a\" not found!" fname)
    (return-from read-gpdata-file))
  (let ((lm (list '(mlist simp))) (mx (list '($MATRIX SIMP))) (ml nil) (is-in-matrix nil) (nl 0))
    ;;      ^ Make sure new cons cells are generated.
    (with-open-file (infile fname :if-does-not-exist nil)
      (loop for line = (read-line infile nil 'eof)
         until (eq line 'eof)
         do
           (tagbody
            s0
              (if is-in-matrix
                  (if (is-line-blank-p line) ; blank line indicates end of a matrix.
                      (progn
                        (when (> nl 0) ; there are valid data in matrix.
                          (push (nreverse mx) lm)
                          (setf mx (list '($MATRIX SIMP))
                                nl 0
                                is-in-matrix nil))) ; not in reading a matrix anymore.
                      (progn ; non-blank like,
                        (unless (first-nonblank-is-commentschar-p line) ; valid matrix data
                          (setf ml (cons '(mlist simp) (parse-string-to-list-of-numbers line)))
                          (push ml mx)
                          (incf nl)))))
              (unless is-in-matrix
                (unless (is-line-blank-p line) ; a new matrix starts.
                  (setf is-in-matrix t)
                  (go s0))))))
    (when (> nl 0)
      (push (nreverse mx) lm))
    (nreverse lm)))

(defun $read_gpdata (fname)
  (read-gpdata-file fname))

;;; Write lists and matrices into gnuplot data file.
(defun write-gpdata-stream (stream val)
  (when (consp (car val))
    (if (atom (cadr val))
        (progn
          (mapcar (lambda (v) (format stream " ~24,16,3,,,,'EE" v)) (cdr val))
          (format stream "~%"))
        (if (eql (caar val) '$matrix)
            (progn
              (format stream "~a Matrix~%" (elt *commentschars* 0))
              (mapcar (lambda (v) (write-gpdata-stream stream v))
                      (cdr val))
              (format stream "~%"))
            (progn
              (format stream "~a List~%" (elt *commentschars* 0))
              (mapcar (lambda (v) (write-gpdata-stream stream v))
                      (cdr val)))))))
(defun write-gpdata-file (fname val)
  (with-open-file (outfile fname :direction :output :if-exists :supersede)
    (write-gpdata-stream outfile val)))

(defun $write_gpdata (fname val)
  (write-gpdata-file fname val)
  fname)

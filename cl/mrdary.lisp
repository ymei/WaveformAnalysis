(in-package #:cl-user)

(defpackage #:mrdary
  (:use #:cffi #:cl)
  (:export #:init-f
           #:read-all
           #:free-f
           #:free-m
           #:value-at
           #:min-max
           #:rows
           #:columns
           #:with-data-file
           #:with-read-data-file
           #:get-array-of-line
           #:get-list-of-line
           #:hist1
           #:hist2
           #:hist2p))
(in-package #:mrdary)

(define-foreign-library libWavAna
  (:darwin "libWavAna.dylib")
  (:unix "libWavAna.so")
  (t (:default "libWavAna")))
(use-foreign-library libWavAna)

;; This is tested to be true on
;; Linux:  32-bit size_t is unsigned int (same as unsigned long) which is 4 bytes in size
;;         64-bit size_t is unsigned long, which is 8 bytes in size
;; MacOSX: 64-bit size_t is unsigned long, which is 8 bytes in size
;; So it is safe to define size_t as unsigned long to make it automatically correct on both 32-bit
;; and 64-bit platforms
(defctype size_t :unsigned-long)
;; here a lisp-side corresponding type specifier (a function) is defined
(deftype size_t () '(unsigned-byte #+x86-64 64
                                   #+x86 32))
(defcstruct mrdary_handle
  "mrdary_handle"
  (fp :pointer)
  (linebuf :pointer)
  (column size_t)
  (row size_t)
  (rowmax size_t)
  (marray (:pointer :double))
  (min (:pointer :double))
  (max (:pointer :double)))

(defcfun "mrdary_init_f" :pointer (fname (:pointer :char)) (rowmax size_t))
(defcfun ("mrdary_read_all" read-all) size_t (hdl :pointer))
(defcfun ("mrdary_free") :int (hdl :pointer))
(defcfun "mrdary_value_mn" :pointer (hdl :pointer) (m size_t) (n size_t))
(defcfun "mrdary_min" :pointer (hdl :pointer) (i size_t))
(defcfun "mrdary_max" :pointer (hdl :pointer) (i size_t))

(defconstant +rowmax-default+ 1000000)

;; Ensure that the handle is valid
;; "The handle is either not a pointer or is a null pointer."
(define-condition invalid-handle-error (error)
  ((handle :initarg :handle :reader handle)))
(declaim (inline ensure-valid-handle))
(defun ensure-valid-handle (hdl)
  "Signal error when the handle is either not a pointer or is a null pointer."
  (if (or (not (pointerp hdl)) (null-pointer-p hdl))
      (error 'invalid-handle-error :handle hdl)))

(defun init-f (fname &optional (rowmax +rowmax-default+))
  "Initialize the mrdary library.
rowmax is default to 1e6.
return: mrdary_hdl"
  (declare (type (size_t) rowmax))
  (let ((ret (with-foreign-string (fn (#+ccl ccl:native-translated-namestring
                                       #-ccl namestring fname))
               (mrdary-init-f fn rowmax))))
    (assert (and (pointerp ret) (not (null-pointer-p ret)))
            (fname)
            "init-f on file \"~S\" failed" fname)
    ret))

(defun free-f (hdl)
  "Free the handle"
  (ensure-valid-handle hdl)
  (mrdary-free hdl))

(defmacro free-m (%hdl)
  "Free the handle when %hdl is a global variable"
  `(progn (free-f ,%hdl)
          (setf ,%hdl (null-pointer))))

;; signal error when row and column indicies are out of bound
(define-condition index-out-of-bound-error (error)
  ((indices :initarg :indices :reader indices)))

(declaim (inline value-at (setf value-at)))
(defun value-at (hdl m n)
  "Returns the value at (m, n)"
  (declare (type (size_t) m n)
           (optimize (debug 1) (safety 0) (speed 3)))
  (ensure-valid-handle hdl)
  (if (or (>= m (rows hdl)) (>= n (columns hdl)))
      (error 'index-out-of-bound-error :indices '(m n)))
  (let ((ptr (mrdary-value-mn hdl m n)))
    ;(or (null-pointer-p ptr)
    (the double-float (mem-ref ptr :double))))

(defun (setf value-at) (val hdl m n)
  "Sets the value at (m, n)."
  (declare (type (size_t) m n) (type double-float val)
           (optimize (debug 1) (safety 0) (speed 3)))
  (ensure-valid-handle hdl)
  (if (or (>= m (rows hdl)) (>= n (columns hdl)))
      (error 'index-out-of-bound-error :indices '(m n)))
  (let ((ptr (mrdary-value-mn hdl m n)))
    ;(or (null-pointer-p ptr)
    (setf (mem-ref ptr :double) val)))

(defun rows (hdl)
  "return the number of rows"
  (ensure-valid-handle hdl)
  (foreign-slot-value hdl 'mrdary_handle 'row))

(defun columns (hdl)
  "return the number of columns"
  (ensure-valid-handle hdl)
  (foreign-slot-value hdl 'mrdary_handle 'column))

(defun min-max (hdl i)
  "return a list (min max) of row i"
  (declare (type (size_t) i))
  (ensure-valid-handle hdl)
  (let ((pmin (mrdary-min hdl i))
        (pmax (mrdary-max hdl i)))
    (unless (or (null-pointer-p pmin) (null-pointer-p pmax))
      (list (mem-ref pmin :double) (mem-ref pmax :double)))))

(defmacro with-data-file ((%handle %fname &optional (%rowmax +rowmax-default+)) &body %body)
  `(let ((,%handle nil))
     (unwind-protect
          (progn (setf ,%handle (init-f ,%fname ,%rowmax))
                 ,@%body)
       (when ,%handle (free-f ,%handle)))))

(defmacro with-read-data-file ((%handle %fname &optional (%rowmax +rowmax-default+)) &body %body)
  `(let ((,%handle nil))
     (unwind-protect
          (progn (setf ,%handle (init-f (namestring ,%fname) ,%rowmax))
                 (read-all ,%handle)
                 ,@%body)
       (when ,%handle (free-f ,%handle)))))

(defun get-array-of-line (handle i)
  (ensure-valid-handle handle)
  (let* ((col (columns handle))
         (ary (make-array col :element-type 'double-float)))
    (dotimes (j col)
      (let ((tmp (value-at handle i j)))
        ; (if (eql tmp t) (assert tmp)
            (setf (aref ary j) tmp)))
    ary))

(defun get-list-of-line (handle i)
  (ensure-valid-handle handle)
  (let* ((col (columns handle))
         (list nil))
    (dotimes (j col)
      (setf list (nconc list (list (value-at handle i j)))))
    list))

(defun hist1 (hdl col binspec &key filter tranfun)
  "1D histogram
   return a 2d array holding the histogram
   binspec has the form '(n (min max)), and '(min max) can be nil
   filter has the form #'(lambda (hdl i) ...) and returns t or nil
   tranfun has the form #'(lambda (hdl i) ...) and returns '(v add)"
  (declare (optimize (debug 3)))
  (let* ((n (car binspec))
         (minmax (cadr binspec))
         (min 0.0d0)(max 1.0d0)
         (h (make-array (list n 2) :element-type 'double-float :initial-element 0.0d0))
         (dx 1.0d0)(c 0)(j 0)(v 0.0d0)(add 1.0d0))
    (declare (type double-float min max dx v add) (type fixnum c j col n))
    (if minmax
        (progn (setf min (coerce (car minmax) 'double-float))
               (setf max (coerce (cadr minmax) 'double-float)))
        (let ((mm (min-max hdl col)))
          (setf min (car mm)) (setf max (cadr mm))))
    ;; fill the first column of the array
    (setf dx (the double-float (/ (- max min) n)))
    (dotimes (i n) (setf (aref h i 0) (+ min (* i dx))))
    ;; fill the counts
    (dotimes (i (rows hdl))
      ;; when index is out of range, fill 0.0d0
      (if tranfun
          (let ((ret (funcall tranfun hdl i)))
            (setf v (coerce (car ret) 'double-float))
            (setf add (coerce (cadr ret) 'double-float)))
          (let ((tmp (value-at hdl i col)))
            (setf v tmp)
            (setf add 1.0d0)))
      (handler-case (setf j (floor (/ (- v min) dx)))
        (arithmetic-error () (setf j -1)))
      (when (and (>= j 0) (< j n) (if filter (funcall filter hdl i) t))
        (incf (aref h j 1) add)
        (incf c)))
    (values h c dx)))

(defun hist2 (hdl cn cm xbinspec ybinspec &key filter tranfun)
  "2D histogram
   cn : column in hdl for col in histogram (x)
   cm : column in hdl for row in histogram (y)
   [x,y]binspec has the form '(n (min max)), and '(min max) can be nil
   filter has the form #'(lambda (hdl i) ...) and returns t or nil
   tranfun has the form #'(lambda (hdl i) ...) and returns '(vx vy add)
   returns the result in a 2D array of m x n"
  (declare (optimize (debug 1) (safety 0) (speed 3)))
  (let* ((nn (car xbinspec)) (mmn (cadr xbinspec))
         (nm (car ybinspec)) (mmm (cadr ybinspec))
         (mmin 0.0d0) (mmax 1.0d0) (nmin 0.0d0) (nmax 1.0d0)
         (h (make-array (list nm nn) :element-type 'double-float :initial-element 0.0d0))
         (dy 1.0d0) (dx 1.0d0) (c 0) (m 0) (n 0) (vy 0.0d0) (vx 0.0d0) (add 1.0d0))
    (declare (type double-float mmin mmax nmin nmax dy dx vy vx)
             (type fixnum cm cn nm nn c m n))
    (if mmm
        (progn (setf mmin (coerce (car mmm) 'double-float))
               (setf mmax (coerce (cadr mmm) 'double-float)))
        (let ((mm (min-max hdl cm)))
          (setf mmin (car mm)) (setf mmax (cadr mm))))
    (if mmn
        (progn (setf nmin (coerce (car mmn) 'double-float))
               (setf nmax (coerce (cadr mmn) 'double-float)))
        (let ((mm (min-max hdl cn)))
          (setf nmin (car mm)) (setf nmax (cadr mm))))
    ;; increment
    (setf dx (/ (- nmax nmin) nn))
    (setf dy (/ (- mmax mmin) nm))
    ;; fill the counts
    (dotimes (i (rows hdl))
      (if tranfun
          (let ((ret (funcall tranfun hdl i)))
            (setf vx (coerce (car ret) 'double-float))
            (setf vy (coerce (cadr ret) 'double-float))
            (setf add (coerce (caddr ret) 'double-float)))
          (let ((tmpn (value-at hdl i cn))
                (tmpm (value-at hdl i cm)))
            (setf vx tmpn)
            (setf vy tmpm)))
      (handler-case (setf m (floor (/ (- vy mmin) dy)))
        (arithmetic-error () (setf m -1)))
      (handler-case (setf n (floor (/ (- vx nmin) dx)))
        (arithmetic-error () (setf n -1)))
      (when (and (>= m 0) (< m nm) (>= n 0) (< n nn) (if filter (funcall filter hdl i) t))
        (incf (aref h m n) add)
        (incf c)))
    (values h (list mmin dy nmin dx) c)))

(defun hist2p (hdl cn cm cz xbinspec ybinspec &key filter tranfun)
  "2D profile
   cn : column in hdl for col (x) in histrogram
   cm : column in hdl for row (y) in histrogram
   cz : column in hdl for z value
   [x,y]binspec has the form '(n (min max)), and '(min max) can be nil
   filter has the form #'(lambda (hdl i) ...) and returns t or nil
   tranfun has the form #'(lambda (hdl i) ...) and returns '(vx vy z)
   returns the result in a 2D array of m x n"
  (declare (optimize (debug 1) (safety 0) (speed 3)))
  (let* ((nn (car xbinspec)) (mmn (cadr xbinspec))
         (nm (car ybinspec)) (mmm (cadr ybinspec))
         (mmin 0.0d0) (mmax 1.0d0) (nmin 0.0d0) (nmax 1.0d0)
         (h (make-array (list nm nn) :element-type 'double-float :initial-element 0.0d0))
         (ch (make-array (list nm nn) :element-type 'double-float :initial-element 0.0d0))
         (dy 1.0d0) (dx 1.0d0) (c 0) (m 0) (n 0) (vy 0.0d0) (vx 0.0d0) (add 1.0d0))
    (declare (type double-float mmin mmax nmin nmax dy dx vy vx add)
             (type fixnum cm cn nm nn c m n))
    (if mmm
        (progn (setf mmin (coerce (car mmm) 'double-float))
               (setf mmax (coerce (cadr mmm) 'double-float)))
        (let ((mm (min-max hdl cm)))
          (setf mmin (car mm)) (setf mmax (cadr mm))))
    (if mmn
        (progn (setf nmin (coerce (car mmn) 'double-float))
               (setf nmax (coerce (cadr mmn) 'double-float)))
        (let ((mm (min-max hdl cn)))
          (setf nmin (car mm)) (setf nmax (cadr mm))))
    ;; increment
    (setf dy (/ (- mmax mmin) nm))
    (setf dx (/ (- nmax nmin) nn))
    ;; fill the counts
    (dotimes (i (rows hdl))
      (if tranfun
          (let ((ret (funcall tranfun hdl i)))
            (setf vx (coerce (car ret) 'double-float))
            (setf vy (coerce (cadr ret) 'double-float))
            (setf add (coerce (caddr ret) 'double-float)))
          (let ((tmpm (value-at hdl i cm))
                (tmpn (value-at hdl i cn))
                (tmpa (value-at hdl i cz)))
            (setf vy tmpm)
            (setf vx tmpn)
            (setf add tmpa)))
      (handler-case (setf m (floor (/ (- vy mmin) dy)))
        (arithmetic-error () (setf m -1)))
      (handler-case (setf n (floor (/ (- vx nmin) dx)))
        (arithmetic-error () (setf n -1)))
      (when (and (>= m 0) (< m nm) (>= n 0) (< n nn)
                 (if filter (funcall filter hdl i) t))
        (handler-case (incf (aref h m n) add)
          (arithmetic-error () nil))
        (incf (aref ch m n) 1.0d0)
        (incf c)))
    (dotimes (i nm)
      (dotimes (j nn)
        (when (> (aref ch i j) 0.0d0)
          (setf (aref h i j) (/ (aref h i j) (aref ch i j))))))
    (values h (list mmin dy nmin dx) c)))

;; Old interface, deprecated
;; (defun hist1g (hdl col n &key minmax filter tranfun)
;;   "general version
;;    returns a 2d array holding the histogram
;;    filter is of form #'(lambda (hdl i) ...) and returns t or nil
;;    tranfun returns '(v add)"
;;   (declare (optimize (debug 3)))
;;   (let ((min 0.0d0)(max 1.0d0)
;;         (h (make-array (list n 2) :element-type 'double-float :initial-element 0.0d0))
;;         (dx 1.0d0)(c 0)(j 0)(v 0.0d0)(add 1.0d0))
;;     (declare (type double-float min max dx v add) (type fixnum c j col n))
;;     (if minmax
;;         (progn (setf min (coerce (car minmax) 'double-float))
;;                (setf max (coerce (cadr minmax) 'double-float)))
;;         (let ((mm (min-max hdl col)))
;;           (setf min (car mm)) (setf max (cadr mm))))
;;     ;; fill the first column of the array
;;     (setf dx (the double-float (/ (- max min) n)))
;;     (dotimes (i n) (setf (aref h i 0) (+ min (* i dx))))
;;     ;; fill the counts
;;     (dotimes (i (rows hdl))
;;       ;; when index is out of range, fill 0.0d0
;;       (if tranfun
;;           (let ((ret (funcall tranfun hdl i)))
;;             (setf v (coerce (car ret) 'double-float))
;;             (setf add (coerce (cadr ret) 'double-float)))
;;           (let ((tmp (value-at hdl i col)))
;;             (setf v tmp)
;;             (setf add 1.0d0)))
;;       (handler-case (setf j (floor (/ (- v min) dx)))
;;         (arithmetic-error () (setf j -1)))
;;       (when (and (>= j 0) (< j n) (if filter (funcall filter hdl i) t))
;;         (incf (aref h j 1) add)
;;         (incf c)))
;;     (values h c dx)))
;;
;; (defun hist1d (hdl col n &key minmax filter)
;;   "returns a 2d array holding the histogram
;;    filter is of form #'(lambda (hdl i) ...)"
;;   (hist1g hdl col n :minmax minmax :filter filter))
;;
;; (defun hist2d (hdl cn cm nn nm &key mmn mmm filter)
;;   "cn : column as col (x) nbins nn
;;    cm : column as row (y) nbins nm
;;    store the result in a 2D array of m x n"
;;   (hist2g hdl cn cm nn nm :mmn mmn :mmm mmm :filter filter))
;;
;; (defun hist2p (hdl cn cm nn nm cz &key mmn mmm filter)
;;   "2D profile"
;;   (hist2pg hdl cn cm nn nm cz :mmn mmn :mmm mmm :filter filter))

;(defvar *hdl* (init-f "/home/ymei/read-double-float/data.dat"))

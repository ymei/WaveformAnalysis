;;; variables starting with % is reserved inside of macros to avoid variable capture

(in-package #:cl-user)

(defpackage #:gnuclot
  (:use #:common-lisp 
        #+sbcl #:sb-ext)
  (:shadow #:set #:load)
  (:export #:*fname*
           #:start-gnuplot
           #:stop-gnuplot
           #:stop-all-gnuplot
           #:list-gnuplot-connection-names
           #:switch-to-gnuplot
           #:gnuplot-command
           #:plot
           #:splot
           #:replot
           #:set
           #:unset
           #:show
           #:reset
           #:format-list-to-gnuplot-fun
           #:with-temp-file
           #:delete-all-temp-files
           #:plotlist
           #:plotarray
           #:plotarray2))
(in-package #:gnuclot)

(defvar *debug* t)
(defvar *gnuplot-path* "gnuplot")
(defvar *fname* nil)
(defconstant +default-gnuplot-connection-name+ 0)
(defvar *current-gnuplot-connection-name* +default-gnuplot-connection-name+)
(defvar *gnuplot-connections* (make-hash-table) "currently open connections")
(defvar *current-gnuplot-proc* nil)
(defvar *current-gnuplot-stream* nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; for temp files
(defvar *temp-file-path* #P"/tmp/gnuclot/")
(defvar *temp-file-index* 0)

;;; adapted from ltk
(defun open-process-stream (program args &optional (wt nil))
  "execute program with args a list containing the arguments passed to the program
   if wt is non-nil, the function will wait for the execution of the program to return.
   returns a two way stream connected to stdin/stdout of the program"
  #+clisp (declare (ignore wt))
  (let ((fullstring program))
    (dolist (a args)
      (setf fullstring (concatenate 'string fullstring " " a)))
    #+(or cmu scl)
    (let ((proc (run-program program args :input :stream :output :stream :wait wt
                             #+scl :external-format #+scl :utf-8)))
      (unless proc
        (error "Cannot create process."))
      (make-two-way-stream
       (ext:process-output proc)
       (ext:process-input proc))
      )
    #+clisp (let ((proc (ext:run-program program :arguments args :input :stream :output :stream :wait t)))
               (unless proc
                 (error "Cannot create process."))
               proc
               )
    #+sbcl (let ((proc (sb-ext:run-program program args :input :stream :output :stream :wait wt :search t)))
              (unless proc
                (error "Cannot create process."))
              #+ext-8859-1
              (make-two-way-stream 
               (sb-sys:make-fd-stream 
                (sb-sys:fd-stream-fd (process-output proc))
                :input t :external-format :iso-8859-1)
               (sb-sys:make-fd-stream 
                (sb-sys:fd-stream-fd (process-input proc))
                :output t  :external-format :iso-8859-1))
              #-ext-8859-1
              (values (make-two-way-stream
                       (process-output proc)
                       (process-input proc))
                      proc)
              )
    #+lispworks (system:open-pipe fullstring :direction :io)
    #+allegro (let ((proc (excl:run-shell-command
			    #+mswindows fullstring
			    #-mswindows (apply #'vector program program args)
			    :input :stream :output :stream :wait wt)))
		 (unless proc
		   (error "Cannot create process."))
		 proc
		 )
    #+ecl(ext:run-program program args :input :stream :output :stream
                           :error :output)
    #+openmcl (let ((proc (ccl:run-program program args :input
					    :stream :output :stream :wait wt :sharing :external)))
		 (unless proc
		   (error "Cannot create process."))
		 (make-two-way-stream
		  (ccl:external-process-output-stream proc)
		  (ccl:external-process-input-stream proc)))
    ))

;;; CMUCL, SCL, and SBCL, use a two-way-stream and the constituent
;;; streams need to be closed.
(defun close-process-stream (stream)
  "Close a 'stream open by 'open-process-stream."
  (when *debug*
    (format t "Closing stream: ~S~%" stream))
  (ignore-errors (close stream))
  #+(or cmu scl sbcl)
  (when (typep stream 'two-way-stream)
    (close (two-way-stream-input-stream stream) :abort t)
    (close (two-way-stream-output-stream stream) :abort t))
  nil)

;;; send "quit" and close the stream
(defun close-gnuplot-stream (stream)
  (format stream "quit~%")
  (force-output stream))
  

;;; open a new gnuplot session, or switch to an already exist one
(defun start-gnuplot (&optional (name +default-gnuplot-connection-name+))
  (setf *current-gnuplot-connection-name* name)
  (multiple-value-bind (val stat) (gethash *current-gnuplot-connection-name* *gnuplot-connections*)
    (if stat
        ;; connection is already there
        (progn (format t "Connection ~s is already open, switching to it...~%"
                            *current-gnuplot-connection-name*)
                    (setf *current-gnuplot-stream* val))
        ;; make a new connection
        (multiple-value-bind (stream proc) (open-process-stream *gnuplot-path* nil)
          (setf *current-gnuplot-stream* stream)
          (when *debug* (setf *current-gnuplot-proc* proc))
          (setf (gethash *current-gnuplot-connection-name* *gnuplot-connections*) stream)
          (format t "Connection ~s started!~%" *current-gnuplot-connection-name*)))))


;;; stop the current gnuplot connection by default, or stop the named one
(defun stop-gnuplot (&optional name)
  (flet ((switch-to-next-available-gnuplot-connection ()
           (with-hash-table-iterator (get-entry *gnuplot-connections*)
             (multiple-value-bind (more? key value) (get-entry)
               (if more? 
                   (progn (setf *current-gnuplot-connection-name* key)
                          (setf *current-gnuplot-stream* value))
                   (progn (format t "All gnuplot connections are closed!~%")
                          (setf *current-gnuplot-connection-name* nil)
                          (setf *current-gnuplot-stream* nil)))))))
    (if name
        ;; close the named gnuplot connection
        (multiple-value-bind (val stat) (gethash name *gnuplot-connections*)
          (if stat (progn (format t "Closing the connection ~s~%" name)
                          (close-gnuplot-stream val) (remhash name *gnuplot-connections*)
                          (switch-to-next-available-gnuplot-connection))
              (format t "Connection ~s does not exist!~%" name)))
        ;; close the current gnuplot connection
        (if *current-gnuplot-stream*
            (progn (format t "Closing the current connection ~s~%"
                           *current-gnuplot-connection-name*)
                   (close-gnuplot-stream *current-gnuplot-stream*)
                   (remhash *current-gnuplot-connection-name* *gnuplot-connections*)
                   (switch-to-next-available-gnuplot-connection))
            (format t "No connection is open right now.~%")))))


;;; stop all the connections
(defun stop-all-gnuplot ()
  (maphash #'(lambda (k v)
               (format t "Closing the connection ~s~%" k)
               (close-gnuplot-stream v)) *gnuplot-connections*)
  (setf *current-gnuplot-connection-name* nil)
  (setf *current-gnuplot-stream* nil)
  (clrhash *gnuplot-connections*)
  (format t "All the connections closed!~%"))
  

;;; list the names of connections
(defun list-gnuplot-connection-names ()
  (maphash #'(lambda (k v) (declare (ignore v)) (format t "~s~%" k)) *gnuplot-connections*))


;;; switch to another open connection
(defun switch-to-gnuplot (&optional name)
  (if name
      (multiple-value-bind (val stat) (gethash name *gnuplot-connections*)
        (declare (ignore val))
        (if stat (progn (format t "Switching to connection ~s~%" name)
                        (setf *current-gnuplot-connection-name* name)
                        (setf *current-gnuplot-stream* 
                              (gethash *current-gnuplot-connection-name* *gnuplot-connections*)))
            (format t "Connection ~s doesn't exist~%" name)))
      (format t "Current connection: ~s, Which connection do you want to switch to?~%"
              *current-gnuplot-connection-name*)))
  

;;; send command string to gnuplot
(defun format-gnuplot (str)
  (format t "~a~%" str)
  (format *current-gnuplot-stream* "~a~%" str)
  (force-output *current-gnuplot-stream*)
  ;; read whatever is ready at the input
  (do ()
      ((not (listen *current-gnuplot-stream*)))
    (format t "~a~%" (read-line *current-gnuplot-stream*)))
  (clear-input *current-gnuplot-stream*))


;;; convert a list of objects to gnuplot command string
(defun list-to-gnuplot-command (list)
  "#P\"xxx\" and (\"xxx\") give additional quotation around xxx"
  (flet ((ttc (v)
           (cond ((pathnamep v) (concatenate 'string "\""
                                             (#+ccl ccl:native-translated-namestring
                                              #-ccl namestring v) "\" "))
                 ((listp v) (concatenate 'string "\"" (car v) "\" "))
                 ((symbolp v) (concatenate 'string (string-downcase (symbol-name v)) " "))
                 ((stringp v) (concatenate 'string v " "))
                 (t (concatenate 'string (format nil "~a" v) " ")))))
    (reduce #'(lambda (a b) (concatenate 'string a b)) (mapcar #'ttc list))))


;;; forming a string of command to be sent to gnuplot
(defmacro gnuplot-command (%command &rest %rest)
  (let ((list (cons %command %rest)))
    `(list-to-gnuplot-command (quote ,list))))


;;; function version
(defun gnuplot-command-fun (%command &rest %rest)
  (let ((list (cons %command %rest)))
    (list-to-gnuplot-command list)))


;;; list directly goes to gnuplot
(defun format-list-to-gnuplot-fun (list)
  (format-gnuplot (list-to-gnuplot-command list)))


;;; the "plot" command
(defmacro plot (&rest %rest)
  "use #P\"filename\" to specify file, use (\"xxx\") to have additional quotation around string."
  `(format-gnuplot (gnuplot-command "plot" ,@%rest)))
;;; the "splot" command
(defmacro splot (&rest %rest) `(format-gnuplot (gnuplot-command "splot" ,@%rest)))
;;; the "replot" command
(defmacro replot (&rest %rest) `(format-gnuplot (gnuplot-command "replot" ,@%rest)))
;;; the "set" command
(defmacro set (&rest %rest) `(format-gnuplot (gnuplot-command "set" ,@%rest)))
;;; the "unset" command
(defmacro unset (&rest %rest) `(format-gnuplot (gnuplot-command "unset" ,@%rest)))
;;; the "show" command
(defmacro show (&rest %rest) `(format-gnuplot (gnuplot-command "show" ,@%rest)))
;;; the "reset" command
(defmacro reset (&rest %rest) `(format-gnuplot (gnuplot-command "reset" ,@%rest)))
;;; the "help" command
(defmacro help (&rest %rest) `(format-gnuplot (gnuplot-command "help" ,@%rest)))
;;; the "load" command
(defmacro load (&rest %rest) `(format-gnuplot (gnuplot-command "load" ,@%rest)))


;;; function version
(defun plot-fun (&rest %rest)
  (format-gnuplot (apply #'gnuplot-command-fun (cons "plot" %rest))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; for temp files
(defun delete-all-temp-files ()
  (when (probe-file *temp-file-path*)
    (let ((dir (directory (make-pathname :name :wild :type :wild :defaults *temp-file-path*))))
      (mapcar #'(lambda (a) (format t "Deleting ~s~%" a) (delete-file a)) dir))))

(defun generate-temp-filename ()
  (setf *fname*
        (do ((fname #P"/")) ((not (probe-file fname)) fname)
          (setf fname (merge-pathnames *temp-file-path* 
                                       (make-pathname :name (format nil "clot~6,'0D" ;; sequential
                                                                    (1- (incf *temp-file-index*)))
                                                      ;; (format nil "clot~36R" ;; random
                                                      ;;           (random (expt 36 6)))
                                                      :type "dat"))))))


(defmacro with-temp-file ((%stream) &rest %rest)
  `(let ((fname (generate-temp-filename)))
     (with-open-file (,%stream (ensure-directories-exist fname) :direction :output)
       (progn ,@%rest))
     fname))

;;; use a temp file to store the list and plot it
(defmacro plotlist (%cmd %list &rest %rest)
  `(let* ((l ,%list)
          (fname (with-temp-file (of)
                   (dolist (row l)
                     (if (listp row)
                         (mapcar #'(lambda (a) (format of "~24,16,3,,,,'eE " a)) row)
                         (format of "~24,16,3,,,,'eE " row))
                     (format of "~%")))))
     (format t "Output to file ~s~%" fname)
     (format-list-to-gnuplot-fun (append (list ',%cmd fname) ',%rest))))


;;; use a temp file to store the array and plot it
(defmacro plotarray (%cmd %ary &rest %rest)
  "Expecting 1-dim or 2-dim arrays
   %cmd is given so you can choose to use plot or splot"
  `(let* ((a ,%ary)
          (dim (array-dimensions a))
          (m (car dim))
          (n (or (cadr dim) 1))
          (fname (with-temp-file (of)
                      (dotimes (i m)
                        (if (= n 1) (format of "~24,16,3,,,,'eE " (aref a i))
                            (dotimes (j n)
                              (format of "~24,16,3,,,,'eE " (aref a i j))))
                        (format of "~%")))))
     (format t "Output to file ~s~%" fname)
     (format-list-to-gnuplot-fun (append (list ',%cmd fname) ',%rest))))


;;; plot a 2D array using the `snake' scan line format
(defmacro plotarray2 (%cmd %ary-l &rest %rest)
  `(multiple-value-bind (h l c) ,%ary-l
     (let* ((dims (array-dimensions h))
            (nm (car dims))
            (nn (cadr dims))
            (mmin (car l)) (dy (cadr l)) (nmin (caddr l)) (dx (cadddr l))
            (fname (with-temp-file (of)
                     (dotimes (m nm)
                       (dotimes (n nn)
                         (format of "~24,16,3,,,,'eE ~24,16,3,,,,'eE ~24,16,3,,,,'eE~%"
                                 (+ nmin (* dx n)) (+ mmin (* dy m)) (aref h m n)))
                       (format of "~%")))))
       (format t "Output to file ~s, count: ~d~%" fname c)
       (format-list-to-gnuplot-fun (append (list ',%cmd fname) ',%rest)))))

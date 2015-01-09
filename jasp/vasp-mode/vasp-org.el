;;; vasp-org.el --- Org utilities for vasp

;;; Commentary:
;; 

;;; Code:

(defun vaspdoc (keyword arg)
  "Lookup vasp documentation for KEYWORD.
Looks at http://cms.mpi.univie.ac.at/wiki/index.php.
With prefix arg, search in google."
  (interactive "sKeyword: \nP")
  (if arg
      (eww (format
		   "http://www.google.com/search?q=%s"		   
		   (url-hexify-string (format "vasp %s" keyword))))
    (eww (format "http://cms.mpi.univie.ac.at/wiki/index.php/%s"  (s-upcase keyword )))))

;; these are just links to python documentation
(org-add-link-type
 "mod"
 (lambda (arg)
   (pydoc arg))
 (lambda (path desc format)
   (cond
    ((eq format 'latex)
     (format "\\texttt{%s}" path)))))

(org-add-link-type
 "func"
 (lambda (arg)
   (pydoc arg))
 (lambda (path desc format)
   (cond
    ((eq format 'latex)
     (format "\\texttt{%s}" path)))))

;; link to incar documentation
(org-add-link-type
 "incar"
  (lambda (keyword)
    (shell-command (format "firefox http://cms.mpi.univie.ac.at/wiki/index.php/%s" keyword) nil))
  ; this function is for formatting
  (lambda (keyword link format)
   (cond
    ((eq format 'html)
     (format "<a href=http://cms.mpi.univie.ac.at/wiki/index.php/%s>%s</a>" keyword keyword))
    ((eq format 'latex)
     (format "\\href{http://cms.mpi.univie.ac.at/wiki/index.php/%s}{%s}"  keyword keyword)))))

(provide 'vasp-org)

;;; vasp-org.el ends here

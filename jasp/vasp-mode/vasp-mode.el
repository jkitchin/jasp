;; vasp-mode for syntax highighting of input files
;; John Kitchin 3/12/2012
;;
;; TODO
;; - fontify output files
;; - figure out how to get help on each keyword, preferrably through info

;; set these files to open in vasp-mode
(add-to-list 'auto-mode-alist '("\\INCAR\\'" . vasp-mode))
(add-to-list 'auto-mode-alist '("\\POSCAR\\'" . vasp-mode))
(add-to-list 'auto-mode-alist '("\\KPOINTS\\'" . vasp-mode))
(add-to-list 'auto-mode-alist '("\\POSCAR\\'" . vasp-mode))

;; these are the input for INCAR
(defvar vasp-keywords
  '("NGX" "NGY" "NGZ" "NGXF" "NGYF" "NGZF"
    "NBANDS" "NBLK"  "SYSTEM" "NWRITE" "ISTART"
    "ICHARG" "ISPIN" "MAGMOM" "INIWAV" "ENCUT"
    "PREC" "NELM" "NELMIN" "NELMDL" "EDIFF" "EDIFFG"
    "NSW" "NBLOCK" "KBLOCK" "IBRION" "ISIF" "IWAVPR"
    "ISYM" "SYMPREC" "LCORR" "POTIM" "TEBEG" "TEEND"
    "SMASS" "NPACO" "APACO" "POMASS" "ZVAL" "RWIGS"
    "NELECT" "NUPDOWN" "EMIN" "EMAX" "ISMEAR" "SIGMA"
    "ALGO" "IALGO" "LREAL" "ROPT" "GGA" "VOSKOWN"
    "DIPOL" "AMIX" "BMIX" "WEIMIN" "EBREAK" "DEPER"
    "TIME" "LWAVE" "LCHARG" "LVTOT" "LVHAR" "LELF"
    "LORBIT" "NPAR" "LSCALAPACK" "LSCALU" "LASYNC"
    ))

(defvar vasp-other-words
  '(; POSCAR
    "Cartesian" "Direct"
    "Selective dynamics" "T" "F"
    ; KPOINTS
    "Tetrahedra" "T"
    "Automatic mesh" "A"
    "Monkhorst-Pack" "M" "m"
    "Gamma" "G" "g"))

;; chemical symbols to be highlighted
(defvar vasp-element-words
  '("H"  "He" "Li" "Be" "B"
    "C"  "N"  "O"  "F"
    "Ne" "Na" "Mg" "Al" "Si"
    "P"  "S"  "Cl" "Ar" "K"
    "Ca" "Sc" "Ti" "V"  "Cr"
    "Mn" "Fe" "Co" "Ni" "Cu"
    "Zn" "Ga" "Ge" "As" "Se"
    "Br" "Kr" "Rb" "Sr" "Y"
    "Zr" "Nb" "Mo" "Tc" "Ru"
    "Rh" "Pd" "Ag" "Cd" "In"
    "Sn" "Sb" "Te" "I"  "Xe"
    "Cs" "Ba" "La" "Ce" "Pr"
    "Nd" "Pm" "Sm" "Eu" "Gd"
    "Tb" "Dy" "Ho" "Er" "Tm"
    "Yb" "Lu" "Hf" "Ta" "W"
    "Re" "Os" "Ir" "Pt" "Au"
    "Hg" "Tl" "Pb" "Bi" "Po"
    "At" "Rn" "Fr" "Ra" "Ac"
    "Th" "Pa" "U"  "Np" "Pu"
    "Am" "Cm" "Bk" "Cf" "Es"
    "Fm" "Md" "No" "Lr"))

(setq all-keywords (append vasp-keywords vasp-other-words))

;; auto-generate the regexp to syntax highlight these words
(defvar vasp-keywords-regexp (regexp-opt all-keywords 'words))
(defvar vasp-elements-regexp (regexp-opt vasp-element-words 'words))

;; list for each font-lock
(setq vasp-font-lock-keywords
      `((,vasp-keywords-regexp . font-lock-keyword-face)
	(,vasp-elements-regexp . font-lock-comment-face)))


;; clear from memory since we dont need these keywords anymore
;; (setq vasp-keywords nil) ;

;; the command to comment/uncomment text
(defun vasp-comment-dwim (arg)
"Comment or uncomment current line or region in a smart way.
For detail, see `comment-dwim'."
   (interactive "*P")
   (require 'newcomment)
   (let ((deactivate-mark nil)
	 (comment-start "#")
	 (comment-end ""))
     (comment-dwim arg)))

;;
(defun vasp-help ()
  "Open a browser window showing documentation for the word under the point.
Defaults to a simple keyword search.
Uses `search-site-url' to do the actual search.
"
  (interactive)
  (require 'url)
  (browse-url
   (apply 'search-site-url
          (thing-at-point 'symbol)
          '("http://cms.mpi.univie.ac.at/vasp/vasp/")
            )))

(defun vasp-ag ()
  "runs ag on the POSCAR file in current directory"
  (interactive)
  (shell-command "ag POSCAR"))

(defun vasp-run ()
  "runs vasp"
  (interactive)
  (start-process-shell-command "vasp" "*scratch*" "vasp"))

(defun vasp-queue ()
  "submits job to a queue system"
  (interactive)
  (shell-command "echo \"pwd; vasp\" | qsub -d `pwd` -l walltime=12:00"))

;; define a menu for vasp-mode
(defvar vasp-mode-map nil "Keymap for vasp-mode")

(when (not vasp-mode-map) ;this is only run when this variable does not exist
  (setq vasp-mode-map (make-sparse-keymap))
  (define-key vasp-mode-map (kbd "C-c C-h") 'vasp-help)
  (define-key vasp-mode-map (kbd "C-c C-v") 'vasp-ag)
  (define-key vasp-mode-map (kbd "C-c C-r") 'vasp-run)
  (define-key vasp-mode-map (kbd "C-c C-q") 'vasp-queue)

;; define the menu
  (define-key vasp-mode-map [menu-bar] (make-sparse-keymap))

  (let ((menuMap (make-sparse-keymap "Vasp")))
    (define-key vasp-mode-map [menu-bar vasp] (cons "Vasp" menuMap))
    (define-key menuMap [Queue] '("Queue VASP" . vasp-queue))
    (define-key menuMap [Run] '("Run VASP" . vasp-run))
    (define-key menuMap [View] '("View POSCAR" . vasp-ag))
    (define-key menuMap [Help] '("Help on VASP keyword" . vasp-help))
    ))

;; define the major mode
(define-derived-mode vasp-mode fundamental-mode
"vasp mode"
"Major mode for editing vasp input files and reading output file"

(setq font-lock-defaults '(vasp-font-lock-keywords))

; setup comment style
(define-key vasp-mode-map [remap comment-dwim] 'vasp-comment-dwim)
;; comment style: #...
  (modify-syntax-entry ?# "< b" vasp-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" vasp-mode-syntax-table)
) ; end major-mode definition

; completion see: http://www.masteringemacs.org/articles/2012/01/16/pcomplete-context-sensitive-completion-emacs/

(provide 'vasp-mode)

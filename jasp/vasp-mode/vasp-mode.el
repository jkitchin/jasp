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

;; http://davidavraamides.net/blog/2008/07/22/mode-aware-google-help-in-emacs/
(defun search-site-url (keyword &optional site inurl lucky)
  "Do a Google search for KEYWORD. Restrict to SITE and INURL, if specified.
Jump to best match (I Feel Lucky) if LUCKY set.
"
  (concat "http://www.google.com/"
          (format "search?q=%s" (url-hexify-string keyword))
          (if site (format "+site:%s" (url-hexify-string site)))
          (if inurl (format "+inurl:%s" (url-hexify-string inurl)))
          (if lucky "&btnI")))

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
  (shell-command "ase-gui POSCAR"))

(defun vasp-run ()
  "runs vasp"
  (interactive)
  (start-process-shell-command "vasp" "*scratch*" "vasp"))

(defun vasp-queue ()
  "submits job to a queue system"
  (interactive)
  (shell-command "echo \"pwd; vasp\" | qsub -d `pwd` -l walltime=12:00"))

(defun vasp-nval (functional chemical-symbol)
  "return number of valence electrons for the POTCAR
functional can be lower case. chemical-symbol must be correct"
  (interactive "sFunctional: \nsChemical Symbol: ")
  (shell-command (format "awk 'NR==2{print;exit}' %s" 
			 (format "%s/potpaw_%s/%s/POTCAR"
				 (getenv "VASP_PP_PATH")
				 (upcase functional)
				 chemical-symbol))))
    
;; define a menu for vasp-mode
(defvar vasp-mode-map
  (let ((map (make-sparse-keymap)))
    (define-key map (kbd "C-c C-h") 'vasp-help)
    (define-key map (kbd "C-c C-v") 'vasp-ag)
    (define-key map (kbd "C-c C-r") 'vasp-run)
    (define-key map (kbd "C-c C-q") 'vasp-queue)
    map)
  "Keymap for vasp-mode")

(defvar vasp-mode-hook nil "*List of functions to call when entering vasp mode.")

(require 'easymenu)

(easy-menu-define my-menu vasp-mode-map "Vasp menu"
  '("Vasp-mode"
    ["Queue VASP" vasp-queue]
    ["run VASP" vasp-run]
    ["View POSCAR" vasp-ag]
    ["Help on VASP keyword" vasp-help]))

;; define the major mode
(define-derived-mode vasp-mode text-mode
  "vasp mode"
  "Major mode for editing vasp input files and reading output file"

  (setq font-lock-defaults '(vasp-font-lock-keywords))

  ;; setup comment style
  (define-key vasp-mode-map [remap comment-dwim] 'vasp-comment-dwim)
  ;; comment style: #...
  (modify-syntax-entry ?# "< b" vasp-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" vasp-mode-syntax-table)
  ) ; end major-mode definition

; completion see: http://www.masteringemacs.org/articles/2012/01/16/pcomplete-context-sensitive-completion-emacs/

(provide 'vasp-mode)

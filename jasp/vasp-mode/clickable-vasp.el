(let* ((keywords (with-current-buffer
		     (find-file-noselect
		      (expand-file-name
		       "vasp.org"
		       (file-name-directory
			(locate-library "vasp-mode"))))
		   (org-map-entries
		    (lambda ()
		      (nth 4 (org-heading-components))))))
       (contents (with-current-buffer
		     (find-file-noselect
		      (expand-file-name
		       "vasp.org"
		       (file-name-directory
			(locate-library "vasp-mode"))))
		   (org-map-entries
		    (lambda ()
		      (save-restriction
			(org-narrow-to-subtree)
			(buffer-string))))))
       (cons-list (loop for key in keywords for cont in contents
			collect (cons key cont)))
       (KEYWORDS (mapcar 'upcase keywords)))

  (button-lock-set-button
   (regexp-opt (append keywords KEYWORDS))
   nil
   :additional-property 'vasp
   :help-echo `(lambda (window object position)
		(goto-char position)
		(let* ((elem (downcase
			      (get-surrounding-text-with-property
			       'vasp))))
		  (with-current-buffer
		      (find-file-noselect
		       (expand-file-name
			"vasp.org"
			(file-name-directory
			 (locate-library "vasp-mode"))))
		    (org-open-link-from-string (format  "[[*%s]]" elem))
		    (save-restriction
		      (org-narrow-to-subtree)
		      (buffer-string)))))))

(provide 'clickable-vasp)

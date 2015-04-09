pro aprates, C, B, A_s, A_b, alpha = alpha, beta = beta, CL = CL, plot =plot, $
	s_bar = s_bar , hi_int = hi_int, lo_int = lo_int
;; This is an IDL version of CIAO task aprates
;; By Rui She

;; Input 
;;     C, source count
;;     B, background count
;;     A_s, area of source aperture   
;;     A_b, area of background aperture
;;	   alpha, psf fraction in source aperture, default is 1
;;	   beta, psf fraction in background aperture, default is 0
;;	   CL, is confidence level

;; Output: 
;;	   s_bar, source intensity expectation
;;	   hi_int, higher confidence interval
;;	   lo_int, lower confidence interval

if not keyword_set(alpha) then alpha = 1d
if not keyword_set(beta) then beta = 0d
if not keyword_set(CL) then cl = 0.68d

C = double(C)
B = double(B)
A_s = double(A_s)
A_b = double(A_b)

r = A_b/A_s
f = alpha
g = beta

;;;; there is equation:
;;;; C = f * s + b
;;;; B = g * s + r * b

s_MLE = (r * C - B) / (r * f - g) ;;  MLE method source intensity

if s_MLE gt 50 then begin;; use Gaussian pdf
	print, "Use Gaussian PDF..."
	sig_s = sqrt((r^2*C+B)/(r*f-g)^2)
	al = (1 - cl)/2d
	nsig = gauss_cvf(al)
	lo_int = s_MLE - nsig * sig_s
	hi_int = s_MLE + nsig * sig_s
	print, format = '("Expectation of source intensity is ",F0.3, ", with confidence interval of [", F0.3,", ", F0.3, "], at a confidence level of ", F0.3)', s_MLE, lo_int, hi_int, CL

endif else begin   ;; use bayesian estimate
	Print, "Use Bayesian estimate..."
	;; range to calculate pdf
	s = indgen(1000)/1000d  * ( 6 * s_MLE ) 

	;; calculation pdf
	p_sum = 0d
	for k = 0L, C do begin
		for j = 0L, B do begin
			ps_jk = f^k * g^j * r^(B-j) * s^(k+j) * exp(-s*(f+g)) * $
				(gamma(C+1)*gamma(B+1)*gamma( C+B -k -j +1) ) / $
				(gamma(k+1)*gamma(C-k+1)*gamma(j+1)*gamma(B-j+1)* (1+r)^(C+B-k-j+1) )
			P_sum += ps_jk
		endfor
	endfor

	;; prabilities density function of s 
	p_cb = ( r * f - g ) / (gamma(C+1)* gamma(B+1))
	p_s = p_sum * p_cb

	;; check for converging
	sump = total(p_s) *(s[1]-s[0])
	if abs(sump-1) ge 5d-3 then message, "Not Converging!"

	;; calculate expectation of s
	s_bar = total(s*p_s) * (s[1]-s[0])

	;; cumulative distribution of p_s
	cdf = total(p_s, /cumulative) * (s[1]-s[0])

	;; confidence interval
	;;; two ways to give confidence interval, and symmetry confidence level is prefer
	a = (1-CL)/2d
	lo_cl = a 
	hi_cl = 1d - a
	lo_int = interpol(s,cdf,lo_cl)
	hi_int = interpol(s,cdf,hi_cl)

	if keyword_set(plot) then begin
		plot, s, p_s, xtitle="Source intensity", ytitle="P"
		oplot, [lo_int,lo_int], [0,1]
		oplot, [hi_int,hi_int], [0,1]
	endif

	print, format = '("Expectation of source intensity is ",F0.3, ", with confidence interval of [", F0.3,", ", F0.3, "], at a confidence level of ", F0.3)', s_bar, lo_int, hi_int, CL

endelse

end

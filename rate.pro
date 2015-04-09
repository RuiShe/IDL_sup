function rate, x, n , cl = cl
;; by Rui She
;; This procedure is a Bayesian approach calculating the ratio (incidence-rate) at some confidence levels
;; To estimate the ratio in fact is to estimate the parameter p of a binomial distribution. Beta distribution is the distribution of parameter p in binomial distribution. 

;; n is size of total sample, can be an array
;; x is size of part sample, can be an array
;; cl is Confidence Level, default is 0.90

ns = n_elements(x) 
if ns ne n_elements(n) then message, 'Dimension of X and N should be the same.'

if not keyword_set(cl) then cl = 0.9d
hicut = (1d + cl) / 2d
locut = (1d - cl) / 2d

para = indgen(1000)/1000d  ;; p in [0,1] 

;;## use uniform prior distribution 
_alpha = 1 & _beta = 1  ;; this is for a uniform prior distribution

res = dblarr(2, ns)

for i =0L, ns -1 do begin
	if n[i] lt x[i] or x[i] lt 0 then message, 'This pro requires: X le N, and X ge 0'
	level = ibeta( x[i] + _alpha , n[i] - x[i] + _beta , para )

	res[0,i] = interpol(para, level, locut)
	res[1,i] = interpol(para, level, hicut)
endfor

return, res
end

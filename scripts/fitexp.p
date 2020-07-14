## gnuplot script for quick exp fittining

set terminal png size 800,800
set output 'massfit.png'

# store errors for fit params in m_err etc
set fit errorvariables 

f(x) = b + a*exp(-m*x)
m=0.2

fit f(x) 'correl_avg' u 1:4 every ::3 via a,b,m

# print the fitted mass in the plot (need to set before plotting)
set label GPFUN_f at graph .05,.95
set label sprintf("a = %g +/- %g", a, a_err) at graph .05,.90 
set label sprintf("a = %g +/- %g", b, b_err) at graph .05,.85
set label sprintf("m = %g +/- %g", m, m_err) at graph .05,.80

plot 'correl_avg' u 1:4:5 every ::2 with errorbars, f(x)

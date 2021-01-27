# trapezoids

Code for trapezoid model of oscillating mRNA signals. See plots/trapezoid_params_schematic.pdf for param details.

main functions:  
/code/functions/make_trapezoid_signal.m : make the trapezoid
/code/functions/integrate_trapezoid_signal.m : integrate trapezoid with a decay term to get accumulated mRNA.

scripts:  
/code/scripts/test_make_trapezoid_signal : for 1 set of params compute and plot trapezoid and accumulated mRNA.
/code/scripts/vary_trapezoids: for a range of parameters compute and plot trapezoids and accumulated mRNA.

To Do:  
protein signal

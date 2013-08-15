
Explanation of the different files:
-----------------------------------

new_main.m:
  the main driver for the matlab code


loadquadrature.m:
  loads angular quadratures (Gauss Legendre) for 1D slab geometry
  this functons calls the Gauss-Lendre quadrature function, GLNodeWt.m

GLNodeWt.m
  creates a Gauss-Legendre quadrature

loadmydata.m:
  loads the cross sections (tot and sca)
        the boundary (qsa) and volumetric (qva) angular sources
        the mesh


build_matrices.m: 
  constructs the linear system matrix by calling
  (1) compute_elem1.m to build the elementary matrices
  (2) compute_T1.m to assemble the global system matrix

compute_elem1.m 
  builds the elementary matrices
  this routine calls comp_gamma_delta.m

comp_gamma_delta.m
  computes the gamma and delta stabilization parameters

compute_T1.m 
  assembles the global system matrix
   
direct_solve_psi.m:
  solves for the angular flux solution 
  and calls compu_phi to compute the scalar flux

compu_phi.m:
  computes the scalar flux  

myplot.m
  customized plotting


-----------------------------------
How to run the code ???
-----------------------------------

--> simply type: new_main in the matlab console

By default, DG1 with standard upwind is called first 
(this is usually a very good estimates of the true answer)
and then DG(porder) with one's choice for stabilization
is called next. Both results are plotted.

-----------------------------------
What options can I play with ???
-----------------------------------
In new_main.m:

line 6: porder: either 0 or 1 for DG0 or DG1

line 8: sn, the angular quadrature order. must be an even number

line 16: dataID, this is the ID of the test case
dataID = 1,...,5 are the 5 problems treated in the 2011 JCP paper
dataID=6 is a new 3-region problem
No other datID's are code but you can easily add one by looking at loadmydata.m

line 20 and 21
the delta_0 and gamma_0 values of the stabilization parameters


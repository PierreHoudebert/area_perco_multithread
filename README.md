# area_perco_multithread
Percolation of the area-interaction model


 1) This program want to approximate the percolation threshold of the area interaction process, see for instance
 Dereudre & Houdebert,  Sharp phase transition for the continuum Widom-Rowlinson model
 for a description of this Gibbs Point process.
 2) The process is sampled using a FK representation with the generalized continuum random cluster model (grcm) , see
 Houdebert,  Phase transition of the non-symmetric Continuum Potts model
 for a description of the grcm
 3) The grcm is sampled using a MCMC birth and death algorithm, see the book
 Moller & Waagepetersen, Statistical Inference and Simulation for Spatial Point Processes, section 7.1.2
 for a description of the algorithm
 
 The program does:
 
 a) sampling of the grcm using MCMC algorith, with a number of step given by the parameter "Number_MCMC" line 588
 b) remove some cluster according to the FK representation to get an area-interaction configuration
 c) test the connectivity of the center of the box to the boundary
 d) do this a number of time to get a Law of Large Number approximation of the probability, given by the parameter Number_LLN line 589
 e) increment the activity parameter and do it again, in order to have beta -> percolation_threshold(beta)
 
 This program uses the openmp parallel for loop.
 
 
 Parameters
 
 - max_p, line 60: maximum size of the sampled configuration. If too small the program crashes, if too large the program is slower
 - max_n, line 62: maximum number of neighbors of a point. same as before
 - A and B, line 583: size of the box ([0,A]x[0,B])
 - beta, line 585: inverse temperature of area
 - zz, line 586: initial activity of area (we are incrementing the activity)
- zz_max, line 588: maximal activity;
- step, line 589: increment of the activity;
- Number_MCMC, line 592: number of MCMC iteration
- Number_LLN, line 593:   number of iteration Law Large Number
 
 Output:
 
 two txt files
- "parameters.txt" containing all the parameters;
- "results.txt" containing on each line, the acticvity, the percolation probability and the experimental intensity
 
To compile the code
 g++ -fopenmp -O1 -march=native grcm_area_perco.cpp  (-O1 may be replace by -O2 or -O3)
 
 To run  ./a.out

#INSTRUCTIONS
NSVD_Var.m calculates the singular values of the Khatri-Rao products of the PARAFAC factor matrices. Also it calculates intermediate estimates of CORCONDIA, MSE and Missing Data MSE.

NSVD_demo.m calculates the final estimates of NSVD, CORCONDIA, MSE and Missing Data MSE, and generates the corresponding plots.

Note that the code currently works only for 3-mode tensors, although it is fairly straightforward to extend to higher dimensions.

----------------------
---- Dependencies ----
----------------------
1)Tensor Toolbox: can be downloaded at https://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html.

2)efficient_corcondia.m: an efficient implementation of the Core Consistency Diagnostic as shown in [1] which can be found at https://www.cs.ucr.edu/~epapalex/src/efficient_corcondia.zip.

3)N-way Toolbox: required for calculating Missing Data MSE, and can be downloaded at https://www.mathworks.com/matlabcentral/fileexchange/1088-the-n-way-toolbox.

4)Matlab R2017a: the Matlab 'filloutliers' method is utilized.





[1] Papalexakis, Evangelos E., and Christos Faloutsos. "Fast efficient and scalable core consistency diagnostic for the parafac decomposition for big sparse tensors." 2015 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP). IEEE, 2015.
 

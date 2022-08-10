#
#  SPCAfunction2.py
#  
#    
#  August 10, 2022
#
#  Performs Sparse PCA of numeric matrix using 
#  the SparsePCA function in sklearn-decomposition
#
#   Modules needed for the SPCA function
#
import numpy as np
import pandas as pd
#
from sklearn.decomposition import SparsePCA
# 
#  **** SPCA function *******
#
def SPCA( PSIdf, *args, nPCs=15, s_alpha=0.8, max_it=400, it_tol=1e-04 ):
#
#   PSIdf is a pandas data frame with events (features) as columns, no missing values.  
#   So, PSI matrix may need to be transposed before calling this ftn.
#
#   Event names are assumed to be the row names, not in a column.
#
#   The following arguments can be set in the call: 
#   nPCs (15), number of PCs to estimate
#   s_alpha (0.8), lasso loading shrinkage parameter (larger leads to sparser)
#   it_tol (1e-04), objective function convergence criterion value.
#   max_iter (400), maximum number of iterations.
#
#   Possible enhancement: check types & values of incoming arguments
#
#
#    import numpy as np
#    import pandas as pd
#
#    from sklearn.decomposition import SparsePCA
#
# -- Make numpy array with the event names/codes
#
    EVnames = np.array(PSIdf.columns.values)
#
# -- Mean-center each column (event)
#
    PSIdf = PSIdf - PSIdf.mean(axis=0)
#   chk=PSIdf.mean(axis=0)   #average by column
#   should be very small values, length=# of events
#
#
# -- SparsePCA execution
#
#   Set up the SPCA object with parameter values
    transformer1 = SparsePCA(n_components=int(nPCs), alpha=s_alpha, ridge_alpha=0, max_iter=int(max_it), tol=it_tol, method='lars',verbose=True, random_state=0)
#
#   Perform the SPCA fitting to the PSI values
    transformer1.fit(PSIdf.values)
#
#   Possible enhancement: wrap .fit call in try() or the like to handle error(s)
#
# -- Extract the matrix of loadings
#
    res1=transformer1.components_
#   res1.shape
#
#  -- Identify events that have non-0 loading in any of the PCs
#
#   Binarize
#
    res1[res1 != 0] = 1
#  
#  Extract corresponding names and return as pandas Series
#
    return pd.Series( EVnames[np.any(res1,axis=0)] )


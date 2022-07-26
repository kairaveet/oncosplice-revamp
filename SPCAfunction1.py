#
#  SPCAfunction1.py
#  
#    
# July 26, 2022
#
#  Performs Sparse PCA of numeric matrix using 
#  the SparsePCA function in sklearn-decomposition
#
#  import:  from SPCAfunction1 import SPCA
# 
# **** SPCA function *******
#
def SPCA( PSIdf=None, SPCAparms=None ):
#
#   PSIdf is a pandas data frame with events (features) as columns, no missing values.  
#   So, PSI matrix may need to be transposed when calling this ftn.
#
#   Event names are assumed to be the row names, not in a column.
#
#   SPCAparms is a pandas series with entries for parameters (default value): 
#   nPCs (15), number of PCs to estimate
#   s_alpha (0.8), lasso loading shrinkage parameter (larger leads to sparser)
#   it_tol (1e-04), objective function convergence criterion value.
#   max_iter (400), maximum number of iterations.
#
#
    import numpy as np
    import pandas as pd
#
    from sklearn.decomposition import SparsePCA
#
# -- Make array with the event names/codes
#
    EVnames = np.array(PSIdf.columns.values)
#
# -- Mean-center each column (event)
#
    PSIdf = PSIdf - PSIdf.mean(axis=0)
#   chk=PSIdf.mean(axis=0)   #average by column
#   should be very small values, length=# of events
#
# -- SparsePCA execution
    transformer1 = SparsePCA(n_components=int(SPCAparms['PCs']), alpha=SPCAparms['s_alpha'], ridge_alpha=0, max_iter=int(SPCAparms['max_it']), tol=SPCAparms['it_tol'], method='lars',verbose=True, random_state=0)
#
    PSInp = pd.DataFrame.to_numpy(PSIdf)
#
    transformer1.fit(PSInp)
#
# -- Further process the matrix of loadings
#
    res1=transformer1.components_
#   res1.shape
#
#   Identify events that have non-0 loading in >0 PCs
#   Binarize, then sum down columns
    res1[res1 != 0]=1
    event_sums = np.sum(res1,axis=0)
#
    res1TF=event_sums !=0
#
# -- Extract event names and return
    keptEvents =  pd.Series( EVnames[res1TF] )
#
    return keptEvents


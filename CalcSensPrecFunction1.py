#
#  CalcSensPrecFunction1.py
# 
#  CalcSensPrec function 
#    
# July 27, 2022
# 
# **** SPCA function *******
#
def CalcSensPrec( target = None, qry = None ):
#
#   target is a pandas Series containing the target (reference) event names  
#   qry is a pandas Series containing the event names to assess vs. target
#
#   Returns a pandas series with the following elements:
#     nt = number of elements in target
#     nq = number of elements in query
#     nOverlap = number of elements in intersection
#     sens = nOverlap / nt  -- proportion, not percent
#     prec = nOverlap / nq  -- proportion, not percent
#
#
# ** Could check that target & qry are type=pd.series & that each has size > 0 ***
#
    import numpy as np
    import pandas as pd
#
#   Set up pandas series for summary stats.
    res={ 'nt' : np.nan, 'nq' : np.nan, 'nOverlap' : np.nan , 'sens' : np.nan, 'prec' : np.nan }
#
    res = pd.Series(res)
#
    nt = target.size
    nq = qry.size
    tq=pd.Series( np.intersect1d(target, qry) )
    nOverlap = tq.size
    sens = float( nOverlap / nt )
    prec = float( nOverlap / nq )
#
    res['nt'] = nt
    res['nq'] = nq
    res['nOverlap'] = nOverlap
    res['sens'] = sens
    res['prec'] = prec
#
    return res






#
#  IDKfunction1.txt
# 
#  IDK function 
#    
# July 25, 2022
#
#  Can identify K using 1 or more of the following methods:
#   MADA, CORRINT, MOM, PCAF0, TLE
# **** IDK function *******
def IDK( PSIdf, Kmethods ):
#   Set up pandas series for estimates of K.  -1 for not done or error.
#   PSIdf is a pandas data frame with events (features) as columns
#   event names are not needed; can be row names but not a column.
#   Kmethods is a pandas series with True/False for each available method.
#   Available methods are: MADA, CORRINT, MOM, PCAF0, TLE.
#
    import skdim
    import numpy as np
    import pandas as pd 

    Kests={ 'MADA' : np.nan, 'CORRINT' : np.nan, 'MOM' : np.nan , 'PCAF0' : np.nan, 'TLE' : np.nan }
#
    Kests = pd.Series(Kests)
#
    if Kmethods['MADA'] :
        mada = skdim.id.MADA().fit(PSIdf)
        Kests['MADA'] = mada.dimension_
# 
    if Kmethods['CORRINT'] :
        corrint = skdim.id.CorrInt().fit(PSIdf)
        Kests['CORRINT'] = corrint.dimension_
#
    if Kmethods['MOM'] :
        Mom = skdim.id.MOM().fit(PSIdf)
        Kests['MOM'] = mom.dimension_
#
    if Kmethods['PCAF0'] :
        pcaf0 = skdim.id.lPCA(ver='FO').fit(PSIdf)
        Kests['PCAF0'] = pcaf0.dimension_
#
    if Kmethods['TLE'] :
        tle=skdim.id.TLE(epsilon=0.0001).fit(PSIdf)
        Kests['TLE'] = tle.dimension_
#
    return Kests

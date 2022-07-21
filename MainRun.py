
# Load required libraries
import pandas as pd

# Read PSI file
psi_file_path = "blah/blah"
formatted_psi_file, metadata = format_psi_matrix(psi_file_path)

# Find the non-redundant splicing events

# remove_redundant_splicing_events_random_event() -- we are not a fan of this for this step but just fyi

corr_threshold = 0.8 #change this as needed
list_of_events = remove_redundant_splicing_events(formatted_psi_file, corr_threshold, metadata)

# Subset the PSI matrix to only the non-redundant events
formatted_psi_file_1 = formatted_psi_file[list_of_events, ]

# K estimation (primarily used for determining PCs for SPCA or PCA-based feature selection; could be used for sample clustering for NMF)

## Here the input for k estimation is non-redundant events

### Now we apply different feature selection approaches ###

## Output for feature selection methods will be lists of events and any text files will bbe added here rather than inside the function

# SPCA

# PCA-based feature selection

# ICA-based feature selection

# Gene shaving

# K estimation (here it is primarily for NMF (sample clustering))

## Here the input for k estimation is feature selection output events

## At the end of this step, we will have inputs ready for the NMF algorithm



fn='/Users/tha8tf/Documents/testsOncosplice/BlockIDFile.txt'
BlockData=np.loadtxt(fn, dtype='float', comments='#', delimiter='\t', converters=None, skiprows=1, usecols=None)
BlockData.shape
# 159, 366

fn='/DS/OncoSplice/May22/SpliceICGSoutputFromNS/FullDataOnly.txt'
FullData=np.loadtxt(fn, dtype='float', comments='#', delimiter='\t', converters=None, skiprows=1, usecols=None)
FullData.shape
# 1051, 366

#   SHOULD HAVE FEATURES AS COLUMNS
BlockDataT=np.transpose(BlockData)
BlockDataT.shape
#. (366, 159)

FullDataT=np.transpose(FullData)
FullDataT.shape
# (366, 1051)



from RoundWrapper import *
from removeRedundantSplicingEvents import *

# Read PSI file
psi_file_path = "/Users/tha8tf/Documents/testsOncosplice/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-75p.txt"
formatted_psi_file, metadata = format_psi_matrix(psi_file_path)

# Find the non-redundant splicing events
corr_threshold = 0.8  # change this as needed ### this step may not have been done in Meenakshi's version
list_of_events = remove_redundant_splicing_events(formatted_psi_file, corr_threshold, metadata)  # 5625 non-unique events with cor =0.8
list_of_events = list(set(list_of_events))  # 5553 unique events with cor = 0.8

# Subset the PSI matrix to only the non-redundant events
formatted_psi_file_1 = formatted_psi_file.loc[list_of_events, :]  # 5553 by 366 matrix now

# Apply different feature selection approaches
# Output for feature selection methods will be lists of events and any text files will be added here rather than inside the function
# Load PCA-based feature selection events
pca_events_path = "/Users/tha8tf/Documents/testsOncosplice/PCA_NEW_30PCS_Events.txt"  # here I am just importing the events that was pre-determined
pca_events = pd.read_table(pca_events_path, names=["events"])
pca_events = pca_events["events"].values.tolist()  # 1086 events

imputed_psi_file_path = "/Users/tha8tf/Documents/testsOncosplice/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-75p_median_imputed.txt"
formatted_psi_file_imp, metadata_imp = format_psi_matrix(imputed_psi_file_path)
formatted_psi_file_imp_1 = formatted_psi_file_imp.loc[list_of_events, :]

# Round 1 OncoSplice
final_clusters_1, depleted_psi_file_after_round_imp_1, depleted_psi_file_after_round_1 = round_wrapper(full_psi_file=formatted_psi_file_1, full_imputed_psi_file=formatted_psi_file_imp_1, highly_variable_events=pca_events, rank=2, strictness="tough")

'''
Running NMF analysis with rank set to 2...
--- 0.854280948638916 seconds ---
Finding top 150 differential splicing events per NMF cluster...
Number of unique differential events found for this round is 240
Determined 240 unique differential events across all clusters
Running Linear SVM on final 2 clusters with...
Depleting events using a correlation threshold of 0.3...
--- 20.7911319732666 seconds ---
'''

# Round 2 OncoSplice
final_clusters_2, depleted_psi_file_after_round_imp_2, depleted_psi_file_after_round_2 = round_wrapper(full_psi_file=formatted_psi_file_1, full_imputed_psi_file=formatted_psi_file_imp_1, highly_variable_events=depleted_psi_file_after_round_1.index.to_list(), strictness="tough")

'''
Running NMF analysis with rank set to 17...
--- 30.588845014572144 seconds ---
Finding top 150 differential splicing events per NMF cluster...
Less than 100 significant diff events for Cluster 9
Less than 100 significant diff events for Cluster 13
Number of unique differential events found for this round is 1598
Determined 1598 unique differential events across all clusters
Running Linear SVM on final 15 clusters with...
Depleting events using a correlation threshold of 0.3...
--- 181.09785795211792 seconds ---

depleted_psi_file_after_round2 has dimensions 0, 366
final_clusters2 has dimensions 366, 15
'''
#################################################### TESTING WITH STRCTNESS = EASY ####################################################

# Round 1 OncoSplice
final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(full_psi_file=formatted_psi_file_1, full_imputed_psi_file=formatted_psi_file_imp_1, highly_variable_events=pca_events, rank=2, strictness="easy")

'''
Running NMF analysis with rank set to 2...
--- 0.7707772254943848 seconds ---
Finding top 150 differential splicing events per NMF cluster...
Number of unique differential events found for this round is 240
Determined 240 unique differential events across all clusters
Running Linear SVM on final 2 clusters with...
Depleting events using a correlation threshold of 0.3...
--- 21.12402081489563 seconds ---
'''

# Round 2 OncoSplice
final_clusters_4, depleted_psi_file_after_round_imp_4, depleted_psi_file_after_round_4 = round_wrapper(full_psi_file=formatted_psi_file_1, full_imputed_psi_file=formatted_psi_file_imp_1, highly_variable_events=depleted_psi_file_after_round_3.index.to_list(), strictness="easy")

'''
Running NMF analysis with rank set to 18...
--- 38.35409998893738 seconds ---
Finding top 150 differential splicing events per NMF cluster...
Less than 100 significant diff events for Cluster 7
Less than 100 significant diff events for Cluster 14
Number of unique differential events found for this round is 1683
Determined 1683 unique differential events across all clusters
Running Linear SVM on final 16 clusters with...
Depleting events using a correlation threshold of 0.3...
--- 193.77845096588135 seconds ---

'''

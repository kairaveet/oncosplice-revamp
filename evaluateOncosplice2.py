from RoundWrapper import *
from removeRedundantSplicingEvents import *
from PCAbasedFeatureSelection import *
from sklearn.metrics import f1_score
import os

# Read PSI file
psi_file_path = "/Users/tha8tf/Documents/testsOncosplice/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-75p.txt"
formatted_psi_file, metadata = format_psi_matrix(psi_file_path)  # 7456 events

# Find the non-redundant splicing events
corr_threshold = 0.8  # change this as needed ### this step may not have been done in Meenakshi's version
list_of_events = remove_redundant_splicing_events(formatted_psi_file, corr_threshold, metadata)  # 5625 non-unique events with cor =0.8
list_of_events = list(set(list_of_events))  # 5553 unique events with cor = 0.8

# Subset the PSI matrix to only the non-redundant events
formatted_psi_file_1 = formatted_psi_file.loc[list_of_events, :]  # 5553 by 366 matrix now

# Same steps/events for the imputed psi matrix
imputed_psi_file_path = "/Users/tha8tf/Documents/testsOncosplice/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-75p_median_imputed.txt"
formatted_psi_file_imp, metadata_imp = format_psi_matrix(imputed_psi_file_path)
formatted_psi_file_imp_1 = formatted_psi_file_imp.loc[list_of_events, :] # 5553 by 366 matrix now

# Read in mutation metadata of samples in this dataset
mutation_metadata = pd.read_excel('/Users/tha8tf/Leucegene-positive-controls-F1.xlsx', index_col=0)

# Set working directory to save evaluations
os.chdir("/Users/tha8tf/Documents/testsOncosplice/force_broad_on_evaluations/")

# feature selection prior to Round 1
pca_events_round1 = get_events_from_pca(formatted_psi_file_imp, formatted_psi_file, list_of_events, corr_threshold=0.4, n_components=30)  # 1242 events

# Round 1 OncoSplice (tough)
final_clusters_1, depleted_psi_file_after_round_imp_1, depleted_psi_file_after_round_1 = round_wrapper(filename="Round1_Tough_forcebroad", full_psi_file=formatted_psi_file_1, full_imputed_psi_file=formatted_psi_file_imp_1, highly_variable_events=pca_events_round1, rank=2, strictness="tough")

'''
Running NMF analysis with rank set to 2...
--- 0.9148058891296387 seconds ---
0
1
Finding top 150 differential splicing events per NMF cluster...
Number of unique differential events found for this round is 236
Determined 236 unique differential events across all clusters
Running Linear SVM on final 2 clusters with...
Depleting events using a correlation threshold of 0.3...
--- 21.24185800552368 seconds ---


4250 events in the depleted psi file after round 1
'''

round1_clusters_tough = final_clusters_1

# Round 2 OncoSplice (tough)
depleted_events_round1 = depleted_psi_file_after_round_1.index.to_list()  # 4254 events
metadata_depleted_events_round1 = metadata.loc[depleted_events_round1, :]

# feature selection prior to Round 2
pca_events_round2 = get_events_from_pca(formatted_psi_file_imp, formatted_psi_file, depleted_events_round1, corr_threshold=0.4, n_components=30)  # 665 events

# Round 2 Oncosplice (tough)
final_clusters_2, depleted_psi_file_after_round_imp_2, depleted_psi_file_after_round_2 = round_wrapper(filename="Round2_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_1, full_imputed_psi_file=depleted_psi_file_after_round_imp_1, highly_variable_events=pca_events_round2, strictness="tough")

'''
Running NMF analysis with rank set to 18...
--- 31.068728923797607 seconds ---
0
Number of sample groups is < 2; Mann-Whitney test not conducted
differential expression not compatible for cluster 0
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
Finding top 150 differential splicing events per NMF cluster...
Less than 100 significant diff events for Cluster 3
Less than 100 significant diff events for Cluster 5
Less than 100 significant diff events for Cluster 9
Less than 100 significant diff events for Cluster 11
Less than 100 significant diff events for Cluster 12
Less than 100 significant diff events for Cluster 14
Number of unique differential events found for this round is 1306
Determined 1306 unique differential events across all clusters
Running Linear SVM on final 11 clusters with...
Depleting events using a correlation threshold of 0.3...
--- 144.0546190738678 seconds ---

2455 events in depleted psi file after round 2
'''

round2_clusters_tough = final_clusters_2

# Round 3 OncoSplice (tough)
depleted_events_round2 = depleted_psi_file_after_round_2.index.to_list()  # 2477 events
metadata_depleted_events_round2 = metadata.loc[depleted_events_round2, :]

# feature selection prior to Round 3
pca_events_round3 = get_events_from_pca(formatted_psi_file_imp, formatted_psi_file, depleted_events_round2, corr_threshold=0.4, n_components=30)  # 56 events

# Round 3 Oncosplice (tough)
final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(filename="Round3_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_2, full_imputed_psi_file=depleted_psi_file_after_round_imp_2, highly_variable_events=pca_events_round3, strictness="tough")

final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(filename="Round3_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_2, full_imputed_psi_file=depleted_psi_file_after_round_imp_2, highly_variable_events=pca_events_round3, strictness="tough", min_differential_events=50)

final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(filename="Round3_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_2, full_imputed_psi_file=depleted_psi_file_after_round_imp_2, highly_variable_events=pca_events_round3, strictness="tough", min_differential_events=20)

'''
Running NMF analysis with rank set to 15...
--- 5.547612190246582 seconds ---
0
Number of sample groups is < 2; Mann-Whitney test not conducted
differential expression not compatible for cluster 0
1
2
3
Mann-Whitney test not conducted: A group has fewer samples than minimum group size:  3  <  5
differential expression not compatible for cluster 3
4
5
6
7
8
9
10
Mann-Whitney test not conducted: A group has fewer samples than minimum group size:  4  <  5
differential expression not compatible for cluster 10
11
12
Mann-Whitney test not conducted: A group has fewer samples than minimum group size:  4  <  5
differential expression not compatible for cluster 12
13
14
Finding top 150 differential splicing events per NMF cluster...
Less than 20 significant diff events for Cluster 1
Less than 20 significant diff events for Cluster 2
Less than 20 significant diff events for Cluster 4
Less than 20 significant diff events for Cluster 5
Less than 20 significant diff events for Cluster 6
Less than 20 significant diff events for Cluster 7
Less than 20 significant diff events for Cluster 8
Less than 20 significant diff events for Cluster 9
Less than 20 significant diff events for Cluster 13
Less than 20 significant diff events for Cluster 14
Number of unique differential events found for this round is 150
Determined 150 unique differential events across all clusters
Running Linear SVM on final 1 clusters with...
Traceback (most recent call last):
  File "/opt/miniconda3/envs/Oncosplice-FeatureSelection/lib/python3.9/site-packages/IPython/core/interactiveshell.py", line 3505, in run_code
    exec(code_obj, self.user_global_ns, self.user_ns)
  File "<ipython-input-9-0c8a05459997>", line 1, in <module>
    final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(filename="Round3_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_2, full_imputed_psi_file=depleted_psi_file_after_round_imp_2, highly_variable_events=pca_events_round3, strictness="tough", min_differential_events=20)
  File "/Users/tha8tf/PycharmProjects/Oncosplice-FeatureSelection/RoundWrapper.py", line 78, in round_wrapper
    final_clusters = classify(train=centroids, imputed_psi_file_with_diff_events=psi_file_with_diff_events_imp, groups=np.arange(np.shape(refined_binarized_output_f)[0]), conservation=conservation)  # need to figure out a way to give more meaningful cluster name instead of an integer after this step (such as C1, etc etc).
  File "/Users/tha8tf/PycharmProjects/Oncosplice-FeatureSelection/linearSVM.py", line 28, in classify
    regression.fit(train.transpose(), groups)  # train should be # of ranks by # of differential splicing events ("samples" as per python documentation are the ranks/clusters) and "features" are splicing events and groups is # of ranks/clusters by # of patients/samples
  File "/opt/miniconda3/envs/Oncosplice-FeatureSelection/lib/python3.9/site-packages/sklearn/svm/_classes.py", line 257, in fit
    self.coef_, self.intercept_, n_iter_ = _fit_liblinear(
  File "/opt/miniconda3/envs/Oncosplice-FeatureSelection/lib/python3.9/site-packages/sklearn/svm/_base.py", line 1162, in _fit_liblinear
    raise ValueError(
ValueError: This solver needs samples of at least 2 classes in the data, but the data contains only one class: 0
'''

depleted_psi_file_after_round_1.to_csv("depleted_psi_file_after_round_1_tough.txt", sep='\t')
depleted_psi_file_after_round_2.to_csv("depleted_psi_file_after_round_2_tough.txt", sep='\t')
pd.DataFrame(pca_events_round1).to_csv("pca_events_round_1_tough.txt", sep='\t')
pd.DataFrame(pca_events_round2).to_csv("pca_events_round_2_tough.txt", sep='\t')
pd.DataFrame(pca_events_round3).to_csv("pca_events_round_3_tough.txt", sep='\t')

######### do not use "easy" for strictness -- we always want the splicing events that are not correlated with ANY of the clusters we have. This means that corr < 0.3 across all clusters, meaning in correlationdepletion script, it should always be np.all(mask,axis=0) or strictness = tough #########

# Round 1 OncoSplice (tough) dep corr 0.4
final_clusters_1, depleted_psi_file_after_round_imp_1, depleted_psi_file_after_round_1 = round_wrapper(filename="Round1_Tough_forcebroad", full_psi_file=formatted_psi_file_1, full_imputed_psi_file=formatted_psi_file_imp_1, highly_variable_events=pca_events_round1, rank=2, strictness="tough", depletion_corr_threshold=0.4)

'''
Running NMF analysis with rank set to 2...
--- 0.9148058891296387 seconds ---
0
1
Finding top 150 differential splicing events per NMF cluster...
Number of unique differential events found for this round is 236
Determined 236 unique differential events across all clusters
Running Linear SVM on final 2 clusters with...
Depleting events using a correlation threshold of 0.3...
--- 21.24185800552368 seconds ---


4965 events in the depleted psi file after round 1
'''

round1_clusters_dep_04 = final_clusters_1

# Round 2 OncoSplice (tough) , depletion_corr_threshold=0.4
depleted_events_round1 = depleted_psi_file_after_round_1.index.to_list()  # 4965 events
metadata_depleted_events_round1 = metadata.loc[depleted_events_round1, :]

# feature selection prior to Round 2
pca_events_round2 = get_events_from_pca(formatted_psi_file_imp, formatted_psi_file, depleted_events_round1, corr_threshold=0.4, n_components=30)  # 863 events

# Round 2 Oncosplice (tough)
final_clusters_2, depleted_psi_file_after_round_imp_2, depleted_psi_file_after_round_2 = round_wrapper(filename="Round2_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_1, full_imputed_psi_file=depleted_psi_file_after_round_imp_1, highly_variable_events=pca_events_round2, strictness="tough", depletion_corr_threshold=0.4)

'''
Running NMF analysis with rank set to 19...
--- 36.6733820438385 seconds ---
0
Number of sample groups is < 2; Mann-Whitney test not conducted
differential expression not compatible for cluster 0
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
Finding top 150 differential splicing events per NMF cluster...
Less than 100 significant diff events for Cluster 11
Less than 100 significant diff events for Cluster 15
Number of unique differential events found for this round is 1523
Determined 1523 unique differential events across all clusters
Running Linear SVM on final 14 clusters with...
Depleting events using a correlation threshold of 0.4...
--- 177.86131978034973 seconds ---

4025 events in depleted psi file after round 2
'''

round2_clusters_dep_04 = final_clusters_2

# Round 3 OncoSplice (tough)
depleted_events_round2 = depleted_psi_file_after_round_2.index.to_list()  # 2477 events
metadata_depleted_events_round2 = metadata.loc[depleted_events_round2, :]

# feature selection prior to Round 3
pca_events_round3 = get_events_from_pca(formatted_psi_file_imp, formatted_psi_file, depleted_events_round2, corr_threshold=0.4, n_components=30)  # 309 events

# Round 3 Oncosplice (tough)
final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(filename="Round3_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_2, full_imputed_psi_file=depleted_psi_file_after_round_imp_2, highly_variable_events=pca_events_round3, strictness="tough")

final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(filename="Round3_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_2, full_imputed_psi_file=depleted_psi_file_after_round_imp_2, highly_variable_events=pca_events_round3, strictness="tough", min_differential_events=50)

final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(filename="Round3_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_2, full_imputed_psi_file=depleted_psi_file_after_round_imp_2, highly_variable_events=pca_events_round3, strictness="tough", min_differential_events=20)

'''
Running NMF analysis with rank set to 19...
--- 17.714022874832153 seconds ---
0
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
Finding top 150 differential splicing events per NMF cluster...
Less than 100 significant diff events for Cluster 2
Less than 100 significant diff events for Cluster 3
Less than 100 significant diff events for Cluster 7
Less than 100 significant diff events for Cluster 10
Less than 100 significant diff events for Cluster 11
Less than 100 significant diff events for Cluster 13
Less than 100 significant diff events for Cluster 14
Less than 100 significant diff events for Cluster 16
Number of unique differential events found for this round is 1141
Determined 1141 unique differential events across all clusters
Running Linear SVM on final 11 clusters with...
Depleting events using a correlation threshold of 0.4...
--- 148.7050461769104 seconds ---
'''

round3_clusters_dep_04 = final_clusters_3

depleted_psi_file_after_round_1.to_csv("depleted_psi_file_after_round_1_tough_dep04.txt", sep='\t')
depleted_psi_file_after_round_2.to_csv("depleted_psi_file_after_round_2_tough_dep04.txt", sep='\t')
pd.DataFrame(pca_events_round1).to_csv("pca_events_round_1_tough_dep04.txt", sep='\t')
pd.DataFrame(pca_events_round2).to_csv("pca_events_round_2_tough_dep04.txt", sep='\t')
pd.DataFrame(pca_events_round3).to_csv("pca_events_round_3_tough_dep04.txt", sep='\t')

# Round 4 OncoSplice (tough)
depleted_events_round3 = depleted_psi_file_after_round_3.index.to_list()  # 2477 events
metadata_depleted_events_round3 = metadata.loc[depleted_events_round3, :]

# feature selection prior to Round 3
pca_events_round4 = get_events_from_pca(formatted_psi_file_imp, formatted_psi_file, depleted_events_round3, corr_threshold=0.4, n_components=30)  # 309 events

# Round 3 Oncosplice (tough)
final_clusters_4, depleted_psi_file_after_round_imp_4, depleted_psi_file_after_round_4 = round_wrapper(filename="Round4_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_3, full_imputed_psi_file=depleted_psi_file_after_round_imp_3, highly_variable_events=pca_events_round4, strictness="tough")

final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(filename="Round3_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_2, full_imputed_psi_file=depleted_psi_file_after_round_imp_2, highly_variable_events=pca_events_round3, strictness="tough", min_differential_events=50)

final_clusters_3, depleted_psi_file_after_round_imp_3, depleted_psi_file_after_round_3 = round_wrapper(filename="Round3_Tough_forcebroad", full_psi_file=depleted_psi_file_after_round_2, full_imputed_psi_file=depleted_psi_file_after_round_imp_2, highly_variable_events=pca_events_round3, strictness="tough", min_differential_events=20)

round4_clusters_dep_04 = final_clusters_4

######### do not use "easy" for strictness -- we always want the splicing events that are not correlated with ANY of the clusters we have. This means that corr < 0.3 across all clusters, meaning in correlationdepletion script, it should always be np.all(mask,axis=0) or strictness = tough #########

# Perform evaluations using f1 score for dep 04

for k in range(9):
    for i in range(np.shape(round1_clusters_dep_04)[1]):
        y_true = mutation_metadata.iloc[:, k]
        print(mutation_metadata.columns.values[k])
        y_pred = round1_clusters_dep_04.iloc[:, i]
        print(f1_score(y_true, y_pred))

for k in range(9):
    for i in range(np.shape(round2_clusters_dep_04)[1]):
        y_true = mutation_metadata.iloc[:, k]
        print(mutation_metadata.columns.values[k])
        y_pred = round2_clusters_dep_04.iloc[:, i]
        print(f1_score(y_true, y_pred))

for k in range(9):
    for i in range(np.shape(round3_clusters_dep_04)[1]):
        y_true = mutation_metadata.iloc[:, k]
        print(mutation_metadata.columns.values[k])
        y_pred = round3_clusters_dep_04.iloc[:, i]
        print(f1_score(y_true, y_pred))

for k in range(9):
    for i in range(np.shape(round4_clusters_dep_04)[1]):
        y_true = mutation_metadata.iloc[:, k]
        print(mutation_metadata.columns.values[k])
        y_pred = round4_clusters_dep_04.iloc[:, i]
        print(f1_score(y_true, y_pred))

round4_clusters_dep_04_min50 = final_clusters_4
for k in range(9):
    for i in range(np.shape(round4_clusters_dep_04_min50)[1]):
        y_true = mutation_metadata.iloc[:, k]
        print(mutation_metadata.columns.values[k])
        y_pred = round4_clusters_dep_04_min50.iloc[:, i]
        print(f1_score(y_true, y_pred))


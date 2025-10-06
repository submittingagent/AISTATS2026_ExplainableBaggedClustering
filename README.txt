Summary of tables and figures:

Table 1: synthetic dataset description, text
Table 2: clustering performances Synthetic and Iris
Table 3: Synthetic MI results with different m number of features
Table 4: Synthetic comparison with literature
Table 5: Iris MI results with different m number of features
Table 6: Iris comparison with literature 

Figure 1: Synthetic dataset kde
Figure 2: Iris dataset kde

------------------------------------------------------------------------------

Scripts:

-functions.R
-datasets.R
-kde.R
-packages_literaturecomparison.R
-tables_results_proposal.R

------------------------------------------------------------------------------

How to run the code:

1. run datasets.R to generate the datasets
2. -run kde.R for Fig.1,2
   -run tables_results_proposal.R for Table 2-6.
    Select the dataset type in the script. Select the desired m (for feature dropout) in the script.

Outputs: 
-iris.csv, iris dataset
-synthetic*.csv, 20 instances of the synthetic dataset
-table_clustering_results.txt, Table 2
-table_MI.txt, results for unweighted mutual information; for different m, the results of Table 3, 5 (synthetic, iris)
-table_WMI.txt, results for weighted mutual information; for different m, part of the results of Table 4, 6(synthetic, iris)
-table_comparison.txt, results for the literature methods, Table 4, 6 (synthetic, iris)


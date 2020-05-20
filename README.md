# LGEA_Cell_Signature

R implementation of our statistical method to predict the LGEA single cell signature genes.

Developed by Shuyang Zhao.

Our method is comprised of a logistic regression model and a ranking system. The logistic regression model uses elastic net regularization to predict each gene's probability being a cell type signature, based on the integration of multiple signature metrics including cell specific p-values of differential expression tests, gene expression effect size, frequency and sensitivity. Then the ranking system was used to define the signature genes for each cell type.

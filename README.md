# MIMIC
MIMIC: Mixed Integer programming Model to Identify Cell type specific marker panel
  Version 0.0.1
  Last updated: Sep. 18, 2020

# Motivation
Cell sorting aims to separate a heterogeneous mixture of cells and further research into the hierarchical topology of these cells according to the intracellular process (DNA, RNA, and protein interaction) or extracellular properties (morphology and surface marker expression). The ENCODE produces large amounts of such data and provides the opportunity and challenge in cell sorting. However, the challenge always goes with the high dimension. Here we aim to propose a novel model to identify the cell type specific marker panel and finally assist in cell sorting.
# Method
We developed a novel mixed integer programming model for identify the cell type specific marker panel. This model directly selected the cell type specific markers with simultaneously maintaining the hierarchical topology among different cell types given the number of markers. This mixed integer programming model allows us to go through all the optimal combinations by varying parameter from 1 to n (the number of markers). Moreover, we can check their accuracy and compare the selected combinations. In particular, an optimal cell type specific marker panel can be selected by balancing the number of selected markers and the classification accuracy.

# Software
CPLEX Optimizer is required in your MATLAB. CPLEX is available on the website: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/.
This version of the program is in very preliminary stage and provided just for testing purpose. The program is still under development.
Dataset
The dataset used in our paper serves as an example to demonstrate the implementation for MIMIC. The raw data(.xlsx) and processed data(.mat) are available as follows.
# References
•	Meng Zou, Zhana Duren, Qiuyue Yuan, Henry Li, Henry Li, Andrew Hutchins, Wing Hung Wong, Yong Wang, An optimization method to identify cell type specific marker panel for cell sorting. In submission.

This is the computer programs for 
"PDC-MAKES: A Conditional Screening Method for Controlling False Discoveries in High-dimensional Multi-response Setting” 
by Wei Xiong, Han Pan, and Tong Shen.

You first need to install R packages "energy","MFSIS","doParallel","mvtnorm" , "knockoff" and "Ball".  Then The "PCSIS.py" Python script should be placed in the appropriate directory within the MFSIS R package folder. Specifically, you need to place the file in the python subfolder under the main MFSIS package directory.

./codes/

The file 'functions.R' includes the primary functions employed in our study, encompassing those for the PDC-Screen, rPDC-Screen, PDC-MAKES, and rPDC-MAKES methods, along with the modified PC-Screen and PC-Knockoff methods, which are designed to handle both univariate response and multi-response settings.

The file 'Example1.R' contains the code that executes Example 1 of the simulation studies presented in the main text.

The file 'Example2.R' contains the code that executes Example 2 of the simulation studies presented in the main text.

The file 'Example3.R' contains the code that executes Example 3 of the simulation studies presented in the main text.

The file 'Example4.R' contains the code that executes Example 4 of the simulation studies presented in the main text.

./data/

The files "dna," "genes," "10 X genes," and "Scheetz2006.rds" contain the datasets referenced in Section 5.3 of the main text, as well as in Sections G.2 and G.3 of the supplementary materials, respectively. The "dna" and "genes" datasets represent the response variable Y and the covariates X, respectively, in the dataset described in Section 5.3.
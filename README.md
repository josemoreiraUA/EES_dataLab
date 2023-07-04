# EES_dataLab
Run the project as a developer (still no User Interface available):

-> Atribute the WKT files path to the variable wktFiles in vcp_fire.py

-> Also in vcp_fire.py define the name of the file where you want the solution to be saved (variable fileName)

-> Run vcp_fire.py

-> When the program finishes running, open the file where you choose to save the solutions

-> Go to the end of the file and copy the correspondences array
    -- The solution is a tridimensional array of correspondences, where the first dimension is the correspondence between a source and target point, the second dimension is the array of solutions between a source and target polygon and the third are all the correspondences for the dataset

-> Paste the array to the variable correspondences in vcp_test.py

-> Run vcp_test.py and the correspondences and polygon interpolation will be rendered

File Organization:
-> The solutions for the 195 polygons are in the folders myImplementation/results/tigas195_Jaccard_FullCorr and myImplementation/results/tigas195_JaccardDivDist_FullCorr
-> The solutions for the 11 polygons are in the folders myImplementation/results/tigas13_Jaccard and myImplementation/results/tigas13_JaccardAndDistDiv
-> All the wkt files (datasets) are in the folder myImplementation/dataset/tigas13 (for the 11 polygons) and myImplementation/dataset/tigas226 (for the 195 polygons)
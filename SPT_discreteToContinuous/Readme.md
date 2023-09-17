# Data and Resources Paper. Creating Moving Region Representations: a Case Study on the Spread of a Forest Fire

This folder is part of the EES DataLab Project repository and contains the support materials for the article submitted to ACM SIGSpatial 2023.

## Abstract
Our research focuses on generating spatiotemporal data from real-world observations to represent the evolution of phenomena of interest as moving regions.
The case study is the creation of a dataset to represent the spread of a controlled forest fire from aerial images captured using a drone.
We present an overview of the data acquisition and preparation steps, and describe the optimization strategy implemented to establish a vertex correspondence between the regions that delimit the burned area at discrete time instants. These data can be used to create a continuous representation of the evolution of the burned area over time.

## Organization

1. *VideoX_PolYYY_ZZZ*. Four videos showing the results obtained using the small (11 snapshots) and large (195 snapshots) datasets using the JaccardIndex (JI) and the Combined Jaccard Index (CJI) as measures to evaluate vertex correspondences between consecutive pairs of regions (snapshots).

2. Folder *Code*. Holds the programs implemented in Python to create correspondences between regions and visualize the results.

3. Folder *Data*.

* Folder *Source*.
    a) Link *SourceVideo*. To watch the video from which the snapshots were taken.

    b) Folder *tigas13_WKT__pol11* (small dataset). Contains 13 files with the coordinates of the regions segmented from the selected snapshots in the video. Data are represented as polygons in WKT format. Snapshots 12 and 13 were not used because they have polygons with holes.

    c) Folder *tigas226_WKT__pol195* (large dataset). Similar to the previous one. Snapshot 195 in the large dataset corresponds to snapshot 11 in the small dataset. Therefore, both datasets represent the same portion of the original video.

* Folder *Results*.
    a) Images folder. Contains the figures prepared for the paper.

    b) Results folder. Contains the data generated during the optimization of the polygon matching, including the local minima. The vertex correspondences are at the end of each file and are represented as a list in Python. The list has 3 levels: the innermost level represents the matches between a pair of vertices. The intermediate level represents the correspondences between a pair of polygons, a source and a target. The outermost level represents correspondences between paired polygons.

* *ODM_OrthophotoMap_XXXXX__1GB* are links to high-resolution images captured before and after the burn. The image before can be used, for example, to make a characterization of the vegetation.

We have included links to materials available elsewhere to save storage space.

## Usage
1.  Run vcp_fire.py to create the vertex correspondences.
    * The variable *wktFiles* holds the path to the WKT files to be processed.
    * The variable *fileName* defines the file where the results are stored. 

2. Run vcp_test.py to view the interpolation.
    * Open the file where the results were saved, go to the end and copy the list of correspondences.
    * Paste the list of correspondences as the value of the variable *correspondences*.
    * Run the program.

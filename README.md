# Biga2020 is a repository for code used in the paper: "A dynamic, spatially periodic, micro-pattern of HES5 underlies neurogenesis in the mouse spinal cord"

It is organised in multiple sections containing custom routines for:

1) Model- A multicellular and stochastic model of HES5 containing self-repression as well as an inter-cellular Hill function with time delay to represent how protein in one cell can affect transcription in a neighbouring cell via Notch-Delta; cells are arrayed in a 2D grid.

2) Data_Analysis - Analysis of spatial periodicity from kymograph data using the auto-correlation and Lomb-Scargle methods with statistical significance testing; also contains quantification of persistence of gene expression.

3) Microcluster_Quantification - Image processing routines for detection and automated counting of cells in a microcluster based on nuclear-segmented images.

The routines in each folder can be downloaded and run independently of one another.

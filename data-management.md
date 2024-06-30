During these two weeks, we'll be collecting many sessions of data from various mice. All of that data needs to be preprocessed before it can be used for visualization and analysis, and keeping track of which session has had what kind of processing done to it quickly becomes messy!

To keep things organized, data on the server (e.g. `Z:\NSB_2024\03_Mouse\Neuropixels-data\`) should live in one of two specific locations (**incoming** and **preprocessed**) using a standardized folder naming scheme and structure. 

Freshly acquired data goes in the **incoming** folder. It's a "reference copy" of the raw data that should not be touched. 

Each mouse has its own top-level folder within **incoming**, such as `M415`. Within each mouse's folder are folders of individual recording sessions, that should be named `Mxxx-2024-07-yy_description`; for instance `M415-2024-07-05_gapcrossing1`.

To preprocess data, make a copy of the incoming folder you want to use to local storage on your computer. Do not try to write or modify the **incoming** data. After preprocessing is done, you copy the preprocessed data to the **preprocessed** folder, again using the same folder structure and naming scheme.

# LowRankDownscalingSMAP
Code for paper "Low Rank Gap-Filling and Downcaling of SMAP Soil Moisture" by Beale, Bras, and Romberg (Ecohydrology 2025). Includes subset of data and learned parameters needed to produce example results.

* `plot_results.m` -- plots superresolution examples (Figures 3 and 4 from the paper), along with example gap-filling results
* `superresolve_LR_image_blockwise.m` -- separate function for performing superresolution
* `gap_fill_with_eofs.m` -- separate function for performing gap-filling

# LowRankDownscalingSMAP
Code for paper "Low Rank Gap-Filling and Downcaling of SMAP Soil Moisture" by Beale, Bras, and Romberg (Ecohydrology 2025). Includes subset of data and learned parameters needed to reproduce certain results.

* `plot_results.m` -- plots superresolution examples (Figures 3 and 4 from the paper), along with example gap-filling results
* `superresolve_LR_image_blockwise.m` -- separate function for performing superresolution
* `gap_fill_with_eofs.m` -- separate function for performing gap-filling
* `plot_global_error.m` -- produces global error plots for specified error type (Figure 6 from paper) and global error histograms (Figure 7 from paper)
* `get_errors_for_36km_ease_cell.m` -- function to return expected error values of the superresolution result at a specified EASE grid cell

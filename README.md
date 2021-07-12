# FragmentR v1.0

FragmentR is an R tool for DNA fragment analysis, developed and maintained by the Walk Lab at MSU Bozeman.

Instructions, documentation, and tutorials can be found at:

* https://github.com/nvpinkham/FragmentR

FragmentR has been successfully installed and run on Mac OS X and Windows



WORK FLOW: 
1. Summerize fsa (Fragment Analysis Data) files: FSA files contain multiple channels with chromatograms generated at a specific wavelength/dye during the run. Peaks are called in for the query channel and the channel with the DNA ladder. The size of DNA fragments in the query channel is interpolated using a linear model generated from the ladder channel. FSA Files are summarized in a summary matrix containing peak heights and the fragement size (basepairs).
2. Find match: The summary matrix of the query file is compaired to each summary matrix in the database. This is done by by calculating the Bray-Curtis distance between each entry in the database and the query. Normailization and peak binning is done for each comparision. When comparing two peak summaries if two peaks are within 2.5 base pairs of eachother they are combinned. Peak heights are normalized to largest peak in the peak summary matrix. 
3.  Interpretation and visual insepection: Matching is done by measureing the bray curtis disance between peak heights. When comparing two peak summaries if two peaks are within 2.5 base pairs of eachother they are combinned. chromatograms are classified as either a good match (<0.10 distance), questionable match (0.10 - 0.20 disttance), or poor match (>.20 distance). If a match is questionable conisder visually inspecting the chromatograms. We have found a log transformation can increase the influence of smaller peaks in the chromotogram and help shed some light on these questionable matches. 




Details:
The default query channel 1 and the default ladder channel is 4
The default ladder is a rox 1000 ladder. 

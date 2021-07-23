# FragmentR v1.0

FragmentR is an R tool for DNA fragment analysis, developed and maintained by the Walk Lab at MSU Bozeman.

Instructions, documentation, and tutorials can be found at:

* https://github.com/nvpinkham/FragmentR

FragmentR has been successfully installed and run on Mac OS X and Windows



WORK FLOW: 
1.  Summarize fsa (Fragment Analysis Data) files: FSA files contain multiple channels with chromatograms generated at a specific wavelength/dye during the run. Peaks are called in for the query channel and the channel with the DNA ladder. The size of DNA fragments in the query channel is interpolated using a linear model generated from the ladder channel. FSA Files are summarized in a summary matrix containing peak heights and the fragment size (base pairs).![image](https://user-images.githubusercontent.com/47755049/126802400-8b75d3a2-d269-4da5-9515-a8d1eab3d06a.png)
2.	Find match: The summary matrix of the query file is compared to each summary matrix in the database. This is done by calculating the Bray-Curtis distance between each entry in the database and the query. Normalization and peak binning are done for each comparison. When comparing two peak summaries if two peaks are within 2.5 base pairs of each other they are combined. Peak heights are normalized to largest peak in the peak summary matrix.
3.	Interpretation and visual inspection: Matching is done by measuring the Bray-Curtis distance between peak heights. When comparing two peak summaries if two peaks are within 2.5 base pairs of each other they are combined. chromatograms are classified as either a good match (<0.10 distance), questionable match (0.10 - 0.20 distance), or poor match (>.20 distance). If a match is questionable consider visually inspecting the chromatograms. We have found a log transformation can increase the influence of smaller peaks in the chromatogram and help shed some light on these questionable matches.

SET UP:

place "Cdiff_DB_list.5.4.rds" and "F-RibotypingFiles" in the same dirrectoy as "Call_FSA.R"
place all query FSA files in "Files_to_analyzed" folder.
run "Call_FSA.R" script from the commandline. 




Details:
The default query channel 1 and the default ladder channel is 4
The default ladder is a rox 1000 ladder. 

if you are on unix system "Rscript Call_FSA.R"

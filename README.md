# Experiment data analyses
Project repository the combined DIVERCE monoculture experiment analyses.

# Aims, goals, etc.
Explore trait-dependent growth and density-dependent trait change across the two model systems and evaluate the frequency of feedbacks between these.

# Repository structure

Most important files are shown here:

```
.
├── data               // where we keep the raw data
│   ├── ciliates
│       └── DIVERCE_TdB_Ciliates_Traits_Cut.csv
│   └── cyanobacteria
│       └── mono_data.csv
├── figures            // where we save the figs
│   ├── fig_name.pdf
│   ├── ...
│   └── fig_name.pdf
├── Tools and data.R   // funcs used for analysis
├── Saving plots.R     // funcs used for plots
├── Analysis.R         // source to make plots
├── main.qmd           // the actual article
├── supplementary.qmd  // the actual article
└── README.md
```

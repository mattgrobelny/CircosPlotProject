# Circos Plot Project

Goal: Generate a circos plot from a three-spined stickleback fish FST and the RNAseq gene expression fold change dataset.

## Circos Plot Construction
![](https://github.com/mattgrobelny/Data-viz-Circle-plot/blob/master/output_plots/plots_contruction_large.gif "Circos Plot Contstruction")  

## What is a Circos Plot?
A visualization of comparative statistics values for each chromosome for a given organism.  

### What data is Plotted?  

### How do I interpret this plot?
**Overview**  
- Each circle of arcs represent the whole genome of a three-spined stickleback split up into individual chromosomes.
- Each chromosome arc is colored according to its corresponding statistics for a specific nucleotide (Fst) or a region average (Div). Each arc has its genomic coordinates plotted with 5mb breaks (dashed grey lines) with labels at every 10mb intervals.
- The color key for each statistic is in the top right corner ranging from 0 to 1 p-values for Fst and normalized (0 to 1) fold change values for gene expression.

**Fst Values (inner circle)**
- Shows which region of the genome has population wide differentiation, with regions with a low number of single nucleotide polymorphism (SNPs) meaning these regions are fixed and have low variability are colored dark green.
- Regions which contain lots of SNP or variability are white.  

**RNAseq gene expression fold change *DIV* (outer circle)**
- Maps gene expression fold change values for a RNAseq dataset.
- Regions of genes which where highly expressed are colored red, while regions which were relatively less expressed are colored blue.

![](https://github.com/mattgrobelny/CircosPlotProject/blob/master/output_plots/jpg/12Grobelny_data_viz-1.jpg "Final Circos Plot")  

### What is the code work flow?

### What does each function do?

## Which Libraries Where Used to Generate the Plot?
- [Cairo](https://www.cairographics.org) - 2d graphics library
- [pycairo](https://www.cairographics.org/documentation/pycairo/3/) - python connector for the cairo library
- matplotlib - color maps
- numpy - some math functions

## Ideas for Future Improvements
- [] Clean up code.  
- [] Organize functions in separate script.  
- [] Rewrite with each chromosome arc of the circos plot as an object.  
- [] Add a web interface for easy data upload and sider based parameter adjustment.  

![](https://upload.wikimedia.org/wikipedia/commons/thumb/9/9e/Culaea_inconstans_1908.jpg/250px-Culaea_inconstans_1908.jpg "three-spined stickleback (Gasterosteus aculeatus)")

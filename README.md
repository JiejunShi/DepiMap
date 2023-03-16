# DepiMap - ***D***ifferential ***Epi***genomic ***Map***

## DepiMap is used to analyze and visulize the dynamic epigenomic profiles among samples/conditions.

## Authors
- Jiejun Shi (shij@tongji.edu.cn)
## Dependencies
- R with following packages
  - getopt
## Usage
```
$ Rscript DepiMap.r -h
Usage: Rscript DepiMap.r [-[-help|h]] [-[-SampleMatrixes|s] <character>] [-[-SampleNames|n] <character>] [-[-Zmax|z] <character>] [-[-Colors|c] <character>] [-[-SmoothWindow|S] <integer>] [-[-Subtrahend|u] <character>] [-[-Minuend|m] <character>] [-[-DeltaZmax|Z] <character>] [-[-DeltaColors|C] <character>] [-[-IDfiles|i] <character>] [-[-IDorder|I]] [-[-WhichOrder|w] <character>] [-[-OrderInterval|v] <character>] [-[-Decreasing|d] <character>] [-[-DeltaOrder|o] <character>] [-[-LabelIDfile|l] <character>] [-[-ColorResolution|r] <integer>] [-[-UseRaster|R]] [-[-OutTiff|t] <character>] [-[-PanelWidth|W] <integer>] [-[-TotalHeight|H] <integer>] [-[-Xtick|x] <character>]
    -h|--help               useage
    -s|--SampleMatrixes     Sample matrixes seperated by ",". REQUIRED.
    -n|--SampleNames        Names list("," seperated) in the same length as SampleMatrixes. REQUIRED.
    -z|--Zmax               Zmax list("," seperated) in the length of either 1 or the same as SampleMatrixes. REQUIRED.
    -c|--Colors             Color list("," seperated) in the length of either 1 or the same as SampleMatrixes. Default is blue.
    -S|--SmoothWindow       Window size(number of matrix columns) used for Smoothing. Default is 1, i.e. no smoothing.
    -u|--Subtrahend         Sample Names("," seperated) which will be used as Subtrahend in Delta-Heatmaps.
    -m|--Minuend            Sample Name which will be used as Minuend in Delta-Heatmaps, only 1 Minuend sample supported.
    -Z|--DeltaZmax          Zmax of Delta-Heatmaps("," seperated) in the length of either 1 or the same as Subtrahend Samples. REQUIRED if Subtrahend and Minuend are detected.
    -C|--DeltaColors        3 colors("," seperated)corresponding to -DeltaZmax,0,DeltaZmax in Delta-Heatmaps. Default is "black,white,orange".
    -i|--IDfiles            ID list files("," seperated). One ID each row. Only the rows in matrixes whose id(4th column) included are presented in heatmap. Multiple horizontal panels will presented if multiple ID files detected. If absent, all rows are presented.
    -I|--IDorder            Whether sort each horizontal panel by IDfiles.
    -w|--WhichOrder         Sample Names which will be used to order each horizontal panel. The length must be the same as IDfiles count. Default is the first Sample Name.
    -v|--OrderInterval      Two integers indicate the columns which are used to order matrix rows. If absent, all columns are used.
    -d|--Decreasing         List of 1 and/or 0 means whether order Decreasingly or not. The length must be the same as IDfiles count. Default is 1.
    -o|--DeltaOrder         List of 1 and/or 0 means whether order by Delta-Heatmap or Raw heatmap. The length must be the same as IDfiles count. The Delta-Heatmap must be plotted if you want to order by it. Default is 0.
    -l|--LabelIDfile        A file of IDs which will be labeled in heatmap. One ID each row.
    -r|--ColorResolution    Divide 1 into this number of colors. Default is 1.
    -R|--UseRaster          Whether useRaster for the heatmaps.
    -t|--OutTiff            Output tiff file name. Default is "Output_Heatmap.tiff".
    -W|--PanelWidth         Width of each Column Panel. Default is 300.
    -H|--TotalHeight        Total Height of the plot. Default is 1200.
    -x|--Xtick              Xtick interval length indicated by either 2(TSS upstream and downstream) or 3(TSS upstream, genebody and TES downstream) integers seperated by ",". Default is 5000,5000.
```
**Note:** *deeptools-computeMatrix* is recommended to generate the input SampleMatrix.

## Citation
DepiMap was introduced and used in papers below. 

Shi, J. *et al*. Drosophila Brahma complex remodels nucleosome organizations in multiple aspects. [***Nucleic Acids Research*** 42(15):9730-9739 (2014).](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gku717)

Hsu, C.C., Shi, J. *et al*. Recognition of histone acetylation by the GAS41 YEATS domain promotes H2A.Z deposition in non-small cell lung cancer. [***Genes & Development*** 32:58-69 (2018).](http://genesdev.cshlp.org/content/32/1/58.long)

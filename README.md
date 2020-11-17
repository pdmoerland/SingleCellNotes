# SingleCellNotes

## General
* [awesome-single-cell](https://github.com/seandavi/awesome-single-cell): List of software packages (and the people developing these methods) for single-cell data analysis, including RNA-seq, ATAC-seq, etc.
* [scRNA-tools](https://www.scrna-tools.org/): A catalogue of tools for analysing single-cell RNA sequencing data.
* A curated [database](https://docs.google.com/spreadsheets/d/1En7-UV0k0laDiIfjFkdn7dggyR7jIk3WH8QgXaMOZF0/edit#gid=0) of of single-cell transcriptomics studies.
* [Single-cell reading list](https://github.com/gtca/single-cell-reading-list): A curated selection of blog posts and papers on single-cell data analysis. 
* [Single Cell Genomics Day](https://satijalab.org/scgd/): Each year the lab of Rahul Satija organizes a one-day workshop highlighting recent developments in the field.

## Courses
* [Analysis of single cell RNA-seq data](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html): A long-running course from the Cambridge Bioinformatics training unit (Martin Hemberg and others). See also their GitHub [repository](https://github.com/hemberg-lab/scRNA.seq.course). This [course](https://broadinstitute.github.io/2020_scWorkshop/) at the Broad Institute is based on it and offers some interesting extensions (on CITE-Seq for example).
* [MGC/BioSB Course - Single Cell Analysis](https://github.com/LeidenCBC/MGC-BioSB-SingleCellAnalysis2020): This course covers the practicalities of single-cell sample prep and analysis with a particular focus on single-cell RNA-seq libraries. 

## Tutorials
* [Orchestrating Single-Cell Analysis with Bioconductor](http://bioconductor.org/books/release/OSCA/): Very comprehensive on-line book that
teaches you how to make use of cutting-edge Bioconductor tools to process, analyze, visualize, and explore scRNA-seq data. The companion paper can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7358058/).
* [Current best practices in single‐cell RNA‐seq analysis: a tutorial](https://www.embopress.org/doi/full/10.15252/msb.20188746): A highly readable tutorial by Malte Luecken and Fabian Theis. There is a companion GitHub [repository](https://github.com/theislab/single-cell-tutorial) with the scripts.
* [Seurat vignettes](https://satijalab.org/seurat/vignettes.html): [Seurat](https://satijalab.org/seurat/) remains the most comprehensive R toolkit for single genomics analysis accompanied by a large collection of vignettes. 
* [Scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html): [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) is a scalable Python toolkit for analyzing single-cell gene expression data and comes with a rich set of tutorials.
[Kallisto|bustools](https://www.kallistobus.tools/tutorials):  Google Colab notebooks on working with kallisto in combination with bustools.

## Trajectory inference
* [single-cell-pseudotime](https://github.com/agitter/single-cell-pseudotime): Overview of single-cell RNA-seq pseudotime estimation algorithms.
* [dynmethods](https://github.com/dynverse/dynmethods): A collection of 55 trajectory inference methods. To run any of these methods, interpret the results and visualise the trajectory, see the [dyno](https://github.com/dynverse/dyno) package.

## Gene regulatory network inference

### Reviews
* [A comprehensive survey of regulatory network inference methods using single cell RNA sequencing data](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa190/5904505): Nguyen et al. Briefings in Bioinformatics, bbaa190 (2020).
* [Integration of single-cell multi-omics for gene regulatory network inference](https://www.sciencedirect.com/science/article/pii/S2001037020303226): Hu et al. Computational and Structural Biotechnology Journal. 2020;18:1925-1938.
* [Network Inference from Single-Cell Transcriptomic Data](https://link.springer.com/protocol/10.1007%2F978-1-4939-8882-2_10): Todorov et al. (2019) Network Inference from Single-Cell Transcriptomic Data. In: Sanguinetti G., Huynh-Thu V. (eds) Gene Regulatory Networks. Methods in Molecular Biology, vol 1883. Humana Press, New York, NY.
* [Mapping gene regulatory networks from single-cell omics data](https://academic.oup.com/bfg/article/17/4/246/4803107): Fiers et al. Briefings in Functional Genomics, Volume 17, Issue 4, July 2018, Pages 246–254.
* [Learning regulatory models for cell development from single cell transcriptomic data](https://www.sciencedirect.com/science/article/pii/S2452310017301257): Babtie et al. Current Opinion in Systems Biology 5 (2017): 72-81.

### Benchmark studies
* [Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7098173/): Pratapa et al. Nat Methods. 2020;17(2):147-154. 
* [Evaluating methods of inferring gene regulatory networks highlights their lack of performance for single cell gene expression data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2217-z): Chenet al. BMC Bioinformatics 19, 232 (2018). 

### Methods
* [scPADGRN](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007471&rev=2): Zheng et al. (2020) scPADGRN: A preconditioned ADMM approach for reconstructing dynamic gene regulatory network using single-cell RNA sequencing data. PLoS Comput Biol 16(7): e1007471. [Code](https://github.com/xzheng-ac/scPADGRN)
* [GRISLI](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa576/5858974): Aubin-Frankowski et al. Gene regulation inference from single-cell RNA-seq data with linear differential equations and velocity inference, Bioinformatics, btaa576, 2020. [Code](https://github.com/PCAubin/GRISLI)
* [Scribe](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30036-3): Qiu et al. Inferring causal gene regulatory networks from coupled single-cell expression dynamics using Scribe. Cell Syst. 2020;10(3):265-274.e11. [Code](https://github.com/cole-trapnell-lab/Scribe)
* [SINGE](https://www.biorxiv.org/content/10.1101/534834v1): Deshpande et al. (2019) Network inference with Granger causality ensembles on single-cell transcriptomic data. bioRxiv 534834. [Code](https://github.com/gitter-lab/SINGE)
* [WASABI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2798-1): Bonnaffoux et al. WASABI: a dynamic iterative framework for gene regulatory network inference. BMC Bioinformatics 20, 220 (2019).
* [M&NEM](https://academic.oup.com/bioinformatics/article/34/17/i964/5093248): Pirkl et al. Single cell network analysis with a mixture of Nested Effects Models, Bioinformatics, Volume 34, Issue 17, 2018, Pages i964–i971. [Code](https://github.com/cbg-ethz/mnem/)
* [AR1MA1 - VBEM](https://academic.oup.com/bioinformatics/article/34/6/964/4222631): Sanchez-Castillo et al. A Bayesian framework for the inference of gene regulatory networks from time and pseudo-time series data, Bioinformatics, Volume 34, Issue 6, 2018, Pages 964–970. [Code](https://github.com/mscastillo/GRNVBEM)
* [SCENIC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5937676/): Aibar et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods. 2017;14(11):1083-1086. [R code](https://github.com/aertslab/SCENIC), [Python code](https://github.com/aertslab/pySCENIC)
* [SCODE](https://academic.oup.com/bioinformatics/article/33/15/2314/3100331: Matsumoto et al. SCODE: an efficient regulatory network inference algorithm from single-cell RNA-Seq during differentiation, Bioinformatics, Volume 33, Issue 15, 2017, Pages 2314–2321. [Code](https://github.com/hmatsu1226/SCODE)
* [SINCERITIES](https://academic.oup.com/bioinformatics/article/34/2/258/4158033): Gao et al. SINCERITIES: inferring gene regulatory networks from time-stamped single cell transcriptional expression profiles, Bioinformatics, Volume 34, Issue 2, 2018, Pages 258–266. [Code](https://github.com/CABSEL/SINCERITIES)
* [inferenceSnapshot](https://academic.oup.com/bioinformatics/article/31/12/i89/216346): Ocone et al. Reconstructing gene regulatory dynamics from high-dimensional single-cell snapshot data, Bioinformatics, Volume 31, Issue 12, 2015, Pages i89–i96. [Code](https://www.helmholtz-muenchen.de/fileadmin/ICB/software/inferenceSnapshot.zip)
* [SCNS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4374163/). Moignard et al. Decoding the regulatory network of early blood development from single-cell gene expression measurements. Nat Biotechnol. 2015;33(3):269-276. [Code](https://github.com/swoodhouse/SCNS-Toolkit)

### Simulators
* [SERGIO](https://www.sciencedirect.com/science/article/pii/S2405471220302878): Dibaeinia et al. SERGIO: A single-cell expression simulator guided by gene regulatory networks. Cell Syst. 2020;11(3):252-271.e11. [Code](https://github.com/PayamDiba/SERGIO)
* [dyngen](https://www.biorxiv.org/content/10.1101/2020.02.06.936971v3): Cannoodt et al. (2020) dyngen: a multi-modal simulator for spearheading new single-cell omics analyses. bioRxiv 2020.02.06.936971. [Code](https://github.com/dynverse/dyngen)

## Spatial transcriptomics

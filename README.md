# SingleCellNotes

## Journal club
**Week 1: April 22, 2021** (Perry). [SCENIC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5937676/): Aibar et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods. 2017;14(11):1083-1086. [R code](https://github.com/aertslab/SCENIC), [Python code](https://github.com/aertslab/pySCENIC). Follow-up reading:
* [A scalable SCENIC workflow for single-cell gene regulatory network analysis](https://www.nature.com/articles/s41596-020-0336-2): Van de Sande et al. Nat Protoc 15, 2247–2276 (2020). [GitHub](https://github.com/aertslab/SCENICprotocol)

**Week 2: May 6, 2021** (Huub). [A mechanistic pan-cancer pathway model informed by multi-omics data interprets stochastic cell fate responses to drugs and mitogens](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005985): Bouhaddou et al. PLoS Comput Biol. 2018 Mar 26;14(3):e1005985.  Follow-up reading:
* [Efficient Parameter Estimation Enables the Prediction of Drug Response Using a Mechanistic Pan-Cancer Pathway Model](https://www.sciencedirect.com/science/article/pii/S2405471218304381): Fröhlich et al. Cell Syst. 2018 Dec 26;7(6):567-579.e6. 

**Week 3: May 20, 2021** (Utkarsh):
* [GENIE3](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0012776): Huynh-Thu et al. Inferring Regulatory Networks from Expression Data Using Tree-Based Methods(2010). PLOS ONE 5(9): e12776. [GitHub](https://github.com/vahuynh/GENIE3)
* [dynGENIE3](https://www.nature.com/articles/s41598-018-21715-0): Huynh-Thu et al. dynGENIE3: dynamical GENIE3 for the inference of gene networks from time series expression data. Sci Rep 8, 3384 (2018). [GitHub](https://github.com/vahuynh/dynGENIE3)

**Week 4: June 3, 2021** (Aldo). [SCODE](https://academic.oup.com/bioinformatics/article/33/15/2314/3100331): Matsumoto et al. SCODE: an efficient regulatory network inference algorithm from single-cell RNA-Seq during differentiation, Bioinformatics, Volume 33, Issue 15, 2017, Pages 2314–2321. [Code](https://github.com/hmatsu1226/SCODE)

**Week 5: June 17, 2021** (Wessel). [Ridge estimation of network models from time‐course omics data](https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201700195): Miok et al. Biom J. 2019 Mar;61(2):391-405. [R code](https://cran.r-project.org/web/packages/ragt2ridges/index.html). Follow-up reading:
* [Ridge estimation of inverse covariance matrices from high-dimensional data](https://www.sciencedirect.com/science/article/abs/pii/S0167947316301141): Van Wieringen et al., 2016. Computational Statistics & Data Analysis, 103, pp.284-303. [R code](https://cran.r-project.org/web/packages/rags2ridges/index.html)
* [Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time‐course omics data](https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201500269): Miok et al. Biom J. 2017 Jan;59(1):172-191. 

**Week 6: October 7, 2021** (Utkarsh). [Scribe](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30036-3): Qiu et al. Inferring causal gene regulatory networks from coupled single-cell expression dynamics using Scribe. Cell Syst. 2020;10(3):265-274.e11. [Code](https://github.com/cole-trapnell-lab/Scribe)

**Week 7**. [scPred](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1862-5): Alquicira-Hernandez et al. scPred: accurate supervised method for cell-type classification from single-cell RNA-seq data. Genome Biology, 20:264 (2019). [Code](https://github.com/powellgenomicslab/scPred/). Follow-up reading:
* [A comparison of automatic cell identification methods for single-cell RNA sequencing data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z): Abdelaal et al. Genome Biology, 20:194 (2019). 

**Week 8**. [Beyond Predictions in Neural ODEs: Identification and Interventions](https://arxiv.org/abs/2106.12430): Aliee et al. arXiv (2021) 2106.12430.

## General
* [awesome-single-cell](https://github.com/seandavi/awesome-single-cell): List of software packages (and the people developing these methods) for single-cell data analysis, including RNA-seq, ATAC-seq, etc.
* [scRNA-tools](https://www.scrna-tools.org/): A catalogue of tools for analysing single-cell RNA sequencing data.
* A curated [database](https://docs.google.com/spreadsheets/d/1En7-UV0k0laDiIfjFkdn7dggyR7jIk3WH8QgXaMOZF0/edit#gid=0) of of single-cell transcriptomics studies.
* [Single-cell reading list](https://github.com/gtca/single-cell-reading-list): A curated selection of blog posts and papers on single-cell data analysis. 
* [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell): Featuring 295 studies and 11,739,593 cells (and counting).
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
* [Kallisto|bustools](https://www.kallistobus.tools/tutorials):  Google Colab notebooks on working with kallisto in combination with bustools.

## Atlases and compendia
* [The Tabula Sapiens: a single cell transcriptomic atlas of multiple organs from individual human donors](https://www.biorxiv.org/content/10.1101/2021.07.19.452956v1.full): The Tabula Sapiens Consortium, Stephen R Quake (2021) bioRxiv 2021.07.19.452956. [Data](https://tabula-sapiens-portal.ds.czbiohub.org/)
* [Single-nucleus cross-tissue molecular reference maps to decipher disease gene function](https://www.biorxiv.org/content/10.1101/2021.07.19.452954v1.full): Eraslan et al. (2021) bioRxiv 2021.07.19.452954.
* [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell)
* [PangloaDB](https://panglaodb.se/index.html)

## Trajectory inference
* [single-cell-pseudotime](https://github.com/agitter/single-cell-pseudotime): Overview of single-cell RNA-seq pseudotime estimation algorithms.
* [dynmethods](https://github.com/dynverse/dynmethods): A collection of 55 trajectory inference methods. To run any of these methods, interpret the results and visualise the trajectory, see the [dyno](https://github.com/dynverse/dyno) package.

### RNA velocity
* [RNA velocity—current challenges and future perspectives](https://www.embopress.org/doi/full/10.15252/msb.202110282): Bergen et al. Mol Syst Biol (2021)17:e10282.
### Real time versus pseudo time
* [Lineage tracing on transcriptional landscapes links state to fate during differentiation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7608074/): Weinreb et al. Science. 2020;367(6479):eaaw3381.

## Visualization
* [Comparison of visualization tools for single-cell RNAseq data ](https://academic.oup.com/nargab/article/2/3/lqaa052/5877814): Cakir et al. NAR Genomics and Bioinformatics, 2(3):lqaa052, 2020.

## Gene regulatory network inference

### Reviews
* [Gene regulatory network inference in single-cell biology](https://www.sciencedirect.com/science/article/abs/pii/S2452310021000184): Akers et al. Current Opinion in Systems Biology, Volume 26, 2021, Pages 87-97.
* [A comprehensive survey of regulatory network inference methods using single cell RNA sequencing data](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa190/5904505): Nguyen et al. Briefings in Bioinformatics, bbaa190 (2020).
* [Integration of single-cell multi-omics for gene regulatory network inference](https://www.sciencedirect.com/science/article/pii/S2001037020303226): Hu et al. Computational and Structural Biotechnology Journal. 2020;18:1925-1938.
* [Network Inference from Single-Cell Transcriptomic Data](https://link.springer.com/protocol/10.1007%2F978-1-4939-8882-2_10): Todorov et al. (2019) Network Inference from Single-Cell Transcriptomic Data. In: Sanguinetti G., Huynh-Thu V. (eds) Gene Regulatory Networks. Methods in Molecular Biology, vol 1883. Humana Press, New York, NY.
* [Mapping gene regulatory networks from single-cell omics data](https://academic.oup.com/bfg/article/17/4/246/4803107): Fiers et al. Briefings in Functional Genomics, Volume 17, Issue 4, July 2018, Pages 246–254.
* [Learning regulatory models for cell development from single cell transcriptomic data](https://www.sciencedirect.com/science/article/pii/S2452310017301257): Babtie et al. Current Opinion in Systems Biology 5 (2017): 72-81.

### Benchmark studies
* [Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7098173/): Pratapa et al. Nat Methods. 2020;17(2):147-154. 
* [Evaluating methods of inferring gene regulatory networks highlights their lack of performance for single cell gene expression data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2217-z): Chen et al. BMC Bioinformatics 19, 232 (2018). 

### Methods
* [locaTE](https://www.biorxiv.org/content/10.1101/2023.01.08.523176v1.full): Zhang et al. (2023) Dynamical information enables inference of gene regulation at single-cell scale. bioRxiv 2023.01.08.523176. [Code](https://github.com/zsteve/locaTE.jl)
* [scGeneRAI](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac1212/6984592): Keyl et al. (2023) Single-cell gene regulatory network prediction by explainable AI. Nucleic Acids Research, gkac1212 [Code](https://github.com/PhGK/scGeneRAI)
* [dynDeepDRIM](https://academic.oup.com/bib/article/23/6/bbac424/6720420): Xu et al. (2022) dynDeepDRIM: a dynamic deep learning model to infer direct regulatory interactions using time-course single-cell gene expression data/ Briefings in Bioinformatics, Volume 23, Issue 6, bbac424.
* [SINGE](https://www.sciencedirect.com/science/article/pii/S2211124722000493): Deshpande et al. (2022) Network inference with Granger causality ensembles on single-cell transcriptomic data. Cell Reports, Volume 38, Issue 6. [Code](https://github.com/gitter-lab/SINGE)
* [locCSN](https://www.pnas.org/doi/full/10.1073/pnas.2113178118): Wang et al. (2021) Constructing local cell-specific networks from single-cell data. PNAS, 118 (51) e2113178118. [Code](https://github.com/xuranw/locCSN)
* [DeepSEM](https://www.nature.com/articles/s43588-021-00099-8): Shu et al. (2021) Modeling gene regulatory networks using neural network architectures. Nature Computational Science, 1:491–501. [Code](https://github.com/HantaoShu/DeepSEM)
* [scPADGRN](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007471&rev=2): Zheng et al. (2020) scPADGRN: A preconditioned ADMM approach for reconstructing dynamic gene regulatory network using single-cell RNA sequencing data. PLoS Comput Biol 16(7): e1007471. [Code](https://github.com/xzheng-ac/scPADGRN)
* [GRISLI](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa576/5858974): Aubin-Frankowski et al. Gene regulation inference from single-cell RNA-seq data with linear differential equations and velocity inference, Bioinformatics, btaa576, 2020. [Code](https://github.com/PCAubin/GRISLI)
* [Scribe](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30036-3): Qiu et al. Inferring causal gene regulatory networks from coupled single-cell expression dynamics using Scribe. Cell Syst. 2020;10(3):265-274.e11. [Code](https://github.com/cole-trapnell-lab/Scribe)
* [WASABI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2798-1): Bonnaffoux et al. WASABI: a dynamic iterative framework for gene regulatory network inference. BMC Bioinformatics 20, 220 (2019).
* [M&NEM](https://academic.oup.com/bioinformatics/article/34/17/i964/5093248): Pirkl et al. Single cell network analysis with a mixture of Nested Effects Models, Bioinformatics, Volume 34, Issue 17, 2018, Pages i964–i971. [Code](https://github.com/cbg-ethz/mnem/)
* [AR1MA1 - VBEM](https://academic.oup.com/bioinformatics/article/34/6/964/4222631): Sanchez-Castillo et al. A Bayesian framework for the inference of gene regulatory networks from time and pseudo-time series data, Bioinformatics, Volume 34, Issue 6, 2018, Pages 964–970. [Code](https://github.com/mscastillo/GRNVBEM)
* [SINCERITIES](https://academic.oup.com/bioinformatics/article/34/2/258/4158033): Gao et al. SINCERITIES: inferring gene regulatory networks from time-stamped single cell transcriptional expression profiles, Bioinformatics, Volume 34, Issue 2, 2018, Pages 258–266. [Code](https://github.com/CABSEL/SINCERITIES)
* [SCENIC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5937676/): Aibar et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods. 2017;14(11):1083-1086. [R code](https://github.com/aertslab/SCENIC), [Python code](https://github.com/aertslab/pySCENIC)
* [SCODE](https://academic.oup.com/bioinformatics/article/33/15/2314/3100331): Matsumoto et al. SCODE: an efficient regulatory network inference algorithm from single-cell RNA-Seq during differentiation, Bioinformatics, Volume 33, Issue 15, 2017, Pages 2314–2321. [Code](https://github.com/hmatsu1226/SCODE)
* [inferenceSnapshot](https://academic.oup.com/bioinformatics/article/31/12/i89/216346): Ocone et al. Reconstructing gene regulatory dynamics from high-dimensional single-cell snapshot data, Bioinformatics, Volume 31, Issue 12, 2015, Pages i89–i96. [Code](https://www.helmholtz-muenchen.de/fileadmin/ICB/software/inferenceSnapshot.zip)
* [SCNS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4374163/). Moignard et al. Decoding the regulatory network of early blood development from single-cell gene expression measurements. Nat Biotechnol. 2015;33(3):269-276. [Code](https://github.com/swoodhouse/SCNS-Toolkit)

### Simulators
* [dyngen](https://www.nature.com/articles/s41467-021-24152-2): Cannoodt et al. (2021) Spearheading future omics analyses using dyngen, a multi-modal simulator of single cells. Nat Commun 12, 3942. [Code](https://github.com/dynverse/dyngen)
* [SERGIO](https://www.sciencedirect.com/science/article/pii/S2405471220302878): Dibaeinia et al. SERGIO: A single-cell expression simulator guided by gene regulatory networks. Cell Syst. 2020;11(3):252-271.e11. [Code](https://github.com/PayamDiba/SERGIO)
* [BEELINE](https://github.com/Murali-group/Beeline): Evaluation framework built on top of dynverse and BoolODE; used in Pratapa et al. (see above)

### Applications
* [Predicting gene regulatory networks from cell atlases](https://www.life-science-alliance.org/content/3/11/e202000658): Fønss Møller et al. Life Science Alliance (2020), 3:(11), e202000658.

## Integration

### Reviews
* [Computational principles and challenges in single-cell data integration](https://www.nature.com/articles/s41587-021-00895-7): Argelaguet et al. Nat Biotechnol (2021).

### Methods
* [scArches: Mapping single-cell data to reference atlases by transfer learning](https://www.nature.com/articles/s41587-021-01001-7): Lotfollahi et al. Nat Biotechnol (2021). [Code](https://github.com/theislab/scarches)
* [Seurat v4: Integrated analysis of multimodal single-cell data](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3): Hao et al. Cell, 184(13): 3573-3587 (2021).
* [Efficient and precise single-cell reference atlas mapping with Symphony](https://www.biorxiv.org/content/10.1101/2020.11.18.389189v2.full): Kang et al. (2021) bioRxiv 2020.11.18.389189. [Code](https://github.com/immunogenomics/symphony)
* [Fast, sensitive and accurate integration of single-cell data with Harmony](https://www.nature.com/articles/s41592-019-0619-0): Korsunsky et al. Nature Methods 16 (2019): 1-8.
[Code](https://github.com/immunogenomics/harmony)
* [Seurat v3: Comprehensive integration of single-cell data](https://www.sciencedirect.com/science/article/pii/S0092867419305598): Stuart et al. Cell 177.7 (2019): 1888-1902. [Code](https://satijalab.org/seurat/)
* [LIGER: Single-cell multi-omic integration compares and contrasts features of brain cell identity](https://www.sciencedirect.com/science/article/pii/S0092867419305045): Welch et al. Cell 177.7 (2019): 1873-1887. [Code](https://macoskolab.github.io/liger/). See also: [Iterative single-cell multi-omic integration using online learning](https://www.nature.com/articles/s41587-021-00867-x): Gao et al. Nat Botechnol 39, 1000–1007 (2021).

## Spatial transcriptomics

### Reviews
* [Advances in spatial transcriptomic data analysis](https://genome.cshlp.org/content/31/10/1706.full): Dries et al. Genome Research (2021).
* [Spatial omics and multiplexed imaging to explore cancer biology](https://www.nature.com/articles/s41592-021-01203-6): Lewis et al. Nat Methods (2021).
* [Exploring tissue architecture using spatial transcriptomics](https://www.nature.com/articles/s41586-021-03634-9): Rao et al. Nature 596, 211–220 (2021).
* [Integrating single-cell and spatial transcriptomics to elucidate intercellular tissue dynamics](https://www.nature.com/articles/s41576-021-00370-8): Longo et al. Nature Review Genetics (2021).
* [Uncovering an organ’s molecular architecture at single-cell resolution by spatially resolved transcriptomics](https://www.cell.com/trends/biotechnology/fulltext/S0167-7799(20)30140-2): Liao et al. Trends Biotechnol. 2020:S0167-7799(20)30140-2.
* [From whole-mount to single-cell spatial assessment of gene expression in 3D](https://www.nature.com/articles/s42003-020-01341-1): Waylen et al. Commun Biol 3, 602 (2020).
* [Deciphering cell–cell interactions and communication from gene expression](https://www.nature.com/articles/s41576-020-00292-x): Armingol et al. Nat Rev Genet (2020).

### Technologies
* [Slide-seq: A scalable technology for measuring genome-wide expression at high spatial resolution](https://science.sciencemag.org/content/363/6434/1463.full): Rodriques et al. Science. 2019;363(6434):1463-1467.

### Spatial reconstruction
* [Perler: Model-based prediction of spatial gene expression via generative linear mapping](https://www.biorxiv.org/content/10.1101/2020.05.21.107847v3.full): Okochi et al. (2020). bioRxiv 2020.05.21.107847.
* [novoSpaRc: Gene expression cartography](https://www.nature.com/articles/s41586-019-1773-3): Nitzan et al. Nature 576.7785 (2019): 132-137. [Code](https://github.com/rajewsky-lab/novosparc)
* [DistMap: The Drosophila embryo at single-cell transcriptome resolution](https://science.sciencemag.org/content/358/6360/194): Karaiskos et al. Science 358.6360 (2017): 194-199. [Code](https://github.com/rajewsky-lab/distmap)
* [High-throughput spatial mapping of single-cell RNA-seq data to tissue of origin](https://www.nature.com/articles/nbt.3209): Achim et al. Nat Biotechnol 33, 503–509 (2015).
* [Seurat: Spatial reconstruction of single-cell gene expression data](https://www.nature.com/articles/nbt.3192): Satija et al. Nat Biotechnol 33, 495–502 (2015). [Code]()

### Deconvolution
* [DestVI: Multi-resolution deconvolution of spatial transcriptomics data reveals continuous patterns of inflammation](https://www.biorxiv.org/content/10.1101/2021.05.10.443517v1.full): Lopez et al. (2021) bioRxiv 2021.05.10.443517. [Code](https://scvi-tools.org/)
* [stereoscope: Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography](https://www.nature.com/articles/s42003-020-01247-y): Andersson et al. Communications Biology 3:565 (2020).
* [SPOTlight: Seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes](https://www.biorxiv.org/content/10.1101/2020.06.03.131334v1):
Elosua et al. (2020) bioRxiv 2020.06.03.131334. [Code](https://github.com/MarcElosua/SPOTlight)

### Cell-cell interaction
* [Integrated intra- and intercellular signaling knowledge for multicellular omics analysis](https://www.embopress.org/doi/full/10.15252/msb.20209923): Türei et al. Mol Syst Biol (2021)17:e9923. [Web](https://omnipathdb.org/) 
* [SpaOTsc: Inferring spatial and signaling relationships between cells from single cell transcriptomic data](https://www.nature.com/articles/s41467-020-15968-5): Cang et al. Nature Communications 11.1 (2020): 1-13. [Code](https://github.com/zcang/SpaOTsc)
* [CSOmap: Reconstruction of cell spatial organization from single-cell RNA sequencing data based on ligand-receptor mediated self-assembly](https://www.nature.com/articles/s41422-020-0353-2): Ren et al. Cell Res 30, 763–778 (2020).
* [CellChat: Inference and analysis of cell-cell communication using CellChat](https://www.biorxiv.org/content/10.1101/2020.07.21.214387v1.full): Jin et al. (2020). bioRxiv 2020.07.21.214387. [Code](https://github.com/sqjin/CellChat), [Web](http://www.cellchat.org/)

## Perturbations
* [Machine learning for perturbational single-cell omics](https://www.cell.com/cell-systems/fulltext/S2405-4712(21)00202-7): Ji et al. Cell Systems 12(6) 522-537 (2021).

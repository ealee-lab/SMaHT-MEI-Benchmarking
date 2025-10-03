# SMaHT-MEI-Benchmarking
The SMaHT (Somatic Mosaicism across Human Tissues) MEI (Mobile Element Insertion) Working Group establishes a multi-platform integrative framework for accurate sMEI detection and source tracing in normal human tissues. 

## Data Availability 
* Sequencing samples from SMaHT: https://data.smaht.org/
* Sequencing samples from CASTLE project: https://github.com/CASTLE-Panel/castle

## Code Availability 
### sMEI detection benchmark
* xTea_mosaic (v0.1.9) : https://github.com/parklab/xTea
* MELT (v2.2.2): https://melt.igs.umaryland.edu/index.php
* RetroSom (v2): https://github.com/Czhuofu/RetroSomV2
* PALMER (v2.0.1): https://github.com/WeichenZhou/PALMER
* xTea_long (v0.1.0): https://github.com/parklab/xTea (--branch xTea_long_release_v0.1.0)
* cuteSV (v2.1.2): https://github.com/tjiangHIT/cuteSV
* Sniffles2 (v2.6.3): https://github.com/fritzsedlazeck/Sniffles
* Runninng commands for above WGS-based methods: ./benchmark/running_command.txt
* TEnCATS: https://dx.doi.org/10.17504/protocols.io.kqdg3q66ev25/v1 (For the molecular protocol), https://github.com/Boyle-Lab/NanoPal-Snakemake (For NanoPal)
* HAT-seq: https://github.com/ShaynaMallett/HATseq-pipeline/
* Call set for HapMap mixture: https://doi.org/10.5281/zenodo.17254344



### Haplotype phasing and DSA
* DSA-specific repeat library preparation in xTea_long: https://github.com/parklab/xTea/tree/master/xtea/rep_lib_prep
* LRPhasing: https://github.com/wjhlang/LRPhasing
* PhaseBlockExtension: https://github.com/wjhlang/PhaseBlockExtension 
* DSA-based analysis: https://github.com/wjhlang/SMaHT-sMEI
* Running commands for DSA-aligned samples: ./benchmark/running_command.txt
* Call set for in silico tumor-normal mixture: https://doi.org/10.5281/zenodo.17254345
* BL2009 haplotyp-resolved DSA: https://doi.org/10.5281/zenodo.17254345

## Citation 
* Wang et al., [Multi-platform framework for mapping somatic retrotransposition in human tissues](), bioxriv, 2025,

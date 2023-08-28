# SATCfinder
SATCfinder is an RNA-Sequencing pipeline designed for detection of **S**plicing **A**t **T**andem **C**AG repeat tracts. 
Although designed for CAG repeats, SATCfinder's scripts may be useful for any sort of tandem repeats. 

## What SATCfinder does
To find splice junctions in RNA-seq, you need reads which span the junction and have at least a few bases in each exon. 
In cases where a tandem CAG repeat is acting as the splice acceptor, the second exon will be highly repetitive and may 
not align well.

To overcome these mapping issues, SATCfinder computationally removes CAG repeats from reads, then aligns the non-repetitive portion to
the genome. The non-repetitive portion could be directly adjacent to a repeat in the genome (left), indicating that the
original RNA wasn't spliced. Or, the non-repetitive portion could be distant from a repeat in the genome (right),
 suggesting that splicing to an expanded CAG repeat occurred.

<img src="SATCfinder_conceptual.png" height="456" width="800">
 
## Requirements + Installation
Python 3
```
# Prepare packages
pip install "numpy>=1.20.3"
pip install "pandas>=1.2.4"
pip install "gtfparse>=1.2.1"
pip install "pysam>=0.16.0.1"
pip install "biopython>=1.79"
```

[BBTools ≥38.86](https://sourceforge.net/projects/bbmap/)

[STAR ≥2.7.1a](https://github.com/alexdobin/STAR)

SATCfinder has not been tested with previous versions of these tools, but may work.
```
# Fetch SATCfinder pipeline and scripts
git clone https://github.com/AnkurJainLab/SATCfinder.git
```

## Example usage

<details>
<summary>Processing an RNA-seq dataset with SATCfinder</summary>
</details>

<details>
<summary>Find trimmed ends</summary>
</details>

<details>
<summary>Annotating repeats in genomes</summary>
</details>



## Citation
If you find this work useful, please cite our publication at TODO

## License


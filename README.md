# repertoire_assembly

[Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for assembling repertoire data.

At the moment it contains rules to:

- Download `TRUST4` IMGT reference
- Down-sample a proportion of the sequences in a paired `fastq.gz` files into another pair of `fastq.gz` files
- Use `TRUST4` to assembly clonotypes based on paired compressed `fastq.gz` files

## References

- [Song L, Cohen D, Ouyang Z, Cao Y, Hu X, Liu XS. _"TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data_". Nature Methods (2021)](https://www.nature.com/articles/s41592-021-01142-2)
- [Greiff V, Miho E, Menzel U, Reddy ST. _"Bioinformatic and Statistical Analysis of Adaptive Immune Repertoires"_. Trend in Immunology (2015)](https://www.sciencedirect.com/science/article/abs/pii/S1471490615002239)

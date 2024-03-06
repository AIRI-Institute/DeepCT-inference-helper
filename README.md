# DeepCT inference helper

A Snakemake-based pipeline to obtain DeepCT predictions
for a VCF file.

## How to run

0. Download and unpack large files (previously stored via git-lfs):
```
wget https://nxt.2a2i.org/index.php/s/CAQYtAkPnyRJ42F/download -O DeepCT-inference-helper-large-files.zip
unzip DeepCT-inference-helper-large-files.zip
mv DeepCT-inference-helper-large-files/cache* cache/
mv DeepCT-inference-helper-large-files/* models/
rm -rf DeepCT-inference-helper-large-files*
```

1. Install necessary software:
```
conda env create -f environment.yml
conda activate deepct-inference-helper
```
(Yes, this is not the proper Snakemake way,
but currently not every used program has
a Snakemake wrapper implemented.
Also DeepCT has very fragile requirements.)

2. Pull DeepCT:
```
git submodule init
git submodule update
```

3. Provide paths to genome files and liftover chains 
in the `config/config.yml` (there's an example). 
You'll need faidx'ed hg38 and hg19, and chains for both 
hg19→hg38 and hg38→hg19.

4. Put your hg38 VCF in the `input/` directory. 

5. Run Snakemake:
```
snakemake --cores all output/file.vcf
```
Note that you put file in the `input/` but request the `output/`.

## What this does

**TLDR:** You put your hg38 VCF in `input/`, 
it converts variants to hg19, runs the model on them, 
adds three new annotation fields under the `INFO/` 
(`DEEPCT_CHANGE`, `DEEPCT_ORGANS`, `DEEPCT_CELLS`) 
to variants with significant hits, converts them back to hg38, 
and puts hg38 VCF in `output/`.

**In detail:**

1. Two-stage liftover with _CrossMap_: one as usual, 
and then variants that failed due to allele mismatch with the hg19 
will have their REF and ALT swapped, 
most of them will be successfully lifted.

2. Filtering. It will remove indels, variants on chrM, 
and variants with non-AGTC alleles.

3. Cache check. Variants are checked against the cache
(simply via _bedtools intersect_) and separated in two parts.

4. Convert uncached to TSV. This is the format (CHROM POS REF ALT)
that our fork of Selene uses.

5. Create a configuration file for Selene. 

6. Run Selene.

7. Extract predictions (from numpy-array) and add them to VCF.

8. Copy old cache and newly predicted variants to a new cache 
and then replace the old one with it.

9. Meanwhile the second part (cached) gets annotated from a cache.

10. Merge both parts back.

11. Two-stage liftover just as in pt.1, 
only this time back to hg38 from hg19.

## Computing resources

* Inference is slow.
One Tesla V100 processes ~6.5 variants per second. 
24 cores of Xeon E5-2690 process 1 variant in 20 seconds.

* That's why there's a cache implemented.
All your processed variants will be added to 
your local cache at `cache/cache.vcf.gz`.
Cache is another potential bottleneck by itself. 
At ~1M cached variants its delays are already in single minutes.

* RAM requirements are ~1G per 2000 variants. 
It's better to split your files in chunks (of like 50K).
Just remember keep the VCF headers in all of them.

* GPU memory usage is stable at 2800 MB.

* You can change CUDA device in the 
`templates/inference.template.yml` at line 27.
Or use just CPU (comment out that line). Please refer to the 
[Selene documentation](https://selene.flatironinstitute.org/master/overview/cli.html).

## General notes

* You might want to clean `tmp/` manually afterwards.

* You might (and probably will) lose some variants
due to lifting to hg19 and back.

* Indels will be forcefully filtered out.
Selene has to be modified further to support these.
Mitochondrial variants are also skipped.

* Cache is in hg19. Your input and output files are in hg38.

## Authors & License

DeepCT by Sindeeva et al. 
[Cell type-specific interpretation of noncoding variants using deep learning-based methods](https://doi.org/10.1101/2021.12.31.474623).

Provided under Apache License 2.0.

© 2022 Autonomous Non-Profit Organization 
"Artificial Intelligence Research Institute" (AIRI). 
All rights reserved.


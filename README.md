# yolotrim

## Heuristic read trimming for PacBio IsoSeq data

`yolotrim` removes poly-A and primer sequences from Iso-Seq FASTQ files when you don't have the HiFi BAM files need to run [lima](https://lima.how/get-started.html).

## Usage

```sh
yolotrim --in weird.fastq.gz --out weird.trimmed.fastq.gz
```

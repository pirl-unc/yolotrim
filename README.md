# yolotrim
## Heuristic read trimming for PacBio IsoSeq data

Supposedly the CCS IsoSeq reads we received should be just the transcript sequence with nothing extraneous like adapters or primers. And yet, they all start and end with suspiciously similar sequence fragments. While we wait for our sequencing CRO to dig deeper into this mystery, let's YOLO trim the reads. 

## Usage

```sh
yolotrim --in weird.fastq.gz --out weird.trimmed.fastq.gz
```

# Illumina PE Quality Control

A pipeline for Illumina read QC using `Trimmomatic` and `Flash`.

## Dependencies

### External programs

- [FixPairs](https://github.com/mooreryan/FixPairs)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [FLASH](https://ccb.jhu.edu/software/FLASH/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) -- included
- [Java](https://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html)
- [Ruby](https://www.ruby-lang.org/)

### Ruby Gems

- [Optimist](https://rubygems.org/gems/optimist)
- [abort_if](https://rubygems.org/gems/abort_if)
- [systemu](https://rubygems.org/gems/systemu)

## Pipeline opts

Some are hardcoded....

```
SEED_MISMATCHES = 2
PALINDROME_CLIP_THRESHOLD = 30
SIMPLE_CLIP_THRESHOLD = 10
WINDOW_SIZE = 4
QUAL = 15

MIN_LEN = 50
MAX_OVERLAP = 250
```
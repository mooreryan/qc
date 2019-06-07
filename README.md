# Illumina PE Quality Control

A pipeline for Illumina read QC using `Trimmomatic` and `Flash`.

## Installation

This assumes you already have a working version of Ruby on your system.


```bash
# Download the code.
wget https://github.com/mooreryan/qc/archive/v0.5.4.tar.gz

# Extract the source code from the archive.
tar xzf v0.5.4.tar.gz

# Enter source directory
cd qc-0.5.4

# Change the permissions for the two pipeline scripts.
chmod 755 qc.rb && qc_multilib_wrapper.rb

# Update PATH variable
export PATH=${PATH}:${PWD}

# Install Ruby dependencies
gem install bundler && bundle install
```

### Dependencies

Make sure you've met all the dependencies!

#### External programs

- [FixPairs](https://github.com/mooreryan/FixPairs)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [FLASH](https://ccb.jhu.edu/software/FLASH/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) -- included
- [Java](https://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html)
- [Ruby](https://www.ruby-lang.org/)

#### Ruby Gems

- [Bundler](https://rubygems.org/gems/bundler) -- to install Ruby dependencies
- [Optimist](https://rubygems.org/gems/optimist)
- [abort_if](https://rubygems.org/gems/abort_if)
- [systemu](https://rubygems.org/gems/systemu)

These can be installed by running `bundle install` in the directory where you downloaded the source code.

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
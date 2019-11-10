#!/usr/bin/env ruby

# Copyright 2015 - 2018 Ryan Moore
# Contact: moorer@udel.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

Signal.trap("PIPE", "EXIT")

require "optimist"
require "abort_if"
require "systemu"
require "fileutils"

require_relative "lib/core_ext/process"
require_relative "lib/qc/utils"
require_relative "lib/qc/version"

include AbortIf
include QC::Utils
Process.extend CoreExt::Process

opts = Optimist.options do
  version QC::Version::VERSION_BANNER

  banner <<-EOS
#{QC::Version::VERSION_BANNER}

  Run QC pipeline on Illumina reads.

  When specifying forward and reverse reads, you could use globs like
  this: --forward my_reads/*.1.fq.gz --reverse my_reads/*.2.fq.gz

  Gives each library its own nice little folder.  It assumes that the
  sample name is in the reads (it gets the sample name from the
  forward reads specifically).  For example, any of the following
  would be considered to have the sample name RI_1A1:

    - RI_1A1.1.fq.gz
    - RI_1A1.1.fq.bz2
    - RI_1A1.1.fastq.gz
    - RI_1A1.1.fastq.bz2
    - RI_1A1.f.fq.gz
    - RI_1A1.f.fq.bz2
    - RI_1A1.f.fastq.gz
    - RI_1A1.f.fastq.bz2
    - RI_1A1.forward.fq.gz
    - RI_1A1.forward.fq.bz2
    - RI_1A1.forward.fastq.gz
    - RI_1A1.forward.fastq.bz2

  Also note that the sample name will be determined if the reads do
  not end in `.gz` or `.bz2` as well.  Finally, the letters can be
  upper or lowercase, and it should turn out okay.

  If the pipeline actually finishes, it will put an empty file in the
  outdirectory called qc_v#{QC::Version::VERSION}_done.

  Don't forget to check the log file and make sure everything looks
  good!!! (Both the log that this program writes to stderr as well as
  the one in the final output directory...they report different
  things.)

  Options:
  EOS

  opt(:forward, "forward reads", type: :strings)
  opt(:reverse, "reverse reads", type: :strings)

  opt(:threads, "Threads", type: :integer, default: 10)

  opt(:outdir,
      "Output directory",
      type: :string,
      default: "qc")
  opt(:gzip_output,
      "Compress output or not.",
      default: true)

  opt(:trimmomatic_headcrop,
      "Remove N number of bases from the start of reads",
      short: "c",
      default: 0)
  opt(:trimmomatic_seed_mismatches,
      "Seed mismatches for trimmomatic",
      default: 2)
  opt(:trimmomatic_palindrome_clip_threshold,
      "Palindrome clip threshold for trimmomatic",
      default: 30)
  opt(:trimmomatic_simple_clip_threshold,
      "Simple clip threshold for trimmomatic",
      default: 10)
  opt(:trimmomatic_window_size,
      "Window size for trimmomatic",
      default: 4)
  opt(:trimmomatic_window_quality,
      "Window quality for trimmomatic",
      default: 15)
  opt(:trimmomatic_min_len,
      "Min. length for trimmomatic",
      default: 50)
  opt(:trimmomatic_jar,
      "Path to trimmomatic jar",
      default: File.join(File.dirname(__FILE__),
                         "bin",
                         "trimmomatic-0.35",
                         "trimmomatic-0.35.jar"))
  opt(:adapters,
      "Path to adapter sequences to trim",
      default: File.join(File.dirname(__FILE__),
                         "bin",
                         "trimmomatic-0.35",
                         "adapters",
                         "all.fa"))

  opt(:flash_max_overlap,
      "Max. overlap before penalty for FLASH",
      default: 250)

  opt(:bowtie_idx,
      "The bowtie2 index to screen reads against " +
      "(can provide more than one)",
      type: :strings)
  opt(:bowtie_seed,
      "The seed for bowtie",
      default: 123123)
end

check_files *opts[:forward]
check_files *opts[:reverse]

gzip_option = opts[:gzip_output] ? "--gzip-output" : "--no-gzip-output"

$stderr.puts opts.inspect

abort_unless opts[:forward].count == opts[:reverse].count,
             "You need to provide equal numbers of forward " +
             "and reverse files"

qc = File.join File.dirname(__FILE__), "qc.rb"

opts[:forward].each_with_index do |for_f, idx|
  rev_f = opts[:reverse][idx]

  AbortIf.logger.info { "Working on lib ##{idx} " +
                        "(#{for_f} and #{rev_f})" }

  basename = remove_fq_ext File.basename(for_f)

  outd = File.join opts[:outdir], basename

  if opts[:bowtie_idx]
    cmd = "ruby #{qc} " +
          "--forward #{for_f} " +
          "--reverse #{rev_f} " +
          "--threads #{opts[:threads]} " +
          "#{gzip_option} " +
          "--trimmomatic-headcrop #{opts[:trimmomatic_headcrop]} " +
          "--trimmomatic-seed-mismatches #{opts[:trimmomatic_seed_mismatches]} " +
          "--trimmomatic-palindrome-clip-threshold #{opts[:trimmomatic_palindrome_clip_threshold]} " +
          "--trimmomatic-simple-clip-threshold #{opts[:trimmomatic_simple_clip_threshold]} " +
          "--trimmomatic-window-size #{opts[:trimmomatic_window_size]} " +
          "--trimmomatic-window-quality #{opts[:trimmomatic_window_quality]} " +
          "--trimmomatic-min-len #{opts[:trimmomatic_min_len]} " +
          "--trimmomatic-jar #{opts[:trimmomatic_jar]} " +
          "--adapters #{opts[:adapters]} " +
          "--flash-max-overlap #{opts[:flash_max_overlap]} " +
          "--outdir #{outd} " +
          "--basename #{basename} " +
          "--bowtie-idx #{opts[:bowtie_idx].join(" ")} " +
          "--bowtie-seed #{opts[:bowtie_seed]}"
  else
    cmd = "ruby #{qc} " +
          "--forward #{for_f} " +
          "--reverse #{rev_f} " +
          "--threads #{opts[:threads]} " +
          "#{gzip_option} " +
          "--trimmomatic-headcrop #{opts[:trimmomatic_headcrop]} " +
          "--trimmomatic-seed-mismatches #{opts[:trimmomatic_seed_mismatches]} " +
          "--trimmomatic-palindrome-clip-threshold #{opts[:trimmomatic_palindrome_clip_threshold]} " +
          "--trimmomatic-simple-clip-threshold #{opts[:trimmomatic_simple_clip_threshold]} " +
          "--trimmomatic-window-size #{opts[:trimmomatic_window_size]} " +
          "--trimmomatic-window-quality #{opts[:trimmomatic_window_quality]} " +
          "--trimmomatic-min-len #{opts[:trimmomatic_min_len]} " +
          "--trimmomatic-jar #{opts[:trimmomatic_jar]} " +
          "--adapters #{opts[:adapters]} " +
          "--flash-max-overlap #{opts[:flash_max_overlap]} " +
          "--outdir #{outd} " +
          "--basename #{basename} " +
          "--bowtie-seed #{opts[:bowtie_seed]}"

  end

  Process.run_it! cmd
end

AbortIf.logger.info { "qc.rb wrapper has finished" }
AbortIf.logger.info { "Final outdir is #{opts[:outdir]}" }

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

require "abort_if"
require "systemu"
require "fileutils"
require "optimist"

require_relative "lib/core_ext/process"
require_relative "lib/qc/utils"
require_relative "lib/qc/version"

include AbortIf
include AbortIf::Assert
include QC::Utils

Process.extend CoreExt::Process

Signal.trap("PIPE", "EXIT")

opts = Optimist.options do
  version QC::Version::VERSION_BANNER

  banner <<-EOS
#{QC::Version::VERSION_BANNER}

  Run QC pipeline on Illumina reads.

  It is HARDCODED to use phred33 right now (new illumina, sanger) as
  sometimes very few reads will be in the 1U and 2U files and
  Trimmomatic can't determine it.

  This version of the pipeline only works with a single forward and a
  single reverse read set.  If you git it a glob like reads/*.1.fq it
  will just use the first one.

  If you think you want to use this program, you probably really want
  to use qc_multilib_wrapper.rb instead.

  Options:
  EOS

  opt(:forward, "forward", type: :string)
  opt(:reverse, "reverse", type: :string)

  opt(:threads, "Threads", type: :integer, default: 10)

  opt(:basename,
      "Basename for read files",
      default: "reads")
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

BOWTIE_SEED = opts[:bowtie_seed]

TRIMMO = opts[:trimmomatic_jar]
abort_unless_file_exists TRIMMO

TRIMSEQS = opts[:adapters]
abort_unless_file_exists TRIMSEQS

FLASH = `which flash`.chomp
abort_if FLASH.empty?, "Missing flash"

if opts[:bowtie_idx]
  BOWTIE = `which bowtie2`.chomp
  abort_if BOWTIE.empty?, "Missing bowtie2"

  FIX_PAIRS = `which FixPairs`.chomp
  abort_if FIX_PAIRS.empty?, "Missing FixPairs"
end

abort_if opts[:basename].empty?,
         "--basename option was empty!"

BASENAME = opts[:basename]

HEADCROP = opts[:trimmomatic_headcrop]

SEED_MISMATCHES = opts[:trimmomatic_seed_mismatches]
PALINDROME_CLIP_THRESHOLD = opts[:trimmomatic_palindrome_clip_threshold]
SIMPLE_CLIP_THRESHOLD = opts[:trimmomatic_simple_clip_threshold]

WINDOW_SIZE = opts[:trimmomatic_window_size]
QUAL = opts[:trimmomatic_window_quality]
THREADS = opts[:threads]

MIN_LEN = opts[:trimmomatic_min_len]

MAX_OVERLAP = opts[:flash_max_overlap]

java = `which java`.chomp
abort_if java.empty?, "Missing java"

now = Time.now.strftime "%Y%m%d%H%M%S%L"
big_log = File.join opts[:outdir], "qc_log.#{now}.txt"
baseout = File.join opts[:outdir], BASENAME

$stderr.puts opts.inspect

check_files opts[:forward]
check_files opts[:reverse]

# TODO check that index files exist

abort_if File.exists?(opts[:outdir]),
         "Outdir #{opts[:outdir]} already exists!"
FileUtils.mkdir_p opts[:outdir]

in_forward = opts[:forward]
in_reverse = opts[:reverse]

baseout += ".adpater_trimmed"

adapter_trim!(in1: in_forward,
              in2: in_reverse,
              baseout: baseout,
              log: big_log)

out_1P = baseout + "_1P"
out_1U = baseout + "_1U"
out_2P = baseout + "_2P"
out_2U = baseout + "_2U"

check_files out_1P, out_1U, out_2P, out_2U

out = out_1U + ".qual_trimmed"
out_1U = qual_trim_se!(inf: out_1U, out: out, log: big_log)
check_files out_1U

out = out_2U + ".qual_trimmed"
out_2U = qual_trim_se!(inf: out_2U, out: out, log: big_log)
check_files out_2U

flash_dir = File.join opts[:outdir], "flash"
flash!(in1: out_1P, in2: out_2P, outdir: flash_dir, log: big_log)

out_flash_single = File.join opts[:outdir],
                             "#{BASENAME}.adapter_trimmed.flash_combined"
out_flash_1P = File.join opts[:outdir],
                         "#{BASENAME}.adapter_trimmed.flash_notcombined_1P"

out_flash_2P = File.join opts[:outdir],
                         "#{BASENAME}.adapter_trimmed.flash_notcombined_2P"

check_files out_flash_single, out_flash_1P, out_flash_2P

out = out_flash_single + ".qual_trimmed"
out_flash_single = qual_trim_se!(inf: out_flash_single,
                                 out: out,
                                 log: big_log)

check_files out_flash_single

baseout = out_flash_1P.sub(/_1P$/, "") + ".qual_trimmed"
qual_trim_pe!(in1: out_flash_1P,
              in2: out_flash_2P,
              baseout: baseout,
              log: big_log)

out_flash_1P = baseout + "_1P"
out_flash_2P = baseout + "_2P"
out_flash_1U = baseout + "_1U"
out_flash_2U = baseout + "_2U"

# these are surviving outfiles
check_files out_flash_1P,
            out_flash_2P,
            out_flash_single,
            out_flash_1U,
            out_flash_2U,
            out_1U,
            out_2U

out_paired_1 = File.join opts[:outdir], "#{BASENAME}.1.fq"
out_paired_2 = File.join opts[:outdir], "#{BASENAME}.2.fq"
out_unpaired = File.join opts[:outdir], "#{BASENAME}.U.fq"


Process.run_it! "mv #{out_flash_1P} #{out_paired_1}"
Process.run_it! "mv #{out_flash_2P} #{out_paired_2}"

Process.run_it! "cat " +
                "#{out_flash_single} " +
                "#{out_flash_1U} " +
                "#{out_flash_2U} " +
                "#{out_1U} " +
                "#{out_2U} " +
                "> #{out_unpaired}"

Process.run_it "rm " +
               "#{out_flash_single} " +
               "#{out_flash_1U} " +
               "#{out_flash_2U} " +
               "#{out_1U} " +
               "#{out_2U}"

######################################################################
# genome screen
###############

# Now the files are
#   out_paired_1
#   out_paired_2
#   out_unpaired

if opts[:bowtie_idx]

  # Each iteration will replace the out_paired_? variables with the
  # reads that make it through the screen.
  opts[:bowtie_idx].each do |index_fname|
    idx_base = File.basename index_fname

    # make screen outdir
    reads_that_hit_genome_dir = File.join opts[:outdir],
                                          "those_that_hit_#{idx_base}"
    FileUtils.mkdir_p reads_that_hit_genome_dir

    # screen the reads
    out_paired_1_good_reads, out_paired_1_bad_reads =
                             screen!(index: index_fname,
                                     reads: out_paired_1,
                                     seed: BOWTIE_SEED,
                                     log: big_log)

    out_paired_2_good_reads, out_paired_2_bad_reads =
                             screen!(index: index_fname,
                                     reads: out_paired_2,
                                     seed: BOWTIE_SEED + 1,
                                     log: big_log)

    out_unpaired_good_reads, out_unpaired_bad_reads =
                             screen!(index: index_fname,
                                     reads: out_unpaired,
                                     seed: BOWTIE_SEED + 2,
                                     log: big_log)

    reads_that_hit_genome = File.join reads_that_hit_genome_dir,
                                      "reads_that_hit_#{idx_base}.U.fq"

    # combine all the reads that did align
    Process.run_it! "cat " +
                    [out_paired_1_bad_reads,
                     out_paired_2_bad_reads,
                     out_unpaired_bad_reads].join(" ") +
                    " > #{reads_that_hit_genome}"

    # Delete the files that were just catted into the new file
    Process.run_it! "rm " +
                    [out_paired_1_bad_reads,
                     out_paired_2_bad_reads,
                     out_unpaired_bad_reads].join(" ")

    # fix the pair info for the reads that did not hit the genome
    out1, out2, outU = fix_pairs!(in1: out_paired_1_good_reads,
                                  in2: out_paired_2_good_reads,
                                  log: big_log,
                                  outdir: opts[:outdir])

    # remove the two files that went into fix pairs command
    Process.run_it! "rm " +
                    [out_paired_1_good_reads,
                     out_paired_2_good_reads].join(" ")

    # Replace the old out_unpaired file with the new clean reads.
    Process.run_it! "cat " +
                    [out_unpaired_good_reads, outU].join(" ") +
                    " > #{out_unpaired}"
    Process.run_it! "rm " + [out_unpaired_good_reads, outU].join(" ")

    # Replece the old out_paired_1 file with the new screened paired 1
    # reads.
    Process.run_it! "cat " +
                    "#{out1} > #{out_paired_1}"
    Process.run_it "rm #{out1}"

    # Replece the old out_paired_1 file with the new screened paired 1
    # reads.
    Process.run_it! "cat " +
                    "#{out2} > #{out_paired_2}"
    Process.run_it "rm #{out2}"
  end
end

###############
# genome screen
######################################################################

if opts[:gzip_output]
  pigz = `which pigz`.chomp
  if pigz.empty?
    gzip = `which gzip`.chomp

    if gzip.empty?
      AbortIf.logger.warn { "Cannot locate pigz or gzip. Output " +
                            "files will not be zipped." }
    else
      Process.run_it "#{gzip} " +
                     "#{out_unpaired} " +
                     "#{out_paired_1} " +
                     "#{out_paired_2}"
    end
  else
    Process.run_it "#{pigz} --processes #{THREADS} " +
                   "#{out_unpaired} " +
                   "#{out_paired_1} " +
                   "#{out_paired_2}"

  end
end

done_file = File.join(opts[:outdir],
                      "qc_v#{QC::Version::VERSION}_done")
Process.run_it "touch #{done_file}"

AbortIf.logger.info { "QC finished" }
AbortIf.logger.info { "Output directory: #{opts[:outdir]}" }

#!/usr/bin/env ruby

# Copyright 2015 - 2016 Ryan Moore
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

require "parse_fasta"
require "abort_if"
require "systemu"
require "fileutils"
require "trollop"

require_relative "lib/core_ext/process"
require_relative "lib/qc/utils"

include AbortIf
include AbortIf::Assert
include QC::Utils

Process.extend CoreExt::Process

Signal.trap("PIPE", "EXIT")

VERSION = "
    Version: 0.2.1
  Copyright: 2015 - 2016 Ryan Moore
    Contact: moorer@udel.edu
    Website: https://github.com/mooreryan/qc
    License: GPLv3

"

opts = Trollop.options do
  version VERSION

  banner <<-EOS
#{VERSION}

  Run QC pipeline on Illumina reads.

  Options:
  EOS

  opt(:forward, "forward", type: :strings)
  opt(:reverse, "reverse", type: :strings)

  opt(:threads, "Threads", type: :integer, default: 10)

  opt(:outdir, "Output directory", type: :string,
      default: "qc")
end

TRIMMO = File.join File.dirname(__FILE__),
                   "bin",
                   "trimmomatic-0.35",
                   "trimmomatic-0.35.jar"

FLASH = File.join File.dirname(__FILE__),
                  "bin",
                  "flash"

SEED_MISMATCHES = 2
PALINDROME_CLIP_THRESHOLD = 30
SIMPLE_CLIP_THRESHOLD = 10

WINDOW_SIZE = 4
QUAL = 15
THREADS = opts[:threads]

MIN_LEN = 50

MAX_OVERLAP = 250

TRIMSEQS = File.join File.dirname(__FILE__),
                     "bin",
                     "trimmomatic-0.35",
                     "adapters",
                     "TruSeq3-PE-both.fa"

java = `which java`.chomp
abort_if java.empty?, "Missing java"

now = Time.now.strftime "%Y%m%d%H%M%S%L"
big_log = File.join opts[:outdir], "qc_log.#{now}.txt"
baseout = File.join opts[:outdir], "reads"

check_files *opts[:forward]
check_files *opts[:reverse]

in_forward = File.join opts[:outdir], "tmp.1.fq"
in_reverse = File.join opts[:outdir], "tmp.2.fq"

abort_if File.exists?(opts[:outdir]),
         "Outdir #{opts[:outdir]} already exists!"

FileUtils.mkdir_p opts[:outdir]

cat_fastq_files in_forward, *opts[:forward]
cat_fastq_files in_reverse, *opts[:reverse]

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
                             "reads.adapter_trimmed.flash_combined"
out_flash_1P = File.join opts[:outdir],
                         "reads.adapter_trimmed.flash_notcombined_1P"

out_flash_2P = File.join opts[:outdir],
                         "reads.adapter_trimmed.flash_notcombined_2P"

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

out_paired_1 = File.join opts[:outdir], "reads.1.fq"
out_paired_2 = File.join opts[:outdir], "reads.2.fq"
out_unpaired = File.join opts[:outdir], "reads.unpaired.fq"

outfasta_d = File.join opts[:outdir], "for_idba"
FileUtils.mkdir_p outfasta_d

out_paired_interleaved_fa = File.join outfasta_d,
                                      "reads.1_and_2_interleaved.fa"
out_unpaired_fa = File.join outfasta_d,
                            "reads.unpaired.fa"

Process.run_it! "mv #{out_flash_1P} #{out_paired_1}"
Process.run_it! "mv #{out_flash_2P} #{out_paired_2}"

Process.run_it! "cat " +
                "#{out_flash_single} " +
                "#{out_flash_1U} " +
                "#{out_flash_2U} " +
                "#{out_1U} " +
                "#{out_2U} " +
                "> #{out_unpaired}"

Process.run_it! "rm " +
                "#{in_forward} " +
                "#{in_reverse} " +
                "#{out_flash_single} " +
                "#{out_flash_1U} " +
                "#{out_flash_2U} " +
                "#{out_1U} " +
                "#{out_2U}"


fq2fa = `which fq2fa`.chomp
if fq2fa.empty?
  AbortIf.logger.warn { "Cannot locate fq2fa program. No IDBA " +
                        "files will be made" }
else
  Process.run_it! "#{fq2fa} " +
                  "--filter " +
                  "--merge " +
                  "#{out_paired_1} " +
                  "#{out_paired_2} " +
                  "#{out_paired_interleaved_fa}"

  Process.run_it! "#{fq2fa} " +
                  "--filter " +
                  "#{out_unpaired} " +
                  "#{out_unpaired_fa}"
end

pigz = `which pigz`.chomp
if pigz.empty?
  gzip = `which gzip`.chomp

  if gzip.empty?
    AbortIf.logger.warn { "Cannot locate pigz or gzip. Output " +
                          "files will not be zipped." }
  else
    Process.run_it "#{gzip} " +
                    "#{out_unpaired} " +
                    "#{out_unpaired_fa} " +
                    "#{out_paired_1} " +
                    "#{out_paired_2} " +
                    "#{out_paired_interleaved_fa}"
  end
else
  Process.run_it "#{pigz} --best --processes #{THREADS} " +
                 "#{out_unpaired} " +
                 "#{out_unpaired_fa} " +
                 "#{out_paired_1} " +
                 "#{out_paired_2} " +
                 "#{out_paired_interleaved_fa}"

end

AbortIf.logger.info { "QC finished" }
AbortIf.logger.info { "Output directory: #{opts[:outdir]}" }

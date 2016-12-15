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

# Changelog ##########################################################
#
# v0.3.0 - 2016-12-14 (RMM) - If provided only a single library, don't
#                             cat the input files.
#
# v0.3.1 - 2016-12-14 (RMM) - Bugfix
#
# v0.3.2 - 2016-12-14 (RMM) - Don't rm #{in_forward} and #{in_reverse}
#                             if there is only one library, cos in
#                             this case, these are the actual input
#                             files.
#
######################################################################

require "parse_fasta"
require "abort_if"
require "systemu"
require "fileutils"
require "trollop"

# Monkey patch to handle the odd gzipped Illumina files from the
# 2016_09 MMGs
class FastqFile < File
  def each_record_fast
    count = 0
    header = ''
    sequence = ''
    description = ''
    quality = ''

    begin
      begin
        f = Zlib::GzipReader.open(self)
        f = IO.popen("gzip -cd #{self.path}")
        AbortIf.logger.debug { "f is #{f.class}, #{f.inspect}" }
      rescue Zlib::GzipFile::Error => e
        f = self
      end

      num = 0
      f.each_line do |line|
        # debug line[0..10]
        num += 1
        line.chomp!

        case count
        when 0
          header = line[1..-1]
        when 1
          sequence = line
        when 2
          description = line[1..-1]
        when 3
          count = -1
          quality = line
          yield(header, sequence, description, quality)
        end

        count += 1
      end

      f.close if f.instance_of?(Zlib::GzipReader)
      STDERR.puts "DEBUG -- Total lines from pf is #{num}"
      return f
    ensure
      f.close if f
    end
  end
end

require_relative "lib/core_ext/process"
require_relative "lib/qc/utils"

include AbortIf
include AbortIf::Assert
include QC::Utils

Process.extend CoreExt::Process

Signal.trap("PIPE", "EXIT")

VERSION = "
    Version: 0.3.2
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

abort_if File.exists?(opts[:outdir]),
         "Outdir #{opts[:outdir]} already exists!"
FileUtils.mkdir_p opts[:outdir]

# only one library don't cat the files
if opts[:forward].length == 1 && opts[:reverse].length == 1
  AbortIf.logger.info { "Only one forward and one reverse file " +
                        "provided. Not catting files." }
  in_forward = opts[:forward].first
  in_reverse = opts[:reverse].first
else
  AbortIf.logger.info { "Multiple libraries provided. Catting files" }

  in_forward = File.join opts[:outdir], "tmp.1.fq"
  in_reverse = File.join opts[:outdir], "tmp.2.fq"

  cat_fastq_files in_forward, *opts[:forward]
  cat_fastq_files in_reverse, *opts[:reverse]
end

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

# only one library delete infiles
if opts[:forward].length == 1 && opts[:reverse].length == 1
  AbortIf.logger.info { "Only one forward and one reverse file " +
                        "provided. Not gonna delete them!" }

  Process.run_it! "rm " +
                  "#{out_flash_single} " +
                  "#{out_flash_1U} " +
                  "#{out_flash_2U} " +
                  "#{out_1U} " +
                  "#{out_2U}"
else
  Process.run_it! "rm " +
                  "#{in_forward} " +
                  "#{in_reverse} " +
                  "#{out_flash_single} " +
                  "#{out_flash_1U} " +
                  "#{out_flash_2U} " +
                  "#{out_1U} " +
                  "#{out_2U}"
end


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

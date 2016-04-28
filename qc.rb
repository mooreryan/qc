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

include AbortIf
include AbortIf::Assert

module CoreExtensions
  module Process
    def run_it *a, &b
      logger.debug { "Running: #{a.join(" and ")}" }

      exit_status, stdout, stderr = systemu *a, &b

      puts stdout unless stdout.empty?
      $stderr.puts stderr unless stderr.empty?

      exit_status.exitstatus
    end

    def run_it! *a, &b
      exit_status = self.run_it *a, &b

      abort_unless exit_status.zero?,
                   "ERROR: non-zero exit status (#{exit_status})"

      exit_status
    end
  end
end
Process.extend CoreExtensions::Process

def cat_fastq_files outf, *fnames
  File.open(outf, "w") do |f|
    fnames.each do |fname|
      FastqFile.open(fname).each_record_fast do |head, seq, desc, qual|
        f.puts "@#{head}\n#{seq}\n+#{desc}\n#{qual}"
      end
    end
  end
end

def check_files *fnames
  fnames.each do |fname|
    abort_unless_file_exists fname
  end
end

def eputs msg=""
  $stderr.puts msg
end

def seqcount fname
  if fname.match(/.gz$/)
    num_seqs = (`gunzip -c #{fname} | wc -l`.to_f / 4).round
  else
    num_seqs = (`wc -l #{fname}`.to_f / 4).round
  end

  num_seqs
end

def qual_trim_se!(inf:, out:, log:)
  count = seqcount inf
  if count >= 1
    cmd = "java -jar #{TRIMMO} SE " +
          "-threads #{THREADS} " +
          "#{inf} " +
          "#{out} " +
          "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
          "MINLEN:#{MIN_LEN} " +
          ">> #{log} 2>&1"

    Process.run_it! cmd
  else
    warn "WARN: no reads in #{inf}"
  end

  Process.run_it! "rm #{inf}"

  out
end

def qual_trim_pe! in1:, in2:, baseout:, log:
                     cmd = "java -jar #{TRIMMO} PE " +
                           "-threads #{THREADS} " +
                           "#{in1} " +
                           "#{in2} " +
                           "-baseout #{baseout} " +
                           "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
                           "MINLEN:#{MIN_LEN} " +
                           ">> #{log} 2>&1"

  Process.run_it! cmd
  Process.run_it! "rm #{in1} #{in2}"
end

def flash! in1:, in2:, outdir:, log:
              cmd = "#{FLASH} " +
                    "--threads #{THREADS} " +
                    "--output-prefix flashed " +
                    "--max-overlap #{MAX_OVERLAP} " +
                    "#{in1} " +
                    "#{in2} " +
                    "--output-directory #{outdir} " +
                    ">> #{log} 2>&1"

  Process.run_it! cmd
  Process.run_it! "rm #{in1} #{in2}"
  Process.run_it!("mv #{outdir}/flashed.extendedFrags.fastq " +
                  "#{outdir}/../reads.adapter_trimmed.flash_combined")
  Process.run_it!("mv #{outdir}/flashed.notCombined_1.fastq " +
                  "#{outdir}/../reads.adapter_trimmed.flash_notcombined_1P")
  Process.run_it!("mv #{outdir}/flashed.notCombined_2.fastq " +
                  "#{outdir}/../reads.adapter_trimmed.flash_notcombined_2P")
end

def adapter_trim!(in1:, in2:, baseout:, log:)
  # Trim the adapters
  cmd = "java -jar #{TRIMMO} PE " +
        "-threads #{THREADS} " +
        "-baseout #{baseout} " +
        "#{in1} " +
        "#{in2} " +
        "ILLUMINACLIP:" +
        "#{TRIMSEQS}:" +
        "#{SEED_MISMATCHES}:" +
        "#{PALINDROME_CLIP_THRESHOLD}:" +
        "#{SIMPLE_CLIP_THRESHOLD} " +
        ">> #{log} 2>&1"

  Process.run_it! cmd
end


Signal.trap("PIPE", "EXIT")

VERSION = "
    Version: 0.1.0
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
      default: "one_lib_with_flash")
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


gzip = `which pigz`.chomp
gzip = `which gzip`.chomp if gzip.empty?

unless gzip.empty?
  Process.run_it! "#{gzip} --best --processes #{THREADS} " +
                  "#{out_unpaired} " +
                  "#{out_paired_1} " +
                  "#{out_paired_2}"
end

#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

require "trollop"
require "abort_if"
require "systemu"
require "fileutils"

require_relative "lib/core_ext/process"
require_relative "lib/qc/utils"


include AbortIf
include QC::Utils
Process.extend CoreExt::Process



VERSION = "
    Version: 0.4.0
    Copyright: 2015 - 2017 Ryan Moore
    Contact: moorer@udel.edu
    Website: https://github.com/mooreryan/qc
    License: GPLv3

"

opts = Trollop.options do
  version VERSION

  banner <<-EOS
#{VERSION}

  When you don't want your libraries combined into one big set of
  reads, use this wrapper instead of the qc.rb script.

  You can also use this if you only have one library as well.

  IDBA takes one file so use the --idba option to try and make it.

  Options:
  EOS

  opt(:forward, "forward", type: :strings)
  opt(:reverse, "reverse", type: :strings)

  opt(:threads, "Threads", type: :integer, default: 10)

  opt(:outdir, "Output directory", type: :string,
      default: "qc")

  opt(:idba, "Make outfiles for IDBA", type: :boolean,
      defaul: false)

  opt(:bowtie_idx, "The bowtie2 index to screen reads against " +
                   "(can provide more than one)",
      type: :strings)

end

check_files *opts[:forward]
check_files *opts[:reverse]

abort_unless opts[:forward].count == opts[:reverse].count,
             "You need to provide equal numbers of forward " +
             "and reverse files"

if opts[:idba]
  fq2fa = `which fq2fa`.chomp
  if fq2fa.empty?
    AbortIf.logger.warn { "Cannot locate fq2fa program. No IDBA " +
                          "files will be made" }
  end
end

if opts[:idba] && !fq2fa.empty?
  idba_dir = File.join opts[:outdir], "for_idba"
  FileUtils.mkdir_p idba_dir

  all_paired_1 = File.join idba_dir,
                           "all_reads.1.fq"
  all_paired_2 = File.join idba_dir,
                           "all_reads.2.fq"
  all_unpaired = File.join idba_dir,
                           "all_reads.U.fq"

  all_paired_interleaved_fa = File.join idba_dir,
                                        "all_reads.1_and_2_interleaved.fa"
  all_unpaired_fa = File.join idba_dir,
                              "all_reads.U.fa"
end

opts[:forward].each_with_index do |for_f, idx|
  rev_f = opts[:reverse][idx]

  AbortIf.logger.info { "Working on lib ##{idx} " +
                        "(#{for_f} and #{rev_f})" }

  basename = File.basename for_f

  outd = File.join opts[:outdir], basename

  cmd = "ruby qc_single_lib.rb " +
        "-f #{for_f} " +
        "-r #{rev_f} " +
        "-t #{opts[:threads]} " +
        "-o #{outd} " +
        "-b #{opts[:bowtie_idx].join(" ")}"

  Process.run_it! cmd

  if opts[:idba] && !fq2fa.empty?
    AbortIf.logger.info { "Adding to the IDBA tmp files" }

    one = File.join outd, "reads.1.fq.gz"
    two = File.join outd, "reads.2.fq.gz"
    unp = File.join outd, "reads.U.fq.gz"

    cmd = "cat #{one} >> #{all_paired_1}"
    Process.run_it! cmd

    cmd = "cat #{two} >> #{all_paired_2}"
    Process.run_it! cmd

    cmd = "cat #{unp} >> #{all_unpaired}"
    Process.run_it! cmd
  end
end

if opts[:idba] && !fq2fa.empty?
  AbortIf.logger.info { "Making IDBA files" }
  Process.run_it! "#{fq2fa} " +
                  "--filter " +
                  "--merge " +
                  "#{all_paired_1} " +
                  "#{all_paired_2} " +
                  "#{all_paired_interleaved_fa}"

  Process.run_it! "#{fq2fa} " +
                  "--filter " +
                  "#{all_unpaired} " +
                  "#{all_unpaired_fa}"

  Process.run_it "rm #{all_paired_1} " +
                 "#{all_paired_2} " +
                 "#{all_unpaired}"

  pigz = `which pigz`.chomp
  if pigz.empty?
    gzip = `which gzip`.chomp

    if gzip.empty?
      AbortIf.logger.warn { "Cannot locate pigz or gzip. IDBA " +
                            "files will not be compressed." }
    else
      Process.run_it "#{gzip} " +
                     "#{all_paired_interleaved_fa} " +
                     "#{all_unpaired_fa}"
    end
  else
    Process.run_it "#{pigz} --best --processes #{opts[:threads]} " +
                   "#{all_paired_interleaved_fa} " +
                   "#{all_unpaired_fa}"
  end
end

AbortIf.logger.info { "qc.rb wrapper has finished" }
AbortIf.logger.info { "Final outdir is #{opts[:outdir]}" }

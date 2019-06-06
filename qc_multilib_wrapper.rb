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



VERSION = "
    Version: #{QC::Version::VERSION}
    Copyright: 2015 - 2019 Ryan Moore
    Contact: moorer@udel.edu
    Website: https://github.com/mooreryan/qc
    License: GPLv3

"

opts = Optimist.options do
  version QC::Version::VERSION

  banner <<-EOS
#{VERSION}

  When you don't want your libraries combined into one big set of
  reads, use this wrapper instead of the qc.rb script.

  You can also use this if you only have one library as well.

  Options:
  EOS

  opt(:forward, "forward", type: :strings)
  opt(:reverse, "reverse", type: :strings)

  opt(:threads, "Threads", type: :integer, default: 10)

  opt(:outdir,
      "Output directory",
      type: :string,
      default: "qc")
  opt(:gzip_output,
      "Compress output or not.",
      default: true)

  opt(:headcrop,
      "Remove N number of bases from the start of reads",
      short: "c",
      default: 0)
  opt(:bowtie_idx,
      "The bowtie2 index to screen reads against " +
      "(can provide more than one)",
      type: :strings)
end

check_files *opts[:forward]
check_files *opts[:reverse]

p opts

abort_unless opts[:forward].count == opts[:reverse].count,
             "You need to provide equal numbers of forward " +
             "and reverse files"

qc = File.join File.dirname(__FILE__), "qc.rb"

opts[:forward].each_with_index do |for_f, idx|
  rev_f = opts[:reverse][idx]

  AbortIf.logger.info { "Working on lib ##{idx} " +
                        "(#{for_f} and #{rev_f})" }

  basename = File.basename for_f

  outd = File.join opts[:outdir], basename

  if opts[:bowtie_idx]
    cmd = "ruby #{qc} " +
          "-f #{for_f} " +
          "-r #{rev_f} " +
          "-t #{opts[:threads]} " +
          "-g #{opts[:gzip_output]} " +
          "-c #{opts[:headcrop]} " +
          "-o #{outd} " +
          "-b #{opts[:bowtie_idx].join(" ")}"
  else
    cmd = "ruby #{qc} " +
          "-f #{for_f} " +
          "-r #{rev_f} " +
          "-t #{opts[:threads]} " +
          "-g #{opts[:gzip_output]} " +
          "-c #{opts[:headcrop]} " +
          "-o #{outd}"
  end

  Process.run_it! cmd
end

AbortIf.logger.info { "qc.rb wrapper has finished" }
AbortIf.logger.info { "Final outdir is #{opts[:outdir]}" }

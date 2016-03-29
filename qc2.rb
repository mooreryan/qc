#!/usr/bin/env ruby

# TODO globbing broken

# TODO non-gzipped in files broken

# TODO need to run in one directory back of where the infiles are
# located

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods

Ryan.req *%w[parse_fasta fileutils]

opts = Trollop.options do
  banner <<-EOS

  Naming convention:
    folder/base_name.crazy_illumina_sample_name.{fq,fastq}.{,gz}

  --forward and --reverse can take comma separated lists of file names

  -f apple.r1.fq,pie.r1.fq -r apple.r2.fq,pie.r2.fq

  Options:
  EOS

  opt(:forward, "Input file", type: :string)
  opt(:reverse, "Input file", type: :string)
  opt(:threads, "Num threads", type: :integer, default: 3)
  opt(:outdir, "Output directory", type: :string, default: ".")
end

Ryan.check_file opts[:forward], :forward
Ryan.check_file opts[:reverse], :reverse
Ryan.try_mkdir(opts[:outdir])

TRIMMO = File.join File.dirname(__FILE__),
                   "bin",
                   "trimmomatic-0.35",
                   "trimmomatic-0.35.jar"

FLASH = File.join File.dirname(__FILE__),
                  "bin",
                  "flash"

zip = `which pigz`.chomp
unless $?.exitstatus.zero?
  zip = `which gzip`.chomp

  unless $?.exitstatus.zero?
    abort "[ERROR] Need gzip or pigz"
  end
end

def get_baseout fname
  "#{fname.split(".").first}.adapters_removed.fq"
end

def get_library fname
  the_basename = File.basename(fname)

  the_basename.split(".").first
end

trimseqs = File.join File.dirname(__FILE__),
                     "bin",
                     "trimmomatic-0.35",
                     "adapters",
                     "TruSeq3-PE-both.fa"

seed_mismatches = 2
palindrome_clip_threshold = 30
simple_clip_threshold = 10

WINDOW_SIZE = 4
QUAL = 15
THREADS = opts[:threads]

MIN_LEN = 50

MAX_OVERLAP = 200

def trim_se fname, dir
  # quality trim
  log = File.join dir, "log.trimmomatic.sliding_window.on_#{fname.gsub("/", "_")}.txt"
  cmd = "java -jar #{TRIMMO} SE " +
        "-threads #{THREADS} " +
        "#{fname} " +
        "#{fname.sub("adapters_removed.flashed", "adapters_removed.flashed.trimmed")} " +
        "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
        "MINLEN:#{MIN_LEN} " +
        ">> #{log} 2>&1"
  Ryan.run_it cmd
end

def trim_pe forward, reverse, dir
  # quality trim
  log = File.join dir, "log.trimmomatic.sliding_window.on_#{forward.gsub("/", "_")}.txt"
  cmd = "java -jar #{TRIMMO} PE " +
        "-threads #{THREADS} " +
        "#{forward} #{reverse} " +
        "-baseout #{forward.sub("adapters_removed.flashed.paired_1", "adapters_removed.flashed.trimmed.paired")} " +
        "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
        "MINLEN:#{MIN_LEN} " +
        ">> #{log} 2>&1"
  Ryan.run_it cmd
end



file_pairs = [forward_files, reverse_files].transpose
# TODO ERROR globbing wont work after the first set of file pairs,
# need to be specific about files or at least use the current basename
file_pairs.each_with_index do |(forw, rev), idx|
  # get basename
  baseout = get_baseout forw

  # TODO base has the dir as well
  base = baseout.sub(".fq", "")
  dir = File.dirname(forw)

  forw_lib_name = get_library forw
  rev_lib_name = get_library rev

  if forw_lib_name == rev_lib_name
    lib_name = forw_lib_name
  else
    abort "ERROR: lib name mismatch: #{forw_lib_name}, #{rev_lib_name}"
  end

  qc_d = File.join dir, "qc"
  Ryan.try_mkdir qc_d

  intermediate_d = File.join qc_d, "intermediate_files"
  Ryan.try_mkdir intermediate_d

  final_d = File.join qc_d, "final"
  Ryan.try_mkdir final_d

  unless  baseout == get_baseout(rev)
    abort "ERROR: basenames incorrect"
  end

  # trim adapters
  log = File.join dir, "log.trimmomatic.trim_adapters.txt"
  cmd = "java -jar #{TRIMMO} PE " +
        "-threads #{opts[:threads]} " +
        "-baseout #{baseout} " +
        "#{forw} " +
        "#{rev} " +
        "ILLUMINACLIP:" +
        "#{trimseqs}:" +
        "#{seed_mismatches}:" +
        "#{palindrome_clip_threshold}:" +
        "#{simple_clip_threshold} " +
        ">> #{log} 2>&1"

  Ryan.run_it cmd

  # flash reads that are still paired
  log = File.join dir, "log.flash.txt"
  cmd = "#{FLASH} " +
        "--threads #{opts[:threads]} " +
        "--output-prefix flashed " +
        "--max-overlap #{MAX_OVERLAP} " +
        "#{base}_1P.fq #{base}_2P.fq " +
        ">> #{log} 2>&1"

  Ryan.run_it cmd

  # clean up flashed stuff
  output = File.join Dir.pwd, "flashed.*"
  cmd = "mv #{output} #{dir}"
  Ryan.run_it cmd

  # combine the trimmed and flashed reads
  ext_frags = File.join dir, "flashed.extendedFrags.fastq"
  unpaired_forw = "#{base}_1U.fq"
  unpaired_rev = "#{base}_2U.fq"

  # TODO base contains the trimmed part
  unpaired_reads = "#{base}.flashed.unpaired.fq"
  cmd = "cat #{ext_frags} #{unpaired_forw} #{unpaired_rev} " +
        "> #{unpaired_reads}"
  Ryan.run_it cmd

  new_forw = "#{base}.flashed.paired_1.fq"
  new_rev = "#{base}.flashed.paired_2.fq"
  cmd = "mv #{File.join(dir, 'flashed.notCombined_1.fastq')} " +
        "#{new_forw}"
  Ryan.run_it cmd

  cmd = "mv #{File.join(dir, 'flashed.notCombined_2.fastq')} " +
        "#{new_rev}"
  Ryan.run_it cmd

  trim_pe new_forw, new_rev, dir
  trim_se unpaired_reads, dir

  # clean up
  cmd = "mv #{File.join(dir, 'flashed.*')} #{intermediate_d}"
  Ryan.run_it cmd

  cmd = "mv #{File.join(dir, "*.adapters_removed.flashed.trimmed.*")} #{final_d}"
  Ryan.run_it cmd

  cmd = "mv #{File.join(dir, '*_{1,2}{U,P}.fq')} #{intermediate_d}"
  Ryan.run_it cmd

  cmd = "mv #{File.join(dir, "log.*")} #{intermediate_d}"
  Ryan.run_it cmd

  cmd = "mv #{File.join(dir, "*.adapters_removed.flashed.*")} #{intermediate_d}"
  Ryan.run_it cmd

  # add in the unpaired reads causing by quality trimming
  cmd = "cat #{File.join(final_d, '*_1U.*')} #{File.join(final_d, '*_2U.*')} #{unpaired_reads.sub('unpaired', 'trimmed.unpaired').sub(dir, final_d)} > #{unpaired_reads.sub('unpaired', 'trimmed.all.unpaired').sub(dir, final_d)}"
  Ryan.run_it cmd

  # remove fq intermediate files
  # remove_these = File.join intermediate_d, "*.{fq,fastq}"
  # cmd = "rm #{remove_these}"
  # Ryan.run_it cmd

  # zip everything
  # cmd = "#{zip} --recursive #{qc_d}/*"
  # Ryan.run_it cmd
end

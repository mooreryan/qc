# very specific

#!/usr/bin/env ruby

def trim_se! infile:, outfile:, log:
  count = linecount infile
  if count >= 4
    cmd = "java -jar #{TRIMMO} SE " +
          "-threads #{THREADS} " +
          "#{infile} " +
          "#{outfile} " +
          "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
          "MINLEN:#{MIN_LEN} " +
          ">> #{log} 2>&1"

    Ryan.run_it! cmd
  else
    warn "WARN: no reads in #{infile}"
  end

  Ryan.run_it! "rm #{infile}"
end

def trim_pe! in1:, in2:, baseout:, log:
  cmd = "java -jar #{TRIMMO} PE " +
        "-threads #{THREADS} " +
        "#{in1} " +
        "#{in2} " +
        "-baseout #{baseout} " +
        "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
        "MINLEN:#{MIN_LEN} " +
        ">> #{log} 2>&1"

  Ryan.run_it! cmd
  Ryan.run_it! "rm #{in1} #{in2}"
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

  Ryan.run_it! cmd
  Ryan.run_it! "rm #{in1} #{in2}"
end

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods

Ryan.req *%w[parse_fasta]

opts = Trollop.options do
  banner <<-EOS

  Info

  Options:
  EOS

  opt(:six_for, "Six forward", type: :string, long: "s1", short: "s")
  opt(:six_rev, "Six reverse", type: :string, long: "s2", short: "i")

  opt(:two_for, "Two forward", type: :string, long: "t1", short: "r")
  opt(:two_rev, "Two reverse", type: :string, long: "t2", short: "w")

  opt(:threads, "Threads", type: :integer, default: 10, short: "t")

  opt(:outdir, "Output directory", type: :string,
      default: "qc3_with_flash")
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

big_log = File.join opts[:outdir], "log.big.txt"

interleave_reads =
  "/home/moorer/vendor/khmer/scripts/" +
  "interleave-reads.py"

trim_adapters_log =
  File.join opts[:outdir],
            "log.trimmomatic.trim_adapters." +
            "on_six_and_two.txt"

qual_trim_log =
  File.join opts[:outdir],
            "log.trimmomatic.sliding_window.on_" +
            "six_and_two.adapters_trimmed_paired.txt"

qual_trim_se_log =
  File.join opts[:outdir],
            "log.trimmomatic.sliding_window.on_" +
            "six_and_two.adapters_trimmed_se.txt"

trim_adapters_baseout =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed.fq"

quality_trim_baseout =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed.flashed_not_extended." +
            "quality_trimmed.fq"

quality_trim_baseout_extended_frags =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed.flashed_extended." +
            "quality_trimmed.fq"

both_adapters_trimmed_for_paired =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_1P.fq"

both_adapters_trimmed_rev_paired =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_2P.fq"

flash_log =
  File.join opts[:outdir],
            "log.flash_on." +
            "both_adapters_trimmed.txt"

both_adapters_trimmed_for_unpaired =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_1U.fq"

both_adapters_trimmed_for_unpaired_qual_trimmed =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_1U.quality_trimmed.fq"

both_adapters_trimmed_rev_unpaired =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_2U.fq"

both_adapters_trimmed_rev_unpaired_qual_trimmed =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_2U.quality_trimmed.fq"

# flash output
flashed_unpaired_extended_frags =
  File.join opts[:outdir],
            "flashed.extendedFrags.fastq"

flashed_notcombined_paired_1 =
  File.join opts[:outdir],
            "flashed.notCombined_1.fastq"

flashed_notcombined_paired_2 =
  File.join opts[:outdir],
            "flashed.notCombined_2.fastq"

flashed_notcombined_quality_trimmed_baseout =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed.flashed_not_extended." +
            "quality_trimmed.fq"

flashed_hist =
  File.join opts[:outdir],
            "flashed.hist"

flashed_histogram =
  File.join opts[:outdir],
            "flashed.histogram"

adapter_trimmed_and_flashed_1 =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed.flashed_not_combined.1.fq"
adapter_trimmed_and_flashed_2 =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed.flashed_not_combined.2.fq"

all_unpaired = File.join opts[:outdir],
                         "six_and_two.all_unpaired.fq"
glob = File.join opts[:outdir], "*_{1,2}U.fq"

interleaved_reads =
  File.join opts[:outdir],
            "six_and_two.paired_interleaved.fq"

final_for =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed." +
            "quality_trimmed_1P.fq"

final_rev =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed." +
            "quality_trimmed_2P.fq"

six_for = Ryan.check_file(opts[:six_for], :six_for)
six_rev = Ryan.check_file(opts[:six_rev], :six_rev)

two_for = Ryan.check_file(opts[:two_for], :two_for)
two_rev = Ryan.check_file(opts[:two_rev], :two_rev)

Ryan.try_mkdir(opts[:outdir])

both_for = File.join opts[:outdir], "six_and_two.1.fq"
both_rev = File.join opts[:outdir], "six_and_two.2.fq"

if opts[:six_for].match(/.gz$/)
  Ryan.run_it "gunzip #{opts[:six_for]}"
  opts[:six_for] = opts[:six_for].sub(/.gz$/, "")
end

if opts[:six_rev].match(/.gz$/)
  Ryan.run_it "gunzip #{opts[:six_rev]}"
  opts[:six_rev] = opts[:six_rev].sub(/.gz$/, "")
end

if opts[:two_for].match(/.gz$/)
  Ryan.run_it "gunzip #{opts[:two_for]}"
  opts[:two_for] = opts[:two_for].sub(/.gz$/, "")
end

if opts[:two_rev].match(/.gz$/)
  Ryan.run_it "gunzip #{opts[:two_rev]}"
  opts[:two_rev] = opts[:two_rev].sub(/.gz$/, "")
end

Ryan.run_it! "cat #{opts[:six_for]} #{opts[:two_for]} > #{both_for}"
Ryan.run_it! "cat #{opts[:six_rev]} #{opts[:two_rev]} > #{both_rev}"

def linecount fname
  `wc -l #{fname}`.split(" ").first.to_i
end

six_for_count = linecount opts[:six_for]
two_for_count = linecount opts[:two_for]
both_for_count = linecount both_for

two_rev_count = linecount opts[:two_rev]
six_rev_count = linecount opts[:six_rev]
both_rev_count = linecount both_rev

warn %w[lib which count].join "\t"
warn ["six", "forward", six_for_count].join "\t"
warn ["two", "forward", two_for_count].join "\t"
warn ["both", "forward", both_for_count].join "\t"
warn ["six", "reverse", six_rev_count].join "\t"
warn ["two", "reverse", two_rev_count].join "\t"
warn ["both", "reverse", both_rev_count].join "\t"
warn ""

if six_for_count + two_for_count == both_for_count
  warn "Forward counts match"
else
  abort "ERROR: line count problem for six forward and two forward"
end

if six_rev_count + two_rev_count == both_rev_count
  warn "Reverse counts match"
else
  abort "ERROR: line count problem for six reverse and two reverse"
end

if both_for_count == both_rev_count
  warn "Combined forward and reverse counts match"
else
  abort "ERROR: different number of forward and reverse reads"
end

# now the input files should all be NOT gzipped
Ryan.run_it! "pigz #{opts[:six_for]} #{opts[:six_rev]} #{opts[:two_for]} #{opts[:two_rev]}"

# Add back the gz extension
opts[:six_for] += ".gz"
opts[:six_rev] += ".gz"
opts[:two_for] += ".gz"
opts[:two_rev] += ".gz"

# remove the previous log of the same name if it exists
if File.exists? trim_adapters_log
  Ryan.run_it! "rm #{trim_adapters_log}"
end

if File.exists? qual_trim_log
  Ryan.run_it! "rm #{qual_trim_log}"
end

# Trim the adapters
cmd = "java -jar #{TRIMMO} PE " +
      "-threads #{opts[:threads]} " +
      "-baseout #{trim_adapters_baseout} " +
      "#{both_for} " +
      "#{both_rev} " +
      "ILLUMINACLIP:" +
      "#{TRIMSEQS}:" +
      "#{SEED_MISMATCHES}:" +
      "#{PALINDROME_CLIP_THRESHOLD}:" +
      "#{SIMPLE_CLIP_THRESHOLD} " +
      ">> #{big_log} 2>&1"

Ryan.run_it! cmd
# remove files just processed
Ryan.run_it! "rm #{both_for} #{both_rev}"

exit

# trim_pe! in1: both_adapters_trimmed_for_paired,
#          in2: both_adapters_trimmed_rev_paired

# # flash reads that are still paired
# flash! in1: both_adapters_trimmed_for_paired,
#        in2: both_adapters_trimmed_rev_paired,
#        outdir: opts[:outdir],
#        log: big_log

# trim_se! infile: flashed_unpaired_extended_frags,
#          outfile: quality_trim_baseout_extended_frags,
#          log: big_log

# trim_pe! in1: flashed_notcombined_paired_1,
#          in2: flashed_notcombined_paired_2,
#          baseout: flashed_notcombined_quality_trimmed_baseout,
#          log: big_log

# trim_se! infile: both_adapters_trimmed_for_unpaired,
#          outfile: both_adapters_trimmed_for_unpaired_qual_trimmed,
#          log: big_log

# trim_se! infile: both_adapters_trimmed_rev_unpaired,
#          outfile: both_adapters_trimmed_rev_unpaired_qual_trimmed,
#          log: big_log

# quality trim the paired reads
cmd = "java -jar #{TRIMMO} PE " +
      "-threads #{THREADS} " +
      "#{adapter_trimmed_and_flashed_1} " +
      "#{adapter_trimmed_and_flashed_2} " +
      "-baseout #{quality_trim_baseout} " +
      "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
      "MINLEN:#{MIN_LEN} " +
      ">> #{qual_trim_log} 2>&1"

Ryan.run_it! cmd
Ryan.run_it! "rm #{adapter_trimmed_and_flashed_1} " +
             "#{adapter_trimmed_and_flashed_2}"

# quality trim the single reads
cmd = "java -jar #{TRIMMO} SE " +
      "-threads #{THREADS} " +
      "#{flashed_unpaired_extended_frags} " +
      "#{quality_trim_baseout_extended_frags} " +
      "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
      "MINLEN:#{MIN_LEN} " +
      ">> #{qual_trim_se_log} 2>&1"

Ryan.run_it! cmd
Ryan.run_it! "rm #{flashed_unpaired_extended_frags}"

exit

Ryan.run_it! "cat #{glob} > #{all_unpaired}"
Ryan.run_it! "rm #{glob}"
Ryan.run_it! "rm #{both_adapters_trimmed_for_paired}"
Ryan.run_it! "rm #{both_adapters_trimmed_rev_paired}"
Ryan.run_it! "rm #{both_for} #{both_rev}"

cmd = "grep -h Surviving " +
      "#{File.join(opts[:outdir], "log.*")}"
stats = `#{cmd}`
warn stats

Ryan.run_it! "#{interleave_reads} -o " +
            "#{interleaved_reads} " +
            "#{final_for} #{final_rev}"

Ryan.run_it! "rm #{final_for} #{final_rev}"

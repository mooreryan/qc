# very specific

#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path("~"), "lib", "ruby", "ryan.rb")
require_relative methods

Ryan.req *%w[parse_fasta]

opts = Trollop.options do
  banner <<-EOS

  Info

  Options:
  EOS

  opt(:six_for, "Input file", type: :string)
  opt(:six_rev, "Input file", type: :string)

  opt(:two_for, "Input file", type: :string)
  opt(:two_rev, "Input file", type: :string)

  opt(:threads, "Threads", type: :integer, default: 10)

  opt(:outdir, "Output directory", type: :string,
      default: ".")
end

TRIMMO = File.join File.dirname(__FILE__),
                   "bin",
                   "trimmomatic-0.35",
                   "trimmomatic-0.35.jar"

SEED_MISMATCHES = 2
PALINDROME_CLIP_THRESHOLD = 30
SIMPLE_CLIP_THRESHOLD = 10

WINDOW_SIZE = 4
QUAL = 15
THREADS = opts[:threads]

MIN_LEN = 50

MAX_OVERLAP = 200

TRIMSEQS = File.join File.dirname(__FILE__),
                     "bin",
                     "trimmomatic-0.35",
                     "adapters",
                     "TruSeq3-PE-both.fa"

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

trim_adapters_baseout =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed.fq"

quality_trim_baseout =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed." +
            "quality_trimmed.fq"

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

# Ryan.run_it! "pigz #{opts[:six_for]} #{opts[:six_rev]} #{opts[:two_for]} #{opts[:two_rev]}"

Ryan.run_it! "rm #{trim_adapters_log}"

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
      ">> #{trim_adapters_log} 2>&1"

Ryan.run_it! cmd

both_adpaters_trimmed_for_paired =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_1P.fq"

both_adpaters_trimmed_rev_paired =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_2P.fq"

both_adpaters_trimmed_for_unpaired =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_1U.fq"

both_adpaters_trimmed_rev_unpaired =
  File.join opts[:outdir],
            "six_and_two.adapters_trimmed_2U.fq"

cmd = "java -jar #{TRIMMO} PE " +
      "-threads #{THREADS} " +
      "#{both_adpaters_trimmed_for_paired} " +
      "#{both_adpaters_trimmed_rev_paired} " +
      "-baseout #{quality_trim_baseout} " +
      "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
      "MINLEN:#{MIN_LEN} " +
      ">> #{qual_trim_log} 2>&1"

Ryan.run_it! cmd

all_unpaired = File.join opts[:outdir],
                         "six_and_two.all_unpaired.fq"
glob = File.join opts[:outdir], "*_{1,2}U.fq"

Ryan.run_it! "cat #{glob} > #{all_unpaired}"
Ryan.run_it! "rm #{glob}"
Ryan.run_it! "rm #{both_adpaters_trimmed_for_paired}"
Ryan.run_it! "rm #{both_adpaters_trimmed_rev_paired}"
Ryan.run_it! "rm #{both_for} #{both_rev}"

cmd = "grep -h Surviving " +
      "#{File.join(opts[:outdir], "log.*")}"
stats = `#{cmd}`
warn stats

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

Ryan.run_it! "#{interleave_reads} -o " +
            "#{interleaved_reads} " +
            "#{final_for} #{final_rev}"

Ryan.run_it! "rm #{final_for} #{final_rev}"

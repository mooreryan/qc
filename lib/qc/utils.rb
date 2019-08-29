require "securerandom"

module QC
  module Utils
    def cat_fastq_files outf, *fnames
      File.open(outf, "w") do |f|
        fnames.each do |fname|
          AbortIf.logger.debug { "Writing #{fname} to #{outf}" }

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

    # Tries to remove some common sequence file extensions from the
    # given string.
    #
    # Any of would output just 'apple'....
    # apple.1.fq.gz
    # apple.1.fq
    # apple.2.fq.gz
    # apple.2.fq
    # apple.1.fastq.gz
    # apple.1.fastq
    # apple.2.fastq.gz
    # apple.2.fastq
    # apple.1.fq.bz2
    # apple.1.fq
    # apple.2.fq.bz2
    # apple.2.fq
    # apple.1.fastq.bz2
    # apple.1.fastq
    # apple.2.fastq.bz2
    # apple.2.fastQ
    # apple.r.fq
    # apple.forward.fq.bz2
    # apple.for.fastQ.gz
    def remove_fq_ext str
      str.sub(/\.(?:[12]|[fr]|forward|reverse|for|rev)\.(?:fq|fastq)\.*(?:gz|bz2)*/i, "")
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
              "-phred33 " +
              "#{inf} " +
              "#{out} " +
              "HEADCROP:#{HEADCROP} " +
              "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
              "MINLEN:#{MIN_LEN} " +
              ">> #{log} 2>&1"

        Process.run_it! cmd
      else
        AbortIf.logger.warn { "No reads in #{inf}. Not running " +
                              "qual_trim_se!()" }

        # Make a empty placeholder file
        Process.run_it! "touch #{out}"

        AbortIf.logger.info { "Made fake file #{out}" }
      end

      Process.run_it "rm #{inf}"

      out
    end

    def qual_trim_pe! in1:, in2:, baseout:, log:
                         cmd = "java -jar #{TRIMMO} PE " +
                               "-threads #{THREADS} " +
                               "-phred33 " +
                               "#{in1} " +
                               "#{in2} " +
                               "-baseout #{baseout} " +
                               "HEADCROP:#{HEADCROP} " +
                               "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
                               "MINLEN:#{MIN_LEN} " +
                               ">> #{log} 2>&1"

      Process.run_it! cmd
      Process.run_it "rm #{in1} #{in2}"
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
      Process.run_it "rm #{in1} #{in2}"
      Process.run_it!("mv #{outdir}/flashed.extendedFrags.fastq " +
                      "#{outdir}/../#{BASENAME}.adapter_trimmed.flash_combined")
      Process.run_it!("mv #{outdir}/flashed.notCombined_1.fastq " +
                      "#{outdir}/../#{BASENAME}.adapter_trimmed.flash_notcombined_1P")
      Process.run_it!("mv #{outdir}/flashed.notCombined_2.fastq " +
                      "#{outdir}/../#{BASENAME}.adapter_trimmed.flash_notcombined_2P")
    end

    def adapter_trim!(in1:, in2:, baseout:, log:)
      # Trim the adapters
      cmd = "java -jar #{TRIMMO} PE " +
            "-threads #{THREADS} " +
            "-baseout #{baseout} " +
            "-phred33 " +
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

    # Treat all reads as unpaired so that --un gives all reads that
    # did not align to the genome.
    def screen!(index:, reads:, log:)
      # use bowtie2 to screen reads against against a genome
      good_reads = reads + ".did_not_align.fq"
      bad_reads  = reads + ".did_align.fq"

      cmd = "#{BOWTIE} -x #{index} " +
            "-U #{reads} " +
            "--sensitive --end-to-end " +
            "--threads #{THREADS} " +
            "--un #{good_reads} " +
            "--al #{bad_reads} " +
            "-S /dev/null " +
            ">> #{log} 2>&1"

      Process.run_it! cmd

      [good_reads, bad_reads]
    end

    def fix_pairs!(in1:, in2:, outdir:, log:)
      basename = File.join outdir,
                           "qc_#{SecureRandom.hex(20)}"

      cmd = "#{FIX_PAIRS} #{in1} #{in2} #{basename} " +
            ">> #{log} 2>&1"

      Process.run_it! cmd

      out1 = basename + ".1.fq"
      out2 = basename + ".2.fq"
      outU = basename + ".U.fq"

      [out1, out2, outU]
    end
  end
end

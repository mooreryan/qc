def trim_se! infile, outfile, log
  cmd = "java -jar #{TRIMMO} SE " +
        "-threads #{THREADS} " +
        "#{infile} " +
        "#{outfile} " +
        "SLIDINGWINDOW:#{WINDOW_SIZE}:#{QUAL} " +
        "MINLEN:#{MIN_LEN} " +
        ">> #{log} 2>&1"

  Ryan.run_it! cmd
  Ryan.run_it! "rm #{infile}"
end

def trim_pe in1, in2, baseout, log
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

def flash in1, in2, outdir, log
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

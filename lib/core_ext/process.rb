module CoreExt
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
                   "ERROR: non-zero exit status (#{exit_status}) " +
                   "while running cmd(s) '#{a.join(" and ")}'"

      exit_status
    end
  end
end

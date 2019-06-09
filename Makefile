test_d = test_files
test_outd = $(test_d)/test_output

bowtie_index_d = bowtie_indices/bos_taurus/2016_06_11
bowtie_index = $(bowtie_index_d)/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna

forward_reads = $(test_d)/*.1.fq.gz
reverse_reads = $(test_d)/*.2.fq.gz

threads = 2

test_qc:
	rm -r $(test_outd); ./qc.rb -f $(forward_reads) -r $(reverse_reads) -t $(threads) -o $(test_outd) --gzip-output

test_qc_bowtie:
	rm -r $(test_outd); ./qc.rb -f $(forward_reads) -r $(reverse_reads) -t $(threads) -o $(test_outd) --gzip-output -b $(bowtie_index)

test_qc_multilib:
	rm -r $(test_outd); ./qc_multilib_wrapper.rb -f $(forward_reads) -r $(reverse_reads) -t $(threads) -o $(test_outd) --gzip-output

test_qc_multilib_bowtie:
	rm -r $(test_outd); ./qc_multilib_wrapper.rb -f $(forward_reads) -r $(reverse_reads) -t $(threads) -o $(test_outd) --gzip-output -b $(bowtie_index)


test: test_qc test_qc_bowtie test_qc_multilib test_qc_multilib_bowtie

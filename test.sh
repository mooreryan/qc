#!/bin/bash

rm -r test_files/test/qc4_no_flash; pigz test_files/test/* && time ruby qc4_no_flash.rb -s test_files/test/six.r1.fq.gz -i test_files/test/six.r2.fq.gz -r test_files/test/two.r1.fq.gz -w test_files/test/two.r2.fq.gz -t 2 -o test_files/test/qc4_no_flash && tree test_files/test/qc4_no_flash

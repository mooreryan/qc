#!/bin/bash

ruby qc.rb -f test_files/test/six.r1.fq.gz,test_files/test/two.r1.fq.gz,test_files/test2/apple.r1.fq.gz,test_files/test2/pie.r1.fq.gz -r test_files/test/six.r2.fq.gz,test_files/test/two.r2.fq.gz,test_files/test2/apple.r2.fq.gz,test_files/test2/pie.r2.fq.gz

tree test_files

#!/bin/sh
#$ -cwd
#$ -S /bin/sh
python Run_FusionScan.py test_data/test.fa test 75 00.00.00.00:8100 ---phred33 -P 1 -ms 2

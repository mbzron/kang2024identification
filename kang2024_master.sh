#!/bin/bash

SECONDS=0

Rscript kang2024_testing_1.r
Rscript kang2024_testing_2.r
Rscript kang2024_testing_3.r
Rscript compute_Aj.r
Rscript kang2024_testing_post_Aj.r
Rscript kang2024_testing_4.r

# display elapsed time (minutes)
elapsed_time=$(($SECONDS / 60))
echo "Elapsed time: $elapsed_time minutes"
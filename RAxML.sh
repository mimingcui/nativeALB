#!/bin/bash

module load RAxML/8.2.9
raxmlHPC-PTHREADS-SSE3 \
-n ALB.output \
-s 359.txt \
-m GTRGAMMA \
-f a \
-p 164655 \
-x 12345 \
-# 1000 \
-T 8
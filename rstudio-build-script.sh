#!/bin/bash
R CMD INSTALL --no-multiarch --with-keep.source ../survTMB
R -e "devtools::test()"

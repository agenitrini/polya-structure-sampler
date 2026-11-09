#!/bin/sh
# Copyright (C) 2026 Antoine Genitrini -- SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE)
# bench.sh -- measurements of time and of the mean number of examined pairs
# (steps / (n ln n)) for the four families. Output: a Markdown table.
# FAMILIES variable (optional): subset of the families, e.g.
#   FAMILIES=binary ./bench.sh
set -e
cd "$(dirname "$0")"
[ -x ./polyagen ] || make
get() { sed -n "s/.* $1=\([^ s]*\).*/\1/p"; }
printf '| family | n | samples | mean time / structure (s) | steps/(n ln n) | precomputation (s) |\n'
printf '|---|---:|---:|---:|---:|---:|\n'
for fam in ${FAMILIES:-polya hier identity binary}; do
  for spec in 1000:2000 2000:1000 10000:500 100000:100 1000000:20 10000000:3; do
    n=${spec%%:*}; c=${spec##*:}
    line=$(./polyagen --family $fam -n $n --count $c --seed 2026 --format none --stats 2>&1 | grep '^# MEAN family')
    t=$(echo "$line" | get time); r=$(echo "$line" | sed -n 's/.*steps\/(n ln n)=\([0-9.]*\).*/\1/p'); p=$(echo "$line" | get precomputation)
    printf '| %s | %s | %s | %s | %s | %s |\n' $fam $n $c $t $r $p
  done
done

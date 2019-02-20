awk  '{ print $2 }' simulMDmk3.out | uniq -c  | awk '{ if ($1 > 2) print $0}' | more

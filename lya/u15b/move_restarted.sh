#!/bin/bash
# After a MontePython restart (-r), move the old chain files aside (not deleted) once every one of them has been copied into
# a new chain that has already accepted new points, so that no point is counted twice by 'info' or analyze.py.
#   move_restarted.sh <chain folder> <number of chains restarted>
# (The first version compared every new chain with the longest old one and could give up before moving anything:
#  ESPEC_U15b, Nota 3.  The check is now by content, in restarted.py.)
d="$1"; n="$2"; t=0
here="$(cd "$(dirname "$0")" && pwd)"
while [ $t -lt 86400 ]; do
  k=$(python "$here/restarted.py" "$d" | head -1 | awk '{print $2}')
  if [ "${k:-0}" -ge "$n" ]; then python "$here/restarted.py" "$d" --move; exit 0; fi
  sleep 60; t=$((t + 60))
done
echo "WARNING: after 24 h only ${k:-0} of $n old chains of $d are copied; run: python lya/u15b/restarted.py $d --move"

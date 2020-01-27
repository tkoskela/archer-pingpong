set logscale x
set logscale y
set xlabel "Message Size (Bytes)"
set ylabel "Bandwidth (MB/s)"
plot \
"standard_bandwidth.plot" using 1:2 \
title "Standard send" with linespoint, \
"sync_bandwidth.plot" using 1:2 \
title "Synchronous send" with linespoint

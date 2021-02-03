set out "progvsiter.pdf";
set encoding iso_8859_1;
set term pdfcairo enhanced color;
set ylabel "Progress coordinate";
set xlabel "Values over the iteration"
#set xrange [550:700]
#set xtics 100
plot "test_scatter" using :1 with lines lw 2 notitle


set out "free_energy_182_it.png"
set terminal png font Helvetica 26 size 750, 500 enhanced background '#FF000000'#
set ylabel "Probability density" tc rgb "#FFFFFFFF";
set xlabel "progress coordinates" tc rgb "#FFFFFFFF" ;
set key top right tc rgb "#FFFFFFFF" font "Helvetica,26"

set border lc rgb 'white'
set xtics nomirror
set ytics nomirror

set style line 1 lw 4 lt 2 pt 7 lc rgb "cyan"   
set style line 2 lw 4 lt 2 pt 7 lc rgb "orange"
set style line 3 lw 6 lt 2 pt 7 lc rgb "cyan" 
set style line 4 lw 6 lt 2 pt 7 lc rgb "orange"
set style line 5 lw 4 lt 2 pt 7 lc rgb "green"

plot "test" using 1:2 with lines ls 1 title "my code",\
"../feng_1_popc_2/hist_182.dat"  using 1:2 with lines ls 2 title "WESTPA analysis";

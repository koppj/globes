#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "th13chi2.eps"

# Layout
set xrange [1e-4:1e-2]       # Plot range on horizontal axis
set yrange [0:50]            # Plot range on vertical axis
set logscale x               # Use log scale on horizontal axis

# Legend
set title "{/Symbol c}^2 vs. sin^2 2{/Symbol q}_{13} for the wrong hierarchy solution"
set xlabel "sin^2 2{/Symbol q}_{13}" 0,0
set ylabel "{/Symbol c}^2" 0,0
set key graph 0.15, 0.25 spacing 1.3   # Position and vertical spacing of legend

# Do the actual plotting
plot "th13chi2.dat" with lines lt -1 lw 2

# Show resulting EPS figure
system "gv th13chi2.eps"


cd 'C:\cygwin\home\eross\CSI780\Project3\4'
set terminal png
set size square
set xrange [-0.5:199.5]
set yrange [-0.5:199.5]

do for [ii=1:9999] {
  set output ''.ii.'.png'
  plot ''.ii.'.csv' matrix with image
}

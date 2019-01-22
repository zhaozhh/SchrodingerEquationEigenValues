####!/usr/bin/gnuplot
unset key
unset terminal
set terminal epslatex standalone color  14  dashlength 4
#set title "服务发现率" font "Times New Roman,宋体,20"
set output 'temple.tex'
#set xrange [-10 to 10]
set yrange [-10 to 90]
#set xzeroaxis lt 1 #设置零点的坐标轴，线性为暗线
#set yzeroaxis lt 1
#set key bottom right  reverse vertical Left #box #spacing 2
#set key right top reverse horiz  # samplen 3 调整图线示例长度 at 10,0.7
set key right bottom reverse Left #Left 文字左对齐，left：位置左
#set xtics 1 #设置坐标间隔
#set ytics 0.5 
#set mytics 5
#set size 0.9, 0.9
#set origin 0.0, 0.0
#set lmargin 4
#set rmargin 1
#set tmargin 1
#set bmargin 4
#unset key
#set multiplot layout 2,2
set xlabel '$z$'
set ylabel '$V_{\mathrm{III}}(z)$'
#set xtics out nomirror #是坐标的刻度朝外
#set ytics out nomirror
# 设置线的类型
set colorsequence default # |podo|classic
set style line 1 lt 8 lw 4 #lc rgb "black"
set style line 2 lt 7 lw 4 #lc rgb "red”
set style line 3 lt 4 lw 4 #lc rgb "orange"
set style line 4 lt 5 lw 4 #lc rgb "yellow"
set style line 5 lt 2 lw 4 #lc rgb "green"
set style line 6 lt 6 lw 4 #black blue
set style line 7 lt 3 lw 4 #lc rgb "blue"
set style line 8 lt 1 lw 4 #purple


#设置图中的文字和箭头
#set label '$n=0,\;m^2=0$' at 1.5,0
#set arrow from -1,0 to 1,0 nohead          ls 2
#set label '$n=1,\;m^2=16.38$' at 1.5,16.38 
#set arrow from -1,16.38 to 1,16.38 nohead  ls 3
#set label '$n=2,\;m^2=29.85$' at 1.5,29.85
#set arrow from -1,29.85 to 1,29.85 nohead  ls 5
#set label '$n=3,\;m^2=40.73$' at 1.5,40.73
#set arrow from -1,40.73 to 1,40.73 nohead  ls 6

#set label '$V_\text{III}(z)$' at -8,50
#set arrow from -6,50.73 to -2,60  ls 1

plot 'Vzk=1C2=20C3=0p40p50p6.txt' using 1:2 ls 6 title 'C$_3$=0.4' with lines,\
     'Vzk=1C2=20C3=0p40p50p6.txt' using 1:3 ls 5 title 'C$_3$=0.5 ' with lines,\
     'Vzk=1C2=20C3=0p40p50p6.txt' using 1:4 ls 3 title 'C$_3$=0.6 ' with lines
     
#set title '$m^2=0.000000$'
#plot 'waveFunctions_rho-1-lambda-0.1.dat'  using 1:2 ls 1 with lines, \
#     'waveFunctions_rho-1-lambda-0.1.dat'  using (-$1):2 lw 3 lt 1 with lines
#unset title
#set title '$m^2=3.972136$'
#plot 'waveFunctions_rho-1-lambda-0.1.dat'  using 1:3 lw 3 lt 2  with lines, \
#     'waveFunctions_rho-1-lambda-0.1.dat'  using (-$1):(-$3) lw 3 lt 2  with lines

#unset multiplot
reset #重新设置
set output
#!pdfatex temple.tex
#!xelatex  Vzk=1C2=20C3=0p40p50p6.tex 
!latex temple.tex
!dvips -E -o Vzk=1C2=20C3=0p40p50p6.eps  temple.dvi 
#!pstopdf Vzk=1C2=20C3=0p40p50p6.eps
#!open -a preview Vzk=1C2=20C3=0p40p50p6.pdf
#!gsview32 xxx.eps
!del temple*

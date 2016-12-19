file = sprintf("test_etna_%.4d.p_2d",i)

set title sprintf("time = %g seconds",i*10)
splot file u ($1/1000):($2/1000):6:3 lt -2

i = i+1
if(i<50) reread

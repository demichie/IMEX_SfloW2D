set terminal gif enhanced optimize size 720,720 animate delay 10
set output "animation.gif"

# parametri dei grafici
unset key
set xl "W-E (kilometers)"
set yl "S-N (kilometers)"
set xr[500:504]
set yr[4176:4179]
set pm3d map
set contour
unset clabel
set cntrparam lev inc 1600,100,3400 
set cbr[0:10]

i = 0
load "script.gnu"

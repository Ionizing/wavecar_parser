awk 'BEGIN{i=1} /total>/,\
                /\/total>/ \
                 {a[i]=$2 ; b[i]=$3 ; i=i+1} \
     END{for (j=10;j<i-4;j++) print a[j],b[j]}' vasprun.xml > dos.dat

ef=`awk '/efermi/ {print $3}' vasprun.xml`

cat > plotfile << !
set term postscript enhanced eps colour lw 2 "Helvetica" 20
set output "optics.eps"
plot "dos.dat" using (\$1-$ef):(\$2) w lp, "" using (\$1-$ef):(\$2) title "cspline" smooth csplines
!

gnuplot -persist plotfile

rm dos.dat plotfile

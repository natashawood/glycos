

if [ $# -eq 2 ]; then
> $2
    numLines=`wc -l < $1`
    for (( i=1; i<=$numLines; i++)); do
        aN1=`sed -n ''$i'p' $1 | cut -c 12-15 | sed 's/ //g'`
        #rN1=`sed -n ''$i'p' $1 | cut -c 17-21 | sed 's/ //g'`
        rn1=`sed -n ''$i'p' $1 | cut -c 23-26 | sed 's/ //g'`
        aN2=`sed -n ''$i'p' $1 | cut -c 41-46 | sed 's/ //g'`
        #rN2=`sed -n ''$i'p' $1 | cut -c 47-51 | sed 's/ //g'`
        rn2=`sed -n ''$i'p' $1 | cut -c 53-58 | sed 's/ //g'`
        echo "bond mol.$rn1.$aN1 mol.$rn2.$aN2" >> $2
    done
elif [ $# -eq 3 ]; then
> $2
    numLines=`wc -l < $1`
    for (( i=1; i<=$numLines; i++)); do
        aN1=`sed -n ''$i'p' $1 | cut -c 12-15 | sed 's/ //g'`
        #rN1=`sed -n ''$i'p' $1 | cut -c 17-21 | sed 's/ //g'`
        rn1=`sed -n ''$i'p' $1 | cut -c 23-26 | sed 's/ //g'`
        aN2=`sed -n ''$i'p' $1 | cut -c 41-46 | sed 's/ //g'`
        #rN2=`sed -n ''$i'p' $1 | cut -c 47-51 | sed 's/ //g'`
        rn2=`sed -n ''$i'p' $1 | cut -c 53-58 | sed 's/ //g'`
        rn1=$(($rn1 + $3 - 1))
        rn2=$(($rn2 + $3 - 1))
        echo "bond mol.$rn1.$aN1 mol.$rn2.$aN2" >> $2
    done
elif [ $# -eq 4 ]; then
> $2
    numLines=`wc -l < $1`
    for (( i=1; i<=$numLines; i++)); do
        aN1=`sed -n ''$i'p' $1 | cut -c 12-15 | sed 's/ //g'`
        #rN1=`sed -n ''$i'p' $1 | cut -c 17-21 | sed 's/ //g'`
        rn1=`sed -n ''$i'p' $1 | cut -c 23-26 | sed 's/ //g'`
        aN2=`sed -n ''$i'p' $1 | cut -c 41-46 | sed 's/ //g'`
        #rN2=`sed -n ''$i'p' $1 | cut -c 47-51 | sed 's/ //g'`
        rn2=`sed -n ''$i'p' $1 | cut -c 53-58 | sed 's/ //g'`
        rn1=$(($rn1 + $3 - 1))
        rn2=$(($rn2 + $3 - 1))
        if [ $i -eq 1 ]; then
            rn1=$4
            aN1=ND2
        fi
        echo "bond mol.$rn1.$aN1 mol.$rn2.$aN2" >> $2
    done
else
echo "Usage ./LINK2BOND.sh input output [last protein resid] [NLN resid]"
fi

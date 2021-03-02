for i in *.jl; do
    [ -f "$i" ] || break
    echo $i 
    julia $i
done
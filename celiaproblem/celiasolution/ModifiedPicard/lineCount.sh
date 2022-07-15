for f in *.py
do
    echo $f
    cat $f | tr -d " " | grep -v "^$" | grep -v "^#" | wc -l
done

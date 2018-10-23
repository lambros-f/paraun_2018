for genus in `cat frequent_species.lst`;do
ord=$(grep $genus lineages-2018-06-13.csv |grep -i fungi  | head -n1 |cut -f4 -d ',')
echo $genus' '$ord 
done

grep Leotiomycetes frequent_species.orders.lst |cut -f1 -d ' ' > frequent_species.orders.leot.lst


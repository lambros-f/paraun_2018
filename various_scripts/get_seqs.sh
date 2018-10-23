#script for getting the sequences based on a gene id using a sequence file


#echo 'Give list prots output'
#read $1 $2 $3


for genez in `cut -f1 $1`;do
  echo "Processing gene $genez "
 sed -n -e "/$genez/,/>/ p" $2 | head -n -1 >> $3
done

#echo 'telos'

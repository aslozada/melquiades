echo 'remove object & modules'

#FILE=$1

rm *.o
rm *.mod
rm Makefile
rm melquiades.x

echo 'build the Makefile'

perl makemake.perl melquiades.x 

echo 'compile'

make all

#./melquiades.x $FILE

date

rm *.o
rm *.mod

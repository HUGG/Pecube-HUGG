mkdir Pecube
cp README Pecube
mkdir Pecube/NA
cp NA/na.in Pecube/NA
mkdir Pecube/source
cp source/*.f90 Pecube/source
cp source/*.f Pecube/source
cp source/*.c Pecube/source
cp source/*.h Pecube/source
cp source/Makefile Pecube/source
cp source/*.inc Pecube/source
cp source/BUGS.txt Pecube/source
mkdir Pecube/input
cp input/* Pecube/input
mkdir Pecube/data
cp -r data/* Pecube/data
mkdir Pecube/doc
cp doc/* Pecube/doc
cp tarngo.sh Pecube
mkdir Pecube/VTK
mkdir Pecube/bin
cp bin/run.sh Pecube/bin
mkdir Pecube/RUN00
mkdir Pecube/cascade_output
cp cascade_output/* Pecube/cascade_output
tar -cvf Pecube.tar Pecube
gzip Pecube.tar
rm -r Pecube

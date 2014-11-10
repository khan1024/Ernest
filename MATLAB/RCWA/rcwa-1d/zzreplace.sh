#/bin/bash
# Simple shell script to replace 

for name in *
do
  sed -e s/grating_relative_dimensions/grating_relative_dimensions_dimensions/g  <$name >tmpfile
  mv tmpfile $name 
done

rm -f tmpfile
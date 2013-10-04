#! /bin/bash

cp IBIsCO hybridPA
cp IBIsCO atomisticPA

#Run IBIsCO on different inputs
cd hybridPA
./IBIsCO
mv s-md.trj ../hybrid.trj
rm restart
rm s-md.tp
rm timestep
rm s-md.out
rm ERROR
rm config.xyz
cd ../

cd atomisticPA
./IBIsCO
mv s-md.trj ../atomistic.trj
rm restart
rm s-md.tp
rm timestep
rm s-md.out
rm ERROR
rm config.xyz
cd ../

#Run nosetests
nosetests tests.py

#Cleanup

rm hybridPA/IBIsCO
rm atomisticPA/IBIsCO

rm hybrid.trj
rm atomistic.trj
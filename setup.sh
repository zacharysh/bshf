#! /bin/bash
#make clean
echo "Done. Removing target executable.."
rm -f ./schrodingerHF
echo "Done. Building..."
bear -- make
echo "Done."
./schrodingerHF

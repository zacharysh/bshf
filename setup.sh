#! /bin/bash
#make clean
#echo "Done. Removing target executable.."
rm -f output/Li_1s.txt
rm -f ./schrodingerHF
echo "Done. Building..."
bear -- make
echo "Done."
./schrodingerHF
gnuplot output/plot_wavefunction.gnu
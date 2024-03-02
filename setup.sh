#! /bin/bash
#make clean
#echo "Done. Removing target executable.."
rm -f output/Li.txt
rm -f ./schrodingerHF
echo "Done. Building..."
bear -- make
echo "Done."
./schrodingerHF -Z 3 -l 0
gnuplot output/plot_wavefunction.gnu
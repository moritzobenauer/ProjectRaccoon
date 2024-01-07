#!/bin/bash
# RACCOON pdb2gmx testing script, MLO@JGU 2023
# Please adjust your input structure and your forcefield
#
# **This test method ignores hydrogen atoms!**
# ** GROMACS will create charged termini **
#

cd ..
### Generating test sequences with different thresholds for geometry generation
python3 seq_gen_gmx.py
cd gmx

echo "0" "0" | gmx pdb2gmx -f system.pdb -o protein -ff oplsaa-mod2023 -his >pdb2gmx.log 2>&1

if [ $? -ne 0 ]; then
    echo "gmx pdb2gmx error. Please check log file."
else
    echo "gmx pdb2gmx finished successfully."
fi

gmx editconf -f protein.gro -o box.gro -box 10 10 10 >editconf.log 2>&1

if [ $? -ne 0 ]; then
    echo "gmx editconf error. Please check log file."
else
    echo "gmx editconf finished successfully."
fi

gmx grompp -f preminim.mdp -p topol.top -c box.gro -maxwarn 1 >gromp.log 2>&1

if [ $? -ne 0 ]; then
    echo "gmx grompp error. Please check log file."
else
    echo "gmx grompp finished successfully."
fi

gmx mdrun -deffnm topol -v -nsteps 100 >mdrun.log 2>&1

if [ $? -ne 0 ]; then
    echo "gmx mdrun error. Please check log file."
else
    echo "gmx mdrun finished successfully."
    echo "Maximum Force:"
    grep "^Maximum force     =" topol.log | awk '{print $4}'
fi

rm *.tpr
rm *.trr
rm *.*#
rm *.itp
rm *.edr
rm *.top

if [ $? -ne 0 ]; then
    echo ""
else
    echo "Check complete."
fi
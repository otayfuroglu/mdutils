#! /usr/bin/env bash

export GMX_MAXBACKUP=-1
export GMX_MAXWARN=1

file_base=$1
cd $file_base

lig=MOL
cp g16_calculation/$file_base.esp $lig.esp

# rm *.gro *.top *.itp *.rtp *.ndx
antechamber -i "$lig".esp -fi gesp -o "$lig".prepin -fo prepi -c resp -s 2 -rn "$lig" -at gaff -nc 0 -pf yes
parmchk2 -i "$lig".prepin -f prepi -o "$lig".frcmod

tleap -f- <<EOF
source leaprc.gaff
loadamberparams $lig.frcmod
loadamberprep $lig.prepin
saveamberparm $lig $lig.prmtop $lig.inpcrd
quit
EOF

amb2gro_top_gro.py -p "$lig".prmtop -c "$lig".inpcrd -t "$lig"_GMX.top -g "$lig"_GMX.gro -b "$lig"_GMX.pdb

# for MD prepration

# gmx editconf -f "$lig"_GMX.pdb -o "$lig".gro
# # rm topol.top
#
# cp ../templates/topol.top topol.top
# sed -n '/atomtypes/,/moleculetype/{/moleculetype/b;p}' "$lig"_GMX.top > "$lig"_GMX_GAFF.rtp
# sed -n '/moleculetype/,/system/{/system/b;p}' "$lig"_GMX.top > "$lig"_GMX.itp
# sed -i s/LIGAND/"$lig"/g topol.top
#
# echo 0 | gmx genrestr -f "$lig".gro -o posre_ligand.itp -fc 1000 1000 1000
# gmx editconf -f "$lig".gro -o newbox.gro -bt dodecahedron -d 1.0 -center 1.5 1.5 1.0
# gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
# gmx make_ndx -f solv.gro -o index.ndx <<EOF
# q
# EOF
#
#
# rm *.log *.com *.chk *.prmtop *.inpcrd *.pdb *.frcmod *.prepin newbox.gro esout punch qout QOUT *.xyz
# cd ..

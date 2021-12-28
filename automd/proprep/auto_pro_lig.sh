exec>log.txt 2>&1

#scr_dir=/home/modellab/scripts
workdir=$(pwd)
protein=7AMZ
ligand=5HF
con=_model
conf=1
g16dir="$workdir"/"$ligand"_G16_files
yn=n

mkdir "$ligand"_G16_files


cd "$ligand"_G16_files

cp -r ../amber99sb-ildn.ff .

sleep 3


cp ../compounds/"$ligand""$con".sdf .

#converting sdf file to xyz format
obabel -i sdf "$ligand""$con".sdf -o xyz -O "$ligand".xyz -m

obabel -i sdf "$ligand""$con".sdf -o pdb -O "$ligand".pdb

###################################G16####################################################

#Generate ESP after optimization
echo "" >empty
echo "%mem=40GB" >link0
echo "%nprocshared=32" >>link0
echo "%nosave" >>link0
echo "#p opt B3LYP/6-31G*" >route_opt
echo "#p HF/6-31G* Pop=(MK) IOp(6/50=1)"> route_scfesp
echo "#p HF/6-31G* guess=read geom=check Pop=(MK) IOp(6/50=1)">route_optesp
echo "title goes here">title
echo "0 1" >charge
for ((i=1;i<=conf;i++));
do
        #Finalizing G16 input for optimization job
        fn="$ligand""$i"
        echo "%chk=""$fn""_opt.chk" >check
        tail -n +3 <"$fn".xyz >coordinate
        echo "$fn".esp >additionalline
        cat link0 check route_opt empty title empty charge coordinate empty empty > "$fn"_opt.com


        #Finalizing G16 input for ESP job
        echo "%chk=""$fn""_ESP.chk" >>"$fn"_ESP.com
        #cp "$fn"_opt.chk "$fn"_ESP.chk
        cat link0 route_optesp empty title empty charge empty additionalline empty >> "$fn"_ESP.com


       #Finalizing G16 input for ESP without optimization
        echo "%chk=""$fn""_ESP2.chk" >"$fn"_ESP2.com
        cat link0 route_scfesp empty title empty charge coordinate empty additionalline empty >> "$fn"_ESP2.com
        rm empty link0 route_opt route_scfesp route_optesp title charge check coordinate additionalline
done

for ((i=1;i<=conf;i++));do
    fn="$ligand""$i"
    case $yn in
        [yY][eE][sS]|[yY])
                g16 <"$fn"_opt.com> "$fn"_opt.log
                 cp "$fn"_opt.chk "$fn"_ESP.chk
                 g16 <"$fn"_ESP.com> "$fn"_ESP.log
                 ;;
       [nN][oO]|[nN])
               g16 <"$fn"_ESP2.com> "$fn"_ESP2.log
               ;;
       *)
       echo "You did not specify, exiting"
    esac
done

#####################################################



antechamber -i "$ligand""$conf".esp -fi gesp -o "$ligand".prepin -fo prepi -c resp -s 2 -rn "$ligand" -at gaff -nc 0 -pf yes
sleep 3
parmchk2 -i "$ligand".prepin -f prepi -o "$ligand".frcmod

tleap -f- <<EOF
source leaprc.gaff
loadamberparams $ligand.frcmod
loadamberprep $ligand.prepin
saveamberparm $ligand $ligand.prmtop $ligand.inpcrd
quit
EOF



amb2gro_top_gro.py -p "$ligand".prmtop -c "$ligand".inpcrd -t "$ligand"_GMX.top -g "$ligand"_GMX.gro -b "$ligand"_GMX.pdb

cp ../template/topol.top .

sed -n '/atomtypes/,/moleculetype/{/moleculetype/b;p}' "$ligand"_GMX.top > "$ligand"_GMX_GAFF.rtp
sed -n '/moleculetype/,/system/{/system/b;p}' "$ligand"_GMX.top > "$ligand"_GMX.itp
sed -i s/LIGAND/"$ligand"/g topol.top




#######################################################################################################

obabel -i pdb "$ligand"_GMX.pdb -o pdb -O "$ligand"_GMX.pdb

python ../align.py -i "$ligand"_GMX.pdb -r "$ligand".pdb  -o "$ligand"_GMX_aligned.pdb

sed -i -e "1d" "$ligand"_GMX_aligned.pdb
########################################################################################################

gmx editconf -f "$ligand"_GMX_aligned.pdb -o ligand.gro

grep "ATOM  " ../compounds/"$protein".pdb > $g16dir/protein.pdb


gmx pdb2gmx -f protein.pdb -o protein.gro -p topol_protein.top -i posre_protein.itp -ff amber99sb-ildn -water tip3p -ignh

sed -n '/moleculetype/,/endif/{;p}' topol_protein.top > topol_protein.itp
sed -i 's/Protein_chain_A/Protein /' topol_protein.itp
sleep 3


filename='protein.gro'
pro_atm=$(cat ${filename} | wc -l)
pro_atm=$((pro_atm-3))
echo "There are " $pro_atm " atoms in the protein"


filename='ligand.gro'
lig_atm=$(cat ${filename} | wc -l)
lig_atm=$((lig_atm-3))
echo "There are " $lig_atm " atoms in the ligand"

total=$((pro_atm+lig_atm))
echo "Total is " $total
cp protein.gro tmp.gro
sed -i '$ d' tmp.gro
grep "$ligand" ligand.gro >> tmp.gro
tail -n 1 protein.gro >>tmp.gro
sed 0,/$pro_atm/s//$total/ tmp.gro > complex.gro
rm tmp.gro
sleep 3


echo 0 | gmx genrestr -f ligand.gro -o posre_ligand.itp -fc 1000 1000 1000
gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

#################################################################
cp -r ../template/MDP .


filename="MDP"

filepath=$(pwd)
search1="energygrps"
replace1=$ligand
search2="tc-grps"
search3="couple-moltype"

i=0;

for file in $(grep -l -R $search1 $filepath"/"$filename)
do
   sed -e "/$search1/c\energygrps = Protein "$replace1 $file > tempfile.tmp
   mv tempfile.tmp $file

  let i++;

  #echo "Modified: " $file
done

for file in $(grep -l -R $search2 $filepath"/"$filename)
do
   sed -e "/$search2/c\tc-grps     = Protein_$replace1 Water_and_ions " $file > tempfile.tmp
   mv tempfile.tmp $file

  let i++;

  #echo "Modified: " $file
done

for file in $(grep -l -R $search3 $filepath"/"$filename)
do
   sed -e "/$search3/c\couple-moltype           = "$replace1 $file > tempfile.tmp
   mv tempfile.tmp $file

  let i++;

  #echo "Modified: " $file
done




###############################################################
gmx grompp -maxwarn 1 -f MDP/EM/em_steep.mdp -c solv.gro -r solv.gro -p topol.top -o ions.tpr
echo 15 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
gmx make_ndx -f solv_ions.gro -o index.ndx <<EOF
1|13
q
EOF

cd ..

mkdir MDRUNS MDRUNS/EM MDRUNS/EM_1 MDRUNS/EM_2 MDRUNS/NVT MDRUNS/NPT MDRUNS/Production_MD

mkdir MD_"$protein"_"$ligand"

cd $g16dir

cp -r {../amber99sb-ildn.ff,../MDRUNS,./MDP,../normalMD_akyacuda.sh,posre_ligand.itp,topol_protein.itp,posre_protein.itp,"$ligand"_GMX.itp,"$ligand"_GMX_GAFF.rtp,../log.txt,index.ndx,solv_ions.gro,topol.top,ions.tpr} ../MDRUNS


mv  ../MDRUNS ../MD_"$protein"_"$ligand"


cd ..

rm -r $g16dir log.txt 


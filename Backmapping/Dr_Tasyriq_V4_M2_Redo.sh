#!/bin/bash

echo "##########################################################"
echo "############# WELCOME TO PROJECT XXX   ###################"
echo "##########################################################"

# Define the list of loop values
loop_values=(70 99 231 247 529 811 885 931)

for input in "${loop_values[@]}"; do
    echo "Processing for $input ns"

    echo "Script will now run for $input ns"

###################################################################################
## cd /home/pjjauh05/Desktop/Hands-on/backward-v5/"$input"ns ##
###################################################################################
cd /home/pjjauh05/Desktop/Hands-on/backward-v5/"$input"ns

gmx trjconv -f dynamic_fitted.xtc -s dynamic.tpr -n index.ndx -o trj_"$input"ns.gro -dump "$input"000 <<EOF
1
EOF

./initram-v5.sh -f trj_"$input"ns.gro -o aa_"$input"ns.gro -to charmm36 -p topol.top

##############################################################################
### copy top strings in trj_"$input"ns.gro paste into *aa_"$input"ns.gro ###
### Grep top string after ".pdb" from a file and paste it into another file. #
##############################################################################
input2=$(grep -Po ".pdb.\K.*" trj_"$input"ns.gro) 
echo $input2
##############################################
### Paste it into selected file aa_****.gro ##
##############################################
sed -i "1s/.*/Protein - TQ $input2/" aa_"$input"ns.gro


gmx trjconv -f aa_"$input"ns.gro -o aa_"$input"ns.xtc

########################################################################################################
## cp aa_"$input"ns.xtc /home/pjjauh05/Desktop/Hands-on/backward-v5/all_xtc ###
########################################################################################################
cp aa_"$input"ns.xtc all_xtc
cp aa_"$input"ns.gro gro

#cp /home/pjjauh05/Desktop/Hands-on/backward-v5/aa_0ns.xtc /home/pjjauh05/Desktop/Hands-on/backward-v5/All_xtc/ <-- reference


echo "####################################################################"
echo "####################################################################"
echo "################# PROJECT XXX FINISHED RUNNING FOR #################"
echo "                             $input" "NS                          "
echo "####################################################################"
echo "####################################################################"
echo "####################################################################"




input=$(($input+1))
done 


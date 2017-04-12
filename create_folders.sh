#create_folders_2

declare -a reg_type=(aff aff_def def rig rig_aff rig_aff_def)

declare -a atl=(0364 0417 0424 0473 0485 0493 0495 0519 0529 0574 0578 0591 0601 0632 0715 0730 0917 0921 1073 1076 1166  1168 1171)

declare -a targ=(0364 0417 0424 0473 0485 0493 0495 0519 0529 0574 0578 0591 0601 0632 0715 0730 0917 0921 1073 1076 1166  1168 1171)


data_f="/home/fusta/Dropbox/CODING/UBUNTU/git/Multi-Atlas-Segmentation/Registration-Pipeline/data/2.1.initial_rigid_norm_iso/"

targ_f="/home/fusta/Dropbox/CODING/UBUNTU/git/Multi-Atlas-Segmentation/Registration-Pipeline/tmp/2.1.initial_rigid_norm_iso/"

cd $targ_f;

#creating target folders and reg_type folders inside
for i in ${targ[@]}; do 
	echo $i; 
	mkdir $i;
	cd $i;
	mkdir aff aff_def def LGEtoMRA MRAtoLGE rig rig_aff rig_aff_def;
	cd LGEtoMRA;
	mkdir aff aff_def def rig rig_aff rig_aff_def;
	cd ..;
	cd MRAtoLGE;
	mkdir aff aff_def def rig rig_aff rig_aff_def;
	cd ../..;	
done


#create targettoatlas folders inside reg_type folders
to="to";
echo $to;


for t in ${targ[@]}; do
	chmod u+rw $PWD
	cd $targ_f;
	echo $targ_f;
	ls;
	chmod u+rw $PWD
	cd $t;
	echo in_target_folder$t; 
	ls;

	
	for a in ${atl[@]}; do 
		for reg_type in "aff" "aff_def" "def" "rig" "rig_aff" "rig_aff_def"; do
			echo atlas$a; 
			concstr=$t$to$a;
			echo $concstr;
			echo $PWD
			cd $reg_type;
			echo $PWD
			mkdir $concstr;
			#rm $concstr;
			cd ..;
		done
	done
	cd ..;#come out to the target folder
		
done



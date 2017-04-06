#try LGEtoMRA with aff+def

# LGEtoMRA_aff + intl_rigid + AFFINE + DEF+JLF + MRAtoLGE_aff_def

#it is normalized and isotropic.
#initial rigid registered with only TRANSFORMATION

#be le 0.30
#alpha 0.5 beta 4
#rp 444 rs 333

#LGEtoMRA means LGE is the floarting image, MRA is the fixed image
#MRAtoLGE is the vice versa
from timeit import default_timer as timer
import os

#reference image = target image
data_root = '/home/fusta/Dropbox/CODING/UBUNTU/git/Multi-Atlas-Segmentation/Registration-Pipeline/data/2.1.initial_rigid_norm_iso/'
proc_dump = '/home/fusta/Dropbox/CODING/UBUNTU/git/Multi-Atlas-Segmentation/Registration-Pipeline/tmp/2.1.initial_rigid_norm_iso/'
nifty_reg_dir = '/home/fusta/git/nifty_reg/install/bin'
jointfusion_dir = '/home/fusta/git/JLF_ITK_based/PICSL_MALF/'

#use nifty_reg nightly build for functional Nearest Neighbour interpolation
reg_resample = os.path.join(nifty_reg_dir, 'reg_resample')
reg_aladin = os.path.join(nifty_reg_dir, 'reg_aladin')
reg_f3d = os.path.join(nifty_reg_dir, 'reg_f3d')
jointfusion = os.path.join(jointfusion_dir, 'jointfusion')

#FLags
enableCmdFlag = True

LGEtoMRA_registerFlag=True
LGEtoMRA_registerResampleFlag=True

affRegFlag=True
affResampleLabelFLag=True
defRegFlag=True
defResampleLabelFlag=True

fuseAllLabelFlag=True
MRAtoLGE_registerFlag=True
MRAtoLGE_registerResampleFlag=True

calcDSCFlag=0
saveDSCResultsFlag=0

aff_reg_type_LGEtoMRA = 'aff'
def_reg_LGEtoMRA = False#<--
aff_reg_type = 'aff'
def_reg = True #def_reg = True
aff_reg_type_MRAtoLGE = 'aff'
def_reg_MRAtoLGE = True#<--


if def_reg_LGEtoMRA == True:
    global_reg_type_LGEtoMRA = aff_reg_type_LGEtoMRA + '_def'
else:
    global_reg_type_LGEtoMRA = aff_reg_type_LGEtoMRA

if def_reg == True:
    global_reg_type = aff_reg_type + '_def'
else:
    global_reg_type = aff_reg_type

if def_reg_MRAtoLGE == True:
    global_reg_type_MRAtoLGE = aff_reg_type_MRAtoLGE + '_def'
else:
    global_reg_type_MRAtoLGE = aff_reg_type_MRAtoLGE

#PARAMETERS
#Registration Parameters

#' -be 0.010' + ' -le 0.010 0.010' <- how you change it
be='0.30' #be : bending energy
le='0.30 0.30' #le: linear elasticitty
le2 = '0.30'
#Registration resampling parameters:
NN='[0]' #[0] indicates Nearest Neighbour

#label-fusion parameters

#alpha1='0.08'
#alpha2='0.05'
alpha='0.5' #alpha = default [0.1] regularization term,
beta='4' #beta = default [2],  exponent for mapping intensity differnece to joint error

rp = '4x4x4 ' #rp radius: [2x2x2]patch radius for similarity measures
rs = '3x3x3 ' #rs radius: [?] localvsearch radius
rp2 = '444'
rs2 = '333'
# '-m Joint[' + alpha + ',' + beta '] -rp 4x4x4 -rs 3x3x3x '
#  '-m Joint[0.1,2] -rp 4x4x4 -rs 3x3x3x ' <- how you change it
#  '-m Joint[ALPHA,BETA] -rp 4x4x4 -rs 3x3x3x ' <- how you change it

#atlas_dataset = ['0485', '0546', '0578', '1168', '0715']
# #['0485', '0495', '0515', '0529', '0546', '0565', '0578',  '0715', '1115', '1168']
atlas_dataset = ['0485', '0495', '0515', '0529', '0565', '0578', '0715', '1115', '1168']# 2 more atlas dataset down there to change
atlas_dataset_backup = ['0485', '0495', '0515', '0529', '0565', '0578', '0715', '1115', '1168']
target_dataset = ['0485', '0495', '0515', '0529', '0565', '0578', '0715', '1115', '1168']#, '0495', '0515', '0529', '0565', '0578',  '0715', '1115', '1168']#['0485', '0546', '0578', '1168', '0715'] #have to use only one target
target_dataset2 = ['0485', '0495', '0515', '0529', '0565', '0578', '0715', '1115', '1168']#, '0495', '0515', '0529', '0565', '0578',  '0715', '1115', '1168']#['0485', '0546', '0578', '1168', '0715'] #for JLF, copy of target_dataset

start = timer()
#reg_type = 'aff'#it means start with affine, (no rigid)
#global_reg_type = 'aff_def'

for target_id in target_dataset: #target dataset for cross_correlation, if no cross-correlation, uses only one target
    #        atlas_dataset = ['P00', 'P01', 'P02', 'P03', 'P04', 'P05' , 'P06', 'P07', 'P08', 'P09', 'P10' ]
    #start_aff = timer()
    #target_img = os.path.join(data_root, target_id, 'MRA0715.nii')#target is LGE, floating is MRA
    atlas_dataset = ['0485', '0495', '0515', '0529', '0565', '0578',  '0715', '1115', '1168']
    target_img = os.path.join(data_root, target_id, 'MRA' + target_id + '.nii')#target is LGE, floating is MRA

    for a in atlas_dataset:
        if LGEtoMRA_registerFlag == True:
            flo_img = os.path.join(data_root, a, 'LGE' + a + '.nii')
            target_img = os.path.join(data_root, a, 'MRA' + a + '.nii')  # target is LGE, floating is MRA

            aff_mat = os.path.join(proc_dump, a, aff_reg_type_LGEtoMRA, 'LGEtoMRA' + aff_reg_type_LGEtoMRA + '_mat' + a)
            aff_img = os.path.join(proc_dump, a, aff_reg_type_LGEtoMRA, 'LGEtoMRA' + aff_reg_type_LGEtoMRA + '_MRA' + a + '.nii')
            # affine registration
            if aff_reg_type_LGEtoMRA == 'rig':
                cmd_aff_LGEtoMRA = reg_aladin + ' -ref ' + target_img + ' -flo ' + flo_img + ' -res ' + aff_img + ' -aff ' + aff_mat + ' -iso' + ' -rigOnly'

            if aff_reg_type_LGEtoMRA == 'aff':
                cmd_aff_LGEtoMRA = reg_aladin + ' -ref ' + target_img + ' -flo ' + flo_img + ' -res ' + aff_img + ' -aff ' + aff_mat + ' -iso' + ' -affDirect'

            if aff_reg_type_LGEtoMRA == 'rig_aff':
                cmd_aff_LGEtoMRA = reg_aladin + ' -ref ' + target_img + ' -flo ' + flo_img + ' -res ' + aff_img + ' -aff ' + aff_mat + ' -iso'
            print(cmd_aff_LGEtoMRA)
            if enableCmdFlag == True:
                os.system(cmd_aff_LGEtoMRA)

            if def_reg_LGEtoMRA == True:
                def_img = os.path.join(proc_dump, a, global_reg_type, 'LGEtoMRA' + global_reg_type + '_LGE' + a + '.nii')
                def_output_cpp = os.path.join(proc_dump, a, global_reg_type, 'LGEtoMRA' + global_reg_type + '_LGE_CPP' + a + '.nii')
                # deformable registration
                cmd_def_LGEtoMRA = reg_f3d + ' -ref ' + target_img + ' -flo ' + aff_img + ' -res ' + def_img + ' -cpp ' + def_output_cpp + ' -be ' + be + ' -le ' + le
                print(cmd_def_LGEtoMRA)
                if enableCmdFlag == True:
                    os.system(cmd_def_LGEtoMRA)

        if LGEtoMRA_registerResampleFlag == True:
            # flo_lbl1 = os.path.join(data_root, target_id, 'LGE' + target_id + '.nii')
            # Affine registered result will be resampled without deformable
            # flo_lbl2 = os.path.join(proc_dump, target_id, aff_reg_type_LGE_MRA, 'myo_result_iso_norm_be' + be + '-le' + le2 + '-alpha' + alpha + '-beta' + beta + '-rp' + rp2 + '-rs' + rs2 + target_id + '.nii ')
            # global could be either affine or def, incase of affine, flo_lbl_2 will be overwritten
            flo_lbl3 = os.path.join(data_root, a, 'myo' + a + '.nii')
            # flo_lbl4 = os.path.join(proc_dump, target_id, global_reg_type, 'myo_result_iso_norm_be' + be + '-le' + le2 + '-alpha' + alpha + '-beta' + beta + '-rp' + rp2 + '-rs' + rs2 + target_id + '.nii ')

            aff_mat = os.path.join(proc_dump, a, aff_reg_type_LGEtoMRA, 'LGEtoMRA' + aff_reg_type_LGEtoMRA + '_mat' + a)
            # affine
            # res_lbl2 = os.path.join(proc_dump, target_id, 'LGEtoMRA', aff_reg_type_LGE_MRA, 'LGEtoMRA_' + aff_reg_type_LGE_MRA + '_LGE' + target_id + '.nii')
            # affine or deformable label
            #        res_lbl3 = os.path.join(proc_dump, target_id, 'MRAtoLGE', global_reg_type, 'MRAtoLGE_' + global_reg_type + '_myo' + target_id + '.nii')
            res_lbl3 = os.path.join(data_root, a, 'myo' + a + 'forMRA.nii')
            # res_lbl4 = os.path.join(proc_dump, target_id, 'LGEtoMRA', global_reg_type, 'LGEtoMRA_' + global_reg_type + '_myo_result_iso_norm_be' + be + '-le' + le2 + '-alpha' + alpha + '-beta' + beta + '-rp' + rp2 + '-rs' + rs2 + target_id + '.nii ')

            # saving both affine and deformable version of LGE-toMRA registrations
            # cmd_aff_lbl2 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl2 + ' -res ' + res_lbl2 + ' -trans ' + aff_mat + ' -inter ' + NN
            # os.system(cmd_aff_lbl2)

            cmd_aff_lbl3 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl3 + ' -res ' + res_lbl3 + ' -trans ' + aff_mat + ' -inter ' + NN
            if enableCmdFlag == True:
                os.system(cmd_aff_lbl3)

                # cmd_aff_lbl4 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl4 + ' -res ' + res_lbl4 + ' -trans ' + aff_mat + ' -inter ' + NN
                # os.system(cmd_aff_lbl4)
    LGEtoMRA_registerFlag = False
    LGEtoMRA_registerResampleFlag = False

    atlas_dataset.remove(target_id)

    if affRegFlag == True:
        #reg_type = 'aff' #aff or rigid
        for a in atlas_dataset:
            flo_img = os.path.join(data_root, a, 'MRA' + a + '.nii')

            aff_mat = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_mat' + a)
            aff_img = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_MRA' + a + '.nii')

            # affine registration
            if aff_reg_type == 'rig':
                cmd_aff = reg_aladin + ' -ref ' + target_img + ' -flo ' + flo_img + ' -res ' + aff_img + ' -aff ' + aff_mat + ' -iso' + ' -rigOnly'
            if aff_reg_type =='aff':
                cmd_aff = reg_aladin + ' -ref ' + target_img + ' -flo ' + flo_img + ' -res ' + aff_img + ' -aff ' + aff_mat + ' -iso' + ' -affDirect'
            if aff_reg_type == 'rig_aff':
                cmd_aff = reg_aladin + ' -ref ' + target_img + ' -flo ' + flo_img + ' -res ' + aff_img + ' -aff ' + aff_mat + ' -iso'

            print(cmd_aff)
            if enableCmdFlag == True:
                os.system(cmd_aff)

    if affResampleLabelFLag == True:
        for a in atlas_dataset:
            aff_mat = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_mat' + a)

            #images to be registered
            flo_lbl1 = os.path.join(data_root, a, 'LGE' + a + '.nii')
            flo_lbl2 = os.path.join(data_root, a, 'myo' + a + 'forMRA.nii')
            #flo_lbl3 = os.path.join(data_root, a, 'myo_corrected' + a + '.nii')
            flo_lbl4 = os.path.join(data_root, a, 'scar' + a + '.nii')

            #resulted images are written
            aff_lbl1 = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_LGE' + a + '.nii')
            aff_lbl2 = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_myo' + a + '.nii')
            #aff_lbl3 = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_myo_corrected' + a + '.nii')
            aff_lbl4 = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_scar' + a + '.nii')

            # resampling labels using affine registration results
            cmd_aff_lbl1 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl1 + ' -res ' + aff_lbl1 + ' -trans ' + aff_mat + ' -inter ' + NN
            print(cmd_aff_lbl1)
            if enableCmdFlag == True:
                os.system(cmd_aff_lbl1)

            cmd_aff_lbl2 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl2 + ' -res ' + aff_lbl2 + ' -trans ' + aff_mat + ' -inter ' + NN
            if enableCmdFlag == True:
                os.system(cmd_aff_lbl2)

            #cmd_aff_lbl3 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl3 + ' -res ' + aff_lbl3 + ' -trans ' + aff_mat + ' -inter ' + NN
            #os.system(cmd_aff_lbl3)

            cmd_aff_lbl4 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl4 + ' -res ' + aff_lbl4 + ' -trans ' + aff_mat + ' -inter ' + NN
            if enableCmdFlag == True:
                os.system(cmd_aff_lbl4)

    if defRegFlag == True:
        for a in atlas_dataset:
            flo_img = os.path.join(data_root, a, 'MRA' + a + '.nii')

            aff_img = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_MRA' + a + '.nii')
            def_img = os.path.join(proc_dump, a, global_reg_type, a + 'to' + target_id, global_reg_type + '_MRA' + a + '.nii')
            def_output_cpp = os.path.join(proc_dump, a, global_reg_type, a + 'to' + target_id, global_reg_type + '_MRA_CPP' + a + '.nii')

            #deformable registration
            cmd_def = reg_f3d + ' -ref ' + target_img + ' -flo ' + aff_img + ' -res ' + def_img + ' -cpp ' + def_output_cpp + ' -be ' + be + ' -le ' + le

            print(cmd_def)
            if enableCmdFlag == True:
                os.system(cmd_def)

    if defResampleLabelFlag == True:
        for a in atlas_dataset:
            #I am not resampling MRA, bcs it is already transformed above
            flo_lbl1 = os.path.join(data_root, a, 'LGE' + a + '.nii')
            flo_lbl2 = os.path.join(data_root, a, 'myo' + a + 'forMRA.nii')
            #flo_lbl3 = os.path.join(data_root, a, 'myo_corrected' + a + '.nii')
            flo_lbl4 = os.path.join(data_root, a, 'scar' + a + '.nii')

            #flo_lbl2 = os.path.join(data_root, a, 'myo_corrected0715.nii')
            def_lbl1 = os.path.join(proc_dump, a, global_reg_type, a + 'to' + target_id, global_reg_type + '_LGE' + a + '.nii')
            def_lbl2 = os.path.join(proc_dump, a, global_reg_type, a + 'to' + target_id, global_reg_type + '_myo' + a + '.nii')
            #def_lbl3 = os.path.join(proc_dump, a, global_reg_type, a + 'to' + target_id, global_reg_type + '_myo_corrected' + a + '.nii')
            def_lbl4 = os.path.join(proc_dump, a, global_reg_type, a + 'to' + target_id, global_reg_type + '_scar' + a + '.nii')

            def_output_cpp = os.path.join(proc_dump, a,  global_reg_type, a + 'to' + target_id, global_reg_type + '_MRA_CPP' + a + '.nii')

            aff_lbl1 = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_LGE' + a + '.nii')
            aff_lbl2 = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_myo' + a + '.nii')
            #aff_lbl3 = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_myo_corrected' + a + '.nii')
            aff_lbl4 = os.path.join(proc_dump, a, aff_reg_type, a + 'to' + target_id, aff_reg_type + '_scar' + a + '.nii')

            #resampling labels using deformable registration results
            cmd_def_lbl1 = reg_resample + ' -ref ' + target_img + ' -flo ' + aff_lbl1 + ' -res ' + def_lbl1 + ' -trans ' + def_output_cpp + ' -inter ' + NN
            print(cmd_def_lbl1)
            if enableCmdFlag == True:
                os.system(cmd_def_lbl1)

            cmd_def_lbl2 = reg_resample + ' -ref ' + target_img + ' -flo ' + aff_lbl2 + ' -res ' + def_lbl2 + ' -trans ' + def_output_cpp + ' -inter ' + NN
            if enableCmdFlag == True:
                os.system(cmd_def_lbl2)

            #cmd_def_lbl3 = reg_resample + ' -ref ' + target_img + ' -flo ' + aff_lbl3 + ' -res ' + def_lbl3 + ' -trans ' + def_output_cpp + ' -inter ' + NN
            #os.system(cmd_def_lbl3)

            cmd_def_lbl4 = reg_resample + ' -ref ' + target_img + ' -flo ' + aff_lbl4 + ' -res ' + def_lbl4 + ' -trans ' + def_output_cpp + ' -inter ' + NN
            if enableCmdFlag == True:
                os.system(cmd_def_lbl4)


            #cmd_def_lbl2 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl2 + ' -res ' + def_lbl2 + ' -trans ' + def_output_cpp + ' -inter ' + NN
            #print(cmd_def_lbl2)
            #os.system(cmd_def_lbl2)

    #atlas_dataset = ['0485', '0546', '0578', '1168', '0715']#, 'P01', 'P04', 'P06', 'P08', 'P10']#, 'P07', 'P08', 'P09', 'P10']
    #atlas_dataset = atlas_dataset_backup

    target_img = os.path.join(data_root, target_id, 'MRA' + target_id + '.nii')
    if fuseAllLabelFlag == True:
        #start_JLF = timer()
        #reg_type = 'rig_aff_def'
        #here it does label_fusion once for the first target, need to put another loop bove and circulate it 10 times, need to arrange the folders as well
        #for target_subj_id in atlas_dataset:

        atlas2 = ''
        label2 = ''
            #PREPARE DATASET
        #atlas_dataset.remove(target_id)

        for a in atlas_dataset:
                    #now run label_fusion with the reamining

            temp_path_atlas = os.path.join(proc_dump, a, global_reg_type, a + 'to' + target_id, global_reg_type + '_MRA' + a + '.nii ')
            temp_path_label = os.path.join(proc_dump, a, global_reg_type, a + 'to' + target_id, global_reg_type + '_myo' + a + '.nii ')

            #atlas.append(temp_path_atlas)
            atlas2 = atlas2 + '' + temp_path_atlas
            #label.append(temp_path_label)
            label2 = label2 + '' + temp_path_label
            #done with preparation, now run label fusion once

        #atlas_dataset = ['0485', '0546', '0578', '1168', '0715']#, 'P01', 'P04', 'P06', 'P08', 'P10']#, 'P07', 'P08', 'P09', 'P10']
        atlas_seg_result = os.path.join(proc_dump, target_id, global_reg_type, 'myo_result_iso_norm_be' + be + '-le' + le2 + '-alpha' + alpha + '-beta' + beta + '-rp' + rp2 + '-rs' + rs2 + target_id + '.nii ')

        #label fusion
        cmdfuse = jointfusion + ' 3 1 -g ' + atlas2 + '-tg ' + target_img + ' -l ' + label2 + '-m Joint[' + alpha + ',' + beta + '] -rp ' + rp + '-rs ' + rs + atlas_seg_result
        #     '-m Joint[0.1,2] -rp 4x4x4 -rs 3x3x3x ' + atlas_seg_result

        print(cmdfuse)
        if enableCmdFlag == True:
            os.system(cmdfuse)

    if MRAtoLGE_registerFlag == True:
        #reg_type = 'rig_aff_def'
        #reg_type_MRAtoLGE = 'aff_def' #aff or rigid
        aff_reg_type_MRAtoLGE = 'aff'

        if def_reg_MRAtoLGE == True:
            global_reg_type_MRAtoLGE = aff_reg_type_MRAtoLGE + '_def'
        else:
            global_reg_type_MRAtoLGE = aff_reg_type_MRAtoLGE

        target_img = os.path.join(data_root, target_id, 'LGE' + target_id + '.nii')#target is LGE, floating is MRA
        flo_img = os.path.join(data_root, target_id, 'MRA' + target_id + '.nii')

        aff_mat = os.path.join(proc_dump, target_id, 'MRAtoLGE', aff_reg_type_MRAtoLGE, 'MRAtoLGE_' + aff_reg_type_MRAtoLGE + '_mat' + target_id)
        aff_img = os.path.join(proc_dump, target_id, 'MRAtoLGE', aff_reg_type_MRAtoLGE, 'MRAtoLGE_' + aff_reg_type_MRAtoLGE + '_MRA' + target_id + '.nii')
        # affine registration
        if aff_reg_type_MRAtoLGE == 'rig':
            cmd_aff_MRAtoLGE = reg_aladin + ' -ref ' + target_img + ' -flo ' + flo_img + ' -res ' + aff_img + ' -aff ' + aff_mat + ' -iso' + ' -rigOnly'

        if aff_reg_type_MRAtoLGE == 'aff':
            cmd_aff_MRAtoLGE = reg_aladin + ' -ref ' + target_img + ' -flo ' + flo_img + ' -res ' + aff_img + ' -aff ' + aff_mat + ' -iso' + ' -affDirect'

        if aff_reg_type_MRAtoLGE == 'rig_aff':
            cmd_aff_MRAtoLGE = reg_aladin + ' -ref ' + target_img + ' -flo ' + flo_img + ' -res ' + aff_img + ' -aff ' + aff_mat + ' -iso'
        print(cmd_aff_MRAtoLGE)
        if enableCmdFlag == True:
            os.system(cmd_aff_MRAtoLGE)

        if def_reg_MRAtoLGE==True:
            def_img = os.path.join(proc_dump, target_id, 'MRAtoLGE', global_reg_type, 'MRAtoLGE_' + global_reg_type + '_MRA' + target_id + '.nii')
            def_output_cpp = os.path.join(proc_dump, target_id, 'MRAtoLGE', global_reg_type, 'MRAtoLGE_' + global_reg_type + '_MRA_CPP' + target_id + '.nii')
            #deformable registration
            cmd_def_MRAtoLGE = reg_f3d + ' -ref ' + target_img + ' -flo ' + aff_img + ' -res ' + def_img + ' -cpp ' + def_output_cpp + ' -be ' + be + ' -le ' + le
            print(cmd_def_MRAtoLGE)
            if enableCmdFlag == True:
                os.system(cmd_def_MRAtoLGE)


    if MRAtoLGE_registerResampleFlag == True:
        #flo_lbl1 = os.path.join(data_root, target_id_2, 'myo' + target_id_2 + '.nii')
        #flo_lbl2 = os.path.join(data_root, target_id_2, 'myo-corrected' + target_id_2 + '.nii')

        if def_reg_MRAtoLGE == True:
            global_reg_type_MRAtoLGE = aff_reg_type_MRAtoLGE + '_def'
        else:
            global_reg_type_MRAtoLGE = aff_reg_type_MRAtoLGE

        #flo_lbl1 = os.path.join(data_root, target_id, 'LGE' + target_id + '.nii')
        #Affine registered result will be resampled without deformable
        flo_lbl2 = os.path.join(proc_dump, target_id, aff_reg_type_MRAtoLGE, 'myo_result_iso_norm_be' + be + '-le' + le2 + '-alpha' + alpha + '-beta' + beta + '-rp' + rp2 + '-rs' + rs2 + target_id + '.nii ')
        #global could be either affine or def, incase of affine, flo_lbl_2 will be overwritten
        flo_lbl3 = os.path.join(data_root, target_id, 'LGE' + target_id + '.nii')
        flo_lbl4 = os.path.join(proc_dump, target_id, global_reg_type, 'myo_result_iso_norm_be' + be + '-le' + le2 + '-alpha' + alpha + '-beta' + beta + '-rp' + rp2 + '-rs' + rs2 + target_id + '.nii ')

        aff_mat = os.path.join(proc_dump, target_id, 'MRAtoLGE', aff_reg_type_MRAtoLGE, 'MRAtoLGE_' + aff_reg_type_MRAtoLGE + '_mat' + target_id)
        #affine
        res_lbl2 = os.path.join(proc_dump, target_id, 'MRAtoLGE', aff_reg_type_MRAtoLGE, 'MRAtoLGE_' + aff_reg_type_MRAtoLGE + '_LGE' + target_id + '.nii')
        #affine or deformable label
        res_lbl3 = os.path.join(proc_dump, target_id, 'MRAtoLGE', global_reg_type, 'MRAtoLGE_' + global_reg_type + '_LGE' + target_id + '.nii')
        res_lbl4 = os.path.join(proc_dump, target_id, 'MRAtoLGE', global_reg_type, 'MRAtoLGE_' + global_reg_type + '_myo_result_iso_norm_be' + be + '-le' + le2 + '-alpha' + alpha + '-beta' + beta + '-rp' + rp2 + '-rs' + rs2 + target_id + '.nii ')

        #saving both affine and deformable version of LGE-toMRA registrations
        cmd_aff_lbl2 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl2 + ' -res ' + res_lbl2 + ' -trans ' + aff_mat + ' -inter ' + NN
        if enableCmdFlag == True:
            os.system(cmd_aff_lbl2)

        cmd_aff_lbl3 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl3 + ' -res ' + res_lbl3 + ' -trans ' + aff_mat + ' -inter ' + NN
        if enableCmdFlag == True:
            os.system(cmd_aff_lbl3)

        cmd_aff_lbl4 = reg_resample + ' -ref ' + target_img + ' -flo ' + flo_lbl4 + ' -res ' + res_lbl4 + ' -trans ' + aff_mat + ' -inter ' + NN
        if enableCmdFlag == True:
            os.system(cmd_aff_lbl4)


all_alg_time = timer() - start
print('all process took ' + str(all_alg_time) + 'seconds')
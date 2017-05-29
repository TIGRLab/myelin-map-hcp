#!/bin/bash

subject='PRE01_SFD_32010_01'
hcp='/archive/data/PRELAPSE/pipelines/hcp'
raw='/archive/data/PRELAPSE/data/nii'
SUBJECTS_DIR='/archive/data/PRELAPSE/pipelines/freesurfer'

# shortcuts
t1="${hcp}/${subject}/T1w"
t1_native="${t1}/Native"
mni="${hcp}/${subject}/MNINonLinear"
mni_native="${mni}/Native"

# file from ciftify... (text)
fs_labels=${hcp}/${subject}/

l_gm_val="3"
l_wm_val="2"
r_rm_val="42"
r_wm_val="41"

low_res_mesh="32" # for fsaverage_LR32k

myelin_sigma=$(echo "$MyelinMappingFWHM / (2*(sqrt(2*l(2))))" | bc -l)
surface_sigma=$(echo "$SurfaceSmoothingFWHM / (2*(sqrt(2*l(2))))" | bc -l)
correction_sigma="5"
bias_sigma="5"

# create native space ribbon
if [ ! -f ${t1}/ribbon.nii.gz ]; then
    for hemi in L R ; do

        if [ ${hemi} = "L" ]; then
          gm_val=${l_gm_val}
          wm_val=${l_wm_val}
        elif [ ${hemi} = "R" ]; then
          gm_val=${r_gm_val}
          wm_val=${r_wm_val}
        fi

        wb_command -create-signed-distance-volume \
            ${t1_native}/${subject}.${hemi}.white.native.surf.gii \
            ${t1}/T1w.nii.gz \
            ${t1_native}/${subject}.${hemi}.white.native.nii.gz
        wb_command -create-signed-distance-volume \
            ${t1_native}/${subject}.${hemi}.pial.native.surf.gii \
            ${t1}/T1w.nii.gz \
            ${t1_native}/${subject}.${hemi}.pial.native.nii.gz

        # find grey matter mask (ribbon)
        fslmaths \
            ${t1_native}/${subject}.${hemi}.white.native.nii.gz -thr 0 -bin -mul 255 \
            ${t1_native}/${subject}.${hemi}.white_thr0.native.nii.gz
        fslmaths \
            ${t1_native}/${subject}.${hemi}.white_thr0.native.nii.gz -bin \
            ${t1_native}/${subject}.${hemi}.white_thr0.native.nii.gz
        fslmaths \
            ${t1_native}/${subject}.${hemi}.pial.native.nii.gz -uthr 0 -abs -bin -mul 255 \
            ${t1_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        fslmaths \
            ${t1_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz -bin \
            ${t1_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        fslmaths \
            ${t1_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz -mas \
            ${t1_native}/${subject}.${hemi}.white_thr0.native.nii.gz -mul 255 \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz

        # add white matter mask to ribbon file
        fslmaths \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz -bin -mul $GreyRibbonValue \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz
        fslmaths \
            ${t1_native}/${subject}.${hemi}.white.native.nii.gz -uthr 0 -abs -bin -mul 255 \
            ${t1_native}/${subject}.${hemi}.white_uthr0.native.nii.gz
        fslmaths \
            ${t1_native}/${subject}.${hemi}.white_uthr0.native.nii.gz -bin \
            ${t1_native}/${subject}.${hemi}.white_uthr0.native.nii.gz
        fslmaths \
            ${t1_native}/${subject}.${hemi}.white_uthr0.native.nii.gz -mul $WhiteMaskValue \
            ${t1_native}/${subject}.${hemi}.white_mask.native.nii.gz
        fslmaths \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz -add \
            ${t1_native}/${subject}.${hemi}.white_mask.native.nii.gz \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz

        # delete temporary files
        rm ${t1_native}/${subject}.${hemi}.white.native.nii.gz
        rm ${t1_native}/${subject}.${hemi}.white_thr0.native.nii.gz
        rm ${t1_native}/${subject}.${hemi}.pial.native.nii.gz
        rm ${t1_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        rm ${t1_native}/${subject}.${hemi}.white_uthr0.native.nii.gz
        rm ${t1_native}/${subject}.${hemi}.white_mask.native.nii.gz
    done

    # make joined native space ribbon file (across hemispheres)
    fslmaths \
        ${t1_native}/${subject}.L.ribbon.nii.gz -add \
        ${t1_native}/${subject}.R.ribbon.nii.gz \
        ${t1}/ribbon.nii.gz

    rm ${t1_native}/${subject}.L.ribbon.nii.gz
    rm ${t1_native}/${subject}.R.ribbon.nii.gz

    wb_command -volume-label-import \
        ${t1}/ribbon.nii.gz ${fs_labels} ${t1}/ribbon.nii.gz -drop-unused-labels
fi

# create mni space ribbon
if [ ! -f ${mni}/ribbon.nii.gz ]; then
    for hemi in L R ; do

        if [ ${hemi} = "L" ]; then
          gm_val=${l_gm_val}
          wm_val=${l_wm_val}
        elif [ ${hemi} = "R" ]; then
          gm_val=${r_gm_val}
          wm_val=${r_wm_val}
        fi

        wb_command -create-signed-distance-volume \
            ${mni_native}/${subject}.${hemi}.white.native.surf.gii \
            ${mni}/T1w.nii.gz \
            ${mni_native}/${subject}.${hemi}.white.native.nii.gz
        wb_command -create-signed-distance-volume \
            ${mni_native}/${subject}.${hemi}.pial.native.surf.gii \
            ${mni}/T1w.nii.gz \
            ${mni_native}/${subject}.${hemi}.pial.native.nii.gz

        # find grey matter mask (ribbon)
        fslmaths \
            ${mni_native}/${subject}.${hemi}.white.native.nii.gz -thr 0 -bin -mul 255 \
            ${mni_native}/${subject}.${hemi}.white_thr0.native.nii.gz
        fslmaths \
            ${mni_native}/${subject}.${hemi}.white_thr0.native.nii.gz -bin \
            ${mni_native}/${subject}.${hemi}.white_thr0.native.nii.gz
        fslmaths \
            ${mni_native}/${subject}.${hemi}.pial.native.nii.gz -uthr 0 -abs -bin -mul 255 \
            ${mni_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        fslmaths \
            ${mni_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz -bin \
            ${mni_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        fslmaths \
            ${mni_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz -mas \
            ${mni_native}/${subject}.${hemi}.white_thr0.native.nii.gz -mul 255 \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz

        # add white matter mask to ribbon file
        fslmaths \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz -bin -mul ${gm_val} \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz
        fslmaths \
            ${mni_native}/${subject}.${hemi}.white.native.nii.gz -uthr 0 -abs -bin -mul 255 \
            ${mni_native}/${subject}.${hemi}.white_uthr0.native.nii.gz
        fslmaths \
            ${mni_native}/${subject}.${hemi}.white_uthr0.native.nii.gz -bin \
            ${mni_native}/${subject}.${hemi}.white_uthr0.native.nii.gz
        fslmaths \
            ${mni_native}/${subject}.${hemi}.white_uthr0.native.nii.gz -mul ${wm_val} \
            ${mni_native}/${subject}.${hemi}.white_mask.native.nii.gz
        fslmaths \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz -add \
            ${mni_native}/${subject}.${hemi}.white_mask.native.nii.gz \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz

        # delete temporary files
        rm ${mni_native}/${subject}.${hemi}.white.native.nii.gz
        rm ${mni_native}/${subject}.${hemi}.white_thr0.native.nii.gz
        rm ${mni_native}/${subject}.${hemi}.pial.native.nii.gz
        rm ${mni_native}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        rm ${mni_native}/${subject}.${hemi}.white_uthr0.native.nii.gz
        rm${mni_native}/${subject}.${hemi}.white_mask.native.nii.gz
    done

    # make joined native space ribbon file (across hemispheres)
    fslmaths \
        ${mni_native}/${subject}.L.ribbon.nii.gz -add \
        ${mni_native}/${subject}.R.ribbon.nii.gz \
        ${mni}/ribbon.nii.gz
    rm ${mni_native}/${subject}.L.ribbon.nii.gz
    rm ${mni_native}/${subject}.R.ribbon.nii.gz

    wb_command -volume-label-import \
        ${mni}/ribbon.nii.gz ${fs_labels} ${mni}/ribbon.nii.gz -drop-unused-labels
fi

if [ ! -f ${t1}/T2_to_fs.mat ]; then
    bbregister \
        --s PRE01_SFD_32010_01 \
        --mov PRE01_SFD_32010_01_01_T2_09_Sag-CUBE-T2-256x256x25.6x1.0.nii.gz \
        --reg ${t1}/T2_to_fs.dat \
        --fslmat ${t1}/T2_to_fs.mat \
        --init-fsl \
        --t2
fi

if [ ! -f ${t1}/T1_to_fs.mat ]; then
    bbregister \
        --s PRE01_SFD_32010_01 \
        --mov PRE01_SFD_32010_01_01_T1_03_Sag-T1-BRAVO-256x256x25.5x1.0.nii.gz \
        --reg ${t1}/T1_to_fs.dat \
        --fslmat ${t1}/T1_to_fs.mat \
        --init-fsl \
        --t1
fi

if [ ! -f ${t1}/T2_to_fs.nii.gz ]; then
    flirt \
        -in PRE01_SFD_32010_01_01_T2_09_Sag-CUBE-T2-256x256x25.6x1.0.nii.gz \
        -ref ${t1}/T1w.nii.gz \
        -applyxfm -init ${t1}/T2_to_fs.mat \
        -interp spline \
        -out ${t1}/T2_to_fs.nii.gz
fi

if [ ! -f ${t1}/T1_to_fs.nii.gz ]; then
    flirt \
        -in PRE01_SFD_32010_01_01_T1_03_Sag-T1-BRAVO-256x256x25.5x1.0.nii.gz \
        -ref ${t1}/T1w.nii.gz \
        -applyxfm -init ${t1}/T1_to_fs.mat \
        -interp spline \
        -out ${t1}/T1_to_fs.nii.gz
fi

# calculate bias field in freesurfer space
if [ ! -f ${t1}/T1_bias.nii.gz ]; then
    # Form sqrt(T1w*T2w), mask this and normalise by the mean
    fslmaths \
        T1_to_fs.nii.gz \
        -mul T2_to_fs.nii.gz \
        -abs -sqrt T1wmulT2w.nii.gz \
        -odt float
    fslmaths \
        T1wmulT2w.nii.gz \
        -mas ${t1}/brainmask_fs.nii.gz \
        T1wmulT2w_brain.nii.gz
    mean=$(fslstats T1wmulT2w_brain.nii.gz -M)
    fslmaths \
        T1wmulT2w_brain.nii.gz \
        -div ${mean} \
        T1wmulT2w_brain_norm.nii.gz

    # Smooth the normalised sqrt image, using within-mask smoothing : s(Mask*X)/s(Mask)
    fslmaths \
        T1wmulT2w_brain_norm.nii.gz \
        -bin -s ${bias_sigma} SmoothNorm_s${bias_sigma}.nii.gz
    fslmaths \
        T1wmulT2w_brain_norm.nii.gz \
        -s ${bias_sigma} -div SmoothNorm_s${bias_sigma}.nii.gz \
        T1wmulT2w_brain_norm_s${bias_sigma}.nii.gz

    # Divide normalised sqrt image by smoothed version (to do simple bias correction)
    fslmaths \
        T1wmulT2w_brain_norm.nii.gz \
        -div T1wmulT2w_brain_norm_s${bias_sigma}.nii.gz \
        T1wmulT2w_brain_norm_modulate.nii.gz

    # Create a mask using a threshold at Mean - 0.5*Stddev, with filling of holes to remove any non-grey/white tissue.
    std=$(fslstats T1wmulT2w_brain_norm_modulate.nii.gz -S)
    mean=$(fslstats T1wmulT2w_brain_norm_modulate.nii.gz -M)
    lower=$(echo "${mean} - (${std} * 0.5)" | bc -l)

    fslmaths \
        T1wmulT2w_brain_norm_modulate \
        -thr ${lower} -bin -ero -mul 255 T1wmulT2w_brain_norm_modulate_mask
    wb_command -volume-remove-islands \
        T1wmulT2w_brain_norm_modulate_mask.nii.gz \
        T1wmulT2w_brain_norm_modulate_mask.nii.gz

    # Extrapolate normalised sqrt image from mask region out to whole FOV
    fslmaths \
        T1wmulT2w_brain_norm.nii.gz \
        -mas T1wmulT2w_brain_norm_modulate_mask.nii.gz \
        -dilall bias_raw.nii.gz -odt float
    fslmaths \
        bias_raw.nii.gz -s ${bias_sigma} ${t1}/T1_bias.nii.gz
fi

# Use bias field output to create corrected images
if [ ! -f ${t1}/T1_to_fs_corrected_masked.nii.gz ]; then
    fslmaths \
        ${t1}/T1_to_fs.nii.gz \
        -div ${t1}/T1_bias.nii.gz \
        -mas ${t1}/brainmask_fs.nii.gz \
        ${t1}/T1_to_fs_corrected_masked.nii.gz -odt float
fi

if [ ! -f ${t1}/T1_to_fs_corrected.nii.gz ]; then
    fslmaths \
        ${t1}/T1_to_fs.nii.gz \
        -div ${t1}/T1_bias.nii.gz \
        ${t1}/T1_to_fs_corrected.nii.gz -odt float
fi

if [ ! -f ${t1}/T2_to_fs_corrected_masked.nii.gz ]; then
    fslmaths \
        ${t1}/T2_to_fs.nii.gz \
        -div ${t1}/T1_bias.nii.gz \
        -mas ${t1}/brainmask_fs.nii.gz \
        ${t1}/T2_to_fs_corrected_masked.nii.gz -odt float
fi

if [ ! -f ${t1}/T2_to_fs_corrected.nii.gz ]; then
    fslmaths \
        ${t1}/T2_to_fs.nii.gz \
        -div ${t1}/T1_bias.nii.gz \
        ${t1}/T2_to_fs_corrected.nii.gz -odt float
fi

# warp bias field to various spaces
convertwarp --relout --rel --ref="$T1wImageBrain" --premat="$InitialT1wTransform" --warp1="$dcT1wTransform" --out="$OutputOrigT1wToT1w"
convertwarp --relout --rel --ref="$T1wImageBrain" --warp1="$OutputOrigT1wToT1w" --warp2="$AtlasTransform" --out="$OutputOrigT1wToStandard"

convertwarp --relout --rel --ref="$T1wImageBrain" --premat="$InitialT2wTransform" --warp1="$dcT2wTransform" --postmat="$FinalT2wTransform" --out="$OutputOrigT2wToT1w"
convertwarp --relout --rel --ref="$T1wImageBrain" --warp1="$OutputOrigT2wToT1w" --warp2="$AtlasTransform" --out="$OutputOrigT2wToStandard"

applywarp --rel --interp=spline -i "$OrginalT1wImage" -r "$T1wImageBrain" -w "$OutputOrigT1wToT1w" -o "$OutputT1wImage"
fslmaths "$OutputT1wImage" -abs "$OutputT1wImage" -odt float
fslmaths "$OutputT1wImage" -div "$BiasField" "$OutputT1wImageRestore"
fslmaths "$OutputT1wImageRestore" -mas "$T1wImageBrain" "$OutputT1wImageRestoreBrain"

applywarp --rel --interp=spline -i "$BiasField" -r "$T1wImageBrain" -w "$AtlasTransform" -o "$BiasFieldOutput"
fslmaths "$BiasFieldOutput" -thr 0.1 "$BiasFieldOutput"

applywarp --rel --interp=spline -i "$OrginalT1wImage" -r "$T1wImageBrain" -w "$OutputOrigT1wToStandard" -o "$OutputMNIT1wImage"
fslmaths "$OutputMNIT1wImage" -abs "$OutputMNIT1wImage" -odt float
fslmaths "$OutputMNIT1wImage" -div "$BiasFieldOutput" "$OutputMNIT1wImageRestore"
fslmaths "$OutputMNIT1wImageRestore" -mas "$T1wMNIImageBrain" "$OutputMNIT1wImageRestoreBrain"

applywarp --rel --interp=spline -i "$OrginalT2wImage" -r "$T1wImageBrain" -w "$OutputOrigT2wToT1w" -o "$OutputT2wImage"
fslmaths "$OutputT2wImage" -abs "$OutputT2wImage" -odt float
fslmaths "$OutputT2wImage" -div "$BiasField" "$OutputT2wImageRestore"
fslmaths "$OutputT2wImageRestore" -mas "$T1wImageBrain" "$OutputT2wImageRestoreBrain"

applywarp --rel --interp=spline -i "$OrginalT2wImage" -r "$T1wImageBrain" -w "$OutputOrigT2wToStandard" -o "$OutputMNIT2wImage"
fslmaths "$OutputMNIT2wImage" -abs "$OutputMNIT2wImage" -odt float
fslmaths "$OutputMNIT2wImage" -div "$BiasFieldOutput" "$OutputMNIT2wImageRestore"
fslmaths "$OutputMNIT2wImageRestore" -mas "$T1wMNIImageBrain" "$OutputMNIT2wImageRestoreBrain"

# create myelin map
if [ ! -f ${t1}/myelin_map.nii.gz ]; then
    # divide images
    wb_command -volume-math \
        "clamp((T1w / T2w), 0, 100)" \
        myelin_map.nii.gz \
        -var T1w ${t1}/T1_to_fs_corrected.nii.gz \
        -var T2w ${t1}/T2_to_fs_corrected.nii.gz \
        -fixnan 0

    # set color-map
    wb_command -volume-palette \
        ${t1}/myelin_map.nii.gz \
        MODE_AUTO_SCALE_PERCENTAGE \
        -pos-percent 4 96 \
        -interpolate true \
        -palette-name videen_style \
        -disp-pos true \
        -disp-neg false \
        -disp-zero false

    # add to spec file
    wb_command -add-to-spec-file \
        ${t1_native}/${subject}.native.wb.spec INVALID ${t1}/myelin_map.nii.gz
fi

# create ribbon myelin map
if [ ! -f ${t1}/myelin_map_ribbon.nii.gz ]; then
    # divide images
    wb_command -volume-math \
        "(T1w / T2w) * (((ribbon > (${l_gm_val}-0.01)) * (ribbon < (${l_gm_val}+0.01))) + ((ribbon > (${r_gm_val}-0.01)) * (ribbon < (${r_gm_val}+0.01))))" \
        ${t1}/myelin_map_ribbon.nii.gz \
        -var T1w ${t1}/T1_to_fs.nii.gz \
        -var T2w ${t1}/T2_to_fs.nii.gz \
        -var ribbon ${t1}/ribbon.nii.gz

    # set color map
    wb_command -volume-palette
        ${t1}/myelin_map_ribbon.nii.gz \
        MODE_AUTO_SCALE_PERCENTAGE \
        -pos-percent 4 96 \
        -interpolate true \
        -palette-name videen_style \
        -disp-pos true \
        -disp-neg false \
        -disp-zero false

    # add to spec file
    wb_command -add-to-spec-file \
        ${t1_native}/${subject}.native.wb.spec \
        INVALID \
        ${t1}/myelin_map_ribbon.nii.gz
fi

# do surface stuff???
    for hemi in L R ; do
    if [ ${hemi} = "L" ] ; then
        structure="CORTEX_LEFT"
        ribbon=${l_gm_val}
    elif [ ${hemi} = "R" ] ; then
        structure="CORTEX_RIGHT"
        ribbon=${r_gm_val}
    fi

    sphere="${mni_native}/${subject}.${hemi}.sphere.reg.reg_LR.native.surf.gii"

    # make ribbon
    wb_command -volume-math \
        "(ribbon > (${ribbon} - 0.01)) * (ribbon < (${ribbon} + 0.01))" \
        ${t1}/temp_ribbon.nii.gz \
        -var ribbon ${t1}/ribbon.nii.gz
    wb_command -volume-to-surface-mapping \
        ${t1}/myelin_map.nii.gz \
        ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
        ${mni_native}/${subject}.${hemi}.MyelinMap.native.func.gii \
        -myelin-style \
        ${t1}/temp_ribbon.nii.gz \
        ${mni_native}/${subject}.${hemi}.thickness.native.shape.gii \
        ${myelin_sigma}

    rm ${t1}/temp_ribbon.nii.gz

    # regress curvature out of thickness from within roi.native ?
    wb_command -metric-regression \
        ${mni_native}/${subject}.${hemi}.thickness.native.shape.gii \
        ${mni_native}/${subject}.${hemi}.corrThickness.native.shape.gii\
        -roi ${mni_native}/${subject}.${hemi}.roi.native.shape.gii \
        -remove ${mni_native}/${subject}.${hemi}.curvature.native.shape.gii

    # smooth myelin map
    wb_command -metric-smoothing \
        ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
        ${mni_native}/${subject}.${hemi}.MyelinMap.native.func.gii \
        ${surface_sigma} \
        ${mni_native}/${subject}.${hemi}.SmoothedMyelinMap.native.func.gii \
        -roi ${mni_native}/${subject}.${hemi}.roi.native.shape.gii

    # !!!
    # no idea where RefMyelinMap comes from, so we're just doing the MyelinMap -- jdv
    # !!!

    # downsample reference mesh (where does this come from?) to fs
    #wb_command -metric-resample "$AtlasSpaceFolder"/${subject}.${hemi}.RefMyelinMap."$HighResMesh"k_fs_LR.func.gii "$AtlasSpaceFolder"/${subject}.${hemi}.sphere."$HighResMesh"k_fs_LR.surf.gii ${sphere} ADAP_BARY_AREA "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.RefMyelinMap.native.func.gii -area-surfs "$AtlasSpaceFolder"/${subject}.${hemi}.midthickness."$HighResMesh"k_fs_LR.surf.gii ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii -current-roi "$AtlasSpaceFolder"/${subject}.${hemi}.atlasroi."$HighResMesh"k_fs_LR.shape.gii

    # dilate, then mask myelin map
    #wb_command -metric-dilate "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.RefMyelinMap.native.func.gii ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii 30 "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.RefMyelinMap.native.func.gii -nearest
    #wb_command -metric-mask "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.RefMyelinMap.native.func.gii "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.roi.native.shape.gii "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.RefMyelinMap.native.func.gii
    
    #LowResMesh=`echo ${LowResMeshes} | cut -d " " -f 1`
    #for Map in MyelinMap RefMyelinMap ; do

    # reduce memory usage by smoothing on downsampled mesh
    map='MyelinMap' # can add RefMyelinMap later if needed
    wb_command -metric-resample \
        ${mni_native}.${subject}.${hemi}.${map}.native.func.gii \
        ${sphere} \
        ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.sphere.${low_res_mesh}k_fs_LR.surf.gii \
        ADAP_BARY_AREA \
        ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}.${low_res_mesh}k_fs_LR.func.gii \
        -area-surfs ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
        ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.midthickness.${low_res_mesh}k_fs_LR.surf.gii \
        -current-roi ${mni_native}.${subject}.${hemi}.roi.native.shape.gii

    wb_command -metric-smoothing \
        ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.midthickness.${low_res_mesh}k_fs_LR.surf.gii \
        ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}.${low_res_mesh}k_fs_LR.func.gii \
        ${correction_sigma} \
        ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}_s${correction_sigma}.${low_res_mesh}k_fs_LR.func.gii \
        -roi ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.atlasroi.${low_res_mesh}k_fs_LR.shape.gii

    wb_command -metric-resample \
        ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}_s${correction_sigma}.${low_res_mesh}k_fs_LR.func.gii \
        ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.sphere.${low_res_mesh}k_fs_LR.surf.gii \
        ${sphere} \
        ADAP_BARY_AREA \
        ${mni_native}.${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii \
        -area-surfs ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.midthickness.${low_res_mesh}k_fs_LR.surf.gii \
        ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
        -current-roi ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.atlasroi.${low_res_mesh}k_fs_LR.shape.gii

    wb_command -metric-dilate \
        ${mni_native}.${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii \
        ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
        30 \
        ${mni_native}.${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii \
        -nearest

    wb_command -metric-mask \
        ${mni_native}.${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii \
        ${mni_native}.${subject}.${hemi}.roi.native.shape.gii \
        ${mni_native}.${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii

    rm ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}.${low_res_mesh}k_fs_LR.func.gii 
    rm ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}_s${correction_sigma}.${low_res_mesh}k_fs_LR.func.gii
    
    #done

    # these are both commented out in the original code -- jdv
    ##wb_command -metric-smoothing ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii ${mni_native}.${hemi}.RefMyelinMap.native.func.gii ${correction_sigma} ${mni_native}.${hemi}.RefMyelinMap_s${correction_sigma}.native.func.gii -roi ${mni_native}.${hemi}.roi.native.shape.gii 
    #wb_command -metric-smoothing \
    #    ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
    #    ${mni_native}.${subject}${hemi}.MyelinMap.native.func.gii \
    #    ${correction_sigma} \
    #    ${mni_native}.${subject}.${hemi}.MyelinMap_s${correction_sigma}.native.func.gii \
    #    -roi ${mni_native}/${subject}.${hemi}.roi.native.shape.gii
    
    # compute difference maps (I don't think we need these).
    #wb_command -metric-math "(Individual - Reference) * Mask" "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.BiasField.native.func.gii -var Individual "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.MyelinMap_s${correction_sigma}.native.func.gii -var Reference "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.RefMyelinMap_s${correction_sigma}.native.func.gii -var Mask "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.roi.native.shape.gii
    #wb_command -metric-math "(Individual - Bias) * Mask" "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.MyelinMap_BC.native.func.gii -var Individual "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.MyelinMap.native.func.gii -var Bias "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.BiasField.native.func.gii -var Mask "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.roi.native.shape.gii
    #wb_command -metric-math "(Individual - Bias) * Mask" "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.SmoothedMyelinMap_BC.native.func.gii -var Individual "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.SmoothedMyelinMap.native.func.gii -var Bias "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.BiasField.native.func.gii -var Mask "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.roi.native.shape.gii
    #rm "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.MyelinMap_s${correction_sigma}.native.func.gii "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.RefMyelinMap_s${correction_sigma}.native.func.gii
    #for STRING in MyelinMap@func SmoothedMyelinMap@func MyelinMap_BC@func SmoothedMyelinMap_BC@func corrThickness@shape ; do
    #    Map=`echo $STRING | cut -d "@" -f 1`
    #    Ext=`echo $STRING | cut -d "@" -f 2`
    #    wb_command -set-map-name "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}."$Map".native."$Ext".gii 1 ${subject}_${hemi}_"$Map"
    #    wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}."$Map".native."$Ext".gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
    #    wb_command -metric-resample "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}."$Map".native."$Ext".gii ${sphere} "$AtlasSpaceFolder"/${subject}.${hemi}.sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/${subject}.${hemi}."$Map"."$HighResMesh"k_fs_LR."$Ext".gii -area-surfs ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii "$AtlasSpaceFolder"/${subject}.${hemi}.midthickness."$HighResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.roi.native.shape.gii
    #    wb_command -metric-mask "$AtlasSpaceFolder"/${subject}.${hemi}."$Map"."$HighResMesh"k_fs_LR."$Ext".gii "$AtlasSpaceFolder"/${subject}.${hemi}.atlasroi."$HighResMesh"k_fs_LR.shape.gii "$AtlasSpaceFolder"/${subject}.${hemi}."$Map"."$HighResMesh"k_fs_LR."$Ext".gii
    #    for LowResMesh in ${LowResMeshes} ; do
    #        wb_command -metric-resample "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}."$Map".native."$Ext".gii ${sphere} "$AtlasSpaceFolder"/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.sphere.${low_res_mesh}k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}."$Map".${low_res_mesh}k_fs_LR."$Ext".gii -area-surfs ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.midthickness.${low_res_mesh}k_fs_LR.surf.gii -current-roi "$AtlasSpaceFolder"/"$NativeFolder"/${subject}.${hemi}.roi.native.shape.gii
    #        wb_command -metric-mask "$AtlasSpaceFolder"/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}."$Map".${low_res_mesh}k_fs_LR."$Ext".gii "$AtlasSpaceFolder"/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.atlasroi.${low_res_mesh}k_fs_LR.shape.gii "$AtlasSpaceFolder"/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}."$Map".${low_res_mesh}k_fs_LR."$Ext".gii
    #    done
    #done

    wb_command -cifti-create-dense-scalar \
        ${mni_native}/${subject}.MyelinMap.native.dscalar.nii \
        -left-metric ${mni_native}/${subject}.L.MyelinMap.native.func.gii \
        -roi-left ${mni_native}/${subject}.L."$ROI"."$Mesh".shape.gii -right-metric "$Folder"/${subject}.R.${Map}."$Mesh"."$Ext".gii -roi-right "$Folder"/${subject}.R."$ROI"."$Mesh".shape.gii
    wb_command -set-map-names "$Folder"/${subject}.${Map}."$Mesh".dscalar.nii -map 1 "${Subject}_${Map}"
    wb_command -cifti-palette "$Folder"/${subject}.${Map}."$Mesh".dscalar.nii MODE_AUTO_SCALE_PERCENTAGE "$Folder"/${subject}.${Map}."$Mesh".dscalar.nii -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false


done




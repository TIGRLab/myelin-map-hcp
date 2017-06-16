#!/bin/bash -ex
#
# for details see The Minimal Preprocessing Pipelines for the Human Connectome Project Glasser et al.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3720813/
#

hcp='/scratch/jviviano/prelapse_hcp'
assets='/projects/jviviano/code/myelin-map-hcp/assets'
SUBJECTS_DIR='/archive/data/PRELAPSE/pipelines/freesurfer'

# inputs
subject="PRE01_SFD_32002_01"
t1_data="/archive/data/PRELAPSE/data/nii/${subject}/${subject}_??_T1_*.nii.gz"
t2_data="/archive/data/PRELAPSE/data/nii/${subject}/${subject}_??_T2_*.nii.gz"

fwhm_myelin=8
fwhm_surface=20

# shortcuts
t1="${hcp}/${subject}/T1w"
t1_native="${t1}/Native"
mni="${hcp}/${subject}/MNINonLinear"
mni_native="${mni}/Native"
tmp="/tmp/myelin_map_${subject}"

# file from ciftify... (text)
fs_labels="${assets}/lut/FreeSurferAllLut.txt"
reference="${assets}/templates/standard_mesh_atlases/Conte69.MyelinMap_BC.164k_fs_LR.dscalar.nii"

# constants
l_gm_val="3"
l_wm_val="2"
r_gm_val="42"
r_wm_val="41"
high_res_mesh="164" # for fsaverage_LR164k
low_res_mesh="32"   # for fsaverage_LR32k
myelin_sigma=$(echo "${fwhm_myelin} / (2*(sqrt(2*l(2))))" | bc -l)
surface_sigma=$(echo "${fwhm_surface} / (2*(sqrt(2*l(2))))" | bc -l)
correction_sigma="5"
bias_sigma="5"

# set up working directory
if [ -d "${tmp}" ]; then
    rm -r ${tmp}
fi
mkdir -p ${tmp}

# get raw data into freesurfer, then HCP space
if [ ! -f ${t1}/T2_to_fs.mat ]; then
    bbregister \
        --s ${subject} \
        --mov ${t2_data} \
        --reg ${t1}/T2_to_fs.dat \
        --fslmat ${t1}/T2_to_fs.mat \
        --init-fsl \
        --t2
fi

if [ ! -f ${t1}/T1_to_fs.mat ]; then
    bbregister \
        --s ${subject} \
        --mov ${t1_data} \
        --reg ${t1}/T1_to_fs.dat \
        --fslmat ${t1}/T1_to_fs.mat \
        --init-fsl \
        --t1
fi

# get temporary reference T1 from freesurfer
if [ ! -f ${t1}/freesurfer_ref.nii.gz ]; then
    mri_convert ${SUBJECTS_DIR}/${subject}/mri/T1.mgz ${t1}/freesurfer_ref.nii.gz
fi


if [ ! -f ${t1}/T2_to_fs.nii.gz ]; then
    flirt \
        -in ${t2_data} \
        -ref ${t1}/freesurfer_ref.nii.gz \
        -applyxfm -init ${t1}/T2_to_fs.mat \
        -interp spline \
        -out ${t1}/T2_to_fs.nii.gz
fi

if [ ! -f ${t1}/T1_to_fs.nii.gz ]; then
    flirt \
        -in ${t1_data} \
        -ref ${t1}/freesurfer_ref.nii.gz \
        -applyxfm -init ${t1}/T1_to_fs.mat \
        -interp spline \
        -out ${t1}/T1_to_fs.nii.gz
fi

# reorients freesurfer data to HCP data
if [ ! -f ${t1}/T1_to_hcp.nii.gz ]; then
    fslreorient2std ${t1}/T1_to_fs.nii.gz ${t1}/T1_to_hcp.nii.gz
fi

if [ ! -f ${t1}/T2_to_hcp.nii.gz ]; then
    fslreorient2std ${t1}/T2_to_fs.nii.gz ${t1}/T2_to_hcp.nii.gz
fi

if [ ! -f ${t1}/brainmask_hcp.nii.gz ]; then
    fslreorient2std ${t1}/brainmask_fs.nii.gz ${t1}/brainmask_hcp.nii.gz
fi

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
            ${tmp}/${subject}.${hemi}.white.native.nii.gz
        wb_command -create-signed-distance-volume \
            ${t1_native}/${subject}.${hemi}.pial.native.surf.gii \
            ${t1}/T1w.nii.gz \
            ${tmp}/${subject}.${hemi}.pial.native.nii.gz

        # find grey matter mask (ribbon)
        fslmaths \
            ${tmp}/${subject}.${hemi}.white.native.nii.gz -thr 0 -bin -mul 255 \
            ${tmp}/${subject}.${hemi}.white_thr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.white_thr0.native.nii.gz -bin \
            ${tmp}/${subject}.${hemi}.white_thr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.pial.native.nii.gz -uthr 0 -abs -bin -mul 255 \
            ${tmp}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.pial_uthr0.native.nii.gz -bin \
            ${tmp}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.pial_uthr0.native.nii.gz -mas \
            ${tmp}/${subject}.${hemi}.white_thr0.native.nii.gz -mul 255 \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz

        # add white matter mask to ribbon file
        fslmaths \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz -bin -mul ${gm_val} \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.white.native.nii.gz -uthr 0 -abs -bin -mul 255 \
            ${tmp}/${subject}.${hemi}.white_uthr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.white_uthr0.native.nii.gz -bin \
            ${tmp}/${subject}.${hemi}.white_uthr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.white_uthr0.native.nii.gz -mul ${wm_val} \
            ${tmp}/${subject}.${hemi}.white_mask.native.nii.gz
        fslmaths \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz -add \
            ${tmp}/${subject}.${hemi}.white_mask.native.nii.gz \
            ${t1_native}/${subject}.${hemi}.ribbon.nii.gz

        rm ${tmp}/${subject}.${hemi}.*.nii.gz
    done

    # make joined native space ribbon file (across hemispheres)
    fslmaths \
        ${t1_native}/${subject}.L.ribbon.nii.gz -add \
        ${t1_native}/${subject}.R.ribbon.nii.gz \
        ${t1}/ribbon.nii.gz

    rm ${t1_native}/${subject}.?.ribbon.nii.gz

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
            ${tmp}/${subject}.${hemi}.white.native.nii.gz
        wb_command -create-signed-distance-volume \
            ${mni_native}/${subject}.${hemi}.pial.native.surf.gii \
            ${mni}/T1w.nii.gz \
            ${tmp}/${subject}.${hemi}.pial.native.nii.gz

        # find grey matter mask (ribbon)
        fslmaths \
            ${tmp}/${subject}.${hemi}.white.native.nii.gz -thr 0 -bin -mul 255 \
            ${tmp}/${subject}.${hemi}.white_thr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.white_thr0.native.nii.gz -bin \
            ${tmp}/${subject}.${hemi}.white_thr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.pial.native.nii.gz -uthr 0 -abs -bin -mul 255 \
            ${tmp}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.pial_uthr0.native.nii.gz -bin \
            ${tmp}/${subject}.${hemi}.pial_uthr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.pial_uthr0.native.nii.gz -mas \
            ${tmp}/${subject}.${hemi}.white_thr0.native.nii.gz -mul 255 \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz

        # add white matter mask to ribbon file
        fslmaths \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz -bin -mul ${gm_val} \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.white.native.nii.gz -uthr 0 -abs -bin -mul 255 \
            ${tmp}/${subject}.${hemi}.white_uthr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.white_uthr0.native.nii.gz -bin \
            ${tmp}/${subject}.${hemi}.white_uthr0.native.nii.gz
        fslmaths \
            ${tmp}/${subject}.${hemi}.white_uthr0.native.nii.gz -mul ${wm_val} \
            ${tmp}/${subject}.${hemi}.white_mask.native.nii.gz
        fslmaths \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz -add \
            ${tmp}/${subject}.${hemi}.white_mask.native.nii.gz \
            ${mni_native}/${subject}.${hemi}.ribbon.nii.gz

        rm ${tmp}/${subject}.${hemi}.*.nii.gz
    done

    # make joined native space ribbon file (across hemispheres)
    fslmaths \
        ${mni_native}/${subject}.L.ribbon.nii.gz -add \
        ${mni_native}/${subject}.R.ribbon.nii.gz \
        ${mni}/ribbon.nii.gz

    rm ${mni_native}/${subject}.?.ribbon.nii.gz

    wb_command -volume-label-import \
        ${mni}/ribbon.nii.gz ${fs_labels} ${mni}/ribbon.nii.gz -drop-unused-labels
fi

# calculate bias field in hcp space
if [ ! -f ${t1}/T1_bias.nii.gz ]; then

    # calc sqrt(T1w*T2w), mask, and normalise by in-brain mean
    fslmaths \
        ${t1}/T1_to_hcp.nii.gz \
        -mul ${t1}/T2_to_hcp.nii.gz \
        -abs -sqrt ${tmp}/T1wmulT2w.nii.gz \
        -odt float

    fslmaths \
        ${tmp}/T1wmulT2w.nii.gz \
        -mas ${t1}/brainmask_hcp.nii.gz \
        ${tmp}/T1wmulT2w_brain.nii.gz

    mean=$(fslstats ${tmp}/T1wmulT2w_brain.nii.gz -M)
    fslmaths \
        ${tmp}/T1wmulT2w_brain.nii.gz \
        -div ${mean} \
        ${tmp}/T1wmulT2w_brain_norm.nii.gz

    # smooth the normalised sqrt image within mask: s(Mask*X)/s(Mask)
    fslmaths \
        ${tmp}/T1wmulT2w_brain_norm.nii.gz \
        -bin -s ${bias_sigma} ${tmp}/SmoothNorm_s${bias_sigma}.nii.gz
    fslmaths \
        ${tmp}/T1wmulT2w_brain_norm.nii.gz \
        -s ${bias_sigma} -div ${tmp}/SmoothNorm_s${bias_sigma}.nii.gz \
        ${tmp}/T1wmulT2w_brain_norm_s${bias_sigma}.nii.gz

    # divide normalised sqrt image by smoothed version (to do simple bias correction)
    fslmaths \
        ${tmp}/T1wmulT2w_brain_norm.nii.gz \
        -div ${tmp}/T1wmulT2w_brain_norm_s${bias_sigma}.nii.gz \
        ${tmp}/T1wmulT2w_brain_norm_modulate.nii.gz

    # create a mask using a threshold at Mean - 0.5*Stddev, with filling of holes to remove any non-grey/white tissue.
    std=$(fslstats ${tmp}/T1wmulT2w_brain_norm_modulate.nii.gz -S)
    mean=$(fslstats ${tmp}/T1wmulT2w_brain_norm_modulate.nii.gz -M)
    lower=$(echo "${mean} - (${std} * 0.5)" | bc -l)

    fslmaths \
        ${tmp}/T1wmulT2w_brain_norm_modulate.nii.gz \
        -thr ${lower} -bin -ero -mul 255 ${tmp}/T1wmulT2w_brain_norm_modulate_mask.nii.gz
    wb_command -volume-remove-islands \
        ${tmp}/T1wmulT2w_brain_norm_modulate_mask.nii.gz \
        ${tmp}/T1wmulT2w_brain_norm_modulate_mask.nii.gz

    # extrapolate normalised sqrt image from mask region out to whole FOV
    fslmaths \
        ${tmp}/T1wmulT2w_brain_norm.nii.gz \
        -mas ${tmp}/T1wmulT2w_brain_norm_modulate_mask.nii.gz \
        -dilall ${tmp}/bias_raw.nii.gz -odt float
    fslmaths \
        ${tmp}/bias_raw.nii.gz -s ${bias_sigma} ${t1}/T1_bias.nii.gz

    rm ${tmp}/*.nii.gz
fi

# Use bias field output to create corrected images
if [ ! -f ${t1}/T1_to_hcp_corrected_masked.nii.gz ]; then
    fslmaths \
        ${t1}/T1_to_hcp.nii.gz \
        -div ${t1}/T1_bias.nii.gz \
        -mas ${t1}/brainmask_hcp.nii.gz \
        ${t1}/T1_to_hcp_corrected_masked.nii.gz -odt float
fi

if [ ! -f ${t1}/T1_to_hcp_corrected.nii.gz ]; then
    fslmaths \
        ${t1}/T1_to_hcp.nii.gz \
        -div ${t1}/T1_bias.nii.gz \
        ${t1}/T1_to_hcp_corrected.nii.gz -odt float
fi

if [ ! -f ${t1}/T2_to_hcp_corrected_masked.nii.gz ]; then
    fslmaths \
        ${t1}/T2_to_hcp.nii.gz \
        -div ${t1}/T1_bias.nii.gz \
        -mas ${t1}/brainmask_hcp.nii.gz \
        ${t1}/T2_to_hcp_corrected_masked.nii.gz -odt float
fi

if [ ! -f ${t1}/T2_to_hcp_corrected.nii.gz ]; then
    fslmaths \
        ${t1}/T2_to_hcp.nii.gz \
        -div ${t1}/T1_bias.nii.gz \
        ${t1}/T2_to_hcp_corrected.nii.gz -odt float
fi

# warp bias field to various spaces
#convertwarp --relout --rel --ref="$T1wImageBrain" --premat="$InitialT1wTransform" --warp1="$dcT1wTransform" --out="$OutputOrigT1wToT1w"
#convertwarp --relout --rel --ref="$T1wImageBrain" --warp1="$OutputOrigT1wToT1w" --warp2="$AtlasTransform" --out="$OutputOrigT1wToStandard"

#convertwarp --relout --rel --ref="$T1wImageBrain" --premat="$InitialT2wTransform" --warp1="$dcT2wTransform" --postmat="$FinalT2wTransform" --out="$OutputOrigT2wToT1w"
#convertwarp --relout --rel --ref="$T1wImageBrain" --warp1="$OutputOrigT2wToT1w" --warp2="$AtlasTransform" --out="$OutputOrigT2wToStandard"

#applywarp --rel --interp=spline -i "$OrginalT1wImage" -r "$T1wImageBrain" -w "$OutputOrigT1wToT1w" -o "$OutputT1wImage"
#fslmaths "$OutputT1wImage" -abs "$OutputT1wImage" -odt float
#fslmaths "$OutputT1wImage" -div "$BiasField" "$OutputT1wImageRestore"
#fslmaths "$OutputT1wImageRestore" -mas "$T1wImageBrain" "$OutputT1wImageRestoreBrain"

#applywarp --rel --interp=spline -i "$BiasField" -r "$T1wImageBrain" -w "$AtlasTransform" -o "$BiasFieldOutput"
#fslmaths "$BiasFieldOutput" -thr 0.1 "$BiasFieldOutput"

#applywarp --rel --interp=spline -i "$OrginalT1wImage" -r "$T1wImageBrain" -w "$OutputOrigT1wToStandard" -o "$OutputMNIT1wImage"
#fslmaths "$OutputMNIT1wImage" -abs "$OutputMNIT1wImage" -odt float
#fslmaths "$OutputMNIT1wImage" -div "$BiasFieldOutput" "$OutputMNIT1wImageRestore"
#fslmaths "$OutputMNIT1wImageRestore" -mas "$T1wMNIImageBrain" "$OutputMNIT1wImageRestoreBrain"

#applywarp --rel --interp=spline -i "$OrginalT2wImage" -r "$T1wImageBrain" -w "$OutputOrigT2wToT1w" -o "$OutputT2wImage"
#fslmaths "$OutputT2wImage" -abs "$OutputT2wImage" -odt float
#fslmaths "$OutputT2wImage" -div "$BiasField" "$OutputT2wImageRestore"
#fslmaths "$OutputT2wImageRestore" -mas "$T1wImageBrain" "$OutputT2wImageRestoreBrain"

#applywarp --rel --interp=spline -i "$OrginalT2wImage" -r "$T1wImageBrain" -w "$OutputOrigT2wToStandard" -o "$OutputMNIT2wImage"
#fslmaths "$OutputMNIT2wImage" -abs "$OutputMNIT2wImage" -odt float
#fslmaths "$OutputMNIT2wImage" -div "$BiasFieldOutput" "$OutputMNIT2wImageRestore"
#fslmaths "$OutputMNIT2wImageRestore" -mas "$T1wMNIImageBrain" "$OutputMNIT2wImageRestoreBrain"

# create myelin map
if [ ! -f ${t1}/myelin_map.nii.gz ]; then
    # divide images
    wb_command -volume-math \
        "clamp((T1w / T2w), 0, 100)" \
        ${t1}/myelin_map.nii.gz \
        -var T1w ${t1}/T1_to_hcp_corrected_masked.nii.gz \
        -var T2w ${t1}/T2_to_hcp_corrected_masked.nii.gz \
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
        -var T1w ${t1}/T1_to_hcp_corrected_masked.nii.gz \
        -var T2w ${t1}/T2_to_hcp_corrected_masked.nii.gz \
        -var ribbon ${t1}/ribbon.nii.gz

    # set color map
    wb_command -volume-palette \
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

# seperate reference myelin map (from assets)
wb_command -cifti-separate-all \
    ${reference} \
    -left ${mni}/${subject}.L.RefMyelinMap.${high_res_mesh}k_fs_LR.func.gii \
    -right ${mni}/${subject}.R.RefMyelinMap.${high_res_mesh}k_fs_LR.func.gii


# do surface stuff
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
    if [ ! -f ${mni_native}/${subject}.${hemi}.MyelinMap.native.func.gii ]; then
        wb_command -volume-math \
            "(ribbon > (${ribbon} - 0.01)) * (ribbon < (${ribbon} + 0.01))" \
            ${tmp}/tmp_ribbon.nii.gz \
            -var ribbon ${t1}/ribbon.nii.gz

        wb_command -volume-to-surface-mapping \
            ${t1}/myelin_map.nii.gz \
            ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
            ${mni_native}/${subject}.${hemi}.MyelinMap.native.func.gii \
            -myelin-style \
            ${tmp}/tmp_ribbon.nii.gz \
            ${mni_native}/${subject}.${hemi}.thickness.native.shape.gii \
            ${myelin_sigma}

        rm ${tmp}/tmp_ribbon.nii.gz
    fi


    # regress curvature out of thickness from within roi.native ?
    if [ ! -f ${mni_native}/${subject}.${hemi}.RefMyelinMap.native.func.gii ]; then
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

        # downsample reference mesh to fs
        wb_command -metric-resample \
            ${mni}/${subject}.${hemi}.RefMyelinMap.${high_res_mesh}k_fs_LR.func.gii \
            ${mni}/${subject}.${hemi}.sphere.${high_res_mesh}k_fs_LR.surf.gii ${sphere} \
            ADAP_BARY_AREA \
            ${mni_native}/${subject}.${hemi}.RefMyelinMap.native.func.gii \
            -area-surfs ${mni}/${subject}.${hemi}.midthickness.${high_res_mesh}k_fs_LR.surf.gii \
            ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
            -current-roi ${mni}/${subject}.${hemi}.atlasroi.${high_res_mesh}k_fs_LR.shape.gii

        # dilate, then mask myelin map
        wb_command -metric-dilate \
            ${mni_native}/${subject}.${hemi}.RefMyelinMap.native.func.gii \
            ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
            30 \
            ${mni_native}/${subject}.${hemi}.RefMyelinMap.native.func.gii \
            -nearest

        wb_command -metric-mask \
            ${mni_native}/${subject}.${hemi}.RefMyelinMap.native.func.gii \
            ${mni_native}/${subject}.${hemi}.roi.native.shape.gii \
            ${mni_native}/${subject}.${hemi}.RefMyelinMap.native.func.gii
    fi

    for map in MyelinMap RefMyelinMap ; do

        # reduce memory usage by smoothing on downsampled mesh
        if [ ! -f ${mni_native}/${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii ]; then
            wb_command -metric-resample \
                ${mni_native}/${subject}.${hemi}.${map}.native.func.gii \
                ${sphere} \
                ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.sphere.${low_res_mesh}k_fs_LR.surf.gii \
                ADAP_BARY_AREA \
                ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}.${low_res_mesh}k_fs_LR.func.gii \
                -area-surfs ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
                ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.midthickness.${low_res_mesh}k_fs_LR.surf.gii \
                -current-roi ${mni_native}/${subject}.${hemi}.roi.native.shape.gii

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
                ${mni_native}/${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii \
                -area-surfs ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.midthickness.${low_res_mesh}k_fs_LR.surf.gii \
                ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
                -current-roi ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.atlasroi.${low_res_mesh}k_fs_LR.shape.gii

            wb_command -metric-dilate \
                ${mni_native}/${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii \
                ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
                30 \
                ${mni_native}/${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii \
                -nearest

            wb_command -metric-mask \
                ${mni_native}/${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii \
                ${mni_native}/${subject}.${hemi}.roi.native.shape.gii \
                ${mni_native}/${subject}.${hemi}.${map}_s${correction_sigma}.native.func.gii
        fi
    done

    # these are both commented out in the original code -- jdv
    ##wb_command -metric-smoothing ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii ${mni_native}/${hemi}.RefMyelinMap.native.func.gii ${correction_sigma} ${mni_native}/${hemi}.RefMyelinMap_s${correction_sigma}.native.func.gii -roi ${mni_native}.${hemi}.roi.native.shape.gii
    #wb_command -metric-smoothing \
    #    ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
    #    ${mni_native}.${subject}${hemi}.MyelinMap.native.func.gii \
    #    ${correction_sigma} \
    #    ${mni_native}.${subject}.${hemi}.MyelinMap_s${correction_sigma}.native.func.gii \
    #    -roi ${mni_native}/${subject}.${hemi}.roi.native.shape.gii

    # compute difference maps
    if [ ! -f ${mni_native}/${subject}.${hemi}.BiasField.native.func.gii ]; then
        wb_command -metric-math \
            "(Individual - Reference) * Mask" \
            ${mni_native}/${subject}.${hemi}.BiasField.native.func.gii \
            -var Individual ${mni_native}/${subject}.${hemi}.MyelinMap_s${correction_sigma}.native.func.gii \
            -var Reference ${mni_native}/${subject}.${hemi}.RefMyelinMap_s${correction_sigma}.native.func.gii \
            -var Mask ${mni_native}/${subject}.${hemi}.roi.native.shape.gii
        rm ${mni_native}/${subject}.${hemi}.*_s${correction_sigma}.native.func.gii
    fi

    if [ ! -f ${mni_native}/${subject}.${hemi}.MyelinMap_BC.native.func.gii ]; then
        wb_command -metric-math \
            "(Individual - Bias) * Mask" \
            ${mni_native}/${subject}.${hemi}.MyelinMap_BC.native.func.gii \
            -var Individual ${mni_native}/${subject}.${hemi}.MyelinMap.native.func.gii \
            -var Bias ${mni_native}/${subject}.${hemi}.BiasField.native.func.gii \
            -var Mask ${mni_native}/${subject}.${hemi}.roi.native.shape.gii
    fi

    if [ ! -f ${mni_native}/${subject}.${hemi}.SmoothedMyelinMap_BC.native.func.gii ]; then
        wb_command -metric-math \
            "(Individual - Bias) * Mask" \
            ${mni_native}/${subject}.${hemi}.SmoothedMyelinMap_BC.native.func.gii \
            -var Individual ${mni_native}/${subject}.${hemi}.SmoothedMyelinMap.native.func.gii \
            -var Bias ${mni_native}/${subject}.${hemi}.BiasField.native.func.gii \
            -var Mask ${mni_native}/${subject}.${hemi}.roi.native.shape.gii
    fi

    for string in MyelinMap@func SmoothedMyelinMap@func MyelinMap_BC@func SmoothedMyelinMap_BC@func corrThickness@shape ; do
        map=$(echo ${string} | cut -d "@" -f 1)
        ext=$(echo ${string} | cut -d "@" -f 2)

        wb_command -set-map-name \
            ${mni_native}/${subject}.${hemi}.${map}.native.${ext}.gii 1 ${subject}_${hemi}_${map}

        wb_command -metric-palette \
            ${mni_native}/${subject}.${hemi}.${map}.native.${ext}.gii \
            MODE_AUTO_SCALE_PERCENTAGE \
            -pos-percent 4 96 \
            -interpolate true \
            -palette-name videen_style \
            -disp-pos true \
            -disp-neg false \
            -disp-zero false

        # resample native data to high res mesh, and mask
        if [ ! -f ${mni}/${subject}.${hemi}.sphere.${high_res_mesh}k_fs_LR.surf.gii ]; then
            wb_command -metric-resample \
                ${mni_native}/${subject}.${hemi}.${map}.native.${ext}.gii \
                ${sphere} \
                ${mni}/${subject}.${hemi}.sphere.${high_res_mesh}k_fs_LR.surf.gii \
                ADAP_BARY_AREA \
                ${mni}/${subject}.${hemi}.${map}.${high_res_mesh}k_fs_LR.${ext}.gii \
                -area-surfs ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
                ${mni}/${subject}.${hemi}.midthickness.${high_res_mesh}k_fs_LR.surf.gii \
                -current-roi ${mni_native}/${subject}.${hemi}.roi.native.shape.gii

            wb_command -metric-mask \
                ${mni}/${subject}.${hemi}.${map}.${high_res_mesh}k_fs_LR.${ext}.gii \
                ${mni}/${subject}.${hemi}.atlasroi.${high_res_mesh}k_fs_LR.shape.gii \
                ${mni}/${subject}.${hemi}.${map}.${high_res_mesh}k_fs_LR.${ext}.gii
        fi

        # resample native data to low res mesh, and mask
        if [ ! -f ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.sphere.${low_res_mesh}k_fs_LR.surf.gii ]; then
            wb_command -metric-resample \
                ${mni_native}/${subject}.${hemi}.${map}.native.${ext}.gii \
                ${sphere} \
                ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.sphere.${low_res_mesh}k_fs_LR.surf.gii \
                ADAP_BARY_AREA \
                ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}.${low_res_mesh}k_fs_LR.${ext}.gii \
                -area-surfs ${t1_native}/${subject}.${hemi}.midthickness.native.surf.gii \
                ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.midthickness.${low_res_mesh}k_fs_LR.surf.gii \
                -current-roi ${mni_native}/${subject}.${hemi}.roi.native.shape.gii

            wb_command -metric-mask \
                ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}.${low_res_mesh}k_fs_LR.${ext}.gii \
                ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.atlasroi.${low_res_mesh}k_fs_LR.shape.gii \
                ${mni}/fsaverage_LR${low_res_mesh}k/${subject}.${hemi}.${map}.${low_res_mesh}k_fs_LR.${ext}.gii
        fi
    done
done

# Create CIFTI Files
## TODO !!!
#  "${mni}@${high_res_mesh}k_fs_LR@atlasroi"   removed from outer loop -- these files are absent?
#  "${mni}/fsaverage_LR${low_res_mesh}k@${low_res_mesh}k_fs_LR@atlasroi"
## TODO !!!
for string in "${mni_native}@native@roi" ; do
    folder=$(echo ${string} | cut -d "@" -f 1)
    mesh=$(echo ${string} | cut -d "@" -f 2)
    roi=$(echo ${string} | cut -d "@" -f 3)

    for stringII in MyelinMap@func SmoothedMyelinMap@func MyelinMap_BC@func SmoothedMyelinMap_BC@func corrThickness@shape ; do
        map=$(echo ${stringII} | cut -d "@" -f 1)
        ext=$(echo ${stringII} | cut -d "@" -f 2)

        if [ ! -f ${folder}/${subject}.${map}.${mesh}.dscalar.nii ]; then
            wb_command -cifti-create-dense-scalar \
                ${folder}/${subject}.${map}.${mesh}.dscalar.nii \
                -left-metric ${folder}/${subject}.L.${map}.${mesh}.${ext}.gii \
                -roi-left ${folder}/${subject}.L.${roi}.${mesh}.shape.gii \
                -right-metric ${folder}/${subject}.R.${map}.${mesh}.${ext}.gii \
                -roi-right ${folder}/${subject}.R.${roi}.${mesh}.shape.gii

            wb_command -set-map-names \
                ${folder}/${subject}.${map}.${mesh}.dscalar.nii -map 1 "${subject}_${map}"

            wb_command -cifti-palette \
                ${folder}/${subject}.${map}.${mesh}.dscalar.nii \
                MODE_AUTO_SCALE_PERCENTAGE \
                ${folder}/${subject}.${map}.${mesh}.dscalar.nii \
                -pos-percent 4 96 \
                -interpolate true \
                -palette-name videen_style \
                -disp-pos true \
                -disp-neg false \
                -disp-zero false
        fi
    done
done

## TODO !!!
#  "${mni}@${mni}@${high_res_mesh}k_fs_LR"   removed from outer loop -- these files are absent?
#  "${t1}/fsaverage_LR${low_res_mesh}k@${mni}/fsaverage_LR${low_res_mesh}k@${low_res_mesh}k_fs_LR"
#  "${mni}/fsaverage_LR${low_res_mesh}k@${mni}/fsaverage_LR${low_res_mesh}k@${low_res_mesh}k_fs_LR"
## TODO !!!
#Add CIFTI Maps to Spec Files
for string in "${t1_native}@${mni_native}@native" "${mni_native}@${mni_native}@native" ; do
    folder=$(echo ${string} | cut -d "@" -f 1)
    folderII=$(echo ${string} | cut -d "@" -f 2)
    mesh=$(echo ${string} | cut -d "@" -f 3)
    for stringII in MyelinMap_BC@dscalar SmoothedMyelinMap_BC@dscalar corrThickness@dscalar ; do
        map=$(echo ${stringII} | cut -d "@" -f 1)
        ext=$(echo ${stringII} | cut -d "@" -f 2)
        wb_command -add-to-spec-file \
            "${folder}/${subject}.${mesh}.wb.spec" \
            INVALID \
            "${folderII}/${subject}.${map}.${mesh}.${ext}.nii"
    done
done





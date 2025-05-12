#!/usr/bin/env bash
set -euo pipefail

# Turn on nullglob so file matching patterns expands to “nothing” when no file matches
shopt -s nullglob

# Open a Docker container and make a home directory
contId=$(docker run -d gcr.io/ris-registry-shared/fsl6 tail -f /dev/null)
docker exec "$contId" bash -c "mkdir /home/fsl/"
printf "Successfully created Docker container.\nContainer ID: ${contId:0:3}\n"

# Copy the acqParams to the Docker conatiner
docker cp ./fmAcqParams "$contId":/home/fsl/ 2> /dev/null

# Loop through all the input SubjectIds
scriptPath=$(pwd)
for sId in "$@"; do

    # Copy all the files to the container
    cd "../../Data/$sId/Fieldmap"
    for file in *.nii; do
        docker cp "$file" "$contId":/home/fsl/ 2> /dev/null
    done

    # FslMerge
    ## Build an *array* whose elements are the matching file names.
    AP_matches=( *FmAP* )
    PA_matches=( *FmPA* )
    ## Pre-pend “/home/fsl/” to each element
    AP_paths=( "${AP_matches[@]/#/\/home\/fsl\/}" )
    PA_paths=( "${PA_matches[@]/#/\/home\/fsl\/}" )
    ## Collapse the array into one space-separated string.
    AP_string="${AP_paths[*]}"
    PA_string="${PA_paths[*]}"
    APPA_string="${AP_string} ${PA_string}"
    ## Run the command
    docker exec "$contId" bash -c "fslmerge -t /home/fsl/_${sId}_FmAPPA ${APPA_string}"

    # Create the off-resonance image
    topupCmd=$(cat <<EOF
topup \
--imain=/home/fsl/_${sId}_FmAPPA.nii.gz \
--datain=/home/fsl/fmAcqParams \
--config=b02b0.cnf \
--out=/home/fsl/_${sId}_FmTopupCoefs \
--fout=/home/fsl/_${sId}_FmORI \
--iout=/home/fsl/_${sId}_SampleUnwarp \
--logout=/home/fsl/_${sId}_FmTopupLog
EOF
)
    docker exec "$contId" bash -c "$topupCmd"

    # Copy the files back
    docker cp "$contId":"/home/fsl/_${sId}_FmORI.nii.gz" ./ 2> /dev/null
    docker cp "$contId":"/home/fsl/_${sId}_FmTopupCoefs_fieldcoef.nii.gz" ./ 2> /dev/null
    docker cp "$contId":"/home/fsl/_${sId}_FmTopupCoefs_movpar.txt" ./ 2> /dev/null
    docker cp "$contId":"/home/fsl/_${sId}_FmTopupLog" ./ 2> /dev/null
    docker cp "$contId":"/home/fsl/_${sId}_SampleUnwarp.nii.gz" ./ 2> /dev/null

    # Rename files
    mv "./_${sId}_FmTopupCoefs_fieldcoef.nii.gz" "_${sId}_FmTopupCoefs-Splines.nii.gz"
    mv "./_${sId}_FmTopupCoefs_movpar.txt" "_${sId}_FmTopupCoefs-Movpar.txt"

    # Gunzip the resulting images and delete the compressed versions
    gunzip *.gz

    # Delete all files in the container (excluding acqParams)
    docker exec "$contId" bash -c "rm /home/fsl/_${sId}_*"

    # CD back to the OG location
    cd "$scriptPath"
done

# Close and remove the docker container
docker stop "$contId"
docker rm "$contId"
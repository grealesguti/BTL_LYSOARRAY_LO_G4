#!/usr/bin/env bash
#BASEDIR=/storage/af/user/greales/simG4/BTL_LYSOARRAY_LO_G4
BASEDIR=/home/greales/Documents/Coding/BTL_LYSOARRAY_LO_G4

# Define the number of divisions for szloc and len
szloc_divisions=10
len_divisions=10

# Calculate the step size for szloc and len
szloc_step=$(echo "scale=2; 0.9 / $szloc_divisions" | bc)
len_step=$(echo "scale=2; (28.5 - 6) / $len_divisions" | bc)

for ((szloc_num=1; szloc_num<=$szloc_divisions; szloc_num++))
do
    szloc=$(echo "scale=2; 0.1 + ($szloc_num - 1) * $szloc_step" | bc)
    
    for ((len_num=1; len_num<=$len_divisions; len_num++))
    do
        len=$(echo "scale=2; 6 + ($len_num - 1) * $len_step" | bc)
        name="TileSweep_${szloc}_len_${len}"

        $BASEDIR/JobFiles/JobActionTile.sh -z $szloc -l $len -a $name

        # Uncomment the following lines if you want to submit condor jobs
        # condor_submit $BASEDIR/SubFiles/SubmissionFileGC13mmResinMuon.sub
        # condor_submit $BASEDIR/SubFiles/SubmissionFileGC3Muon.sub
        # condor_submit $BASEDIR/SubFiles/SubmissionFileGC1FLResin.sub
        # condor_submit $BASEDIR/SubFiles/SubmissionFileGC13mmResin.sub
        # condor_submit $BASEDIR/SubFiles/SubmissionFileGC3.sub
    done
done


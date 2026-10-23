#!/bin/bash

export GraphPep_root=/home/g/GraphPep_v1.3
export PATH=$GraphPep_root/bin:$PATH
export esm2_path=/home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt
export GraphPep_env=graphpep_env

source /home/g/miniconda3/etc/profile.d/conda.sh
conda activate $GraphPep_env

# -----------------------------------------------------------------------------------------------------------------------------
# command args
if [ $# -lt 2 ];then
    echo ""
    echo "USAGE: `basename $0` protein.pdb decoys.pdb [options]"
    echo ""
    echo "Descriptions:"
    echo "    protein.pdb   : input, the file of protein strutcture(pdb file)"
    echo "    decoys.pdb    : input, the file of peptide decoys structure(pdb file)"
    echo "    -out          : the output filename of predicted scores, default --> score.csv"
    echo "    -each         : flag to output each residue contact confidence score, default --> False"
    echo "    -outeach      : the output filename of predicted each confidence score, default --> score_each.csv"
    echo "    -PB           : perform plausibility checks for models, and output valid number in output file."
    echo ""
        exit 1
fi


pro=$1
pep=$2
fout=score.csv
outeach=False
fouteach=score_each.csv
PB=False
ckpt=$GraphPep_root/ckpt/model_param.ckpt

cdir=$(pwd)
tmpdir=$(mktemp -d)
cd $tmpdir

id=$$

ln -s $cdir/$1 pro${id}.pdb
ln -s $cdir/$2 pep${id}.pdb

while [ $# -gt 2 ];do
    case $3 in

    -out)
        shift
        fout=$3;;
    -ckpt)
        shift
        ckpt=$3;;
    -each)
        outeach=True;;
    -outeach)
        shift
        fouteach=$3;;
    -PB)
        PB=True;;
    *)
        echo " ERROR: wrong command argument \"$3\" !!"
        echo " Type \"$0\" for help !!"
        exit 2;;
    esac
    shift
done

start_time=$(date +%s)

id=$$
$GraphPep_root/bin/renumm.awk pro${id}.pdb | $GraphPep_root/bin/cleanpdb.awk > pro_renum${id}.pdb
$GraphPep_root/bin/renumm.awk pep${id}.pdb | $GraphPep_root/bin/cleanpdb.awk > pep_renum${id}.pdb
pro_ori=pro${id}.pdb
pep_ori=pep${id}.pdb
pro=pro_renum${id}.pdb
pep=pep_renum${id}.pdb

if [[ $outeach == True && $PB == False ]];then
    python $GraphPep_root/bin/predict.py --data_path . --ckpt_path $ckpt --profile $pro --decoysfile $pep --testout_file $fout --outputres --testresout_file tmp_each.csv --esm2_path $esm2_path
    # change residue number to the origin residue number and chain id
    python /mnt/e/Users/graphpep与pepmlm/output/bh3_bclw_six_complex_pilot_20261020/graphpep_patches/resmap_compat.py $pro_ori $pro $pep_ori $pep tmp_each.csv $fouteach

    mv -f $fout $cdir 2>/dev/null
    mv -f $fouteach $cdir 2>/dev/null

elif [[ $outeach == False && $PB == False ]];then
    python $GraphPep_root/bin/predict.py --data_path . --ckpt_path $ckpt --profile $pro --decoysfile $pep --testout_file $fout --esm2_path $esm2_path
    mv -f $fout $cdir 2>/dev/null

elif [[ $outeach == False && $PB == True ]];then
    python $GraphPep_root/bin/predict.py --data_path . --ckpt_path $ckpt --profile $pro --decoysfile $pep --testout_file score${id}.csv --esm2_path $esm2_path
    nmodel=`grep MODEL $pep | wc -l`
    if [ $nmodel == 0 ];then
#      babel $pep pep_renum${id}.sdf &>/dev/null
      python $GraphPep_root/bin/pdb2sdf.py $pep pep_renum${id}.sdf &>/dev/null
      bust pep_renum${id}.sdf -p $pro > pb-valid${id}.out
    else
      > pb-valid${id}.out
      for((i=1;i<=$nmodel;i++))
      do
        getmodel.sh $pep $i
#        python $GraphPep_root/bin/getmodel.py $pep $i
#        babel pep_renum${id}.$i.pdb pep_renum${id}.$i.sdf &>/dev/null
        python $GraphPep_root/bin/pdb2sdf.py pep_renum${id}.$i.pdb pep_renum${id}.$i.sdf &>/dev/null
        bust pep_renum${id}.$i.sdf -p $pro >> pb-valid${id}.out
      done
    fi
    awk 'BEGIN{print "PB-valid"}{print $(NF-2),$(NF-1),$NF}' pb-valid${id}.out | paste -d"," score${id}.csv - > $fout
    mv -f $fout $cdir 2>/dev/null
else 
    python $GraphPep_root/bin/predict.py --data_path . --ckpt_path $ckpt --profile $pro --decoysfile $pep --testout_file score${id}.csv --outputres --testresout_file tmp_each.csv --esm2_path $esm2_path
    # change residue number to the origin residue number and chain id
    python /mnt/e/Users/graphpep与pepmlm/output/bh3_bclw_six_complex_pilot_20261020/graphpep_patches/resmap_compat.py $pro_ori $pro $pep_ori $pep tmp_each.csv $fouteach
    nmodel=`grep MODEL $pep | wc -l`
    if [ $nmodel == 0 ];then
#      babel $pep pep_renum${id}.sdf &>/dev/null
      python $GraphPep_root/bin/pdb2sdf.py $pep pep_renum${id}.sdf &>/dev/null
      bust pep_renum${id}.sdf -p $pro > pb-valid${id}.out
    else
      > pb-valid${id}.out
      for((i=1;i<=$nmodel;i++))
      do
        getmodel.sh $pep $i
#        python $GraphPep_root/bin/getmodel.py $pep $i
#        babel pep_renum${id}.$i.pdb pep_renum${id}.$i.sdf &>/dev/null
        python $GraphPep_root/bin/pdb2sdf.py pep_renum${id}.$i.pdb pep_renum${id}.$i.sdf &>/dev/null
        bust pep_renum${id}.$i.sdf -p $pro >> pb-valid${id}.out
      done
    fi
    awk 'BEGIN{print "PB-valid"}{print $(NF-2),$(NF-1),$NF}' pb-valid${id}.out | paste -d"," score${id}.csv - > $fout

    mv -f $fout $cdir 2>/dev/null
    mv -f $fouteach $cdir 2>/dev/null

fi

cd - >/dev/null
rm -rf $tmpdir

end_time=$(date +%s)
echo predicting time $(($end_time - $start_time))

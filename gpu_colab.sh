#!/bin/bash
. /edisk3/af2/localcolabfold/conda/etc/profile.d/conda.sh
conda activate /edisk3/af2/localcolabfold/colabfold-conda
CURRENTPATH=`pwd`
colabfold_batch "${CURRENTPATH}/msas" "${CURRENTPATH}/predictions"
echo "完成"

###超算中心####
 ###  #!/bin/bash
 ###  #SBATCH --nodes=1
 ###  #SBATCH --job-name="af2"
 ###  #SBATCH --ntasks-per-node=10
 ###  #SBATCH --partition=gpu1
 ###  #SBATCH --gres=gpu:1
 ###  . ~/program/af2/localcolabfold/conda/etc/profile.d/conda.sh
 ###  conda activate  ~/program/af2/localcolabfold/colabfold-conda
 ###  CURRENTPATH=`pwd`
 ###  colabfold_batch "${CURRENTPATH}/msas" "${CURRENTPATH}/predictions"
 ###  echo "完成"

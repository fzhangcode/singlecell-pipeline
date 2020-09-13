#!/bin/bash

function do_mkfastq()
{
    local library=$1
    
cat << EOF | bsub 

#!/bin/bash
#BSUB -J ${library}
#BSUB -o /data/srlab/bwh10x/10X-Core-Pipeline${library}/log/%J.out
#BSUB -e /data/srlab/bwh10x/10X-Core-Pipeline${library}/log/%J.err
#BSUB -q big-multi
#BSUB -M 64000
#BSUB -n 16
#BSUB -N

module load bcl2fastq2/2.19.1
module load casava/1.8.3
cd /data/srlab/bwh10x/10X-Core-Pipeline/$library

/data/srlab/cellranger/cellranger-3.0.2/cellranger mkfastq --id=FASTQS --run=/data/srlab/bwh10x/10X-Core-Pipeline/$library --csv=/data/srlab/bwh10x/10X-Core-Pipeline/$library/sample_sheet.csv --jobmode=local --localcores=16 --localmem=64

EOF
}

RAWBCL=$1

cat ${RAWBCL}/lsf_params_mkfastq | while read library
do
#echo $library
do_mkfastq $library
done


# Configuration file
import os
if len(config) == 0:
  if os.path.isfile("./config.yaml"):
    configfile: "./config.yaml"
  else:
    sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Make sure that all expected variables from the config file are in the config dictionary
configvars = ['rawbcl']
for k in configvars:
        if k not in config:
                config[k] = None

## If any of the file paths is missing, replace it with ""
def sanitizefile(str):
        if str is None:
                str = ''
        return str

config['rawbcl'] = sanitizefile(config['rawbcl'])

import pandas as pd
data = '/data/srlab/bwh10x/rawbcltarlist.txt'
sample_sheet = '/data/srlab/bwh10x/Master-10Xsamplelog-RNA-SeqSamplelog.csv'


rule all: 
  input: 
    expand('{rawbcl}/sample_sheet.csv', rawbcl = config['rawbcl']),
    expand('{rawbcl}/lsf_params_mkfastq', rawbcl = config['rawbcl']),
    expand('{rawbcl}/made_fastqs.txt', rawbcl = config['rawbcl']),
#    expand('{rawbcl}/log/', rawbcl = config['rawbcl']),
#    expand('{rawbcl}/cellranger-3.0.2/', rawbcl = config['rawbcl']),
    expand('{rawbcl}/lsf_params_count', rawbcl = config["rawbcl"])


#rule make_mapping:
#  output:
#    'mapping.txt'
#  run:
#     with open(output[0], 'w') as out:
#       for i in shell("cat {data}", iterable = True):
#        sample = i.split('.')[0]
#        workorder = sample.split('_')[1]
#        if workorder[0:2] == "KW":
#          out.write(sample + ' ' + workorder[0:6] + '\n')
#        else:
#          next


#rule already_processed:
#  input:
#    'mapping.txt'
#  output:
#    'clean_mapping.txt'
#  run:
#    with open(output[0], 'w') as out:
#      for i in shell("cat {input[0]}", iterable = True):
#        rawbcl = i.split()[0]
#        rawbcl_path = '/data/srlab/bwh10x/{}/RTARead1Complete.txt'.format(rawbcl)
#        rawbcl_path_exists = os.path.isfile(rawbcl_path)
#        path = '/data/srlab/bwh10x/{}/sample_sheet.csv'.format(rawbcl)
#        sample_sheet_exists = os.path.isfile(path)
#        if sample_sheet_exists:
#          next
#        elif not rawbcl_path_exists:
#          next
#        else:
#          out.write(i + '\n')


rule mk_sample_sheet:
  input:
    rawbcl = expand('{rawbcl}', rawbcl = config["rawbcl"])
  output:
    '{rawbcl}/sample_sheet.csv'
  run:
    for i in input:
      open('{}/sample_sheet.csv'.format(i), 'w').write('Lane,Sample,Index\n')
      work_order = i.split('_')[1]
      sample_log = pd.read_csv(sample_sheet, sep = ",", index_col = 'work order')
      sample_log_small = sample_log.loc[[work_order]]
      samples = sample_log_small['Library ID'].tolist()
      i7 = sample_log_small['i7 Well'].tolist()
      for j in range(len(samples)):
        open('{}/sample_sheet.csv'.format(i), 'a+').write('1-4,' + samples[j] + ',SI-GA-' + i7[j] + '\n')


rule mkfastq_params:
  input:
    rawbcl = '{rawbcl}'
  output:
    '{rawbcl}/lsf_params_mkfastq'
  shell:
    'echo {input.rawbcl} > {output}'  
      
#rule mkfastq:
#  input:
#    rawbcl_params = '{rawbcl}/lsf_params_mkfastq'
#  output:
#    '{rawbcl}/made_fastqs.txt'
#  shell:
#    'cat {input.rawbcl_params} | bash run_dummy.sh; '
#    'echo {input.rawbcl_params} > {output}'
#    'cat input.rawbcl_params} | bash run_mkfastq.sh'

rule mkfastq:
  input:
    sample_sheet = '{rawbcl}/sample_sheet.csv',
    rawbcl = '{rawbcl}',
    rawbcl_params = '{rawbcl}/lsf_params_mkfastq'
  output:
    '{rawbcl}/made_fastqs.txt'
  shell:  
    'module load bcl2fastq2/2.19.1; '
    'module load casava/1.8.3; '
    'cd /data/srlab/bwh10x/10X-Core-Pipeline/{input.rawbcl}; '
    '/data/srlab/cellranger/cellranger-3.0.2/cellranger mkfastq --id=FASTQS --run=/data/srlab/bwh10x/10X-Core-Pipeline/{input.rawbcl} --csv=/data/srlab/bwh10x/10X-Core-Pipeline/{input.rawbcl}/sample_sheet.csv --jobmode=local --localcores=16 --localmem=64; '
    'cat control_done.txt > {input.rawbcl}/made_fastqs.txt; '

rule control_count:
  input:
    rawbcl = expand('{rawbcl}', rawbcl = config["rawbcl"])
  output:
    expand('{rawbcl}/lsf_params_count', rawbcl = config["rawbcl"])
  run:
# --wait-for-files
    for i in input.rawbcl:
      work_order = i.split('_')[1]
      sample_log = pd.read_csv(sample_sheet, sep = ",", index_col = 'work order')
      sample_log_small = sample_log.loc[[work_order]]
      samples = sample_log_small['Library ID'].tolist()
      species = sample_log_small['Species'].tolist()
      ADT = sample_log_small['ADT (cite-seq) library'].tolist()
      HTO = sample_log_small['HTO (cell-hashing) library'].tolist()
      TCR = sample_log_small['TCR Library'].tolist()
      for j in range(len(samples)):
        if species[j] == "Human":
          genome = "GRCh38"
        elif species[j] == "Mouse":
          genome = "mm10"  
        if ADT[j] != 'No':
          print('ADT is {}'.format(ADT[j]))
          open('{}/lsf_params_count'.format(i), 'w').write(i + '\t' + samples[j] + '\t' + '/data/srlab/external-data/10xgenomics/refdata-cellranger-' + genome + '-3.0.0' + '\t' + 'cellranger-3.0.2' + '\t' + genome + '\t' + 'ADT')
        if HTO[j] != 'No':
          print('HTO is {}'.format(HTO[j]))
          open('{}/lsf_params_count'.format(i), 'w').write(i + '\t' + samples[j] + '\t' + '/data/srlab/external-data/10xgenomics/refdata-cellranger-' + genome + '-3.0.0' + '\t' + 'cellranger-3.0.2' + '\t' + genome + '\t' + 'HTO')
        if TCR[j] != 'No':
          print('TCR is {}'.format(TCR[j]))
          open('{}/lsf_params_count'.format(i), 'w').write(i + '\t' + samples[j] + '\t' + '/data/srlab/external-data/10xgenomics/refdata-cellranger-' + genome + '-3.0.0' + '\t' + 'cellranger-3.0.2' + '\t' + genome + '\t' + 'TCR')
        else:
          print('mRNA')
          open('{}/lsf_params_count'.format(i), 'w').write(i + '\t' + samples[j] + '\t' + '/data/srlab/external-data/10xgenomics/refdata-cellranger-' + genome + '-3.0.0' + '\t' + 'cellranger-3.0.2' + '\t' + genome + '\t' + 'mRNA')



   

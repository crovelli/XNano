import CRABClient
from CRABClient.UserUtilities import config, ClientException#, getUsernameFromCRIC
import yaml
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser

production_tag = datetime.date.today().strftime('%Y%b%d')

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'Xnani_%s' % production_tag

config.section_('Data')
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_bphys/crovelli/nanoaod_X/%s' % (config.General.workArea)
config.Data.inputDBS = 'global'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/run_nano_2022_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

if __name__ == '__main__':

  from CRABAPI.RawCommand import crabCommand
  from CRABClient.ClientExceptions import ClientException
  from http.client import HTTPException
  from multiprocessing import Process

  def submit(config):
      try:
          crabCommand('submit', config = config)
      except HTTPException as hte:
          print("Failed submitting task:",hte.headers)
      except ClientException as cle:
          print("Failed submitting task:",cle)

  parser = ArgumentParser()
  parser.add_argument('-y', '--yaml', default = 'MCsamples_2022.yml', help = 'File with dataset descriptions')
  parser.add_argument('-f', '--filter', default='*', help = 'filter samples, POSIX regular expressions allowed')
  args = parser.parse_args()

  with open(args.yaml) as f:
    doc = yaml.load(f,Loader=yaml.FullLoader) # Parse YAML file
    common = doc['common'] if 'common' in doc else {'data' : {}, 'mc' : {}}
    
    # loop over samples
    for sample, info in doc['samples'].items():
      # Input DBS
      input_dbs = info['dbs'] if 'dbs' in info else None
      isLocal = info['isLocal'] if 'isLocal' in info else False
      # Given we have repeated datasets check for different parts
      parts = info['parts'] if 'parts' in info else [None]
      for part in parts:
        name = sample % part if part is not None else sample
        
        # filter names according to what we need
        if not fnmatch(name, args.filter): continue
        print('submitting', name)

        isMC = info['isMC']

        if not isLocal :
            config.Data.inputDataset = info['dataset'] % part \
                                    if part is not None else \
                                        info['dataset']
            if input_dbs : config.Data.inputDBS = input_dbs 
        else:
            config.Data.userInputFiles = open(info['dataset']).readlines()
            config.Data.outputPrimaryDataset = sample + '_MGv5NLO_pythia8_NanoAOD'

        config.General.requestName = name
        common_branch = 'mc' if isMC else 'data'
        config.Data.splitting = 'FileBased' if isMC else 'LumiBased'
        if not isMC:
            config.Data.lumiMask = info.get(
                'lumimask', 
                common[common_branch].get('lumimask', None)
            )
        else:
            config.Data.lumiMask = ''

        config.Data.unitsPerJob = info.get(
            'splitting',
            common[common_branch].get('splitting', None)
        )
        globaltag = info.get(
            'globaltag',
            common[common_branch].get('globaltag', None)
        )
        
        config.JobType.pyCfgParams = [
            'isMC=%s' % isMC, 
            'reportEvery=1000',
            'tag=%s' % production_tag,
            'globalTag=%s' % globaltag,
        ]
        
        config.JobType.outputFiles = ['_'.join(['xNANO', 'mc' if isMC else 'data', production_tag])+'.root']
        
        print()
        print(config)
        submit(config)


import os
import fnmatch
import subprocess
import datetime

def merge_ROOT_files():
    output_file = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H%M%SZ_merged.root")

    root_files = [
        os.path.join(root, filename)
        for root, dirs, files in os.walk('/eos/cms/store/group/phys_bphys/crovelli/nanoaod_X/Xdata2017_2022Jul08/Charmonium/crab_data_Run2017B/220708_140547/0001/')
        for filename in files
        if fnmatch.fnmatch(filename, '*.root')
    ]

    # Remove 'echo' when you want to go live.
    subprocess.check_call(['echo', 'python haddnano.py', output_file]+root_files)

if __name__ == "__main__":
    merge_ROOT_files()

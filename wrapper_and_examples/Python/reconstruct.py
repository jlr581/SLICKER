import sys
import os
import numpy as np
import pandas as pd

from slicker import SLICKER
welcome_msg="Welcome to SLICKER"

def run_test(testinputpath,testid):
    test_files_path=testinputpath
    prefix=testid
    reco_outpath=test_files_path+"_output"
    print(reco_outpath)
    print(prefix)
    timebase=np.genfromtxt(os.path.join(test_files_path,prefix+'_input_timebase.txt'))
    proxy1=np.genfromtxt(os.path.join(test_files_path,prefix+'_input_proxy1.txt'))
    print(proxy1.shape)
    proxy2=np.genfromtxt(os.path.join(test_files_path,prefix+'_input_proxy2.txt'))
    print(proxy2.shape)
    target=np.genfromtxt(os.path.join(test_files_path,prefix+'_input_target.txt'))
    print("target....")
    print(target)
    reco_engine=SLICKER.Slicker(prefix=testid,output_path=reco_outpath,timebase=timebase,proxies=[proxy1,proxy2],target=target,nmembers=2048)
    reco_engine.prepare_inputs()
    print(reco_engine.filestring)
    reco_engine.run(testid)
    output=reco_engine.parse_output(fmt='pd')
    print(output)


def run_test_read(testpath,testid):
    prefix=testid
    reco_engine=SLICKER.Slicker(output_path=testpath,prefix=prefix,existing=True)
    reco_engine.read_existing(testpath,prefix)
    results=reco_engine.parse_output()




if __name__ == '__main__':
    print(welcome_msg)

    run_test('../Python_examples', 'test5')
    #run_test_read('../Python_examples_output','test5')

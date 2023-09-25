import numpy as np
import pandas as pd
import os
import subprocess
from pathlib import Path
import glob

class Slicker():

    def __init__(self, prefix='test11',output_path='./output_reco', linear=True,nmembers=4096,tol=1e-05,max_time=1200,width=[0.4,1.6],ens_sub=0.5,proxy_inv=False,slicker_path='..',**kwargs):
        '''
        output_path - a string containing the path to write the files used by SLICKER to, used in the
        prefix - a prefix used to identify individual
        timebase - a 1xn array of times strictly increasing for reconstruction output
        target - a 2xm array of form time,value of the target calibration data
        proxies - a list of N 2xp arrays of form time,value for use in the reconstruction
        width - A tuple of two width parameters optional
        linear - Boolean to turn non-linear treatment of proxies on or off
        proxy_inv - A switch for inverting the proxies - suggested set to True if linear is False
        ens_sub - Fraction of ensemble to subset based on stationarity (0-1)
        slicker_path = a string or path to set the location of the compiled SLICKER binary.
        '''
        self.output_path=Path(output_path).absolute()
        self.prefix=prefix
        self.proxies=[]
        self.linear=linear
        self.nmembers=nmembers
        self.tol=tol
        self.max_time=max_time
        self.width=width
        self.ens_sub=ens_sub
        self.proxy_inv=proxy_inv
        self.set_SLICKER_bin(slicker_path) #need to think about this to make sure we have the fortran bin in the right place
        self.input_filename = ""
        self.filestring = ""
        if 'existing' in kwargs:
            if kwargs['existing']:
                self.read_existing(output_path,prefix)
        if 'timebase' in kwargs:
            self.set_timebase(kwargs['timebase'])
        if 'proxies' in kwargs:
            self.set_proxies(kwargs['proxies'])
            self.nproxies=len(self.proxies)
        if 'target' in kwargs:
            self.set_target(kwargs['target'])

    def set_SLICKER_bin(self,p):
        if p:
            self.path_to_slicker = Path(p).absolute()
        else:
            self.path_to_slicker = Path('../fortan').absolute()

    def _check_increasing(self,times):
        '''
        wrapper method to check that all times are increasing
        :return: True or False
        '''
        return np.all(times[:-1]<times[1:])

    def set_timebase(self,tb):
        '''
        Sets the values for the proxy  numpy array in the SLICKER class
        :param proxyies: a list of 2-D array containing times and target values for the reconstruction
        :return:
        '''
        # check timebase shape, needs to be 1xntimes
        if tb.ndim != 1:
            print(tb.shape)
            print('timebase has too many dimensions')
            exit('bye')
            #throw an error
        if self._check_increasing(tb):
            self.timebase=tb
        else:
            print('timebase not strictily increasing values, exiting')
            exit()
            #thro som

    def set_proxies(self,proxies):
        '''
        Sets the values for the proxy  numpy array in the SLICKER class
        :param proxyies: a list of 2-D array containing times and target values for the reconstruction
        :return:
        '''
        for i,p in enumerate(proxies):
            if p.ndim != 2:
                exit('Proxy {} has wrong dimensions, shape must be (2,n), exiting'.format(i))
            if self._check_increasing(p[:,0]):
                self.proxies.append(p)
            else:
                exit('Proxy {} times not strictly increasing values, exiting'.format(i))

    def set_target(self,target):
        '''
        Sets the target numpy array in the SLICKER class. Checks that the times are strictly increasing and returns an
        error if
        :param target: a 2-D array containing times and target values for the reconstruction
        :return:
        '''
        if target.ndim != 2:
            exit('Wrong dimensions, shape must be (2,n), exiting')
        if self._check_increasing(target[:,0]):
            self.target=target
        else:
            exit('Times not strictly increasing values, exiting')

    def prepare_inputs(self):
        '''
        This method setups up the input files to be used by the SLICKER binary. All the information of the path to
        write files to and where the output will be written to, as well as the reconstruction ID and the numpy arrays
        containing the proxies, target and timebase are set in the SLICKER class constructor
        :return:
        '''
        print(self.output_path)
        outpath=self.output_path
        print(self.prefix)
        if not outpath.exists():
            os.mkdir(outpath)
        os.chdir(outpath)
        input_filestring=""
        np.savetxt(os.path.join(outpath,'{}_input_target.txt'.format(self.prefix)),self.target,fmt='%1.6f %1.6f')
        input_filestring+=self.prefix+'_input_target.txt\n'
        input_filestring+='{}\n'.format(self.nproxies)
        for i,p in enumerate(self.proxies):
            np.savetxt(os.path.join(outpath,'{}_input_proxy{}.txt'.format(self.prefix,i+1)),p, fmt='%1.6f %1.6f')
            input_filestring+="{}_input_proxy{}.txt\n".format(self.prefix,i+1)
        np.savetxt(os.path.join(outpath,'{}_input_timebase.txt'.format(self.prefix)),self.timebase,fmt='%1.6f')
        input_filestring+=self.prefix+'_input_timebase.txt\n'
        input_filestring+="{}\n{},{}\n{}\n{}\n".format(self.nmembers,self.width[0],self.width[1],self.max_time,self.tol)
        if self.linear:
            self.output_filename=self.prefix+"_output_SLICKER_linear.txt"
            input_filestring+=self.prefix+"_output_SLICKER_linear.txt\nn\n"
            self.input_filename="input_"+self.prefix+"_linear.txt"
        else:
            self.output_filename=self.prefix+"_output_SLICKER_nonlinear.txt"
            input_filestring+=self.prefix+"_output_SLICKER_nonlinear.txt\ny\n"
            self.input_filename="input_"+self.prefix+"_nonlinear.txt"
        input_filestring+="{}\n".format(self.ens_sub)
        if self.proxy_inv:
            input_filestring+="y\n"
        else:
            input_filestring+="n\n"
        with open(self.input_filename,'w') as f:
            f.write(input_filestring)
        self.filestring=input_filestring
        return

    def run(self,prefix):
        '''
        The main method to perform the reconstruction after the setup of all the input files and makes a system call
        to the Fortran SLICKER binary
        :param prefix: ID of the reconstruction to run
        :return: Return code of the external process called
        '''
        print([os.path.join(self.path_to_slicker,"reconstruct")])
        proc=subprocess.run([os.path.join(self.path_to_slicker,"reconstruct")],input=self.filestring,
            text=True)
        output = proc.stdout
        result = proc.returncode
        return result


    def read_existing(self,output_path,prefix,nonlinear=False):
        '''
        Read an input files that have been previously used for a reconstruction. Does no sanity checking of paths or
        inputs.
        :param output_path: path to find the output from a previous reconstruction, where the existing input files are kept
        :param prefix: prefix for the input files, the ID of the reconstruction
        :return:
        '''
        self.prefix=prefix
        self.output_path=output_path
        print(os.getcwd())
        print(output_path)
        os.chdir(output_path)
        if nonlinear:
            inputfile=glob.glob(os.path.join(output_path,'input*_nonlinear.txt'))[0]
        else:
            inputfile=glob.glob(os.path.join(output_path,'input*_linear.txt'))[0]
        print(inputfile)
        with open(inputfile) as f:
            targname=f.readline().rstrip('\n')
            self.set_target(np.genfromtxt(targname))
            self.nproxies=int(f.readline().rstrip('\n'))
            for n in range(self.nproxies):
                proxyname=f.readline().rstrip('\n')
                self.proxies.append(np.genfromtxt(proxyname))
                print(self.proxies)
            timebasename=f.readline().rstrip('\n')
            self.set_timebase(np.genfromtxt(timebasename))
            self.nmembers=int(f.readline().rstrip('\n'))
            params=f.readline().rstrip('\n').split(',')
            self.width=[float(params[0]),float(params[1])]
            self.max_time=int(f.readline().rstrip('\n'))
            self.tol=float(f.readline().rstrip('\n'))
            self.output_filename=f.readline().rstrip('\n')
            lin=f.readline().rstrip('\n')
            if lin=="y":
                self.linear=True
            else:
                self.linear=False
            self.ens_sub=f.readline().rstrip('\n')
            proxy_inv=f.readline().rstrip('\n')
            if proxy_inv=="y":
                self.proxy_inv=True
            else:
                self.proxy_inv=False


    def parse_output(self,fmt='np'):
        '''
        Reads the output file and returns either an array or pandas dataframe
        :param fmt A string to denote format of return, either 'np' for numpy array or 'pd' for pandas dataframe
        :return:
        '''
        if fmt=='pd':
            output_colnames=['time','M-estimator','95_CI','Qn']
            fullname=os.path.join(self.output_path,self.output_filename)
            self.results=pd.read_csv(fullname, sep='\s+', names=output_colnames, comment='!')
        else:
            self.results=np.genfromtxt(os.path.join(self.output_path,self.output_filename),comments='!')
        return self.results



import numpy as np
import pylab as plt
import sys
import scipy.optimize as opt
import argparse as ar
import os 

parser = ar.ArgumentParser(description='Fit of the RM by a gaussian function')
parser.add_argument('-in',dest='input',help='input file')
parser.add_argument('-pdn',dest='pldir',help='Plot dir name')
parser.add_argument('-s',dest='deprm',help='Start value of the RM.')
parser.add_argument('-m',dest='fen',nargs=2,help='Fit window : RM min RM max.')
parser.add_argument('-sinc',dest='sinc',action='store_true',help='For fit with a cardinal sinus rather than gaussian.')
parser.add_argument('-nbins',dest='nbins',nargs=1,default=1000,help='Bins number for fitted curve computing (defualt = 1000).')
parser.add_argument('-sig',dest='depsig',default=0.1,help='Start value of the sigma (default = 0.1).')
parser.add_argument('-amp',dest='depamp',default=100.,help='Start value of the amplitude (default = 100).')
parser.add_argument('-base',dest='depbase',default=10.,help='Start value of the baseline (default = 10).')
parser.add_argument('-G',dest='save',default='False',help='Name of the saved figure.')
parser.add_argument('-L',dest='log',default='False',help='Name of the output log files.')
parser.add_argument('-rev',dest='rev',action='store_true',help='Reverse abscissa for the plot.')
args = parser.parse_args()

file_name = args.input

data = np.loadtxt(file_name)

if (args.fen) :
	data = np.compress(data[:,0]>float(args.fen[0]),data,axis=0)
	data = np.compress(data[:,0]<float(args.fen[1]),data,axis=0)


def gauss (x,x0,sig,A,B) :

	if args.sinc :
	       return abs(A * np.sinc( (x-x0) * sig )) + B
	else :
		return abs(A) / (np.sqrt(2*np.pi * sig**2)) * np.exp(-(x-x0)**2/(2*sig**2)) + B


fit = opt.curve_fit(gauss,data[:,0],data[:,1])#,p0=(float(args.deprm),float(args.depsig),float(args.depamp),float(args.depbase)),sigma=np.std(data[:,1]))

g1 = gauss(data[:,0],*fit[0])
print '\nRM = {:.4f}'.format(fit[0][0])
print 'Sigma = {:.4f}'.format(fit[0][1])
print 'Amplitude = {:.4f}'.format(fit[0][2])
print 'Baseline = {:.4f}'.format(fit[0][3])
print '\nCovariance matrix =\n',fit[1]
print '\nError on RM position = {:.4f}'.format(np.sqrt(fit[1][0,0]))
print 'Error on sigma = {:.4f}'.format(np.sqrt(fit[1][1,1]))
print '\nchi2 = {:.4f}\n'.format(np.sum( (g1 - data[:,1])**2 ) / len(data[:,1]) )

if (args.log != 'False') :
	_file0 = open(str(args.log),'w')
	_file0.write('Scipy curve_fit results :')
	_file0.write('\n\nRM = {:.4f}'.format(fit[0][0]))
        _file0.write('\n\nRM_err = {:.4f}'.format(fit[0][0]).format(np.sqrt(fit[1][0,0])))
	_file0.write('\nSigma = {:.4f}'.format(fit[0][1]))
	_file0.write('\nAmplitude = {:.4f}'.format(fit[0][2]))
	_file0.write('\nBaseline = {:.4f}'.format(fit[0][3]))
	_file0.write('\n\nCovariance matrix =\n'+str(fit[1]))
	_file0.write('\n\nError on the RM = {:.4f}'.format(np.sqrt(fit[1][0,0])))
	_file0.write('\nError on the sigma = {:.4f}'.format(np.sqrt(fit[1][1,1])))
	_file0.write('\n\nchi2 = {:.4f}\n'.format(np.sum( (g1 - data[:,1])**2 ) / len(data[:,1]) ))
	_file0.close()

if (args.fen) :
	absc = np.linspace(float(args.fen[0]),float(args.fen[1]),int(args.nbins))
else :
	absc = np.linspace(data[0,0],data[-1,0],int(args.nbins))


g1 = gauss(absc,*fit[0])

if args.rev :
	absc *= -1
	data[:,0] *= -1


plt.figure('Fit gauss RM',(20,10))
plt.plot(data[:,0],data[:,1],'Black',marker='s',linewidth=2,label='rmfit')
plt.plot(absc,g1,'Red',linewidth=2.5,label='curve_fit')
plt.xlabel(r'$RM (rad.m^{-2})$',fontsize=18.)
plt.ylabel(r'$Amp (A.U.)$',fontsize=18.)
plt.grid()
plt.tick_params(labelsize=18.,width=2,length=2)
plt.legend(loc=0,fontsize=18,labelspacing=1)
print "RM_err",'{:.4f}'.format(fit[0][0]),'{:.4f}'.format(np.sqrt(fit[1][0,0]))

main_dir = '/data/ezahraoui/J1022+1001/RMfit/'
plot_dir = os.path.join(args.pldir,'plots/')

if (args.save != 'False') :
	plt.savefig( main_dir + plot_dir + str(args.save),format='png')


#plt.show()

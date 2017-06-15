from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import sys,warnings,os

warnings.simplefilter('ignore')

def uv_weighted_avg(values):
	"""Return the weighted average and standard deviation.

	values, weights -- Numpy ndarrays with the same shape.

	samplemean : include the correction factor for the "sample"
	weighted mean (as opposed to the mean of the full population): (see
	http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf)

	"""
	samplemean=True
	weights=1/(uvdata['error'][np.where(uvdata['mag']==values[0])[0][0]:np.where(uvdata['mag']==values[-1])[0][0]+1])**2
	N = len(np.nonzero(weights)[0])
	if N == 1:
		return (values[0])
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)
	if samplemean:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection = 1
	return (average)

def getUVErrors(values):

	samplemean=True
	weights=1/(uvdata['error'][np.where(uvdata['mag']==values[0])[0][0]:np.where(uvdata['mag']==values[-1])[0][0]+1])**2
	N = len(np.nonzero(weights)[0])
	if N == 1:
		return (1/np.sqrt(weights[0]))
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)
	if samplemean:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection = 1
	return (np.sqrt(variance)/np.sqrt(samplecorrection))


def ir_weighted_avg(values):
	"""Return the weighted average and standard deviation.

	values, weights -- Numpy ndarrays with the same shape.

	samplemean : include the correction factor for the "sample"
	weighted mean (as opposed to the mean of the full population): (see
	http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf)

	"""
	samplemean=True
	weights=1/(irdata['error'][np.where(irdata['mag']==values[0])[0][0]:np.where(irdata['mag']==values[-1])[0][0]+1])**2
	N = len(np.nonzero(weights)[0])
	if N == 1:
		return (np.array(values[0]))
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)
	if samplemean:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection = 1
	return (average)

def getIRErrors(values):
	samplemean=True	
	weights=1/(irdata['error'][np.where(irdata['mag']==values[0])[0][0]:np.where(irdata['mag']==values[-1])[0][0]+1])**2
	N = len(np.nonzero(weights)[0])
	if N == 1:
		return (np.array(1/np.sqrt(weights[0])))
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)
	if samplemean:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection = 1
	return (np.sqrt(variance)/np.sqrt(samplecorrection))


def uvweightedMean(x):
	top=0
	bot=0
	errors=uvdata['error'][np.where(uvdata['mag']==x[0])[0][0]:np.where(uvdata['mag']==x[-1])[0][0]+1]
	
	for i in range(len(x)):
		top=top+(x[i]/(errors[i])**2)
		bot=bot+(1/(errors[i])**2)
	return(np.array(top/bot))

def irweightedMean(x):
	top=0
	bot=0
	errors=irdata['error'][np.where(irdata['mag']==x[0])[0][0]:np.where(irdata['mag']==x[-1])[0][0]+1]
	for i in range(len(x)):
		top=top+(x[i]/(errors[i])**2)
		bot=bot+(1/(errors[i])**2)
	return(np.array(top/bot))

def getErrors(x):
	err=0
	N = len(np.nonzero(x)[0])
	for i in range(len(x)):
		err=err+(1/(x[i])**2)
	if N>1:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection=1
	return((1/(err)**.5)/np.sqrt(samplecorrection))
		#rr=
	#return(np.mean(x))

t='Ic'
uvcol={
	'Ib':['U-B','U-V','u-V','U-r'],
	'Ic':['U-V']
}
ircol={
	'Ib':['r-','V-'],
	'Ic':['r-','V-']
}
UVcolorOrder=uvcol[t]
IRcolorOrder=ircol[t]
#UVcolorOrder=['U-B','U-V'] #typeIc
#IRcolorOrder=['r-','V-']	#typeIc
#UVcolorOrder=['U-B','U-V','u-V','U-r'] #typeIb
#IRcolorOrder=['r-','V-'] #typeIb
colors=ascii.read(os.path.join('type'+t,'tables','all'+t+'Colors.dat'))

allUVColors=[]
allIRColors={
	'J':[],
	'H':[],
	'K':[]
}
allUVerr=[]
allIRerr={
	'J':[],
	'H':[],
	'K':[]
}
uvcolors=[]
ircolors={
	'J':[],
	'H':[],
	'K':[]
}
uvtime=[]
irtime={
	'J':[],
	'H':[],
	'K':[]
}
for row in colors:
	for col in UVcolorOrder:
		if row[col]:
			allUVColors.append(row[col])
			allUVerr.append(row[col[0]+col[-1]+'_err'])
			uvcolors.append(col)
			uvtime.append(row['time'])
			break
	for ir in ['J','H','K']:
		for col in IRcolorOrder:
			if row[col+ir]:
				allIRColors[ir].append(row[col+ir])
				allIRerr[ir].append(row[col[0]+ir+'_err'])
				ircolors[ir].append(col+ir)
				irtime[ir].append(row['time'])
				break


uvdata=Table([uvtime,uvcolors,allUVColors,allUVerr],names=['time','color','mag','error'],masked=True)
bin=np.array(np.trunc(uvdata['time']/4))
uvgrouped=uvdata.group_by(bin)
uvtime=uvgrouped['time'].groups.aggregate(np.mean)
#uvmag=uvgrouped['mag'].groups.aggregate(uv_weighted_avg)
#uvmagerr=uvgrouped['mag'].groups.aggregate(getUVErrors)
uvmag=uvgrouped['mag'].groups.aggregate(uvweightedMean)
uvmagerr=uvgrouped['error'].groups.aggregate(getErrors)

temp=Table([np.array(uvtime),np.array(uvmag),np.array(uvmagerr)],names=('time','mag','magerr'))
ascii.write(temp,os.path.join('type'+t,'tables','uv.dat'))
fig=plt.figure()
ax=plt.gca()
ax.errorbar(np.array(uvtime),np.array(uvmag),yerr=np.array(uvmagerr),fmt='x')
ax.invert_yaxis()
plt.title('UV Average Color, Type '+t+' (Bin Size=4 Days)',size=15)
fig.text(0.5, 0.02, 'Time (Since Peak)', ha='center',size=20)
fig.text(0.02, 0.5, 'Color Magnitude (Vega)', va='center', rotation='vertical',size=20)
plt.savefig(os.path.join("type"+t,'plots','UVAverageColor.pdf'),format='pdf')
#plt.show()

for ir in ['J','H','K']:
	irdata=Table([irtime[ir],ircolors[ir],allIRColors[ir],allIRerr[ir]],names=['time','color','mag','error'],masked=True)
	bin=np.array(np.trunc(irdata['time']/2))
	irgrouped=irdata.group_by(bin)
	time=irgrouped['time'].groups.aggregate(np.mean)
	#mag=irgrouped['mag'].groups.aggregate(ir_weighted_avg)
	#magerr=irgrouped['mag'].groups.aggregate(getIRErrors)
	mag=irgrouped['mag'].groups.aggregate(irweightedMean)
	magerr=irgrouped['error'].groups.aggregate(getErrors)
	
	temp=Table([np.array(time),np.array(mag),np.array(magerr)],names=('time','mag','magerr'))
	ascii.write(temp,os.path.join('type'+t,'tables',ir+'.dat'))
	fig=plt.figure()
	ax=plt.gca()
	ax.errorbar(np.array(time),np.array(mag),yerr=np.array(magerr),fmt='x')
	ax.invert_yaxis()
	plt.title('IR Average Color, Type '+t+' (Bin Size=2 Days, Band= '+ir+')',size=15)
	fig.text(0.5, 0.02, 'Time (Since Peak)', ha='center',size=18)
	fig.text(0.02, 0.5, 'Color Magnitude (Vega)', va='center', rotation='vertical',size=18)
	plt.savefig(os.path.join("type"+t,'plots','IRAverage'+ir+"Color.pdf"),format='pdf')

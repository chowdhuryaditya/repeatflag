import numpy as np
from __casac__ import *
from taskinit import casalog
def attachms(vis,nomodify,tblock,field):
	msobj=ms.ms()
	msobj.open(vis,nomodify=nomodify)
	msobj.msselect({'field':field})
	msobj.iterinit(interval=tblock,adddefaultsortcolumns=False)
	endflag=msobj.iterorigin()
	return msobj,endflag

def getchanfreq(vis):
	tbobj=table.table()
	tbobj.open(vis+'/SPECTRAL_WINDOW')
	freq=tbobj.getcol('CHAN_FREQ')
	tbobj.close()
	return freq[:,0]
def getchanmap(visfrom,visto):
	nchanfrom=getchanfreq(visfrom).shape[0]
	nchanto=getchanfreq(visto).shape[0]
	ratio=nchanto//nchanfrom
	excess=nchanto-ratio*nchanfrom
	return nchanfrom,ratio,excess
def transfer(visfrom,visto,fieldfrom='',fieldto='',tblock=300):
	casalog.post("Transfering flags from %s to %s"%(visfrom,visto))
	nchanfrom,ratio,excess=getchanmap(visfrom,visto)
	if(nchanfrom==1):
		singlechan=True
	else:
		singlechan=False
	casalog.post("Number of channels in source ms: %d"%nchanfrom)
	casalog.post("Ratio of number of channels: %d"%ratio)
	casalog.post("Excess channels at the edge (to be flagged): %d"%excess)
	msfrom,endflag=attachms(visfrom,nomodify=True,tblock=tblock,field=fieldfrom)
	msto,endflag=attachms(visto,nomodify=False,tblock=tblock,field=fieldto)
	while(endflag):
		flagsfrom=msfrom.getdata(['FLAG','FLAG_ROW'])
		scan=msfrom.getdata(['SCAN_NUMBER'])['scan_number']
		casalog.post("Transfering frags from scan(s): %s"%np.array2string(np.unique(scan)))
		flagsto=msto.getdata(['FLAG','FLAG_ROW']) 
		
		if(not singlechan):
			flagsto['flag_row']=flagsfrom['flag_row']
			flagsto['flag'][:,:ratio*nchanfrom,:]=np.logical_or(flagsto['flag'][:,:ratio*nchanfrom,:],np.repeat(flagsfrom['flag'],repeats=ratio,axis=1))
			flagsto['flag'][:,ratio*nchanfrom:,:]=True
		else:
			flagsto['flag_row']=flagsfrom['flag_row']
			for ipol in range(flagsto['flag'].shape[0]):
				flagsto['flag'][ipol,:,:]=flagsfrom['flag'][ipol,:,:]
		msto.putdata(flagsto)
		endflagto=msto.iternext()
		endflagfrm=msfrom.iternext()
		endflag=np.logical_and(endflagto,endflagfrm)
		if(not endflag and (endflagfrm and endflagto)):
			print("ERROR! Mismatch in size")
 	msto.iterend()
	msfrom.iterend()
	msto.close()
	msfrom.close()


def repeatflag(visfrom='',visto='',fieldfrom='',fieldto=''):
	transfer(visfrom,visto,fieldfrom,fieldto)
	return None
	

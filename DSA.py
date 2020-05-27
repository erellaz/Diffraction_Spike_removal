# Apply a mute in the linear radon domain to clean linear noise on an image
# Visit erellaz.com for details
# 2017-03-05
#______________________________________________________________________________
#Parameters, change only this block:
#input fit file to process
fitin='C:\Users\User\Desktop\M.fit'
#fitin='C:\Users\User\Desktop\RigelProject\WH3.fit'
maskin='C:\Users\User\Desktop\M_Mask5.fit'
fitout='C:\Users\User\Desktop\M_Radon.fit'
#Splitting in sub image
splitx=4
splity=4
#overlap
overx=200
overy=200
#radiance source coordinates in pixels
radiance=78
# Width in degree of the mute zone
arange=7
#______________________________________________________________________________
# Import astro for fit image format handeling, skimage for Radon transform
from astropy.io import fits
import numpy as np
import math,timeit,os.path,sys
from skimage.transform import radon, iradon
import matplotlib.pyplot as plt
#______________________________________________________________________________
#______________________________________________________________________________
#Plot in regular and radon domain
def disptrans(im,radt):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    ax1.set_title("Space domain")
    ax1.imshow(im)
    ax2.set_title("Radon domain\n")
    ax2.set_xlabel("Projection angle (deg)")
    ax2.set_ylabel("Projection position (pixels)")
    ax2.imshow(radt, aspect='auto')
    fig.tight_layout()
    plt.show()

#______________________________________________________________________________
# See wikipedia Gaussian function for this apodization window
# returns a scale factor between on [0,1]
def gauss(angle,mina,maxa):
    width=float(maxa-mina)
    center=float(mina+width/2)
    x=float(angle)
    return (1-np.exp(-((x-center)**2 / ((.2*width)**2))))
#______________________________________________________________________________
#  See wikipedia Hann function for this apodization window
# returns a scale factor between on [0,1]
def hann(angle,mina,maxa):
    width=float(maxa-mina)
    w=width/(3.0*math.pi)
    x=float(angle)
    if (x-mina)<(width/3):     
        return .5*(1+math.cos((x-mina)/w) )   
    elif (x-mina)>(2*width/3):
        return .5*(1-math.cos((x+(2*width/3)-mina)/w) )    
    else:
        return 0.0
#______________________________________________________________________________
#Apply radon filtering on a 2d array given a radiance angle to filter 
def radfil(im,radiance,arange=10,mutemet='hann'):
    (lx,ly)=im.shape
    minang=max(radiance-arange,0)
    maxang=min(radiance+arange,180)    
    print "Muting=",minang,"to",maxang,"\nSize of chunk to process",lx,ly
    
    #Forward Radon Transform
    print "\nStarting Radon transform"
    #The transform is "squared"
    theta = np.linspace(0., 180., max(im.shape), endpoint=False)
    radt=radon(im,theta=theta)
    print "...radon done.\n"
    
    #__________________________________________________________________
    disptrans(im,radt)
    #__________________________________________________________________
    #Mute in the radon domain
    print "\nMuting in the radon domain...\n"
    
    sf=float(radt.shape[1]/180.0) #scale factor for the angles
    mina=int(minang*sf)
    maxa=int(maxang*sf)
    ortho=int(90*sf)
    for angle in range(mina,maxa):
        for ori in range(0,int(radt.shape[0])-1):
            #choose amongst 3 different ways to mute            
            if mutemet=='door':            
                radt[ori,angle]=0
            elif mutemet=='gauss':
                radt[ori,angle]=radt[ori,angle]*gauss(angle,mina,maxa)
                radt[ori,angle+ortho]=radt[ori,angle+ortho]*gauss(angle,mina,maxa)
            else:
                radt[ori,angle]=radt[ori,angle]*hann(angle,mina,maxa)
                radt[ori,angle+ortho]=radt[ori,angle+ortho]*hann(angle,mina,maxa)
           
           
    #__________________________________________________________________
    print "Starting Inverse Radon transform...\n"
    iradt=iradon(radt,theta=theta)
    print "...Inverse radon done.\n"
    #__________________________________________________________________
    disptrans(iradt,radt)
    #__________________________________________________________________
    if(lx<ly):
        pad=int(ly-lx)/2
        return iradt[pad:lx+pad,:]
    else:
        pad=int(lx-ly)/2
        return iradt[:,pad:ly+pad]
#______________________________________________________________________________
#______________________________________________________________________________
# Here we apodise the edges of the 2D windows, so the overlapping edges, 
# including corners can be stacked to unity.
# 4 linear tapers, one for each edge
def apodise(imf,x1,x2,y1,y2):
    print "Apodise dump:",x1,x2,y1,y2,imf.shape[0],imf.shape[1]   
    x1=2*x1    
    x2=2*x2
    y1=2*y1
    y2=2*y2    
    xx1=float(x1)
    xx2=float(x2)
    yy1=float(y1)
    yy2=float(y2)
      
    for i in range(0,x1):
        imf[i,:]=imf[i,:]*float(i)/(xx1-1.0)
    
    maxx=imf.shape[0]
    for i in range(0,x2):
        imf[maxx-i-1,:]=imf[maxx-i-1,:]*(1.0+(float(i)-(xx2-1.0))/(xx2-1.0))
    
    for i in range(0,y1):
        imf[:,i]=imf[:,i]*float(i)/(yy1-1.0)
    
    maxy=imf.shape[1]
    for i in range(0,y2):
        imf[:,maxy-i-1]=imf[:,maxy-i-1]*(1.0+(float(i)-(yy2-1.0))/(yy2-1.0))
    
    return imf
#______________________________________________________________________________
#______________________________________________________________________________
# Read the fit
start = timeit.timeit()
print "Reading Data Fit file:",fitin
hdulist = fits.open(fitin)
print "\nBasic Fit file information:",
hdulist.info()
print "\nData size:"

if (os.path.isfile(fitout) == True):
    sys.exit("Output file already exists! Aborting.")    


#get the array representing the image
scidata = hdulist[0].data
scidataori=np.copy(scidata) #create a true copy
hdulist.close()
print "Size of the image",scidata.shape,"\n"
sx=int(scidata.shape[1]/splitx)
sy=int(scidata.shape[2]/splity)
res=np.zeros(scidata.shape)
#______________________________________________________________________________
# Read the fit
print "Reading Mask Fit file:",maskin
mhdulist = fits.open(maskin)
print "\nBasic Fit file information:",
mhdulist.info()
print "\nData size:"
#get the array representing the image
maskdata = mhdulist[0].data
hdulist.close()
print "Size of the mask",maskdata.shape,"\n"

#______________________________________________________________________________
# Do star removal with the mask
max_mask=float(max(map(max,maskdata)))
print "Maximum value found in the mask array", max_mask
for channel in range(0,scidata.shape[0]):#loop on all the channels 1 (mono) or 3(RGB)   
    print "Masking Channel=",channel
    for ix in range(0,scidata.shape[1]):
        for iy in range(0,scidata.shape[2]):
            scidata[channel,ix,iy]=scidata[channel,ix,iy]*(1.0-float(maskdata[ix,iy])/max_mask)

fig1 = plt.figure(figsize=(30,15))
ax1 = fig1.add_subplot(111)
ax1.set_title("Star mask")
ax1.imshow(maskdata)
plt.show()
#______________________________________________________________________________
#Loop on the images chunks
for channel in range(0,scidata.shape[0]):#loop on all the channels 1 (mono) or 3(RGB)   
    print "processing Channel=",channel
    for ix in range(0,splitx):
        for iy in range(0,splity):
            #assign the overlaps in the windows, managing edges and corners
            x1=0 if ix==0 else overx
            x2=0 if ix==(splitx-1) else overx
            y1=0 if iy==0 else overy
            y2=0 if iy==(splity-1) else overy
            im=scidata[channel,sx*ix-x1:sx*(ix+1)+x2,sy*iy-y1:sy*(iy+1)+y2]
            #Coordinates of chunck center            
            cx=sx*(ix+.5)
            cy=sy*(iy+.5)          
           
            
            #Radon Filter the chunk            
            imf=radfil(im,radiance,arange,'gauss')
            
           
            fig3 = plt.figure(figsize=(30,15))
            ax3 = fig3.add_subplot(111)
            ax3.set_title("Return from filter and cut pad")
            ax3.imshow(imf)
            plt.show()
            
            #apodize overlap zone of windows for seamless stack
            aimf=apodise(imf,x1,x2,y1,y2)
            
            print "return from Apodize with cut"
            fig4 = plt.figure(figsize=(30,15))
            ax4 = fig4.add_subplot(111)
            ax4.set_title('Return from Apodise')
            ax4.imshow(aimf)
            plt.show()
            
            #stack the windows 
            res[channel,sx*ix-x1:sx*(ix+1)+x2,sy*iy-y1:sy*(iy+1)+y2]=res[channel,sx*ix-x1:sx*(ix+1)+x2,sy*iy-y1:sy*(iy+1)+y2]+aimf

#______________________________________________________________________________
#applying the mask
for channel in range(0,scidata.shape[0]):#loop on all the channels 1 (mono) or 3(RGB)   
    print "DeMasking Channel=",channel
    for ix in range(0,scidata.shape[1]):
        for iy in range(0,scidata.shape[2]):
         res[channel,ix,iy]=(scidataori[channel,ix,iy]*float(maskdata[ix,iy])/max_mask) + (res[channel,ix,iy]*(1.0-float(maskdata[ix,iy])/max_mask))
         
#______________________________________________________________________________
#write to a new fit
print "Writing fit out..." 
hdu = fits.PrimaryHDU(res)
hdulist = fits.HDUList([hdu])
hdulist.writeto(fitout)
end = timeit.timeit()
print "Total time elapsed:", end - start
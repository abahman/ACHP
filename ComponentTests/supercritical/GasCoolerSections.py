import pylab, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.style.use('classic')
mpl.style.use('Elsevier.mplstyle')
#mpl.rcParams['mathtext.fontset'] = 'custom'

fig=pylab.figure()
#ax = fig.add_subplot(111)

pylab.plot(np.r_[-0.01,4.01,4.01,-0.01,-0.01],np.r_[0,0,1,1,0],'k')

pylab.gca().text(2,0.5,'Supercritical',ha='center',va='center',color="w")
pylab.gca().text(2,1.1,'Supercritical Section only',ha='center',va='bottom')

pylab.gca().text(2,-0.1,'$w_{supct}$',ha='center',va='center')

#plot color gradient 1st plot
gradient = np.linspace(0, 1, 2)
gradient = np.vstack((gradient, gradient))
img = pylab.imshow(gradient, extent=[0,4.01,-0.005,1], alpha=0.9,aspect=4, cmap=plt.get_cmap('bwr'))
img.set_clim(-1,0)

y0=1.7 #shift of second axis
pylab.plot(np.r_[-0.01,4.01,4.01,-0.01,-0.01],np.r_[-y0,-y0,-y0+1,-y0+1,-y0],'k')
pylab.plot(np.r_[1.4,1.4],np.r_[-y0,-y0+1],'k-.',color="w")
pylab.gca().text(0.65,-y0+0.5,'Supercritical\n liquid',ha='center',va='center',color="w")
pylab.gca().text(2.65,-y0+0.5,'Supercritical',ha='center',va='center',color="w")
pylab.gca().text(2,-y0+1.1,'Supecritical and Supecritical Liquid Sections',ha='center',va='bottom')

pylab.gca().text(2.65,-y0-0.1,'$w_{supct}$',ha='center',va='center')
pylab.gca().text(0.65,-y0-0.1,'$w_{supliq}$',ha='center',va='center')

#plot color gradient 2nd plot
gradient2 = np.linspace(0, 1, 2)
gradient2 = np.vstack((gradient2, gradient2))
imgplot = pylab.imshow(gradient2, extent=[0,4.01,-1.7,-0.71], alpha=0.9,aspect=4, cmap=plt.get_cmap('brg'))
imgplot.set_clim(-0.01,2.0)

#pylab.xlim(-0.1,4.3)
pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.savefig('GasCoolerSections2.pdf')
pylab.show()
pylab.close()
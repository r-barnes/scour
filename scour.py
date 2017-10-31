#!/usr/bin/env python3
from scipy.interpolate import griddata
import numpy             as np
import pandas            as pd
import matplotlib        as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib.path import Path
import matplotlib.patches as patches
import polyedit
import sys
from scipy.spatial import ConvexHull

#Code to read with Numpy, but I suspect you'll want pandas so you can filter by
#rows easier
# treedat   = np.loadtxt('913aadj.txt')
# seddat    = np.loadtxt('913adjfill.txt')
# treeflags = treedat[:,0]
# treedat   = treedat[:,1:]
# sedflags  = seddat[:,0]
# seddat    = seddat[:,1:]

class ScourInator:
  def __init__(self, sedfile, treefile, step_size=0.1):
    seddat          = pd.read_table(sedfile,    names=['flag','x','y','z'])
    seddat['type']  = 'sed'
    treedat         = pd.read_table(treefile, names=['flag','x','y','z'])
    treedat['type'] = 'tree'
    #Combine data into same data frame
    self.dat = pd.concat((treedat,seddat))

    #Get bounds for interpolated data
    gxmin = min(self.dat['x'])
    gymin = min(self.dat['y'])
    gxmax = max(self.dat['x'])
    gymax = max(self.dat['y'])
    gxran = abs(gxmax-gxmin)
    gyran = abs(gymax-gymin)
    gxmin -= 0.05*gxran
    gxmax += 0.05*gxran
    gymin -= 0.05*gyran
    gymax += 0.05*gyran

    #No region of interest, initially
    self.roi        = None
    self.roi_path   = None
    self.use_roi    = True
    self.roi_points = None
    self.step_size  = step_size

    #Separate out the sediment data
    self.sed_data     = self.dat.loc[self.dat['type'] == 'sed'].copy()

    #Mash points into array
    #Take the columns x,y,z and zip them together into a list of (x,y,z) tuples
    #Turn this list of tuples into a numpy array    
    self.sed_data_arr = np.array(list(zip(self.sed_data['x'], self.sed_data['y'], self.sed_data['z'])))

    #Get z-values. Useful for plotting
    self.zmin = min(self.sed_data['z'])
    self.zmax = max(self.sed_data['z'])
    self.zran = self.zmax-self.zmin

    #Get the outline of the tree
    ttemp         = self.dat.loc[self.dat['type'] == 'tree'] #Extract tree data frame data frame
    verts         = list(zip(ttemp['x'],ttemp['y']))         #Build list of (x,y) pairs of the points of the tree
    codes         = [Path.LINETO]*len(verts)                 #Build list of codes equal to length of list of (x,y) pairs
    codes[0]      = Path.MOVETO                              #Set first code to MOVETO
    codes[-1]     = Path.CLOSEPOLY                           #Set last code to CLOSEPOLY
    self.treepath = Path(verts, codes)                       #Create a path

    #Generate plotting grid with 100 steps along both axes
    self.grid_x, self.grid_y = np.mgrid[gxmin:gxmax+1e-6:self.step_size, gymin:gymax+1e-6:self.step_size]

  def useROI(self, use_roi):
    self.use_roi = use_roi

  def getROI(self):
    fig     = plt.figure()
    ax      = fig.add_subplot(111)
    surfdat = self.getSedSurface(mask=False,zcutoff=None)
    surf    = self.plotSurface(ax, surfdat)
    mc      = polyedit.MaskCreator(ax, poly_xy=self.roi_points)

    #Add the tree to the picture
    patch    = patches.PathPatch(self.treepath, facecolor='#FFA400', lw=1) #Create a patch from the path
    ax.add_patch(patch)                                                    #Draw the patch

    #Plot sediment points we interpolated from
    ax.plot(self.sed_data['x'], self.sed_data['y'], 'k.', ms=3)                      

    #Show the MaskCreator
    plt.show()

    self.roi         = mc.get_mask(self.grid_x, self.grid_y)
    #self.roi_points = mc.get_pts_in_path(list(zip(si.sed_data['x'],si.sed_data['y'])))
    self.roi_points  = mc.verts
    #self.roi_points = mc.get_pts_in_path(self.sed_data_arr[:,0:2])

  #Returns <Included, Excluded>
  def getPoints(self, zcutoff=None):
    allpts = self.sed_data_arr.copy()
    hull   = ConvexHull(allpts)
    hull   = allpts[hull.vertices,]
    if zcutoff is None:
      return allpts, []
    else:
      included = allpts[allpts[:,2]>=zcutoff,:]
      included = np.concatenate((included,hull), axis=0)
      excluded = allpts[allpts[:,2]< zcutoff,:]
      return included, excluded

  def getSedSurface(self, mask=True, zcutoff=None):
    inc, ex = self.getPoints(zcutoff)
    #Interpolate data
    surf = griddata(inc[:,0:2], inc[:,2] , (self.grid_x, self.grid_y), method='linear') #Also: 'linear'
    surf = surf.T

    # print('surf',surf.shape)
    # print('roi',self.roi.shape)
    if mask and not self.roi is None:
     surf[~self.roi] = np.nan

    return surf

  def plotSurface(self,axis,surf,alpha=1):
    surfimg = axis.imshow(surf, origin='lower', extent=(self.grid_x.min(),self.grid_x.max(),self.grid_y.min(),self.grid_y.max()), alpha=alpha)
    return surfimg

  def explore(self, zcutoff=None):
    fig = plt.figure()  #New figure

    inc, ex = self.getPoints(zcutoff)

    this_surf = self.getSedSurface(zcutoff=zcutoff)

    outer_grid = gridspec.GridSpec(6, 6, wspace=0, hspace=0.15, left=0.05, right=1, top=1, bottom=0)

    #Plot the interpolated sediment surface
    ax1 = plt.Subplot(fig, outer_grid[0:3,0:2])
    fig.add_subplot(ax1)
    surf = self.plotSurface(ax1,this_surf)
    surf.set_clim(self.zmin,self.zmax)

    #Add the tree to the picture
    patch = patches.PathPatch(self.treepath, facecolor='#FFA400', lw=1) #Create a patch from the path
    ax1.add_patch(patch)                                                #Draw the patch

    #Plot sediment points we interpolated from
    ax1.plot(inc[:,0], inc[:,1], 'k.', ms=3)

    ax2 = plt.Subplot(fig, outer_grid[3,0:2])
    fig.add_subplot(ax2)    
    ax2.plot(self.sed_data['x'], self.sed_data['z'], '.')
    ax2.set_title('X cross-section')
    if zcutoff:
     ax2.axhline(zcutoff, color='r')

    ax3 = plt.Subplot(fig, outer_grid[4,0:2])
    fig.add_subplot(ax3)    
    ax3.plot(self.sed_data['y'], self.sed_data['z'], '.')
    ax3.set_title('Y cross-section')
    if zcutoff:
      ax3.axhline(zcutoff, color='r')

    inner_grid = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=outer_grid[0:3,3:6], wspace=0.0, hspace=0.0)
    for i, cell in enumerate(inner_grid):
      zcutoff = self.zmin+i*self.zran/8
      ax = plt.Subplot(fig, cell)
      ax.set_title('cutoff = {0:.2}'.format(zcutoff))
      ax.axis('off')
      fig.add_subplot(ax)
      try:
        surf    = self.getSedSurface(zcutoff=zcutoff)
      except:
        print("Failed to get a surface: no bounding points at this zcutoff level!")
        continue
      surfimg = self.plotSurface(ax, surf)
      surfimg.set_clim(self.zmin,self.zmax)

    cutoff_comp = self.areaVzcutoff()
    axcut = plt.Subplot(fig, outer_grid[3:5,3:6])
    fig.add_subplot(axcut)
    axcut.plot(*list(zip(*cutoff_comp)), '.')
    axcut.set_title('Volume vs cut-off depth')

    #fig.tight_layout()
    plt.show()

  def showCutout(self, zcutoff=None):
    fig = plt.figure()  #New figure

    this_surf = self.getSedSurface(mask=False)

    inc, ex = self.getPoints(zcutoff)

    outer_grid = gridspec.GridSpec(2, 2, wspace=0.1, hspace=0.3, left=0.05, right=0.95, top=0.95, bottom=0.05)

    #Plot the interpolated sediment surface
    ax1 = plt.Subplot(fig, outer_grid[0,0])
    ax1.set_title('Overview')
    fig.add_subplot(ax1)
    surf = self.plotSurface(ax1,this_surf)
    surf.set_clim(self.zmin,self.zmax)

    #Add the tree to the picture
    patch = patches.PathPatch(self.treepath, facecolor='#FFA400', lw=1) #Create a patch from the path
    ax1.add_patch(patch)                                                #Draw the patch

    #Plot sediment points we interpolated from
    ax1.plot(inc[:,0], inc[:,1], 'k.', ms=3)    
    ax1.plot(ex[:,0], ex[:,1], 'r.', ms=3)

    ax1col = plt.Subplot(fig, outer_grid[0,0])
    plt.colorbar(surf,ax=ax1col)


    gxmin = min(si.roi_points[:,0])
    gymin = min(si.roi_points[:,1])
    gxmax = max(si.roi_points[:,0])
    gymax = max(si.roi_points[:,1])

    #Plot interpolated surface
    this_surf_zcut = self.getSedSurface(zcutoff=zcutoff)
    ax2 = plt.Subplot(fig, outer_grid[0,1])
    ax2.set_title('Depression-Filled')
    fig.add_subplot(ax2)
    surf_df = self.plotSurface(ax2,this_surf_zcut)
    surf_df.set_clim(self.zmin,self.zmax)

    #Add the tree to the picture
    patch = patches.PathPatch(self.treepath, facecolor='none', ec='black', lw=1, alpha=1) #Create a patch from the path
    ax2.add_patch(patch)                                                #Draw the patch

    #Plot sediment points we interpolated from
    ax2.plot(inc[:,0], inc[:,1], 'k.', ms=3)    
    # for i in range(len(inc[:,0])):
      # ax2.annotate('{0:.2}'.format(inc[i,2]), xy=(inc[i,0], inc[i,1]), textcoords='data')    
    ax2.set_xlim(gxmin,gxmax)
    ax2.set_ylim(gymin,gymax)
    #Plot excluded data points
    ax2.plot(ex[:,0], ex[:,1], 'r.', ms=3)
    # for i in range(len(ex[:,0])):
      # ax2.annotate('{0:.2}'.format(ex[i,2]), xy=(ex[i,0], ex[i,1]), textcoords='data')

    ax2col = plt.Subplot(fig, outer_grid[0,1])
    plt.colorbar(surf,ax=ax2col)

    #Plot original surface
    ax3 = plt.Subplot(fig, outer_grid[1,0])
    ax3.set_title('Original')
    fig.add_subplot(ax3)
    surf = self.plotSurface(ax3,this_surf)
    surf.set_clim(self.zmin,self.zmax)

    #Add the tree to the picture
    patch = patches.PathPatch(self.treepath, facecolor='none', ec='black', lw=1) #Create a patch from the path
    ax3.add_patch(patch)                                                         #Draw the patch

    #Plot sediment points we interpolated from
    ax3.plot(inc[:,0], inc[:,1], 'k.', ms=3)    
    ax3.plot(ex[:,0], ex[:,1], 'r.', ms=3)
    ax3.set_xlim(gxmin,gxmax)
    ax3.set_ylim(gymin,gymax)

    ax3col = plt.Subplot(fig, outer_grid[1,0])
    plt.colorbar(surf,ax=ax3col)

    #Plot difference surface
    ax4 = plt.Subplot(fig, outer_grid[1,1])
    fig.add_subplot(ax4)
    diffsurf = np.abs(this_surf_zcut-this_surf)
    diffsurf[diffsurf<0.001] = np.nan
    surf = self.plotSurface(ax4,diffsurf)
    diffvol = np.nansum((self.step_size**2)*diffsurf)
    ax4.set_title('Difference {0:.3} m³'.format(diffvol))
    #surf.set_clim(self.zmin,self.zmax)

    #Add the tree to the picture
    patch = patches.PathPatch(self.treepath, facecolor='#FFA400', lw=1) #Create a patch from the path
    ax4.add_patch(patch)                                                #Draw the patch

    #Plot sediment points we interpolated from
    ax4.plot(self.sed_data['x'], self.sed_data['y'], 'k.', ms=3)    
    ax4.plot(ex[:,0], ex[:,1], 'r.', ms=3)
    ax4.set_xlim(gxmin,gxmax)
    ax4.set_ylim(gymin,gymax)
    surf.set_cmap('Reds')

    ax4col = plt.Subplot(fig, outer_grid[1,1])
    plt.colorbar(surf,ax=ax4col)

    plt.show()              


  def areaVzcutoff(self, count=10):
    vals  = []
    orig_surf = self.getSedSurface()
    for i in range(count):
      zcutoff = self.zmin+i*self.zran/(count-1)
      try:
        cursurf = self.getSedSurface(zcutoff=zcutoff)
      except:
        print("Failed to find points at zcutoff={0}".format(zcutoff))
        continue
      diff    = cursurf-orig_surf
      vals.append((zcutoff,np.nansum(diff)))

    return vals

  def plot3D(self, zcutoff=None):
    fig   = plt.figure()
    ax    = fig.gca(projection='3d', aspect='equal')
    surf  = ax.plot_surface(self.grid_x, self.grid_y, self.sed_surf)



if len(sys.argv)!=3:
  print("Syntax: {0} <Sediment file> <Tree file>".format(sys.argv[0]))
  sys.exit(-1)


sedfile  = sys.argv[1]
treefile = sys.argv[2]

si = ScourInator(sedfile, treefile)

si.getROI()
#si.showCutout(-0.1)
#si.explore()

#si.plotAll(-0.5)




    #fig.add_subplot(ax2)    
    #ax2.imshow(self.roi, origin='lower', extent=(self.grid_x.min(),self.grid_x.max(),self.grid_y.min(),self.grid_y.max()))
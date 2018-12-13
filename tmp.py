
def visualize_average_run(sourceRoot, targetRoot, animalID, sessID, runList, smooth_fwhm=None, magRatio_thresh=None,\
    analysisDir=None,motionCorrection=False, flip = False, modify_range=True, mask =None):
        anatSource=os.path.join(targetRoot,'Surface')
        motionDir=os.path.join(targetRoot,'Motion')
        motionFileDir=os.path.join(motionDir, 'Registration')


        fileInDir=os.path.join(analysisDir,'phaseDecoding','SingleConditionData','Files')
        figOutDirRoot=os.path.join(analysisDir,'phaseDecoding','SingleConditionData','Figures')
        fileOutDirRoot=os.path.join(analysisDir,'phaseDecoding','SingleConditionData','Files')


        smoothString=''
        threshString=''

        condList = get_condition_list(sourceRoot,sessID,runList)


        for condCount,cond in enumerate(condList):
        print('Current Condition: %s'%(cond))
        cond = condList[condCount]
        figOutDir=os.path.join(figOutDirRoot,'cond%s'%(str(int(cond))))
        if not os.path.exists(figOutDir):
            os.makedirs(figOutDir)

        #LOAD MAPS
        fileName = '%s_cond%i_maps.npz'%(sessID, cond)
        f=np.load(os.path.join(fileInDir,fileName))
        phaseMap=f['phaseMap']
        magRatioMap=f['magRatioMap']

        if smooth_fwhm is not None:
            phaseMap=smooth_array(phaseMap,smooth_fwhm,phaseArray=True)
            magRatioMap=smooth_array(magRatioMap,smooth_fwhm)
            smoothString='_fwhm_'+str(smooth_fwhm)


        #set phase map range for visualization
        if modify_range:
            phaseMapDisplay=np.copy(phaseMap)
            phaseMapDisplay[phaseMap<0]=-phaseMap[phaseMap<0]
            phaseMapDisplay[phaseMap>0]=(2*np.pi)-phaseMap[phaseMap>0]

            rangeMin=0
            rangeMax=2*np.pi
        else:
            phaseMapDisplay=np.copy(phaseMap)
            rangeMin=-np.pi
            rangeMax=np.pi


        #apply threshhold
        if magRatio_thresh is not None:
            phaseMapDisplay[magRatioMap<magRatio_thresh]=np.nan
            threshString='_thresh_'+str(magRatio_thresh)
        else:
            magRatiothresh = np.max(magRatioMap)
            phaseMapDisplay[magRatioMap<magRatio_thresh]=np.nan
            threshString='_thresh_'+str(magRatio_thresh)


        #load surface for overlay
        #READ IN SURFACE


        imFile=anatSource+'/frame0_registered.tiff'
        if not os.path.isfile(imFile):
            imFile=anatSource+'/frame0.tiff'

        imSurf=cv2.imread(imFile,-1)
        szY,szX=imSurf.shape
        imSurf=np.true_divide(imSurf,2**12)*2**8

        if flip:
            print('Flipping Images')
            imSurf = np.fliplr(imSurf)
            phaseMapDisplay = np.fliplr(phaseMapDisplay)

        if motionCorrection:
            #LOAD MOTION CORRECTED BOUNDARIES
            inFile=motionFileDir+'/'+sessID+'_motionCorrectedBoundaries.npz'
            f=np.load(inFile)
            boundaries=f['boundaries']
            padDown=int(boundaries[0])
            padUp=int(szY-boundaries[1])
            padLeft=int(boundaries[2])
            padRight=int(szX-boundaries[3])

            phaseMapDisplay=np.pad(phaseMapDisplay,((padDown,padUp),(padLeft,padRight)),'constant',constant_values=((np.nan, np.nan),(np.nan,np.nan)))
        #plot
        fileName = 'overlay_images_%s_cond%s%s%s.png'%(sessID,str(int(cond)),smoothString,threshString)


        dpi = 80
        szY,szX = imSurf.shape
        # What size does the figure need to be in inches to fit the image?
        figsize = szX / float(dpi), szY / float(dpi)

        # Create a figure of the right size with one axes that takes up the full figure
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0, 0, 1, 1])

        # Hide spines, ticks, etc.
        ax.axis('off')

        ax.imshow(imSurf, 'gray')
        ax.imshow(phaseMapDisplay,'nipy_spectral',alpha=.5,vmin=rangeMin,vmax=rangeMax)

        fig.savefig(os.path.join(figOutDir,fileName), dpi=dpi, transparent=True)
        plt.close()

        #output masked image as well, if indicated
        if mask is not None:
            #load mask
            maskFile=targetRoot+'/masks/Files/'+mask+'.npz'
            f=np.load(maskFile)
            maskM=f['maskM']

            #apply mask
            phaseMapDisplay[maskM==0]=np.nan

            #plot
            outFile=outFile = '%s_cond%s%s%s_phaseMap_mask_%s_image.png'%\
            (figOutDir+sessID,str(int(cond)),smoothString,threshString,mask)

            #Create a figure of the right size with one axes that takes up the full figure
            fig = plt.figure(figsize=figsize)
            ax = fig.add_axes([0, 0, 1, 1])

            # Hide spines, ticks, etc.
            ax.axis('off')
            ax.imshow(imSurf, 'gray')
            ax.imshow(phaseMapDisplay,'nipy_spectral',alpha=.5,vmin=rangeMin,vmax=rangeMax)

            fig.savefig(outFile, dpi=dpi, transparent=True)
            plt.close()


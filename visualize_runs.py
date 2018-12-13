import os
cwd = os.getcwd()
import sys
sys.path.insert(0, cwd)
print('here')
from retino_functions import * 

import optparse

def get_comma_separated_args(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))



def visualize_single_run(sourceRoot, targetRoot, animalID, sessID, runList, smooth_fwhm=None, magRatio_thresh=None,\
    analysisDir=None,motionCorrection=False, flip = False, modify_range=True, mask =None):

    anatSource=os.path.join(targetRoot,'Surface')
    motionDir=os.path.join(targetRoot,'Motion')
    motionFileDir=os.path.join(motionDir, 'Registration')


    fileInDir=os.path.join(analysisDir,'SingleRunData','Files')
    figOutDirRoot=os.path.join(analysisDir,'SingleRunData','Figures')
    fileOutDirRoot=os.path.join(analysisDir,'SingleRunData','Files')
    #for file name
    smoothString=''
    threshString=''

    condList = get_condition_list(sourceRoot,animalID,sessID,runList)

    # runCount = 0
    # run = runList[runCount]

    for runCount,run in enumerate(runList):
        print('Current Run: %s'%(run))
        cond = condList[runCount]
        figOutDir=os.path.join(figOutDirRoot,'cond%s'%(str(int(cond))))
        if not os.path.exists(figOutDir):
            os.makedirs(figOutDir)

        #LOAD MAPS
        fileName = '%s_%s_map.npz'%(sessID, run)
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
        fileName = 'overlay_images_%s_%s%s%s.png'%(sessID,run,smoothString,threshString)


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
            maskFile=targetRoot+'/Sessions/'+sessID+'/masks/Files/'+mask+'.npz'
            f=np.load(maskFile)
            maskM=f['maskM']

            #apply mask
            phaseMapDisplay[maskM==0]=np.nan

            #plot
            outFile=outFile = '%s_%s%s%s_phaseMap_mask_%s_image.png'%\
            (figOutDir+sessID,run,smoothString,threshString,mask)

            #Create a figure of the right size with one axes that takes up the full figure
            fig = plt.figure(figsize=figsize)
            ax = fig.add_axes([0, 0, 1, 1])

            # Hide spines, ticks, etc.
            ax.axis('off')
            ax.imshow(imSurf, 'gray')
            ax.imshow(phaseMapDisplay,'nipy_spectral',alpha=.5,vmin=rangeMin,vmax=rangeMax)

            fig.savefig(outFile, dpi=dpi, transparent=True)
            plt.close()




def visualize_run(options):
    root = options.rootdir
    animalid = options.animalid
    project_dir = options.project_dir
    session = options.session
    run_list = options.run_list

    source_root = os.path.join(root,'raw_data',project_dir,animalid,session)
    target_root = os.path.join(root,'analyzed_data',project_dir,animalid,session)

    framerate, stimfreq = get_run_parameters(source_root, animalid, session, run_list[0])

    interp = options.interpolate
    exclude_edges= options.exclude_edges
    rolling_mean= options.rolling_mean
    time_average = options.time_average
    if time_average is not None:
        time_average = int(time_average)
    motion = options.motion

    ratio_thresh = options.ratio_thresh
    if ratio_thresh is not None:
        ratio_thresh = float(ratio_thresh)
    smooth_fwhm = options.smooth_fwhm
    if smooth_fwhm is not None:
        smooth_fwhm = int(smooth_fwhm)
    flip = options.flip

    analysis_root = os.path.join(target_root,'Analyses')
    analysis_dir=get_analysis_path_phase(analysis_root, stimfreq, interp, exclude_edges, rolling_mean, \
    motion, time_average)

    visualize_single_run(source_root, target_root, animalid, session, run_list, smooth_fwhm, ratio_thresh, analysis_dir, motion,flip)

def extract_options(options):

      parser = optparse.OptionParser()

      # PATH opts:
      parser.add_option('-R', '--root', action='store', dest='rootdir', default='/n/coxfs01/widefield-data', help='data root dir (root project dir containing all animalids) [default: /widefield]')
      parser.add_option('-p', '--project', action = 'store', dest = 'project_dir', default = 'Retinotopy/phase_encoding/Images_Cartesian_Constant', help = 'project directoy [default: Retinotopy/phase_enconding/Images_Cartesian_Constant]')
      parser.add_option('-i', '--animalid', action='store', dest='animalid', default='', help='Animal ID')
      parser.add_option('-S', '--session', action='store', dest='session', default='', help='session dir (format: YYYMMDD')
      parser.add_option('-r', '--run_list', action='callback', dest='run_list', default='',type='string',callback=get_comma_separated_args, help='comma-separated names of run dirs containing tiffs to be processed (ex: run1, run2, run3)')


      #specifications of analysis to visualize
      parser.add_option('-m', action='store_true', dest='motion', default=False, help="use motion corrected data")
      parser.add_option('-n', action='store_true', dest='interpolate', default=True, help="interpolate to an assumed steady frame rate")
      parser.add_option('-e', action='store_true', dest='exclude_edges', default=True, help="exclude first and last cycle of run")
      parser.add_option('-g', '--rollingmean', action='store_true', dest='rolling_mean', default=True, help='Boolean to indicate whether to subtract rolling mean from signal')
      parser.add_option('-w', '--timeaverage', action='store', dest='time_average', default=None, help='Size of time window with which to average frames (integer)')

      #visualization options
      parser.add_option('-f', '--fwhm', action='store', dest='smooth_fwhm', default=None, help='full-width at half-max size of kernel for smoothing')
      parser.add_option('-t', '--thresh', action='store', dest='ratio_thresh', default=None, help='magnitude ratio cut-off threshold')
      parser.add_option('-l', '--flip', action='store_true', dest='flip', default=False, help='boolean to indicate whether to perform horizontal flip on phase map images (to match actual orientation of FOV)')

      parser.add_option('--default', action='store_true', dest='default', default='store_false', help="Use all DEFAULT params, for params not specified by user (prevent interactive)")

      (options, args) = parser.parse_args(options)


      return options


def main(options):
      options = extract_options(options)
      visualize_run(options)


if __name__ == '__main__':
    main(sys.argv[1:])


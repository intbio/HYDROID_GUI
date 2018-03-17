from multiprocessing import Process
from hydroid.HYDROIDexp import assign_peaks_interactive
import tempfile,os

class hydroidConfig(object):
     def __init__(self,laneList):
        self.laneList=laneList
        #we use it like that as regular temfile can not be opend multiple times in windows
        fd,self.configFile=tempfile.mkstemp()
        with os.fdopen(fd,'w') as file:
            file.write('''#
#column - name of column with an array of profile values in lane_profile.xls
#lanme - arbitrary name of the gel lane (e.g. experiment label).
#leftlim, rightlim - indexes of values in profile array that specify the datarange that will be analyzed, in NaN dataset is not truncated
#peakthresh - parameter for automatic peak identification algorithm, use lower values to identify more subtle peaks. However, this may lead to false positives.
#min_dist_left - minimal distance between the individual peaks allowed at the left side of the data range.
#min_dist_right - minimal distance between the individual peaks allowed at the right side of the data range.
#segments - number of segment the data range is split into for linear interpolation of the min_dist value for the automatic peak identification algorithm.
#base - try to substract linear baseline from the profile before peak identification attempt.
#interpolate - try to guess the location of unidentified peaks from the local average spacing betwee the nearby peaks.
#alignpos - the index of position in the profile array used to align all the profiles if they are to be plotted simultaneously (i.e. to overlay Maxam-Gilbert profiles with OH-radical profiles), end - means use the last datapoint in the profile.
#seqpeak - approximate position of the peak that corresponds to a particular postion on DNA sequence specified by seqpos
#seqpos - position on DNA sequence that is attributed to a particular peak identified by seqpeak value.
#addpeaks - a list of peak positions to add manually to those identified automatically by the algorithm, format - values separated by spaces.
#delpeaks - a list of peak positions to be deleted manually to those identified automatically by the algorithm, format - values separaed by spaces.
## if any of the parameter values may be omitted or set to NaN, default gess values will be used.
##If multiple lines with the same lname are present the last one will be read.
column,	lname,				leftlim,	rightlim,	peakthresh,	min_dist_left,	min_dist_right,	segments,	base,	interpolate,	alignpos,	seqpeak,	seqpos,	addpeaks,	delpeaks \n''')
            for lane in self.laneList:
                file.write('%s,%s,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, , \n'%(lane,lane))

        
        
        
class PlotProcess(object):
    def __init__(self,**kwargs):
        FUNC=kwargs.pop('FUNC')
        if FUNC == 'assign_peaks_interactive':
            pf=kwargs.pop('lane_profile_file')
            cf=kwargs.pop('lane_config_file')
            self.plot_process = Process(
                target=assign_peaks_interactive,
                args=(pf,cf,kwargs['lane_name']))
                
        elif FUNC == 'call_peaks_interactive':
            pf=kwargs.pop('lane_profile_file')
            cf=kwargs.pop('lane_config_file')
            self.plot_process = Process(
                target=call_peaks_interactive,
                args=(pf,cf),
                kwargs=kwargs)
        else:
            print 'No action assigned'
            return
        self.plot_process.daemon = True
        self.plot_process.start()
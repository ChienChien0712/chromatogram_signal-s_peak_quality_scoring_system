import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import itertools
import pingouin as pg
import scipy
from multiprocessing import Process, Pool
import time
import math

from queue import Queue

class TargetedMSData():
  def __init__(self,  chromatogramFile, boundaryFile=None):
    start_time = time.time()
    self.chrom = pd.read_csv(chromatogramFile, sep='\t')
    chromatogram_data_load_time = time.time()
    print('Chromatogram data loaded in %.2f seconds' % (chromatogram_data_load_time - start_time))
    self.peak_boundary = None
    if boundaryFile is not None:
      self.peak_boundary = pd.read_csv(boundaryFile)
    peak_boundary_load_time = time.time()
    print('Peak boundary data loaed in %.2f seconds' % (peak_boundary_load_time - chromatogram_data_load_time))
  def getChromData(self, fileName, pepModSeq, nFragmentIon, start=None, end=None, transitions=None):
    sample = self.chrom[(self.chrom['FileName'] == fileName) & (self.chrom['PeptideModifiedSequence'] == pepModSeq)]
    if transitions is not None:
      sample = sample[(sample['FragmentIon'] + '.' +  sample['ProductCharge'].astype(str)).isin(transitions)]
    sample = sample.sort_values(by=['IsotopeLabelType','FragmentIon'], ascending=[False,True])
    
    # Check if all Times and Intensities have the same length, if not, pad 0s in front of the list
    max_len = 0
    for eachTime in list(sample['Times']):
      a = np.array(list(map(float,eachTime.split(','))))
      if len(a) > max_len:
        max_len = len(a)

    all_time = np.array([], dtype='float')
    for i in list(sample['Times']):
      a = np.array(list(map(float,i.split(','))))
      if len(a) == max_len:
        all_time = np.array(a, dtype='float')
        break

    # times = list(map(lambda x: x.split(','), sample['Times']))
    # times = list(map(lambda x: np.array(x, dtype=np.float), times))
    # sample['Times'] = times
    # intensities = list(map(lambda x: x.split(','), sample['Intensities']))
    # intensities = list(map(lambda x: np.array(x, dtype=np.float), intensities))
    # sample['Intensities'] = intensities
    #time
    #1. start
    #2. end
    #3. all_time
    #4. peak_time
    
    if start is None:
      if self.peak_boundary is None:
        start = np.random.choice(all_time, 1)[0]
      else:
        start = list(self.peak_boundary[(self.peak_boundary['File Name']==fileName)&(self.peak_boundary['Peptide Modified Sequence']==pepModSeq)]['Min Start Time'])[0]
    if end is None:
      if self.peak_boundary is None:
        temp_end = np.random.choice(all_time, 1)[0]
        while temp_end == start:
          temp_end = np.random.choice(all_time, 1)[0]
        end = temp_end
      else:
        end = list(self.peak_boundary[(self.peak_boundary['File Name']==fileName)&(self.peak_boundary['Peptide Modified Sequence']==pepModSeq)]['Max End Time'])[0]
    
    if start > end:
      temp = end
      end = start
      start = temp

    #all_time = np.array(list(map(float,list(sample['Times'])[0].split(','))))
    peak_filter = (all_time >= start) & (all_time <= end)
    peak_time = all_time[peak_filter]
    #get columns name
    total_col_name = ['.'.join(i) for i in list(zip(sample['FragmentIon'],map(str,sample['ProductCharge']),sample['IsotopeLabelType']))]

    #create peak_intensity
    intensity = []
    peak_intensity = []
    for i in list(sample['Intensities']):
      a = np.array(list(map(float,i.split(','))))
      if len(a) < max_len:
        diff = max_len - len(a)
        a = np.pad(a, (diff, 0), 'constant', constant_values=0)
      intensity.append(a)
      peak_intensity.append(a[peak_filter])
    intensity = pd.DataFrame(np.matrix(np.array(intensity, dtype='float')).T)
    peak_intensity = pd.DataFrame(np.matrix(np.array(peak_intensity, dtype='float')).T)
    #peak_intensity = pd.DataFrame(np.matrix([np.array(list(map(float,i.split(','))))[peak_filter] for i in list(sample['Intensities'])]).T)
    #intensity = pd.DataFrame(np.matrix([np.array(list(map(float,i.split(',')))) for i in list(sample['Intensities'])]).T)
    peak_intensity.columns = total_col_name
    intensity.columns = total_col_name

    #get column name for endogenous & sntandard
    r = re.compile(".*light")
    endogenous_cols = list(filter(r.match,peak_intensity.columns))
    r = re.compile(".*heavy")
    standard_cols = list(filter(r.match,peak_intensity.columns))
    light_transitions = set(map(lambda x: '.'.join(x.split('.')[0:-1]), endogenous_cols))
    heavy_transitions = set(map(lambda x: '.'.join(x.split('.')[0:-1]), standard_cols))
    intersect_transitions = list(light_transitions & heavy_transitions)
    endogenous_cols = list(map(lambda x: x + '.light', intersect_transitions))
    standard_cols = list(map(lambda x: x + '.heavy', intersect_transitions))
    

    #get transition ion's order for light & heavy
    lightTransitionIonOrder = list(np.amax(peak_intensity[endogenous_cols], axis=0).sort_values(ascending=False).index)
    lightTransitionIonOrder = list(map(lambda x: '.'.join(x.split('.')[0:-1]), lightTransitionIonOrder))[0:nFragmentIon]
    heavyTransitionIonOrder = list(np.amax(peak_intensity[standard_cols], axis=0).sort_values(ascending=False).index)
    heavyTransitionIonOrder = list(map(lambda x: '.'.join(x.split('.')[0:-1]), heavyTransitionIonOrder))[0:nFragmentIon]
    
    #nFragmentIon
    maxIntensity_rank = peak_intensity[standard_cols].max().rank(ascending=False)
    maxIntensity_rank_dict = dict(zip(maxIntensity_rank.index,maxIntensity_rank.values))
    standard_cols_remained = [i for i in standard_cols if maxIntensity_rank_dict[i] <= nFragmentIon]
    endogenous_cols_remained = ['.'.join(i.split('.')[:-1])+'.light' for i in standard_cols_remained]
    


    #1. peak_intensity
    total_col_name = endogenous_cols_remained + standard_cols_remained
    peak_intensity = peak_intensity[total_col_name]
    intensity = intensity[total_col_name]
    #col_name = peak_intensity.columns
    #col_name = heavyTransitionIonOrder
    #2. peak_sum_intensity
    peak_sum_intensity = pd.DataFrame(zip(peak_intensity[endogenous_cols_remained].sum(1),peak_intensity[standard_cols_remained].sum(1)),columns=['sum.light','sum.heavy'])
    
    #3. Area2SumRatio (Note: this value is calculated for all transitions)
    area = []
    for c in total_col_name:
      area.append(scipy.integrate.trapz(peak_intensity[c],peak_time))
    area = pd.DataFrame(area).transpose()
    area.columns = peak_intensity.columns
    
    Area2SumRatio_endo = area[endogenous_cols_remained]/sum(area[endogenous_cols_remained].loc[0]) #Area2SumRatio要分開light Heavy
    Area2SumRatio_stand = area[standard_cols_remained]/sum(area[standard_cols_remained].loc[0]) #Area2SumRatio要分開light Heavy
    Area2SumRatio = pd.concat([Area2SumRatio_endo, Area2SumRatio_stand], axis=1)
    return dict(
      fileName = fileName,
      peptideModifiedSequence = pepModSeq,
      start = start, #left boundary RT
      end = end, #right boundary RT
      peak_time = peak_time, #peak signal RT
      peak_intensity = peak_intensity, #peak signal (stripped)
      peak_sum_intensity = peak_sum_intensity, #sum of same isotope(heavy, light) peak signal
      Area2SumRatio = Area2SumRatio, #AUC/sum(all transition's AUC) of peak signal for each transition
      lightTransitionIonOrder = lightTransitionIonOrder,
      heavyTransitionIonOrder = heavyTransitionIonOrder,
      endogenous_cols = endogenous_cols_remained,
      standard_cols = standard_cols_remained,
      intensity= intensity,
      time = all_time,
      transitions = intersect_transitions
    )
  
  
  

class TargetedMSQC():
  def __init__(self, chromatogramFile, boundaryFile, core_num=1, nFragmentIon=5):
    self.ms_data = TargetedMSData(chromatogramFile, boundaryFile)
    self.filename_list = np.unique(self.ms_data.chrom['FileName'].astype(str))
    self.peptide_list = np.unique(self.ms_data.chrom['PeptideModifiedSequence'].astype(str))
    self.chrom_dict = dict()
    self.core_num = core_num
    self.quality_indexes = dict()
    self.nFragmentIon = nFragmentIon

  def selectChromData(self, fileName, pepModSeq, transitions=None):
    if transitions is not None:
      try:
        return self.ms_data.getChromData(fileName, pepModSeq, self.nFragmentIon, transitions = transitions)
      except Exception as e:
        print(e)
        return None
    if fileName in self.chrom_dict and pepModSeq in self.chrom_dict[fileName]:
      # Using cached data
      chromData = self.chrom_dict[fileName][pepModSeq]
    else:
      try:
        chromData = self.ms_data.getChromData(fileName, pepModSeq, self.nFragmentIon)
        if fileName not in self.chrom_dict:
          self.chrom_dict[fileName] = dict()
        self.chrom_dict[fileName][pepModSeq] = chromData
      except:
        chromData = None
    return chromData
  def updateBoundary(self, chrom=None, fn=None, ps=None, start=None, end=None):
    try:
      if chrom is None:
        chrom_data = self.ms_data.getChromData(fn, ps, self.nFragmentIon, start=start, end=end)
      else:
        chrom_data = chrom
        chrom_data['start'] = start
        chrom_data['end'] = end
        all_time = chrom_data['time']
        all_intensity = chrom_data['intensity']
        peak_filter = (all_time >= start) & (all_time <= end)
        peak_intensity = chrom_data['intensity'][peak_filter].reset_index(drop=True)
        endogenous_cols = chrom_data['endogenous_cols']
        standard_cols = chrom_data['standard_cols']
        peak_time = all_time[peak_filter]
        chrom_data['peak_time'] = peak_time
        chrom_data['peak_intensity'] = peak_intensity
        chrom_data['peak_sum_intensity'] = pd.DataFrame(zip(peak_intensity[endogenous_cols].sum(1),peak_intensity[standard_cols].sum(1)),columns=['sum.light','sum.heavy'])
        
        area = []
        for c in peak_intensity.columns:
          area.append(scipy.integrate.trapz(peak_intensity[c],peak_time))
        area = pd.DataFrame(area).transpose()
        area.columns = peak_intensity.columns

        Area2SumRatio_endo = area[endogenous_cols]/sum(area[endogenous_cols].loc[0]) #Area2SumRatio要分開light Heavy
        Area2SumRatio_stand = area[standard_cols]/sum(area[standard_cols].loc[0]) #Area2SumRatio要分開light Heavy
        chrom_data['Area2SumRatio'] = pd.concat([Area2SumRatio_endo, Area2SumRatio_stand], axis=1)
      if chrom_data is None:
        return None
      fn = chrom_data['fileName']
      ps = chrom_data['peptideModifiedSequence']
      if fn not in self.chrom_dict:
        self.chrom_dict[fn] = dict()
      self.chrom_dict[fn][ps] = chrom_data
    except Exception as e:
      print(e)
      print('Error')
      chrom_data = None
    return chrom_data
  
  def plotChromData(self, chrom = None, fn=None, ps=None, transitions=None):
    plt.rcParams['figure.figsize'] = (12.0, 6.0)
    if chrom is None:
      chrom_data = self.selectChromData(fn, ps, transitions=transitions)
    else:
      chrom_data = chrom
    if chrom_data is None:
      return None
    fn = chrom_data['fileName']
    ps = chrom_data['peptideModifiedSequence']
    fig, (ax1, ax2) = plt.subplots(2, 1)
    for ilt in ['light', 'heavy']:  
      key = 'heavyTransitionIonOrder'
      time = chrom_data['time']
      max_intensity = []
      intensities = []
      for fragment in chrom_data[key]:
        intensity = list(chrom_data['intensity'][fragment + '.' + ilt])
        intensities.append(intensity)
        max_intensity.append(int(max(intensity)))
      exp = len(str(max(max_intensity))) - 2    
      for i in range(len(chrom_data[key])):
        if ilt == 'light':
          ax1.plot(time, intensities[i],  linewidth=1, label=chrom_data[key][i] + '.' + ilt)
        else:
          ax2.plot(time, intensities[i],  linewidth=1, label=chrom_data[key][i] + '.' + ilt)
      if ilt == 'light':
        ax1.set_title(fn + ': ' + ps, fontsize=8)
        ax1.set_ylabel('Intensity(10^%s)'%exp)
        ax1.axvline(chrom_data['start'], linestyle= '--', linewidth=1,color='black')
        ax1.axvline(chrom_data['end'], linestyle= '--', linewidth=1,color='black')
        ax1.legend(fontsize=8,loc='best')
      elif ilt == 'heavy':
        ax2.set_ylabel('Intensity(10^%s)'%exp)
        ax2.set_xlabel('Retention Time')
        ax2.axvline(chrom_data['start'], linestyle= '--', linewidth=1,color='black')
        ax2.axvline(chrom_data['end'], linestyle= '--', linewidth=1,color='black')
        ax2.legend(fontsize=8,loc='best')
    fig.show()
    return fig

  # Intensity(3)
  def calculateTransitionMaxIntensity(self, chromatogram_data, sum_transition=False):
    key = 'peak_sum_intensity' if sum_transition==True else 'peak_intensity'
    peak_intensity_mx = chromatogram_data[key]
    return pd.DataFrame(peak_intensity_mx.max()).transpose()
  def calculateTransitionMaxBoundaryIntensity(self, chromatogram_data, sum_transition=False):
    key = 'peak_sum_intensity' if sum_transition==True else 'peak_intensity'
    peak_intensity_mx = chromatogram_data[key]
    result = peak_intensity_mx.loc[(0, peak_intensity_mx.index[-1]),:].max()
    return pd.DataFrame(result).transpose()
  def calculateTransitionMaxBoundaryIntensityNormalized(self, chromatogram_data, sum_transition=False):
    return self.calculateTransitionMaxBoundaryIntensity(chromatogram_data, sum_transition)/self.calculateTransitionMaxIntensity(chromatogram_data, sum_transition)

  # Shift(4)
  def calculatePairShift(self, chromatogram_data, sum_transition=False): #4
    key = 'peak_sum_intensity' if sum_transition==True else 'peak_intensity'
    peak_intensity_mx = chromatogram_data[key]
    peak_time = chromatogram_data['peak_time']
      
    fragmentIon_list, idx = np.unique(['.'.join(i.split('.')[:-1]) for i in peak_intensity_mx.columns],return_index=True)
    fragmentIon_list = fragmentIon_list[np.argsort(idx)]
    PairShift = []
    for fi in fragmentIon_list:
      light_max_intensity = max(peak_intensity_mx[fi+'.light'])
      light_max_time = peak_time[peak_intensity_mx[fi+'.light'][peak_intensity_mx[fi+'.light'] == light_max_intensity].index[0]]
      heavy_max_intensity = max(peak_intensity_mx[fi+'.heavy'])
      heavy_max_time = peak_time[peak_intensity_mx[fi+'.heavy'][peak_intensity_mx[fi+'.heavy'] == heavy_max_intensity].index[0]]  

      PairShift.append(round(abs(light_max_time-heavy_max_time) / (peak_time[-1]-peak_time[0]),4))
    PairShift = pd.DataFrame(PairShift).transpose()
    PairShift.columns = fragmentIon_list
    return PairShift
  def calculateTransitionShift(self, chromatogram_data, sum_transition=False): #matrix
    key = 'peak_sum_intensity' if sum_transition==True else 'peak_intensity'
    peak_intensity_mx = chromatogram_data[key]
    peak_time = chromatogram_data['peak_time']
      
    shift_mx = np.empty((len(peak_intensity_mx.columns),len(peak_intensity_mx.columns)))
    shift_mx[:] = np.nan
    for a,b in list(itertools.combinations(peak_intensity_mx.columns, 2)):
      a_max_intensity = max(peak_intensity_mx[a])
      a_max_time = peak_time[peak_intensity_mx[a][peak_intensity_mx[a] == a_max_intensity].index[0]]
      b_max_intensity = max(peak_intensity_mx[b])
      b_max_time = peak_time[peak_intensity_mx[b][peak_intensity_mx[b] == b_max_intensity].index[0]]    
      shift = round(abs(a_max_time-b_max_time) / (peak_time[-1]-peak_time[0]),4)
      shift_mx[list(peak_intensity_mx.columns).index(a)][list(peak_intensity_mx.columns).index(b)] = shift
      shift_mx[list(peak_intensity_mx.columns).index(b)][list(peak_intensity_mx.columns).index(a)] = shift

    all_max_time = []
    for a in peak_intensity_mx.columns:
      a_max_intensity = max(peak_intensity_mx[a])
      a_max_time = peak_time[peak_intensity_mx[a][peak_intensity_mx[a] == a_max_intensity].index[0]]    
      all_max_time.append(a_max_time)
    all_max_time = np.array(all_max_time)
    median_max_time = np.median(all_max_time)
    shift_diag = np.round(abs(all_max_time-median_max_time) / (peak_time[-1]-peak_time[0]),4)
    np.fill_diagonal(shift_mx, shift_diag)
    shift_mx = pd.DataFrame(shift_mx)
    shift_mx.columns = peak_intensity_mx.columns
    shift_mx.index = peak_intensity_mx.columns
    return shift_mx
  def calculateTransitionShift_diagonal(self, chromatogram_data, sum_transition=False): #8
    TransitionShift = self.calculateTransitionShift(chromatogram_data, sum_transition)
    diag = np.matrix(TransitionShift).diagonal().tolist()[0]
    transition_name = TransitionShift.columns
    TransitionShift_diag = pd.DataFrame(diag).T
    TransitionShift_diag.columns = transition_name
    return TransitionShift_diag
  def calculateIsotopeShift(self, chromatogram_data, sum_transition=False): #2
    diag = np.matrix(self.calculateTransitionShift(chromatogram_data, sum_transition)).diagonal().tolist()[0]
    light_shift = np.mean(diag[:int(len(diag)/2)])
    heavy_shift = np.mean(diag[int(len(diag)/2):])
    shift_mx = pd.DataFrame(zip([light_shift],[heavy_shift]))
    shift_mx.columns = ['light','heavy']
    return shift_mx
  def calculatePeakGroupShift(self, chromatogram_data, sum_transition=False): #1
    diag = np.matrix(self.calculateTransitionShift(chromatogram_data, sum_transition)).diagonal().tolist()[0]
    return np.mean(diag)
  
  # Jaggedness(3)
  def calculateTransitionJaggedness(self, chromatogram_data, sum_transition=False, flatness_factor = 0.05): #8
    key = 'peak_sum_intensity' if sum_transition==True else 'peak_intensity'
    peak_intensity_mx = chromatogram_data[key]
    peak_diff_mx = pd.DataFrame(np.diff(peak_intensity_mx,axis=0),columns=peak_intensity_mx.columns)
    peak_diff_mx[abs(peak_diff_mx) < flatness_factor*abs(peak_intensity_mx).max()] = 0
    jaggedness = ((abs(pd.DataFrame(np.diff(np.sign(peak_diff_mx),axis=0),columns=peak_intensity_mx.columns)) > 1).sum() - 1) / (len(peak_diff_mx)-1)
    jaggedness = pd.DataFrame(jaggedness).transpose()
    jaggedness.loc[jaggedness.shape[0]] = [0]*jaggedness.shape[1]
    jaggedness = pd.DataFrame(jaggedness.max()).transpose()
    jaggedness = np.round(jaggedness,4)    
    return jaggedness

  def calculateIsotopeJaggedness(self, chromatogram_data, sum_transition=False, flatness_factor = 0.05): #2
    jaggedness = self.calculateTransitionJaggedness(chromatogram_data, sum_transition, flatness_factor)
    r = re.compile(".*light")
    col_name_light = list(filter(r.match,jaggedness.columns))
    jaggedness_light = np.mean(jaggedness[col_name_light].loc[0])
    r = re.compile(".*heavy")
    col_name_heavy = list(filter(r.match,jaggedness.columns))
    jaggedness_heavy = np.mean(jaggedness[col_name_heavy].loc[0])
    jaggedness_mx = pd.DataFrame(zip([jaggedness_light],[jaggedness_heavy]))
    jaggedness_mx.columns = ['light','heavy']
    return jaggedness_mx
  def calculatePeakGroupJaggedness(self, chromatogram_data, sum_transition=False, flatness_factor = 0.05): #1
    jaggedness = self.calculateTransitionJaggedness(chromatogram_data, sum_transition, flatness_factor)
    return  np.mean(jaggedness.loc[0])

  # Similarity(3)
  def calculatePairSimilarity(self, chromatogram_data, sum_transition=False): #matrix
    key = 'peak_sum_intensity' if sum_transition==True else 'peak_intensity'
    peak_intensity_mx = chromatogram_data[key]    
    similarity_mx = peak_intensity_mx.rcorr(method='pearson', decimals=7, padjust='holm',stars=False)
    col_name = similarity_mx.columns
    similarity_mx = np.matrix(similarity_mx)
    np.fill_diagonal(similarity_mx,'0')
    similarity_mx = similarity_mx.astype('float')
    similarity_mx = similarity_mx + similarity_mx.T
    np.fill_diagonal(similarity_mx,1)
    similarity_mx = pd.DataFrame(similarity_mx)
    similarity_mx.columns = col_name
    similarity_mx.index = col_name    
    return similarity_mx
  def calculatePairSimilarity_IsotopePair(self, chromatogram_data, sum_transition=False): #4
    PairSimilarity = self.calculatePairSimilarity(chromatogram_data, sum_transition)
    col_name = PairSimilarity.columns
    pair_name, idx  = np.unique(['.'.join(i.split('.')[:-1]) for i in col_name],return_index=True)
    pair_name = pair_name[np.argsort(idx)]
    IsotopePairSimilarity = []
    for i in pair_name:
      IsotopePairSimilarity.append(PairSimilarity[i+'.light'].loc[i+'.heavy'])
    IsotopePairSimilarity = pd.DataFrame(IsotopePairSimilarity).T
    IsotopePairSimilarity.columns = pair_name    
    return IsotopePairSimilarity
  def calculateIsotopeSimilarity(self, chromatogram_data, sum_transition=False): #2
    similarity_mx = self.calculatePairSimilarity(chromatogram_data, sum_transition)
    similarity_light = np.mean(np.matrix(similarity_mx.filter(like='light',axis=0).filter(like='light',axis=1)))
    similarity_heavy = np.mean(np.matrix(similarity_mx.filter(like='heavy',axis=0).filter(like='heavy',axis=1)))
    similarity_mx = pd.DataFrame(zip([similarity_light],[similarity_heavy]))
    similarity_mx.columns = ['light','heavy']
    return similarity_mx
  def calculatePeakGroupSimilarity(self, chromatogram_data, sum_transition=False): #1
    return np.mean(np.matrix(self.calculatePairSimilarity(chromatogram_data, sum_transition)))
  
  # Symmetry(3)
  def calculateTransitionSymmetry(self, chromatogram_data, sum_transition=False): #8
    key = 'peak_sum_intensity' if sum_transition==True else 'peak_intensity'
    peak_intensity_mx = chromatogram_data[key]       
    left = peak_intensity_mx.loc[:int(np.floor(len(peak_intensity_mx)/2))-1].reset_index(drop=True)
    right = peak_intensity_mx.loc[:int(np.floor(len(peak_intensity_mx)/2))+1:-1].reset_index(drop=True)
    return pd.DataFrame(left.corrwith(right, axis = 0)).transpose()    
  def calculateIsotopeSymmetry(self, chromatogram_data, sum_transition=False): #2
    symmetry = self.calculateTransitionSymmetry(chromatogram_data, sum_transition)
    r = re.compile(".*light")
    col_name_light = list(filter(r.match,symmetry.columns))
    symmetry_light = np.mean(symmetry[col_name_light].loc[0])
    r = re.compile(".*heavy")
    col_name_heavy = list(filter(r.match,symmetry.columns))
    symmetry_heavy = np.mean(symmetry[col_name_heavy].loc[0])
    symmetry_mx = pd.DataFrame(zip([symmetry_light],[symmetry_heavy]))
    symmetry_mx.columns = ['light','heavy']
    return symmetry_mx
  def calculatePeakGroupSymmetry(self, chromatogram_data, sum_transition=False): #1
    symmetry = self.calculateTransitionSymmetry(chromatogram_data, sum_transition)
    return  np.mean(symmetry.loc[0])

  # FWHM(8)
  def calc_fwhm(self, sig,time):
    peakmax = max(sig)
    try:
      left_index = (sig[(sig - peakmax/2) > 0].index[0]-1, sig[(sig - peakmax/2) > 0].index[0])
    except:
      left_index = (np.nan,np.nan)
    try:
      right_index = (sig[(sig - peakmax/2) > 0].index[-1], sig[(sig - peakmax/2) > 0].index[-1]+1)
    except:
      right_index = (np.nan,np.nan)
        
    if (left_index[0] == -1) or (np.isnan(left_index[0])):
      t_left = time[0]
    else:
      t_left = (time[left_index[1]] - time[left_index[0]])/(sig[left_index[1]] - sig[left_index[0]])*(peakmax/2 - sig[left_index[0]]) + time[left_index[0]]

    if (right_index[1] > (len(time)-1)) or (np.isnan(right_index[1])):
      t_right = time[-1]
    else:
      t_right = (time[right_index[1]] - time[right_index[0]])/(sig[right_index[1]] - sig[right_index[0]])*(peakmax/2 - sig[right_index[0]]) + time[right_index[0]]
    fwhm = t_right - t_left
    return fwhm
      
  def calculateTransitionFWHM(self, chromatogram_data, sum_transition=False): #8
    key = 'peak_sum_intensity' if sum_transition==True else 'peak_intensity'
    peak_intensity_mx = chromatogram_data.get(key)
    time = chromatogram_data.get('peak_time')
    TransitionFWHM = []
    for col in peak_intensity_mx.columns:
      TransitionFWHM.append(self.calc_fwhm(peak_intensity_mx[col],time))
    TransitionFWHM = pd.DataFrame(TransitionFWHM).transpose()
    TransitionFWHM.columns = peak_intensity_mx.columns   
    return TransitionFWHM

  def calculateTransitionFWHM2base(self, chromatogram_data, sum_transition=False): #8
    TransitionFWHM = self.calculateTransitionFWHM(chromatogram_data, sum_transition)
    time = chromatogram_data.get('peak_time')
    return TransitionFWHM/(time[-1] - time[0])

  def calculateIsotopeFWHM(self, chromatogram_data, sum_transition=False): #2
    FWHM = self.calculateTransitionFWHM(chromatogram_data, sum_transition)
    r = re.compile(".*light")
    col_name_light = list(filter(r.match,FWHM.columns))
    FWHM_light = np.mean(FWHM[col_name_light].loc[0])
    r = re.compile(".*heavy")
    col_name_heavy = list(filter(r.match,FWHM.columns))
    FWHM_heavy = np.mean(FWHM[col_name_heavy].loc[0])
    FWHM_mx = pd.DataFrame(zip([FWHM_light],[FWHM_heavy]))
    FWHM_mx.columns = ['light','heavy']
    return FWHM_mx

  def calculateIsotopeFWHM2base(self, chromatogram_data, sum_transition=False): #2
    FWHM = self.calculateTransitionFWHM2base(chromatogram_data, sum_transition)
    r = re.compile(".*light")
    col_name_light = list(filter(r.match,FWHM.columns))
    FWHM_light = np.mean(FWHM[col_name_light].loc[0])
    r = re.compile(".*heavy")
    col_name_heavy = list(filter(r.match,FWHM.columns))
    FWHM_heavy = np.mean(FWHM[col_name_heavy].loc[0])
    FWHM_mx = pd.DataFrame(zip([FWHM_light],[FWHM_heavy]))
    FWHM_mx.columns = ['light','heavy']
    return FWHM_mx

  def calculatePeakGroupFWHM(self, chromatogram_data, sum_transition=False): #1
    FWHM = self.calculateTransitionFWHM(chromatogram_data, sum_transition)
    return  np.mean(FWHM.loc[0])

  def calculatePeakGroupFWHM2base(self, chromatogram_data, sum_transition=False): #1
    FWHM = self.calculateTransitionFWHM2base(chromatogram_data, sum_transition)
    return  np.mean(FWHM.loc[0])

  def calculatePairFWHMConsistency(self, chromatogram_data, sum_transition=False): #4
    FWHM = self.calculateTransitionFWHM(chromatogram_data, sum_transition)
    cols = sorted(list(set(['.'.join(i.split('.')[:-1]) for i in FWHM.columns])))
    FWHMConsistency = []
    for c in cols:
      FWHMConsistency.append(abs(FWHM[c+'.light'][0]-FWHM[c+'.heavy'][0])/FWHM[c+'.heavy'][0])
    FWHMConsistency = pd.DataFrame(FWHMConsistency).transpose()
    FWHMConsistency.columns = cols
    return FWHMConsistency

  def calculateMeanIsotopeFWHMConsistency(self, chromatogram_data, sum_transition=False): #8
    TransitionFWHM = self.calculateTransitionFWHM(chromatogram_data, sum_transition)
    TransitionFWHM_crossAll = TransitionFWHM.copy()
    pepModSeq = chromatogram_data['peptideModifiedSequence']
    #FileName_ls = np.unique(self.ms_data.chrom['FileName']) 
    #other_fn = list(FileName_ls)
    #other_fn.remove(chromatogram_data['fileName'])
    other_fn = list(self.filename_list)
    other_fn.remove((chromatogram_data['fileName']))
    for i in other_fn:
      #tmp = self.ms_data.getChromData(fileName=i,pepModSeq=pepModSeq,nFragmentIon=99)
      tmp = self.selectChromData(fileName=i, pepModSeq=pepModSeq)
      if tmp is None:
        continue
      tmp['peak_intensity'] = tmp['peak_intensity'][chromatogram_data['peak_intensity'].columns]
      tmp['peak_sum_intensity'] = tmp['peak_sum_intensity'][chromatogram_data['peak_sum_intensity'].columns]
      TransitionFWHM_crossAll = TransitionFWHM_crossAll.append(self.calculateTransitionFWHM(tmp, sum_transition))
    return abs(TransitionFWHM-np.mean(TransitionFWHM_crossAll))/np.mean(TransitionFWHM_crossAll)   
  
  #Modality(3)
  def calculateTransitionModality(self, chromatogram_data, sum_transition=False, flatness_factor=0.05): #8
    key = 'peak_sum_intensity' if sum_transition==True else 'peak_intensity'
    peak_intensity_mx = chromatogram_data.get(key)
    # find the differential of the peak    
    peak_diff_mx = pd.DataFrame(np.diff(peak_intensity_mx,axis=0),columns=peak_intensity_mx.columns)
    # any differences that are below the flatnessfactor of the maximum peak height are flattened.
    peak_diff_mx[abs(peak_diff_mx) < flatness_factor*abs(peak_intensity_mx).max()] = 0   
    # find the first and last timepoint where the differential changes sign
    first_fall = []
    last_rise = []
    for col_name in peak_diff_mx.columns:
      if len(peak_diff_mx[col_name][peak_diff_mx[col_name]<0]) > 0:
        first_fall.append(peak_diff_mx[col_name][peak_diff_mx[col_name]<0].index[0])
      else:
        first_fall.append(len(peak_intensity_mx[col_name]))
      if len(peak_diff_mx[col_name][peak_diff_mx[col_name]>0]) > 0:    
        last_rise.append(peak_diff_mx[col_name][peak_diff_mx[col_name]>0].index[-1])
      else:
        last_rise.append(-1)

      # if first fall is after last rise, peak cannot be bi or multi-modal, so max.dip is set to 0. Otherwise it is set to the largest fall or rise between the first fall and last rise
    modality_mx = []
    for f,l,col_name in zip(first_fall,last_rise,peak_diff_mx.columns):
      max_dip = 0
      if f < l:
        max_dip = max(abs(peak_diff_mx[col_name][list(range(f,l+1,1))]))

      # The output is the maximum dip normalized by the peak height
      if max(peak_intensity_mx[col_name])==0:
        modality = 0
      else:
        modality = max_dip/max(peak_intensity_mx[col_name])
      modality_mx.append(modality)  
    modality_mx = pd.DataFrame(modality_mx).transpose()
    modality_mx.columns = peak_intensity_mx.columns        
    return modality_mx
  def calculateIsotopeModality(self, chromatogram_data, sum_transition=False, flatness_factor=0.05): #2
    Modality = self.calculateTransitionModality(chromatogram_data, sum_transition)
    r = re.compile(".*light")
    col_name_light = list(filter(r.match,Modality.columns))
    Modality_light = np.mean(Modality[col_name_light].loc[0])
    r = re.compile(".*heavy")
    col_name_heavy = list(filter(r.match,Modality.columns))
    Modality_heavy = np.mean(Modality[col_name_heavy].loc[0])
    Modality_mx = pd.DataFrame(zip([Modality_light],[Modality_heavy]))
    Modality_mx.columns = ['light','heavy']
    return Modality_mx
  def calculatePeakGroupModality(self, chromatogram_data, sum_transition=False, flatness_factor=0.05): #1
    Modality = self.calculateTransitionModality(chromatogram_data, sum_transition)
    return  np.mean(Modality.loc[0])  
  
  # Area Ratio(4)
  def calculateArea2SumRatioCV(self, chromatogram_data): #8(sum補0)
    Area2SumRatio = chromatogram_data['Area2SumRatio']
    pepModSeq = chromatogram_data['peptideModifiedSequence']
    #FileName_ls = np.unique(self.ms_data.chrom['FileName'])
    #other_fn = list(FileName_ls)
    other_fn = list(self.filename_list)
    other_fn.remove(chromatogram_data['fileName'])
      #tmp = self.ms_data.getChromData(fileName=i,pepModSeq=pepModSeq,nFragmentIon=99)
    for i in other_fn:
      tmp = self.selectChromData(fileName = i, pepModSeq=pepModSeq)
      if tmp is None:
        continue
      tmp['Area2SumRatio'] = tmp['Area2SumRatio'][chromatogram_data['Area2SumRatio'].columns]
      Area2SumRatio = Area2SumRatio.append(tmp['Area2SumRatio'])
    return pd.DataFrame(np.std(Area2SumRatio,ddof=1)/np.mean(Area2SumRatio)).T
  def calculatePeakGroupRatioCorr(self, chromatogram_data): #1(sum補相同數值)
    Area2SumRatio = chromatogram_data['Area2SumRatio']
    FragmentIon = sorted(list(set(['.'.join(i.split('.')[:-1]) for i in Area2SumRatio.columns])))
    light_col = [fi+'.light' for fi in FragmentIon]
    heavy_col =[fi+'.heavy' for fi in FragmentIon]
    corr = np.corrcoef(Area2SumRatio[light_col].loc[0], Area2SumRatio[heavy_col].loc[0])[0][1]
    return corr
  def calculatePairRatioConsistency(self, chromatogram_data): #4(sum補0)
    Area2SumRatio = chromatogram_data['Area2SumRatio']
    cols = sorted(list(set(['.'.join(i.split('.')[:-1]) for i in Area2SumRatio.columns])))
    PairRatioConsistency = []
    for c in cols:
      PairRatioConsistency.append(abs(Area2SumRatio[c+'.light'][0]-Area2SumRatio[c+'.heavy'][0])/Area2SumRatio[c+'.heavy'][0])
    PairRatioConsistency = pd.DataFrame(PairRatioConsistency).transpose()
    PairRatioConsistency.columns = cols
    return PairRatioConsistency
  def calculateMeanIsotopeRatioConsistency(self, chromatogram_data): #8(sum補0)
    Area2SumRatio = chromatogram_data['Area2SumRatio']
    pepModSeq = chromatogram_data['peptideModifiedSequence']
    #FileName_ls = np.unique(self.ms_data.chrom['FileName'])
    #other_fn = list(FileName_ls)
    #other_fn.remove(chromatogram_data['fileName'])
    other_fn = list(self.filename_list)
    other_fn.remove(chromatogram_data['fileName'])
    for i in other_fn:
      #tmp = self.ms_data.getChromData(fileName=i,pepModSeq=pepModSeq,nFragmentIon=99)
      tmp = self.selectChromData(fileName = i, pepModSeq=pepModSeq)
      if tmp is None:
        continue
      tmp['Area2SumRatio'] = tmp['Area2SumRatio'][chromatogram_data['Area2SumRatio'].columns]
      Area2SumRatio = Area2SumRatio.append(tmp['Area2SumRatio'])
    return abs(chromatogram_data['Area2SumRatio']-np.mean(Area2SumRatio))/np.mean(Area2SumRatio)

  # RT(1)
  # PeakCenter = (MaxEndTime + MinStartTime)/2
  # MeanIsotopeRTConsistency = abs(PeakCenter - MeanPeakCenter)/MeanPeakCenter
  # MeanPeakCenter: mean of PeakCenter across all samples is calculated for each transition and isotope label.
  def calculateMeanIsotopeRTConsistency(self, chromatogram_data): #1(sum補相同數值)
    # FileName_ls = np.unique(self.ms_data.chrom['FileName'])
    other_fn = list(self.filename_list)
    other_fn.remove(chromatogram_data['fileName'])
    peakRT = [(chromatogram_data['start']+chromatogram_data['end'])/2]
    peakRT_crossAll = peakRT.copy()
    for i in other_fn:
      #tmp = self.ms_data.getChromData(fileName=i, pepModSeq=chromatogram_data['peptideModifiedSequence'], nFragmentIon=99)
      tmp = self.selectChromData(fileName=i, pepModSeq=chromatogram_data['peptideModifiedSequence'])
      if tmp is None:
        continue
      peakRT_crossAll.append((tmp['start']+tmp['end'])/2)
    return abs(peakRT[0]-np.mean(peakRT_crossAll))/np.mean(peakRT_crossAll)
# ----------------Export quality indexes-------------------------------------------------
  def getQualityIndexesByFileName(self, fn):
    peak_group_level = []
    isotope_level = []
    transition_pair_level = []
    transition_level = []
    for ps in self.peptide_list:
      try:
        chrom_data = self.selectChromData(fn, ps)
      except:
        chrom_data = None
      if chrom_data is not None:
        # Peak group level
        peakgroup_result = dict(
          PeakGroupRatioCorr = self.calculatePeakGroupRatioCorr(chrom_data),  
          PeakGroupJaggedness_mean = self.calculatePeakGroupJaggedness(chrom_data, flatness_factor=0.05),
          PeakGroupSymmetry_mean = self.calculatePeakGroupSymmetry(chrom_data),
          PeakGroupSimilarity_mean = self.calculatePeakGroupSimilarity(chrom_data), 
          PeakGroupModality_mean = self.calculatePeakGroupModality(chrom_data, flatness_factor=0.05), 
          PeakGroupShift_mean = self.calculatePeakGroupShift(chrom_data), 
          PeakGroupFWHM_mean = self.calculatePeakGroupFWHM(chrom_data), 
          PeakGroupFWHM2base_mean = self.calculatePeakGroupFWHM2base(chrom_data), 
          MeanIsotopeRTConsistency_acrossAll = self.calculateMeanIsotopeRTConsistency(chrom_data),
          PeakGroupJaggedness_sum = self.calculatePeakGroupJaggedness(chrom_data, sum_transition=True, flatness_factor=0.05),
          PeakGroupSymmetry_sum = self.calculatePeakGroupSymmetry(chrom_data, sum_transition=True), 
          PeakGroupSimilarity_sum = self.calculatePeakGroupSimilarity(chrom_data, sum_transition=True), 
          PeakGroupModality_sum = self.calculatePeakGroupModality(chrom_data, sum_transition=True, flatness_factor=0.05), 
          PeakGroupShift_sum = self.calculatePeakGroupShift(chrom_data, sum_transition=True),
          PeakGroupFWHM_sum = self.calculatePeakGroupFWHM(chrom_data, sum_transition=True),
          PeakGroupFWHM2base_sum = self.calculatePeakGroupFWHM2base(chrom_data, sum_transition=True)
        )
        peak_group_level.append(dict(
          FileName = fn,
          PeptideModifiedSequence = ps,
          PeakGroupRatioCorr = peakgroup_result['PeakGroupRatioCorr'],
          PeakGroupJaggedness_mean = peakgroup_result['PeakGroupJaggedness_mean'],
          PeakGroupSymmetry_mean = peakgroup_result['PeakGroupSymmetry_mean'],
          PeakGroupSimilarity_mean = peakgroup_result['PeakGroupSimilarity_mean'],
          PeakGroupModality_mean = peakgroup_result['PeakGroupModality_mean'],
          PeakGroupShift_mean = peakgroup_result['PeakGroupShift_mean'],
          PeakGroupFWHM_mean = peakgroup_result['PeakGroupFWHM_mean'],
          PeakGroupFWHM2base_mean = peakgroup_result['PeakGroupFWHM2base_mean'],
          MeanIsotopeRTConsistency_acrossAll = peakgroup_result['MeanIsotopeRTConsistency_acrossAll'],
          PeakGroupJaggedness_sum = peakgroup_result['PeakGroupJaggedness_sum'],
          PeakGroupSymmetry_sum = peakgroup_result['PeakGroupSymmetry_sum'],
          PeakGroupSimilarity_sum = peakgroup_result['PeakGroupSimilarity_sum'],
          PeakGroupModality_sum = peakgroup_result['PeakGroupModality_sum'],
          PeakGroupShift_sum = peakgroup_result['PeakGroupShift_sum'],
          PeakGroupFWHM_sum = peakgroup_result['PeakGroupFWHM_sum'],
          PeakGroupFWHM2base_sum = peakgroup_result['PeakGroupFWHM2base_sum']
        ))
        # Isotope level
        Jaggedness_mean = self.calculateIsotopeJaggedness(chrom_data, flatness_factor=0.05) 
        Symmetry_mean = self.calculateIsotopeSymmetry(chrom_data) 
        Similarity_mean = self.calculateIsotopeSimilarity(chrom_data) 
        Modality_mean = self.calculateIsotopeModality(chrom_data, flatness_factor=0.05) 
        Shift_mean = self.calculateIsotopeShift(chrom_data) 
        FWHN_mean = self.calculateIsotopeFWHM(chrom_data)
        FWHM2base_mean = self.calculateIsotopeFWHM2base(chrom_data) 
        Jaggedness_sum = self.calculateIsotopeJaggedness(chrom_data, sum_transition=True, flatness_factor=0.05) 
        Symmetry_sum = self.calculateIsotopeSymmetry(chrom_data, sum_transition=True) 
        Similarity_sum = self.calculateIsotopeSimilarity(chrom_data, sum_transition=True)
        Modality_sum = self.calculateIsotopeModality(chrom_data, sum_transition=True, flatness_factor=0.05) 
        Shift_sum = self.calculateIsotopeShift(chrom_data, sum_transition=True)
        FWHN_sum = self.calculateIsotopeFWHM(chrom_data, sum_transition=True) 
        FWHM2base_sum = self.calculateIsotopeFWHM2base(chrom_data, sum_transition=True)
        for isotope in ['light','heavy']:
          isotope_level.append(dict(
            FileName = fn,
            PeptideModifiedSequence = ps,
            Isotope = isotope,
            IsotopeJaggedness_mean = Jaggedness_mean[isotope][0],
            IsotopeSymmetry_mean = Symmetry_mean[isotope][0],
            IsotopeSimilarity_mean = Similarity_mean[isotope][0],
            IsotopeModality_mean = Modality_mean[isotope][0],
            IsotopeShift_mean = Shift_mean[isotope][0],
            IsotopeFWHM_mean = FWHN_mean[isotope][0],
            IsotopeFWHM2base_mean = FWHM2base_mean[isotope][0],
            IsotopeJaggedness_sum = Jaggedness_sum[isotope][0],
            IsotopeSymmetry_sum = Symmetry_sum[isotope][0],
            IsotopeSimilarity_sum = Similarity_sum[isotope][0],
            IsotopeModality_sum = Modality_sum[isotope][0],
            IsotopeShift_sum = Shift_sum[isotope][0],
            IsotopeFWHM_sum = FWHN_sum[isotope][0],
            IsotopeFWHM2base_sum = FWHM2base_sum[isotope][0]            
          ))

        for sum_transition_per_isotope in [False, True]:
          # Transition pair level
          TransitionPair_Name = np.unique(['.'.join(i.split('.')[:-1]) for i in chrom_data.get('peak_intensity').columns]) if sum_transition_per_isotope==False else ['sum']
          PairRatioConsistency = self.calculatePairRatioConsistency(chrom_data) if sum_transition_per_isotope==False else -1 #sum補-1 (-1)
          PairSimilarity = self.calculatePairSimilarity_IsotopePair(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum)
          PairShift = self.calculatePairShift(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum)
          PairFWHMConsistency = self.calculatePairFWHMConsistency(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum)
          for tp in TransitionPair_Name:
            transition_pair_level.append(dict(
              FileName = fn,
              PeptideModifiedSequence = ps,
              TransitionPair = tp,
              PairRatioConsistency = PairRatioConsistency[tp][0] if tp!='sum' else PairRatioConsistency,
              PairSimilarity = PairSimilarity[tp][0],
              PairShift = PairShift[tp][0],
              PairFWHMConsistency = PairFWHMConsistency[tp][0]
            ))
          
          #Transition level
          Transition_Name = chrom_data.get('peak_intensity').columns if sum_transition_per_isotope==False else ['sum.light', 'sum.heavy']
          TransitionJaggedness = self.calculateTransitionJaggedness(chrom_data, sum_transition=sum_transition_per_isotope, flatness_factor=0.05) #sum (sum.light | sum.heavy)
          TransitionSymmetry = self.calculateTransitionSymmetry(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
          TransitionModality = self.calculateTransitionModality(chrom_data, sum_transition=sum_transition_per_isotope, flatness_factor=0.05) #sum (sum.light | sum.heavy)
          TransitionShift = self.calculateTransitionShift_diagonal(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
          TransitionFWHM = self.calculateTransitionFWHM(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
          TransitionFWHM2base = self.calculateTransitionFWHM2base(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
          TransitionMaxIntensity = self.calculateTransitionMaxIntensity(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
          TransitionMaxBoundaryIntensity = self.calculateTransitionMaxBoundaryIntensity(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
          TransitionMaxBoundaryIntensityNormalized = self.calculateTransitionMaxBoundaryIntensityNormalized(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
          Area2SumRatioCV_acrossAll = self.calculateArea2SumRatioCV(chrom_data) if sum_transition_per_isotope==False else -1 #sum補-1 (-1 | -1)
          MeanIsotopeRatioConsistency_acrossAll = self.calculateMeanIsotopeRatioConsistency(chrom_data) if sum_transition_per_isotope==False else -1 #sum補-1 (-1 | -10) 
          MeanIsotopeFWHMConsistency_acrossAll = self.calculateMeanIsotopeFWHMConsistency(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
          for t in Transition_Name:
            isotope = t.split('.')[-1] # heavy | light
            transition = '.'.join(t.split('.')[0:-1])
            transition_level.append(dict(
              FileName = fn,
              PeptideModifiedSequence = ps,
              Isotope = isotope,
              Transition = transition,
              TransitionJaggedness = TransitionJaggedness[t][0],
              TransitionSymmetry = TransitionSymmetry[t][0],
              TransitionModality = TransitionModality[t][0],
              TransitionShift = TransitionShift[t][0],
              TransitionFWHM = TransitionFWHM[t][0],
              TransitionFWHM2base = TransitionFWHM2base[t][0],
              TransitionMaxIntensity = TransitionMaxIntensity[t][0],
              TransitionMaxBoundaryIntensity = TransitionMaxBoundaryIntensity[t][0],
              TransitionMaxBoundaryIntensityNormalized = TransitionMaxBoundaryIntensityNormalized[t][0],
              Area2SumRatioCV_acrossAll = Area2SumRatioCV_acrossAll[t][0] if transition!='sum' else Area2SumRatioCV_acrossAll,
              MeanIsotopeRatioConsistency_acrossAll = MeanIsotopeRatioConsistency_acrossAll[t][0] if transition!='sum' else MeanIsotopeRatioConsistency_acrossAll,
              MeanIsotopeFWHMConsistency_acrossAll = MeanIsotopeFWHMConsistency_acrossAll[t][0]
            ))
    return dict(
      peak_group_level = peak_group_level,
      isotope_level = isotope_level,
      transition_pair_level = transition_pair_level,
      transition_level = transition_level
    )
  
  def getQualityIndexesByModPepSeq(self, ps):
    peak_group_level = []
    isotope_level = []
    transition_pair_level = []
    transition_level = []
    for fn in self.filename_list:
      try:
        chrom_data = self.selectChromData(fn, ps)
      except:
        chrom_data = None
      if chrom_data is None:
        continue
      # Peak group level
      peakgroup_result = dict(
          PeakGroupRatioCorr = self.calculatePeakGroupRatioCorr(chrom_data),  
          PeakGroupJaggedness_mean = self.calculatePeakGroupJaggedness(chrom_data, flatness_factor=0.05),
          PeakGroupSymmetry_mean = self.calculatePeakGroupSymmetry(chrom_data),
          PeakGroupSimilarity_mean = self.calculatePeakGroupSimilarity(chrom_data), 
          PeakGroupModality_mean = self.calculatePeakGroupModality(chrom_data, flatness_factor=0.05), 
          PeakGroupShift_mean = self.calculatePeakGroupShift(chrom_data), 
          PeakGroupFWHM_mean = self.calculatePeakGroupFWHM(chrom_data), 
          PeakGroupFWHM2base_mean = self.calculatePeakGroupFWHM2base(chrom_data), 
          MeanIsotopeRTConsistency_acrossAll = self.calculateMeanIsotopeRTConsistency(chrom_data),
          PeakGroupJaggedness_sum = self.calculatePeakGroupJaggedness(chrom_data, sum_transition=True, flatness_factor=0.05),
          PeakGroupSymmetry_sum = self.calculatePeakGroupSymmetry(chrom_data, sum_transition=True), 
          PeakGroupSimilarity_sum = self.calculatePeakGroupSimilarity(chrom_data, sum_transition=True), 
          PeakGroupModality_sum = self.calculatePeakGroupModality(chrom_data, sum_transition=True, flatness_factor=0.05), 
          PeakGroupShift_sum = self.calculatePeakGroupShift(chrom_data, sum_transition=True),
          PeakGroupFWHM_sum = self.calculatePeakGroupFWHM(chrom_data, sum_transition=True),
          PeakGroupFWHM2base_sum = self.calculatePeakGroupFWHM2base(chrom_data, sum_transition=True)
      )
      peak_group_level.append(dict(
          FileName = fn,
          PeptideModifiedSequence = ps,
          PeakGroupRatioCorr = peakgroup_result['PeakGroupRatioCorr'],
          PeakGroupJaggedness_mean = peakgroup_result['PeakGroupJaggedness_mean'],
          PeakGroupSymmetry_mean = peakgroup_result['PeakGroupSymmetry_mean'],
          PeakGroupSimilarity_mean = peakgroup_result['PeakGroupSimilarity_mean'],
          PeakGroupModality_mean = peakgroup_result['PeakGroupModality_mean'],
          PeakGroupShift_mean = peakgroup_result['PeakGroupShift_mean'],
          PeakGroupFWHM_mean = peakgroup_result['PeakGroupFWHM_mean'],
          PeakGroupFWHM2base_mean = peakgroup_result['PeakGroupFWHM2base_mean'],
          MeanIsotopeRTConsistency_acrossAll = peakgroup_result['MeanIsotopeRTConsistency_acrossAll'],
          PeakGroupJaggedness_sum = peakgroup_result['PeakGroupJaggedness_sum'],
          PeakGroupSymmetry_sum = peakgroup_result['PeakGroupSymmetry_sum'],
          PeakGroupSimilarity_sum = peakgroup_result['PeakGroupSimilarity_sum'],
          PeakGroupModality_sum = peakgroup_result['PeakGroupModality_sum'],
          PeakGroupShift_sum = peakgroup_result['PeakGroupShift_sum'],
          PeakGroupFWHM_sum = peakgroup_result['PeakGroupFWHM_sum'],
          PeakGroupFWHM2base_sum = peakgroup_result['PeakGroupFWHM2base_sum']
      ))
      # Isotope level
      Jaggedness_mean = self.calculateIsotopeJaggedness(chrom_data, flatness_factor=0.05) 
      Symmetry_mean = self.calculateIsotopeSymmetry(chrom_data) 
      Similarity_mean = self.calculateIsotopeSimilarity(chrom_data) 
      Modality_mean = self.calculateIsotopeModality(chrom_data, flatness_factor=0.05) 
      Shift_mean = self.calculateIsotopeShift(chrom_data) 
      FWHN_mean = self.calculateIsotopeFWHM(chrom_data)
      FWHM2base_mean = self.calculateIsotopeFWHM2base(chrom_data) 
      Jaggedness_sum = self.calculateIsotopeJaggedness(chrom_data, sum_transition=True, flatness_factor=0.05) 
      Symmetry_sum = self.calculateIsotopeSymmetry(chrom_data, sum_transition=True) 
      Similarity_sum = self.calculateIsotopeSimilarity(chrom_data, sum_transition=True)
      Modality_sum = self.calculateIsotopeModality(chrom_data, sum_transition=True, flatness_factor=0.05) 
      Shift_sum = self.calculateIsotopeShift(chrom_data, sum_transition=True)
      FWHN_sum = self.calculateIsotopeFWHM(chrom_data, sum_transition=True) 
      FWHM2base_sum = self.calculateIsotopeFWHM2base(chrom_data, sum_transition=True)
      for isotope in ['light','heavy']:
        isotope_level.append(dict(
          FileName = fn,
          PeptideModifiedSequence = ps,
          Isotope = isotope,
          IsotopeJaggedness_mean = Jaggedness_mean[isotope][0],
          IsotopeSymmetry_mean = Symmetry_mean[isotope][0],
          IsotopeSimilarity_mean = Similarity_mean[isotope][0],
          IsotopeModality_mean = Modality_mean[isotope][0],
          IsotopeShift_mean = Shift_mean[isotope][0],
          IsotopeFWHM_mean = FWHN_mean[isotope][0],
          IsotopeFWHM2base_mean = FWHM2base_mean[isotope][0],
          IsotopeJaggedness_sum = Jaggedness_sum[isotope][0],
          IsotopeSymmetry_sum = Symmetry_sum[isotope][0],
          IsotopeSimilarity_sum = Similarity_sum[isotope][0],
          IsotopeModality_sum = Modality_sum[isotope][0],
          IsotopeShift_sum = Shift_sum[isotope][0],
          IsotopeFWHM_sum = FWHN_sum[isotope][0],
          IsotopeFWHM2base_sum = FWHM2base_sum[isotope][0]
        ))

      for sum_transition_per_isotope in [False, True]:
        # Transition pair level
        TransitionPair_Name = np.unique(['.'.join(i.split('.')[:-1]) for i in chrom_data.get('peak_intensity').columns]) if sum_transition_per_isotope==False else ['sum']
        PairRatioConsistency = self.calculatePairRatioConsistency(chrom_data) if sum_transition_per_isotope==False else 0 #sum補-1 (-1)
        PairSimilarity = self.calculatePairSimilarity_IsotopePair(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum)
        PairShift = self.calculatePairShift(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum)
        PairFWHMConsistency = self.calculatePairFWHMConsistency(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum)
        for tp in TransitionPair_Name:
          transition_pair_level.append(dict(
            FileName = fn,
            PeptideModifiedSequence = ps,
            TransitionPair = tp,
            PairRatioConsistency = PairRatioConsistency[tp][0] if tp!='sum' else PairRatioConsistency,
            PairSimilarity = PairSimilarity[tp][0],
            PairShift = PairShift[tp][0],
            PairFWHMConsistency = PairFWHMConsistency[tp][0]
          ))
        
        #Transition level
        Transition_Name = chrom_data.get('peak_intensity').columns if sum_transition_per_isotope==False else ['sum.light', 'sum.heavy']
        TransitionJaggedness = self.calculateTransitionJaggedness(chrom_data, sum_transition=sum_transition_per_isotope, flatness_factor=0.05) #sum (sum.light | sum.heavy)
        TransitionSymmetry = self.calculateTransitionSymmetry(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
        TransitionModality = self.calculateTransitionModality(chrom_data, sum_transition=sum_transition_per_isotope, flatness_factor=0.05) #sum (sum.light | sum.heavy)
        TransitionShift = self.calculateTransitionShift_diagonal(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
        TransitionFWHM = self.calculateTransitionFWHM(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
        TransitionFWHM2base = self.calculateTransitionFWHM2base(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
        TransitionMaxIntensity = self.calculateTransitionMaxIntensity(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
        TransitionMaxBoundaryIntensity = self.calculateTransitionMaxBoundaryIntensity(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
        TransitionMaxBoundaryIntensityNormalized = self.calculateTransitionMaxBoundaryIntensityNormalized(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
        Area2SumRatioCV_acrossAll = self.calculateArea2SumRatioCV(chrom_data) if sum_transition_per_isotope==False else -1 #sum補-1 (-1 | -1)
        MeanIsotopeRatioConsistency_acrossAll = self.calculateMeanIsotopeRatioConsistency(chrom_data) if sum_transition_per_isotope==False else -1 #sum補-1 (-1 | -10) 
        MeanIsotopeFWHMConsistency_acrossAll = self.calculateMeanIsotopeFWHMConsistency(chrom_data, sum_transition=sum_transition_per_isotope) #sum (sum.light | sum.heavy)
        for t in Transition_Name:
          isotope = t.split('.')[-1] # heavy | light
          transition = '.'.join(t.split('.')[0:-1])
          transition_level.append(dict(
            FileName = fn,
            PeptideModifiedSequence = ps,
            Isotope = isotope,
            Transition = transition,
            TransitionJaggedness = TransitionJaggedness[t][0],
            TransitionSymmetry = TransitionSymmetry[t][0],
            TransitionModality = TransitionModality[t][0],
            TransitionShift = TransitionShift[t][0],
            TransitionFWHM = TransitionFWHM[t][0],
            TransitionFWHM2base = TransitionFWHM2base[t][0],
            TransitionMaxIntensity = TransitionMaxIntensity[t][0],
            TransitionMaxBoundaryIntensity = TransitionMaxBoundaryIntensity[t][0],
            TransitionMaxBoundaryIntensityNormalized = TransitionMaxBoundaryIntensityNormalized[t][0],
            Area2SumRatioCV_acrossAll = Area2SumRatioCV_acrossAll[t][0] if transition!='sum' else Area2SumRatioCV_acrossAll,
            MeanIsotopeRatioConsistency_acrossAll = MeanIsotopeRatioConsistency_acrossAll[t][0] if transition!='sum' else MeanIsotopeRatioConsistency_acrossAll,
            MeanIsotopeFWHMConsistency_acrossAll = MeanIsotopeFWHMConsistency_acrossAll[t][0]
          ))
    return dict(
      peak_group_level = peak_group_level,
      isotope_level = isotope_level,
      transition_pair_level = transition_pair_level,
      transition_level = transition_level
    )
  
  def getQualityIndexesByModPepSeqLevel1(self, ps):
    transition_pair_level = []
    transition_level = []
    for fn in self.filename_list:
      try:
        chrom_data = self.selectChromData(fn, ps)
      except:
        chrom_data = None
      if chrom_data is None:
        continue
      try:
        #for sum_transition_per_isotope in [False, True]:
        # Transition pair level
        TransitionPair_Name = np.unique(['.'.join(i.split('.')[:-1]) for i in chrom_data.get('peak_intensity').columns])
        PairRatioConsistency = self.calculatePairRatioConsistency(chrom_data)
        PairSimilarity = self.calculatePairSimilarity_IsotopePair(chrom_data, sum_transition=False)
        PairShift = self.calculatePairShift(chrom_data, sum_transition=False)
        PairFWHMConsistency = self.calculatePairFWHMConsistency(chrom_data, sum_transition=False)
        for tp in TransitionPair_Name:
          transition_pair_level.append(dict(
            FileName = fn,
            PeptideModifiedSequence = ps,
            TransitionPair = tp,
            PairRatioConsistency = PairRatioConsistency[tp][0],
            PairSimilarity = PairSimilarity[tp][0],
            PairShift = PairShift[tp][0],
            PairFWHMConsistency = PairFWHMConsistency[tp][0]
          ))
        #Transition level
        Transition_Name = chrom_data.get('peak_intensity').columns
        TransitionJaggedness = self.calculateTransitionJaggedness(chrom_data, sum_transition=False, flatness_factor=0.05) #sum (sum.light | sum.heavy)
        TransitionSymmetry = self.calculateTransitionSymmetry(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
        TransitionModality = self.calculateTransitionModality(chrom_data, sum_transition=False, flatness_factor=0.05) #sum (sum.light | sum.heavy)
        TransitionShift = self.calculateTransitionShift_diagonal(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
        TransitionFWHM = self.calculateTransitionFWHM(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
        TransitionFWHM2base = self.calculateTransitionFWHM2base(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
        #TransitionMaxIntensity = self.calculateTransitionMaxIntensity(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
        #TransitionMaxBoundaryIntensity = self.calculateTransitionMaxBoundaryIntensity(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
        TransitionMaxBoundaryIntensityNormalized = self.calculateTransitionMaxBoundaryIntensityNormalized(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
        # Area2SumRatioCV_acrossAll = self.calculateArea2SumRatioCV(chrom_data) if False==False else -1 #sum補-1 (-1 | -1)
        # MeanIsotopeRatioConsistency_acrossAll = self.calculateMeanIsotopeRatioConsistency(chrom_data) if False==False else -1 #sum補-1 (-1 | -10) 
        # MeanIsotopeFWHMConsistency_acrossAll = self.calculateMeanIsotopeFWHMConsistency(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
        for t in Transition_Name:
          isotope = t.split('.')[-1] # heavy | light
          transition = '.'.join(t.split('.')[0:-1])
          transition_level.append(dict(
            FileName = fn,
            PeptideModifiedSequence = ps,
            Isotope = isotope,
            Transition = transition,
            TransitionJaggedness = TransitionJaggedness[t][0],
            TransitionSymmetry = TransitionSymmetry[t][0],
            TransitionModality = TransitionModality[t][0],
            TransitionShift = TransitionShift[t][0],
            TransitionFWHM = TransitionFWHM[t][0],
            TransitionFWHM2base = TransitionFWHM2base[t][0],
            #TransitionMaxIntensity = TransitionMaxIntensity[t][0],
            #TransitionMaxBoundaryIntensity = TransitionMaxBoundaryIntensity[t][0],
            TransitionMaxBoundaryIntensityNormalized = TransitionMaxBoundaryIntensityNormalized[t][0],
            # Area2SumRatioCV_acrossAll = Area2SumRatioCV_acrossAll[t][0] if transition!='sum' else Area2SumRatioCV_acrossAll,
            # MeanIsotopeRatioConsistency_acrossAll = MeanIsotopeRatioConsistency_acrossAll[t][0] if transition!='sum' else MeanIsotopeRatioConsistency_acrossAll,
            # MeanIsotopeFWHMConsistency_acrossAll = MeanIsotopeFWHMConsistency_acrossAll[t][0]
          ))
      except Exception as e:
        print(e)
        continue
    return dict(
      transition_pair_level = transition_pair_level,
      transition_level = transition_level
    )

  def getQualityIndexesByModPepSeqLevel2(self, ps):
    peak_group_level = []
    isotope_level = []
    for fn in self.filename_list:
      try:
        chrom_data = self.selectChromData(fn, ps)
      except:
        chrom_data = None
      if chrom_data is None:
        continue
      try:
        peak_group_level.append(dict(
          FileName = fn,
          PeptideModifiedSequence = ps,
          PeakGroupRatioCorr = self.calculatePeakGroupRatioCorr(chrom_data),
          PeakGroupJaggedness = self.calculatePeakGroupJaggedness(chrom_data, flatness_factor=0.05),
          PeakGroupSymmetry = self.calculatePeakGroupSymmetry(chrom_data),
          PeakGroupSimilarity = self.calculatePeakGroupSimilarity(chrom_data), 
          PeakGroupModality = self.calculatePeakGroupModality(chrom_data, flatness_factor=0.05), 
          PeakGroupShift = self.calculatePeakGroupShift(chrom_data), 
          PeakGroupFWHM = self.calculatePeakGroupFWHM(chrom_data), 
          PeakGroupFWHM2base = self.calculatePeakGroupFWHM2base(chrom_data)
        ))
        IsotopeJaggedness = self.calculateIsotopeJaggedness(chrom_data, flatness_factor=0.05) 
        IsotopeSymmetry = self.calculateIsotopeSymmetry(chrom_data) 
        IsotopeSimilarity = self.calculateIsotopeSimilarity(chrom_data) 
        IsotopeModality = self.calculateIsotopeModality(chrom_data, flatness_factor=0.05) 
        IsotopeShift = self.calculateIsotopeShift(chrom_data) 
        IsotopeFWHM = self.calculateIsotopeFWHM(chrom_data)
        IsotopeFWHM2base = self.calculateIsotopeFWHM2base(chrom_data)
        for isotope in ['light','heavy']:
          isotope_level.append(dict(
            FileName = fn,
            PeptideModifiedSequence = ps,
            Isotope = isotope,
            IsotopeJaggedness = IsotopeJaggedness[isotope][0],
            IsotopeSymmetry = IsotopeSymmetry[isotope][0],
            IsotopeSimilarity = IsotopeSimilarity[isotope][0],
            IsotopeModality = IsotopeModality[isotope][0],
            IsotopeShift = IsotopeShift[isotope][0],
            IsotopeFWHM = IsotopeFWHM[isotope][0],
            IsotopeFWHM2base = IsotopeFWHM2base[isotope][0]
          ))
      except Exception as e:
        print(e)
        continue
    return dict(
      peak_group_level = peak_group_level,
      isotope_level = isotope_level
    )

  def calculateQualityIndexesLevel1(self):
    start = time.time()
    pool = Pool(self.core_num)
    results = pool.map(self.getQualityIndexesByModPepSeqLevel1, self.peptide_list)
    pool.close()
    pool.join()
    end = time.time()
    print('Total calculation time for level 1 (transition and transtion pair): %.2f seconds' % (end - start))
    #print(end - start)
    transition_pair_level = []
    transition_level = []
    serializedData = []
    self.quality_indexes = dict()
    for eachResult in results:
      for eachTransitionPair in eachResult['transition_pair_level']:
        transition_pair_level.append(eachTransitionPair)
      for eachTransition in eachResult['transition_level']:
        transition_level.append(eachTransition)
    end2 = time.time()
    print('Processing results in %.2f seconds' % (end2 - end))

    #transition_level
    for eachItem in transition_level:
      fn = eachItem['FileName']
      ps = eachItem['PeptideModifiedSequence']
      isotope = eachItem['Isotope']
      transition = eachItem['Transition']
      if fn not in self.quality_indexes:
        self.quality_indexes[fn] = dict()
      if ps not in self.quality_indexes[fn]:
        self.quality_indexes[fn][ps] = dict()
      if 'transition_level' not in self.quality_indexes[fn][ps]:
        self.quality_indexes[fn][ps]['transition_level'] = dict()
      if isotope not in self.quality_indexes[fn][ps]['transition_level']:
        self.quality_indexes[fn][ps]['transition_level'][isotope] = dict()
      self.quality_indexes[fn][ps]['transition_level'][isotope][transition] = dict(
        TransitionJaggedness = eachItem['TransitionJaggedness'],
        TransitionSymmetry = eachItem['TransitionSymmetry'],
        TransitionModality = eachItem['TransitionModality'],
        TransitionShift = eachItem['TransitionShift'],
        TransitionFWHM = eachItem['TransitionFWHM'],
        TransitionFWHM2base = eachItem['TransitionFWHM2base'],
        #TransitionMaxIntensity = eachItem['TransitionMaxIntensity'],
        #TransitionMaxBoundaryIntensity = eachItem['TransitionMaxBoundaryIntensity'],
        TransitionMaxBoundaryIntensityNormalized = eachItem['TransitionMaxBoundaryIntensityNormalized'],
        # Area2SumRatioCV_acrossAll = eachItem['Area2SumRatioCV_acrossAll'],
        # MeanIsotopeRatioConsistency_acrossAll = eachItem['MeanIsotopeRatioConsistency_acrossAll'],
        # MeanIsotopeFWHMConsistency_acrossAll = eachItem['MeanIsotopeFWHMConsistency_acrossAll']
      )
    end3 = time.time()
    print('Storing transition level data in %.2f seconds' % (end3 - end2))

    #transition_pair_level
    for eachItem in transition_pair_level:
      fn = eachItem['FileName']
      ps = eachItem['PeptideModifiedSequence']
      tp = eachItem['TransitionPair']
      if fn not in self.quality_indexes:
        self.quality_indexes[fn] = dict()
      if ps not in self.quality_indexes[fn]:
        self.quality_indexes[fn][ps] = dict()
      if 'transition_pair_level' not in self.quality_indexes[fn][ps]:
        self.quality_indexes[fn][ps]['transition_pair_level'] = dict()
      self.quality_indexes[fn][ps]['transition_pair_level'][tp] = dict(
        PairRatioConsistency = eachItem['PairRatioConsistency'],
        PairSimilarity = eachItem['PairSimilarity'],
        PairShift = eachItem['PairShift'],
        PairFWHMConsistency = eachItem['PairFWHMConsistency']
      )
      if tp == 'sum':
        continue
      transition_level = self.quality_indexes[fn][ps]['transition_level']
      eachLightTransitionQuality = transition_level['light'][tp]
      eachHeavyTransitionQuality = transition_level['heavy'][tp]
      serializedData.append({
        'file': fn,
        'peptide': ps,
        'transition': tp,
        'data': [
          # light transition
          eachLightTransitionQuality['TransitionJaggedness'],
          eachLightTransitionQuality['TransitionSymmetry'],
          eachLightTransitionQuality['TransitionModality'],
          eachLightTransitionQuality['TransitionShift'],
          eachLightTransitionQuality['TransitionFWHM'],
          eachLightTransitionQuality['TransitionFWHM2base'],
          #eachLightTransitionQuality['TransitionMaxIntensity'],
          #eachLightTransitionQuality['TransitionMaxBoundaryIntensity'],
          eachLightTransitionQuality['TransitionMaxBoundaryIntensityNormalized'],
          # heavy transition
          eachHeavyTransitionQuality['TransitionJaggedness'],
          eachHeavyTransitionQuality['TransitionSymmetry'],
          eachHeavyTransitionQuality['TransitionModality'],
          eachHeavyTransitionQuality['TransitionShift'],
          eachHeavyTransitionQuality['TransitionFWHM'],
          eachHeavyTransitionQuality['TransitionFWHM2base'],
          #eachHeavyTransitionQuality['TransitionMaxIntensity'],
          #eachHeavyTransitionQuality['TransitionMaxBoundaryIntensity'],
          eachHeavyTransitionQuality['TransitionMaxBoundaryIntensityNormalized'],

          #eachItem['PairRatioConsistency'], # Used all transitions including heavy and light to calculate 
          eachItem['PairSimilarity'],
          eachItem['PairShift'],
          eachItem['PairFWHMConsistency']
        ]
      })
    end4 = time.time()
    print('Storing transition pair level data in %.2f seconds' % (end4 - end3))
    return serializedData
    
  def calculateQualityIndexesLevel2(self):
    start = time.time()
    pool = Pool(self.core_num)
    results = pool.map(self.getQualityIndexesByModPepSeqLevel2, self.peptide_list)
    pool.close()
    pool.join()
    end = time.time()
    print('Total calculation time for level 2 (isotope and peak group): %.2f seconds' % (end - start))
    #print(end - start)
    peak_group_level = []
    isotope_level = []
    serializedData = []
    self.quality_indexes = dict()
    for eachResult in results:
      for eachPeakGroup in eachResult['peak_group_level']:
        peak_group_level.append(eachPeakGroup)
      for eachIsotope in eachResult['isotope_level']:
        isotope_level.append(eachIsotope)
    end2 = time.time()
    print('Processing results in %.2f seconds' % (end2 - end))
    
    #Isotope level
    for eachItem in isotope_level:
      fn = eachItem['FileName']
      ps = eachItem['PeptideModifiedSequence']
      isotope = eachItem['Isotope']
      if fn not in self.quality_indexes:
        self.quality_indexes[fn] = dict()
      if ps not in self.quality_indexes[fn]:
        self.quality_indexes[fn][ps] = dict()
      if 'isotope_level' not in self.quality_indexes[fn][ps]:
        self.quality_indexes[fn][ps]['isotope_level'] = dict()
      self.quality_indexes[fn][ps]['isotope_level'][isotope] = dict(
        IsotopeJaggedness = eachItem['IsotopeJaggedness'],
        IsotopeSymmetry = eachItem['IsotopeSymmetry'],
        IsotopeSimilarity = eachItem['IsotopeSimilarity'],
        IsotopeModality = eachItem['IsotopeModality'],
        IsotopeShift = eachItem['IsotopeShift'],
        IsotopeFWHM = eachItem['IsotopeFWHM'],
        IsotopeFWHM2base = eachItem['IsotopeFWHM2base'],
      )
    end3 = time.time()
    print('Storing isotope level data in %.2f seconds' % (end3 - end2))
    #Peak group level
    for eachItem in peak_group_level:
      fn = eachItem['FileName']
      ps = eachItem['PeptideModifiedSequence']
      if fn not in self.quality_indexes:
        self.quality_indexes[fn] = dict()
      if ps not in self.quality_indexes[fn]:
        self.quality_indexes[fn][ps] = dict()
      self.quality_indexes[fn][ps]['peak_group_level'] = dict(
        PeakGroupRatioCorr = eachItem['PeakGroupRatioCorr'],
        PeakGroupJaggedness = eachItem['PeakGroupJaggedness'],
        PeakGroupSymmetry = eachItem['PeakGroupSymmetry'],
        PeakGroupSimilarity = eachItem['PeakGroupSimilarity'],
        PeakGroupModality = eachItem['PeakGroupModality'],
        PeakGroupShift = eachItem['PeakGroupShift'],
        PeakGroupFWHM = eachItem['PeakGroupFWHM'],
        PeakGroupFWHM2base = eachItem['PeakGroupFWHM2base']
      )
      isotopeLevelLight = self.quality_indexes[fn][ps]['isotope_level']['light']
      isotopeLevelHeavy = self.quality_indexes[fn][ps]['isotope_level']['heavy']
      serializedData.append({
        'file': fn,
        'peptide': ps,
        'data': [
          # Isotope level: light
          isotopeLevelLight['IsotopeJaggedness'],
          isotopeLevelLight['IsotopeSymmetry'],
          isotopeLevelLight['IsotopeSimilarity'],
          isotopeLevelLight['IsotopeModality'],
          isotopeLevelLight['IsotopeShift'],
          isotopeLevelLight['IsotopeFWHM'],
          isotopeLevelLight['IsotopeFWHM2base'],

          # Isotope level: heavy
          isotopeLevelHeavy['IsotopeJaggedness'],
          isotopeLevelHeavy['IsotopeSymmetry'],
          isotopeLevelHeavy['IsotopeSimilarity'],
          isotopeLevelHeavy['IsotopeModality'],
          isotopeLevelHeavy['IsotopeShift'],
          isotopeLevelHeavy['IsotopeFWHM'],
          isotopeLevelHeavy['IsotopeFWHM2base'],

          # Peak group level
          eachItem['PeakGroupRatioCorr'],
          eachItem['PeakGroupJaggedness'],
          eachItem['PeakGroupSymmetry'],
          eachItem['PeakGroupSimilarity'],
          eachItem['PeakGroupModality'],
          eachItem['PeakGroupShift'],
          eachItem['PeakGroupFWHM'],
          eachItem['PeakGroupFWHM2base']
        ]
      })
    end4 = time.time()
    print('Storing peak group level data in %.2f seconds' % (end4 - end3))
    return serializedData

  def serializeLevel1(self, chrom = None, transition=None, fn=None, ps=None, transitions=None):
    if chrom is None:
      chrom_data = self.selectChromData(fn, ps, transitions)
    else:
      chrom_data = chrom
    
    if chrom_data is None:
      return None
    
    #Transition Pair Level
    #PairRatioConsistency = self.calculatePairRatioConsistency(chrom_data)
    PairSimilarity = self.calculatePairSimilarity_IsotopePair(chrom_data, sum_transition=False) #sum (sum)
    PairShift = self.calculatePairShift(chrom_data, sum_transition=False) #sum (sum)
    PairFWHMConsistency = self.calculatePairFWHMConsistency(chrom_data, sum_transition=False) #sum (sum)
    #Transition Level
    TransitionJaggedness = self.calculateTransitionJaggedness(chrom_data, sum_transition=False, flatness_factor=0.05) #sum (sum.light | sum.heavy)
    TransitionSymmetry = self.calculateTransitionSymmetry(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
    TransitionModality = self.calculateTransitionModality(chrom_data, sum_transition=False, flatness_factor=0.05) #sum (sum.light | sum.heavy)
    TransitionShift = self.calculateTransitionShift_diagonal(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
    TransitionFWHM = self.calculateTransitionFWHM(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
    TransitionFWHM2base = self.calculateTransitionFWHM2base(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
    TransitionMaxBoundaryIntensityNormalized = self.calculateTransitionMaxBoundaryIntensityNormalized(chrom_data, sum_transition=False) #sum (sum.light | sum.heavy)
    return [
      # light transition
      TransitionJaggedness[transition + '.light'][0],
      TransitionSymmetry[transition + '.light'][0],
      TransitionModality[transition + '.light'][0],
      TransitionShift[transition + '.light'][0],
      TransitionFWHM[transition + '.light'][0],
      TransitionFWHM2base[transition + '.light'][0],
      TransitionMaxBoundaryIntensityNormalized[transition + '.light'][0],

      TransitionJaggedness[transition + '.heavy'][0],
      TransitionSymmetry[transition + '.heavy'][0],
      TransitionModality[transition + '.heavy'][0],
      TransitionShift[transition + '.heavy'][0],
      TransitionFWHM[transition + '.heavy'][0],
      TransitionFWHM2base[transition + '.heavy'][0],
      TransitionMaxBoundaryIntensityNormalized[transition + '.heavy'][0],
      #PairRatioConsistency[transition][0],
      PairSimilarity[transition][0],
      PairShift[transition][0],
      PairFWHMConsistency[transition][0]
    ]
  def serializeLevel2(self, chrom = None, fn=None, ps=None, transitions=None):
    if chrom is None:
      chrom_data = self.selectChromData(fn, ps, transitions)
    else:
      chrom_data = chrom
    
    if chrom_data is None:
      return None
    PeakGroupRatioCorr = self.calculatePeakGroupRatioCorr(chrom_data),
    PeakGroupJaggedness = self.calculatePeakGroupJaggedness(chrom_data, flatness_factor=0.05),
    PeakGroupSymmetry = self.calculatePeakGroupSymmetry(chrom_data),
    PeakGroupSimilarity = self.calculatePeakGroupSimilarity(chrom_data), 
    PeakGroupModality = self.calculatePeakGroupModality(chrom_data, flatness_factor=0.05), 
    PeakGroupShift = self.calculatePeakGroupShift(chrom_data), 
    PeakGroupFWHM = self.calculatePeakGroupFWHM(chrom_data), 
    PeakGroupFWHM2base = self.calculatePeakGroupFWHM2base(chrom_data)

    IsotopeJaggedness = self.calculateIsotopeJaggedness(chrom_data, flatness_factor=0.05) 
    IsotopeSymmetry = self.calculateIsotopeSymmetry(chrom_data) 
    IsotopeSimilarity = self.calculateIsotopeSimilarity(chrom_data) 
    IsotopeModality = self.calculateIsotopeModality(chrom_data, flatness_factor=0.05) 
    IsotopeShift = self.calculateIsotopeShift(chrom_data) 
    IsotopeFWHM = self.calculateIsotopeFWHM(chrom_data)
    IsotopeFWHM2base = self.calculateIsotopeFWHM2base(chrom_data)
    return [
      IsotopeJaggedness['light'][0],
      IsotopeSymmetry['light'][0],
      IsotopeSimilarity['light'][0],
      IsotopeModality['light'][0],
      IsotopeShift['light'][0],
      IsotopeFWHM['light'][0],
      IsotopeFWHM2base ['light'][0],

      IsotopeJaggedness['heavy'][0],
      IsotopeSymmetry['heavy'][0],
      IsotopeSimilarity['heavy'][0],
      IsotopeModality['heavy'][0],
      IsotopeShift['heavy'][0],
      IsotopeFWHM['heavy'][0],
      IsotopeFWHM2base ['heavy'][0],

      PeakGroupRatioCorr[0],
      PeakGroupJaggedness[0],
      PeakGroupSymmetry[0],
      PeakGroupSimilarity[0],
      PeakGroupModality[0],
      PeakGroupShift[0],
      PeakGroupFWHM[0],
      PeakGroupFWHM2base
    ]

  def serializeLevel1Title(self):
    return [
      'TransitionJaggedness_light',
      'TransitionSymmetry_light',
      'TransitionModality_light',
      'TransitionShift_light',
      'TransitionFWHM_light',
      'TransitionFWHM2base_light',
      #'TransitionMaxIntensity_light',
      #'TransitionMaxBoundaryIntensity_light',
      'TransitionMaxBoundaryIntensityNormalized_light',
      'TransitionJaggedness_heavy',
      'TransitionSymmetry_heavy',
      'TransitionModality_heavy',
      'TransitionShift_heavy',
      'TransitionFWHM_heavy',
      'TransitionFWHM2base_heavy',
      #'TransitionMaxIntensity_heavy',
      #'TransitionMaxBoundaryIntensity_heavy',
      'TransitionMaxBoundaryIntensityNormalized_heavy',
      #'PairRatioConsistency',
      'PairSimilarity',
      'PairShift',
      'PairFWHMConsistency'
    ]
  
  def serializeLevel2Title(self):
    return [
      'IsotopeJaggedness_light',
      'IsotopeSymmetry_light',
      'IsotopeSimilarity_light',
      'IsotopeModality_light',
      'IsotopeShift_light',
      'IsotopeFWHM_light',
      'IsotopeFWHM2base_light',

      'IsotopeJaggedness_heavy',
      'IsotopeSymmetry_heavy',
      'IsotopeSimilarity_heavy',
      'IsotopeModality_heavy',
      'IsotopeShift_heavy',
      'IsotopeFWHM_heavy',
      'IsotopeFWHM2base_heavy',

      'PeakGroupRatioCorr',
      'PeakGroupJaggedness',
      'PeakGroupSymmetry',
      'PeakGroupSimilarity',
      'PeakGroupModality',
      'PeakGroupShift',
      'PeakGroupFWHM',
      'PeakGroupFWHM2base'
    ]

  def calculateQualityIndexes(self):
    start = time.time()
    pool = Pool(self.core_num)
    
    results = pool.map(self.getQualityIndexesByModPepSeq, self.peptide_list)
    pool.close()
    pool.join()
    end = time.time()
    print('Total calculation time: %.2f seconds' % (end - start))
    #print(end - start)
    peak_group_level = []
    isotope_level = []
    transition_pair_level = []
    transition_level = []
    self.quality_indexes = dict()
    for eachResult in results:
      for eachPeak in eachResult['peak_group_level']:
        peak_group_level.append(eachPeak)
      for eachIsotope in eachResult['isotope_level']:
        isotope_level.append(eachIsotope)
      for eachTransitionPair in eachResult['transition_pair_level']:
        transition_pair_level.append(eachTransitionPair)
      for eachTransition in eachResult['transition_level']:
        transition_level.append(eachTransition)

    #peak_group_level
    for eachItem in peak_group_level:
      fn = eachItem['FileName']
      ps = eachItem['PeptideModifiedSequence']
      if fn not in self.quality_indexes:
        self.quality_indexes[fn] = dict()
      if ps not in self.quality_indexes[fn]:
        self.quality_indexes[fn][ps] = dict()
      self.quality_indexes[fn][ps]['peak_group_level'] = dict(
        PeakGroupRatioCorr = eachItem['PeakGroupRatioCorr'],
        PeakGroupJaggedness_mean = eachItem['PeakGroupJaggedness_mean'],
        PeakGroupSymmetry_mean = eachItem['PeakGroupSymmetry_mean'],
        PeakGroupSimilarity_mean = eachItem['PeakGroupSimilarity_mean'],
        PeakGroupModality_mean = eachItem['PeakGroupModality_mean'],
        PeakGroupShift_mean = eachItem['PeakGroupShift_mean'],
        PeakGroupFWHM_mean = eachItem['PeakGroupFWHM_mean'],
        PeakGroupFWHM2base_mean = eachItem['PeakGroupFWHM2base_mean'],
        MeanIsotopeRTConsistency_acrossAll = eachItem['MeanIsotopeRTConsistency_acrossAll'],
        PeakGroupJaggedness_sum = eachItem['PeakGroupJaggedness_sum'],
        PeakGroupSymmetry_sum = eachItem['PeakGroupSymmetry_sum'],
        PeakGroupSimilarity_sum = eachItem['PeakGroupSimilarity_sum'],
        PeakGroupModality_sum = eachItem['PeakGroupModality_sum'],
        PeakGroupShift_sum = eachItem['PeakGroupShift_sum'],
        PeakGroupFWHM_sum = eachItem['PeakGroupFWHM_sum'],
        PeakGroupFWHM2base_sum = eachItem['PeakGroupFWHM2base_sum']
      )

    #isotope_level
    for eachItem in isotope_level:
      fn = eachItem['FileName']
      ps = eachItem['PeptideModifiedSequence']
      isotope = eachItem['Isotope']
      if fn not in self.quality_indexes:
        self.quality_indexes[fn] = dict()
      if ps not in self.quality_indexes[fn]:
        self.quality_indexes[fn][ps] = dict()
      if 'isotope_level' not in self.quality_indexes[fn][ps]:
        self.quality_indexes[fn][ps]['isotope_level'] = dict()      
      self.quality_indexes[fn][ps]['isotope_level'][isotope] = dict(
        IsotopeJaggedness_mean = eachItem['IsotopeJaggedness_mean'],
        IsotopeSymmetry_mean = eachItem['IsotopeSymmetry_mean'],
        IsotopeSimilarity_mean = eachItem['IsotopeSimilarity_mean'],
        IsotopeModality_mean = eachItem['IsotopeModality_mean'],
        IsotopeShift_mean = eachItem['IsotopeShift_mean'],
        IsotopeFWHM_mean = eachItem['IsotopeFWHM_mean'],
        IsotopeFWHM2base_mean = eachItem['IsotopeFWHM2base_mean'],
        IsotopeJaggedness_sum = eachItem['IsotopeJaggedness_sum'],
        IsotopeSymmetry_sum = eachItem['IsotopeSymmetry_sum'],
        IsotopeSimilarity_sum = eachItem['IsotopeSimilarity_sum'],
        IsotopeModality_sum = eachItem['IsotopeModality_sum'],
        IsotopeShift_sum = eachItem['IsotopeShift_sum'],
        IsotopeFWHM_sum = eachItem['IsotopeFWHM_sum'],
        IsotopeFWHM2base_sum = eachItem['IsotopeFWHM2base_sum']
      )

    #transition_pair_level
    for eachItem in transition_pair_level:
      fn = eachItem['FileName']
      ps = eachItem['PeptideModifiedSequence']
      tp = eachItem['TransitionPair']
      if fn not in self.quality_indexes:
        self.quality_indexes[fn] = dict()
      if ps not in self.quality_indexes[fn]:
        self.quality_indexes[fn][ps] = dict()
      if 'transition_pair_level' not in self.quality_indexes[fn][ps]:
        self.quality_indexes[fn][ps]['transition_pair_level'] = dict()
      self.quality_indexes[fn][ps]['transition_pair_level'][tp] = dict(
        PairRatioConsistency = eachItem['PairRatioConsistency'],
        PairSimilarity = eachItem['PairSimilarity'],
        PairShift = eachItem['PairShift'],
        PairFWHMConsistency = eachItem['PairFWHMConsistency']
      )

    #transition_level
    for eachItem in transition_level:
      fn = eachItem['FileName']
      ps = eachItem['PeptideModifiedSequence']
      isotope = eachItem['Isotope']
      transition = eachItem['Transition']
      if fn not in self.quality_indexes:
        self.quality_indexes[fn] = dict()
      if ps not in self.quality_indexes[fn]:
        self.quality_indexes[fn][ps] = dict()
      if 'transition_level' not in self.quality_indexes[fn][ps]:
        self.quality_indexes[fn][ps]['transition_level'] = dict()
      if isotope not in self.quality_indexes[fn][ps]['transition_level']:
        self.quality_indexes[fn][ps]['transition_level'][isotope] = dict()
      self.quality_indexes[fn][ps]['transition_level'][isotope][transition] = dict(
        TransitionJaggedness = eachItem['TransitionJaggedness'],
        TransitionSymmetry = eachItem['TransitionSymmetry'],
        TransitionModality = eachItem['TransitionModality'],
        TransitionShift = eachItem['TransitionShift'],
        TransitionFWHM = eachItem['TransitionFWHM'],
        TransitionFWHM2base = eachItem['TransitionFWHM2base'],
        TransitionMaxIntensity = eachItem['TransitionMaxIntensity'],
        TransitionMaxBoundaryIntensity = eachItem['TransitionMaxBoundaryIntensity'],
        TransitionMaxBoundaryIntensityNormalized = eachItem['TransitionMaxBoundaryIntensityNormalized'],
        Area2SumRatioCV_acrossAll = eachItem['Area2SumRatioCV_acrossAll'],
        MeanIsotopeRatioConsistency_acrossAll = eachItem['MeanIsotopeRatioConsistency_acrossAll'],
        MeanIsotopeFWHMConsistency_acrossAll = eachItem['MeanIsotopeFWHMConsistency_acrossAll']
      )

    return dict(
      peak_group_level = pd.DataFrame(peak_group_level),
      isotope_level = pd.DataFrame(isotope_level),
      transition_pair_level = pd.DataFrame(transition_pair_level),
      transition_level = pd.DataFrame(transition_level)
    )
    #getQualityIndexesByFileName
  def getSerializedTitle(self):
    returnList = [
      'PeakGroupRatioCorr',
      'PeakGroupJaggedness_mean',
      'PeakGroupSymmetry_mean',
      'PeakGroupSimilarity_mean',
      'PeakGroupModality_mean',
      'PeakGroupShift_mean',
      'PeakGroupFWHM_mean',
      'PeakGroupFWHM2base_mean',
      'PeakGroupJaggedness_sum',
      'PeakGroupSymmetry_sum',
      'PeakGroupSimilarity_sum',
      'PeakGroupModality_sum',
      'PeakGroupShift_sum',
      'PeakGroupFWHM_sum',
      'PeakGroupFWHM2base_sum',
      'IsotopeJaggedness_mean',
      'IsotopeSymmetry_mean',
      'IsotopeSimilarity_mean',
      'IsotopeModality_mean',
      'IsotopeShift_mean',
      'IsotopeFWHM_mean',
      'IsotopeFWHM2base_mean',
      'IsotopeJaggedness_sum',
      'IsotopeSymmetry_sum',
      'IsotopeSimilarity_sum',
      'IsotopeModality_sum',
      'IsotopeShift_sum',
      'IsotopeFWHM_sum',
      'IsotopeFWHM2base_sum',
      'IsotopeJaggedness_mean',
      'IsotopeSymmetry_mean',
      'IsotopeSimilarity_mean',
      'IsotopeModality_mean',
      'IsotopeShift_mean',
      'IsotopeFWHM_mean',
      'IsotopeFWHM2base_mean',
      'IsotopeJaggedness_sum',
      'IsotopeSymmetry_sum',
      'IsotopeSimilarity_sum',
      'IsotopeModality_sum',
      'IsotopeShift_sum',
      'IsotopeFWHM_sum',
      'IsotopeFWHM2base_sum',
    ] #43
    for index in range(self.nFragmentIon): #4
      returnList.append('PairRatioConsistency-' + str(index))
      returnList.append('PairSimilarity-' + str(index))
      returnList.append('PairShift-' + str(index))
      returnList.append('PairFWHMConsistency-' + str(index))
      
    for index in range(self.nFragmentIon): #9
      returnList.append('TransitionJaggedness_light-' + str(index))
      returnList.append('TransitionSymmetry_light-' + str(index))
      returnList.append('TransitionModality_light-' + str(index))
      returnList.append('TransitionShift_light-' + str(index))
      returnList.append('TransitionFWHM_light-' + str(index))
      returnList.append('TransitionFWHM2base_light-' + str(index))
      returnList.append('TransitionMaxIntensity_light-' + str(index))
      returnList.append('TransitionMaxBoundaryIntensity_light-' + str(index))
      returnList.append('TransitionMaxBoundaryIntensityNormalized_light-' + str(index))

    for index in range(self.nFragmentIon): #9
      returnList.append('TransitionJaggedness_heavy-' + str(index))
      returnList.append('TransitionSymmetry_heavy-' + str(index))
      returnList.append('TransitionModality_heavy-' + str(index))
      returnList.append('TransitionShift_heavy-' + str(index))
      returnList.append('TransitionFWHM_heavy-' + str(index))
      returnList.append('TransitionFWHM2base_heavy-' + str(index))
      returnList.append('TransitionMaxIntensity_heavy-' + str(index))
      returnList.append('TransitionMaxBoundaryIntensity_heavy-' + str(index))
      returnList.append('TransitionMaxBoundaryIntensityNormalized_heavy-' + str(index))
    returnList.append('MeanIsotopeRTConsistency_acrossAll') #44
    for index in range(self.nFragmentIon): #3
      returnList.append('Area2SumRatioCV_acrossAll_light-' + str(index))
      returnList.append('MeanIsotopeRatioConsistency_acrossAll_light-' + str(index))
      returnList.append('MeanIsotopeFWHMConsistency_acrossAll_light-' + str(index))
    for index in range(self.nFragmentIon): #3
      returnList.append('Area2SumRatioCV_acrossAll_heavy-' + str(index))
      returnList.append('MeanIsotopeRatioConsistency_acrossAll_heavy-' + str(index))
      returnList.append('MeanIsotopeFWHMConsistency_acrossAll_heavy-' + str(index))
    return returnList #28*n + 44
  def getSerializedTitleWithTransition(self):
    returnList = [
      'PeakGroupRatioCorr',
      'PeakGroupJaggedness_mean',
      'PeakGroupSymmetry_mean',
      'PeakGroupSimilarity_mean',
      'PeakGroupModality_mean',
      'PeakGroupShift_mean',
      'PeakGroupFWHM_mean',
      'PeakGroupFWHM2base_mean',

      'PeakGroupJaggedness_sum',
      'PeakGroupSymmetry_sum',
      'PeakGroupSimilarity_sum',
      'PeakGroupModality_sum',
      'PeakGroupShift_sum',
      'PeakGroupFWHM_sum',
      'PeakGroupFWHM2base_sum',

      'IsotopeJaggedness_light_mean',
      'IsotopeSymmetry_light_mean',
      'IsotopeSimilarity_light_mean',
      'IsotopeModality_light_mean',
      'IsotopeShift_light_mean',
      'IsotopeFWHM_light_mean',
      'IsotopeFWHM2base_light_mean',

      'PeakGroupJaggedness_light_sum',
      'PeakGroupSymmetry_light_sum',
      'PeakGroupSimilarity_light_sum',
      'PeakGroupModality_light_sum',
      'PeakGroupShift_light_sum',
      'PeakGroupFWHM_light_sum',
      'PeakGroupFWHM2base_light_sum',

      'IsotopeJaggedness_heavy_mean',
      'IsotopeSymmetry_heavy_mean',
      'IsotopeSimilarity_heavy_mean',
      'IsotopeModality_heavy_mean',
      'IsotopeShift_heavy_mean',
      'IsotopeFWHM_heavy_mean',
      'IsotopeFWHM2base_heavy_mean',

      'IsotopeJaggedness_heavy_sum',
      'IsotopeSymmetry_heavy_sum',
      'IsotopeSimilarity_heavy_sum',
      'IsotopeModality_heavy_sum',
      'IsotopeShift_heavy_sum',
      'IsotopeFWHM_heavy_sum',
      'IsotopeFWHM2base_heavy_sum',

    ] #43
    returnList.append('PairRatioConsistency')
    returnList.append('PairSimilarity')
    returnList.append('PairShift')
    returnList.append('PairFWHMConsistency')
      
    returnList.append('TransitionJaggedness_light')
    returnList.append('TransitionSymmetry_light')
    returnList.append('TransitionModality_light')
    returnList.append('TransitionShift_light')
    returnList.append('TransitionFWHM_light')
    returnList.append('TransitionFWHM2base_light')
    returnList.append('TransitionMaxIntensity_light')
    returnList.append('TransitionMaxBoundaryIntensity_light')
    returnList.append('TransitionMaxBoundaryIntensityNormalized_light')

    returnList.append('TransitionJaggedness_heavy')
    returnList.append('TransitionSymmetry_heavy')
    returnList.append('TransitionModality_heavy')
    returnList.append('TransitionShift_heavy')
    returnList.append('TransitionFWHM_heavy')
    returnList.append('TransitionFWHM2base_heavy')
    returnList.append('TransitionMaxIntensity_heavy')
    returnList.append('TransitionMaxBoundaryIntensity_heavy')
    returnList.append('TransitionMaxBoundaryIntensityNormalized_heavy') #65
     
    
    #summation of all transitions
    returnList.append('PairSimilarity_sum')
    returnList.append('PairShift_sum')
    returnList.append('PairFWHMConsistency_sum')

    returnList.append('TransitionJaggedness_light_sum')
    returnList.append('TransitionSymmetry_light_sum')
    returnList.append('TransitionModality_light_sum')
    returnList.append('TransitionShift_light_sum')
    returnList.append('TransitionFWHM_light_sum')
    returnList.append('TransitionFWHM2base_light_sum')
    returnList.append('TransitionMaxIntensity_light_sum')
    returnList.append('TransitionMaxBoundaryIntensity_light_sum')
    returnList.append('TransitionMaxBoundaryIntensityNormalized_light_sum')

    returnList.append('TransitionJaggedness_heavy_sum')
    returnList.append('TransitionSymmetry_heavy_sum')
    returnList.append('TransitionModality_heavy_sum')
    returnList.append('TransitionShift_heavy_sum')
    returnList.append('TransitionFWHM_heavy_sum')
    returnList.append('TransitionFWHM2base_heavy_sum')
    returnList.append('TransitionMaxIntensity_heavy_sum')
    returnList.append('TransitionMaxBoundaryIntensity_heavy_sum')
    returnList.append('TransitionMaxBoundaryIntensityNormalized_heavy_sum')

    

    # acrossAll
    returnList.append('MeanIsotopeRTConsistency_acrossAll') #44
    returnList.append('Area2SumRatioCV_acrossAll_light')
    returnList.append('MeanIsotopeRatioConsistency_acrossAll_light')
    returnList.append('MeanIsotopeFWHMConsistency_acrossAll_light')
    returnList.append('Area2SumRatioCV_acrossAll_heavy')
    returnList.append('MeanIsotopeRatioConsistency_acrossAll_heavy')
    returnList.append('MeanIsotopeFWHMConsistency_acrossAll_heavy')
    return returnList #70

  def serializeWithTransition(self, fn, ps, transition):
    chrom_data = self.selectChromData(fn, ps)
    if chrom_data is None:
      return None
    peak_level = self.quality_indexes[fn][ps]['peak_group_level']
    isotope_light_level = self.quality_indexes[fn][ps]['isotope_level']['light']
    isotope_heavy_level = self.quality_indexes[fn][ps]['isotope_level']['heavy']
    transition_pair_level = self.quality_indexes[fn][ps]['transition_pair_level']
    transition_level = self.quality_indexes[fn][ps]['transition_level']
    returnList = [
      peak_level['PeakGroupRatioCorr'],
      peak_level['PeakGroupJaggedness_mean'],
      peak_level['PeakGroupSymmetry_mean'],
      peak_level['PeakGroupSimilarity_mean'],
      peak_level['PeakGroupModality_mean'],
      peak_level['PeakGroupShift_mean'],
      peak_level['PeakGroupFWHM_mean'],
      peak_level['PeakGroupFWHM2base_mean'],

      peak_level['PeakGroupJaggedness_sum'],
      peak_level['PeakGroupSymmetry_sum'],
      peak_level['PeakGroupSimilarity_sum'],
      peak_level['PeakGroupModality_sum'],
      peak_level['PeakGroupShift_sum'],
      peak_level['PeakGroupFWHM_sum'],
      peak_level['PeakGroupFWHM2base_sum'],

      isotope_light_level['IsotopeJaggedness_mean'],
      isotope_light_level['IsotopeSymmetry_mean'],
      isotope_light_level['IsotopeSimilarity_mean'],
      isotope_light_level['IsotopeModality_mean'],
      isotope_light_level['IsotopeShift_mean'],
      isotope_light_level['IsotopeFWHM_mean'],
      isotope_light_level['IsotopeFWHM2base_mean'],

      isotope_light_level['IsotopeJaggedness_sum'],
      isotope_light_level['IsotopeSymmetry_sum'],
      isotope_light_level['IsotopeSimilarity_sum'],
      isotope_light_level['IsotopeModality_sum'],
      isotope_light_level['IsotopeShift_sum'],
      isotope_light_level['IsotopeFWHM_sum'],
      isotope_light_level['IsotopeFWHM2base_sum'],

      isotope_heavy_level['IsotopeJaggedness_mean'],
      isotope_heavy_level['IsotopeSymmetry_mean'],
      isotope_heavy_level['IsotopeSimilarity_mean'],
      isotope_heavy_level['IsotopeModality_mean'],
      isotope_heavy_level['IsotopeShift_mean'],
      isotope_heavy_level['IsotopeFWHM_mean'],
      isotope_heavy_level['IsotopeFWHM2base_mean'],

      isotope_heavy_level['IsotopeJaggedness_sum'],
      isotope_heavy_level['IsotopeSymmetry_sum'],
      isotope_heavy_level['IsotopeSimilarity_sum'],
      isotope_heavy_level['IsotopeModality_sum'],
      isotope_heavy_level['IsotopeShift_sum'],
      isotope_heavy_level['IsotopeFWHM_sum'],
      isotope_heavy_level['IsotopeFWHM2base_sum']
    ]
    if transition not in transition_pair_level:
      return None
    eachTransitionPairQuality = transition_pair_level[transition]
    returnList.append(eachTransitionPairQuality['PairRatioConsistency'])
    returnList.append(eachTransitionPairQuality['PairSimilarity'])
    returnList.append(eachTransitionPairQuality['PairShift'])
    returnList.append(eachTransitionPairQuality['PairFWHMConsistency'])

    if (transition not in transition_level['light']) or (transition not in transition_level['heavy']):
      return None
    eachLightTransitionQuality = transition_level['light'][transition]
    returnList.append(eachLightTransitionQuality['TransitionJaggedness'])
    returnList.append(eachLightTransitionQuality['TransitionSymmetry'])
    returnList.append(eachLightTransitionQuality['TransitionModality'])
    returnList.append(eachLightTransitionQuality['TransitionShift'])
    returnList.append(eachLightTransitionQuality['TransitionFWHM'])
    returnList.append(eachLightTransitionQuality['TransitionFWHM2base'])
    returnList.append(eachLightTransitionQuality['TransitionMaxIntensity'])
    returnList.append(eachLightTransitionQuality['TransitionMaxBoundaryIntensity'])
    returnList.append(eachLightTransitionQuality['TransitionMaxBoundaryIntensityNormalized'])

    eachHeavyTransitionQuality = transition_level['heavy'][transition]
    returnList.append(eachHeavyTransitionQuality['TransitionJaggedness'])
    returnList.append(eachHeavyTransitionQuality['TransitionSymmetry'])
    returnList.append(eachHeavyTransitionQuality['TransitionModality'])
    returnList.append(eachHeavyTransitionQuality['TransitionShift'])
    returnList.append(eachHeavyTransitionQuality['TransitionFWHM'])
    returnList.append(eachHeavyTransitionQuality['TransitionFWHM2base'])
    returnList.append(eachHeavyTransitionQuality['TransitionMaxIntensity'])
    returnList.append(eachHeavyTransitionQuality['TransitionMaxBoundaryIntensity'])
    returnList.append(eachHeavyTransitionQuality['TransitionMaxBoundaryIntensityNormalized'])
    
    #summation of all transitions of heavy and light, respectively
    eachTransitionPairQuality = transition_pair_level['sum']
    returnList.append(eachTransitionPairQuality['PairSimilarity'])
    returnList.append(eachTransitionPairQuality['PairShift'])
    returnList.append(eachTransitionPairQuality['PairFWHMConsistency'])

    eachLightTransitionQuality = transition_level['light']['sum']
    returnList.append(eachLightTransitionQuality['TransitionJaggedness'])
    returnList.append(eachLightTransitionQuality['TransitionSymmetry'])
    returnList.append(eachLightTransitionQuality['TransitionModality'])
    returnList.append(eachLightTransitionQuality['TransitionShift'])
    returnList.append(eachLightTransitionQuality['TransitionFWHM'])
    returnList.append(eachLightTransitionQuality['TransitionFWHM2base'])
    returnList.append(eachLightTransitionQuality['TransitionMaxIntensity'])
    returnList.append(eachLightTransitionQuality['TransitionMaxBoundaryIntensity'])
    returnList.append(eachLightTransitionQuality['TransitionMaxBoundaryIntensityNormalized'])
    
    eachHeavyTransitionQuality = transition_level['heavy']['sum']
    returnList.append(eachHeavyTransitionQuality['TransitionJaggedness'])
    returnList.append(eachHeavyTransitionQuality['TransitionSymmetry'])
    returnList.append(eachHeavyTransitionQuality['TransitionModality'])
    returnList.append(eachHeavyTransitionQuality['TransitionShift'])
    returnList.append(eachHeavyTransitionQuality['TransitionFWHM'])
    returnList.append(eachHeavyTransitionQuality['TransitionFWHM2base'])
    returnList.append(eachHeavyTransitionQuality['TransitionMaxIntensity'])
    returnList.append(eachHeavyTransitionQuality['TransitionMaxBoundaryIntensity'])
    returnList.append(eachHeavyTransitionQuality['TransitionMaxBoundaryIntensityNormalized'])

    #acrossAll
    returnList.append(peak_level['MeanIsotopeRTConsistency_acrossAll'])
    eachLightTransitionQuality = transition_level['light'][transition]
    returnList.append(eachLightTransitionQuality['Area2SumRatioCV_acrossAll'])
    returnList.append(eachLightTransitionQuality['MeanIsotopeRatioConsistency_acrossAll'])
    returnList.append(eachLightTransitionQuality['MeanIsotopeFWHMConsistency_acrossAll'])
    
    eachHeavyTransitionQuality = transition_level['heavy'][transition]
    returnList.append(eachHeavyTransitionQuality['Area2SumRatioCV_acrossAll'])
    returnList.append(eachHeavyTransitionQuality['MeanIsotopeRatioConsistency_acrossAll'])
    returnList.append(eachHeavyTransitionQuality['MeanIsotopeFWHMConsistency_acrossAll'])

    return returnList
  
  def serialize(self, fn, ps):
    chrom_data = self.selectChromData(fn, ps)
    heavyTransitionIonOrder = chrom_data['heavyTransitionIonOrder']
    lightTransitionIonOrder = chrom_data['lightTransitionIonOrder']
    peak_level = self.quality_indexes[fn][ps]['peak_group_level']
    isotope_light_level = self.quality_indexes[fn][ps]['isotope_level']['light']
    isotope_heavy_level = self.quality_indexes[fn][ps]['isotope_level']['heavy']
    transition_pair_level = self.quality_indexes[fn][ps]['transition_pair_level']
    transition_level = self.quality_indexes[fn][ps]['transition_level']
    returnList = [
      peak_level['PeakGroupRatioCorr'],
      peak_level['PeakGroupJaggedness_mean'],
      peak_level['PeakGroupSymmetry_mean'],
      peak_level['PeakGroupSimilarity_mean'],
      peak_level['PeakGroupModality_mean'],
      peak_level['PeakGroupShift_mean'],
      peak_level['PeakGroupFWHM_mean'],
      peak_level['PeakGroupFWHM2base_mean'],
      peak_level['PeakGroupJaggedness_sum'],
      peak_level['PeakGroupSymmetry_sum'],
      peak_level['PeakGroupSimilarity_sum'],
      peak_level['PeakGroupModality_sum'],
      peak_level['PeakGroupShift_sum'],
      peak_level['PeakGroupFWHM_sum'],
      peak_level['PeakGroupFWHM2base_sum'],
      isotope_light_level['IsotopeJaggedness_mean'],
      isotope_light_level['IsotopeSymmetry_mean'],
      isotope_light_level['IsotopeSimilarity_mean'],
      isotope_light_level['IsotopeModality_mean'],
      isotope_light_level['IsotopeShift_mean'],
      isotope_light_level['IsotopeFWHM_mean'],
      isotope_light_level['IsotopeFWHM2base_mean'],
      isotope_light_level['IsotopeJaggedness_sum'],
      isotope_light_level['IsotopeSymmetry_sum'],
      isotope_light_level['IsotopeSimilarity_sum'],
      isotope_light_level['IsotopeModality_sum'],
      isotope_light_level['IsotopeShift_sum'],
      isotope_light_level['IsotopeFWHM_sum'],
      isotope_light_level['IsotopeFWHM2base_sum'],
      isotope_heavy_level['IsotopeJaggedness_mean'],
      isotope_heavy_level['IsotopeSymmetry_mean'],
      isotope_heavy_level['IsotopeSimilarity_mean'],
      isotope_heavy_level['IsotopeModality_mean'],
      isotope_heavy_level['IsotopeShift_mean'],
      isotope_heavy_level['IsotopeFWHM_mean'],
      isotope_heavy_level['IsotopeFWHM2base_mean'],
      isotope_heavy_level['IsotopeJaggedness_sum'],
      isotope_heavy_level['IsotopeSymmetry_sum'],
      isotope_heavy_level['IsotopeSimilarity_sum'],
      isotope_heavy_level['IsotopeModality_sum'],
      isotope_heavy_level['IsotopeShift_sum'],
      isotope_heavy_level['IsotopeFWHM_sum'],
      isotope_heavy_level['IsotopeFWHM2base_sum']
    ]
    for index in range(self.nFragmentIon):
      if index < len(heavyTransitionIonOrder):
        transition = heavyTransitionIonOrder[index]
        eachTransitionPairQuality = transition_pair_level[transition]
        #returnList.append(eachTransitionPairQuality['PairRatioConsistency'])
        returnList.append(eachTransitionPairQuality['PairSimilarity'])
        returnList.append(eachTransitionPairQuality['PairShift'])
        returnList.append(eachTransitionPairQuality['PairFWHMConsistency'])
      else: # padding (the padding value needs to be determined)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)

    for index in range(self.nFragmentIon):
      if index < len(heavyTransitionIonOrder):
        transition = heavyTransitionIonOrder[index]
        eachLightTransitionQuality = transition_level['light'][transition]
        returnList.append(eachLightTransitionQuality['TransitionJaggedness'])
        returnList.append(eachLightTransitionQuality['TransitionSymmetry'])
        returnList.append(eachLightTransitionQuality['TransitionModality'])
        returnList.append(eachLightTransitionQuality['TransitionShift'])
        returnList.append(eachLightTransitionQuality['TransitionFWHM'])
        returnList.append(eachLightTransitionQuality['TransitionFWHM2base'])
        returnList.append(eachLightTransitionQuality['TransitionMaxIntensity'])
        returnList.append(eachLightTransitionQuality['TransitionMaxBoundaryIntensity'])
        returnList.append(eachLightTransitionQuality['TransitionMaxBoundaryIntensityNormalized'])
      else:
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)

    for index in range(self.nFragmentIon):
      if index < len(heavyTransitionIonOrder):
        transition = heavyTransitionIonOrder[index]
        eachHeavyTransitionQuality = transition_level['heavy'][transition]
        returnList.append(eachHeavyTransitionQuality['TransitionJaggedness'])
        returnList.append(eachHeavyTransitionQuality['TransitionSymmetry'])
        returnList.append(eachHeavyTransitionQuality['TransitionModality'])
        returnList.append(eachHeavyTransitionQuality['TransitionShift'])
        returnList.append(eachHeavyTransitionQuality['TransitionFWHM'])
        returnList.append(eachHeavyTransitionQuality['TransitionFWHM2base'])
        returnList.append(eachHeavyTransitionQuality['TransitionMaxIntensity'])
        returnList.append(eachHeavyTransitionQuality['TransitionMaxBoundaryIntensity'])
        returnList.append(eachHeavyTransitionQuality['TransitionMaxBoundaryIntensityNormalized'])
      else:
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
    # summation of all transitions (heavy and light, respectively)
    sumTransitionPairQuality = transition_pair_level['sum']
    returnList.append(sumTransitionPairQuality['PairRatioConsistency'])
    returnList.append(sumTransitionPairQuality['PairSimilarity'])
    returnList.append(sumTransitionPairQuality['PairShift'])
    returnList.append(sumTransitionPairQuality['PairFWHMConsistency'])
    sumLightTransitionQuality = transition_level['light']['sum']
    returnList.append(sumLightTransitionQuality['TransitionJaggedness'])
    returnList.append(sumLightTransitionQuality['TransitionSymmetry'])
    returnList.append(sumLightTransitionQuality['TransitionModality'])
    returnList.append(sumLightTransitionQuality['TransitionShift'])
    returnList.append(sumLightTransitionQuality['TransitionFWHM'])
    returnList.append(sumLightTransitionQuality['TransitionFWHM2base'])
    returnList.append(sumLightTransitionQuality['TransitionMaxIntensity'])
    returnList.append(sumLightTransitionQuality['TransitionMaxBoundaryIntensity'])
    returnList.append(sumLightTransitionQuality['TransitionMaxBoundaryIntensityNormalized'])
    sumHeavyTransitionQuality = transition_level['heavy']['sum']
    returnList.append(sumHeavyTransitionQuality['TransitionJaggedness'])
    returnList.append(sumHeavyTransitionQuality['TransitionSymmetry'])
    returnList.append(sumHeavyTransitionQuality['TransitionModality'])
    returnList.append(sumHeavyTransitionQuality['TransitionShift'])
    returnList.append(sumHeavyTransitionQuality['TransitionFWHM'])
    returnList.append(sumHeavyTransitionQuality['TransitionFWHM2base'])
    returnList.append(sumHeavyTransitionQuality['TransitionMaxIntensity'])
    returnList.append(sumHeavyTransitionQuality['TransitionMaxBoundaryIntensity'])
    returnList.append(sumHeavyTransitionQuality['TransitionMaxBoundaryIntensityNormalized'])
    # acrossAll
    returnList.append(peak_level['MeanIsotopeRTConsistency_acrossAll'])
    for index in range(self.nFragmentIon):
      if index < len(heavyTransitionIonOrder): 
        transition = heavyTransitionIonOrder[index]
        eachLightTransitionQuality = transition_level['light'][transition]
        returnList.append(eachLightTransitionQuality['Area2SumRatioCV_acrossAll'])
        returnList.append(eachLightTransitionQuality['MeanIsotopeRatioConsistency_acrossAll'])
        returnList.append(eachLightTransitionQuality['MeanIsotopeFWHMConsistency_acrossAll'])
      else:
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)

    for index in range(self.nFragmentIon):
      if index < len(heavyTransitionIonOrder):
        transition = heavyTransitionIonOrder[index]
        eachHeavyTransitionQuality = transition_level['heavy'][transition]
        returnList.append(eachHeavyTransitionQuality['Area2SumRatioCV_acrossAll'])
        returnList.append(eachHeavyTransitionQuality['MeanIsotopeRatioConsistency_acrossAll'])
        returnList.append(eachHeavyTransitionQuality['MeanIsotopeFWHMConsistency_acrossAll'])
      else:
        returnList.append(-1)
        returnList.append(-1)
        returnList.append(-1)
    return returnList

  def exportAllQualityIndexes(self, outputExcelPath):
    start_time = time.time()
    qi = self.calculateQualityIndexes()
    end_time = time.time()
    print('Export all quality indexes in: %.2f seconds' % (end_time - start_time))
    print('Write to xlsx file...')
    xlsx_writer = pd.ExcelWriter(outputExcelPath, engine='xlsxwriter')
    qi['peak_group_level'].to_excel(xlsx_writer, sheet_name='PeakGroup', index=False)
    qi['isotope_level'].to_excel(xlsx_writer, sheet_name='Isotope', index=False)
    qi['transition_pair_level'].to_excel(xlsx_writer, sheet_name='TransitionPair', index=False)
    qi['transition_level'].to_excel(xlsx_writer, sheet_name='Transition', index=False)
    xlsx_writer.save()
    print('finished exporting to an Excel file %.2f seconds' % (time.time() - end_time))
    print('Total time: %.2f seconds' % (time.time() - start_time))
import doctest
import numpy as np
import obspy
from scipy.signal import hilbert,correlate
import scipy.fft as sy_fft
from seispy import mccc

from obspy.core.stream import Stream

def prepare_sac2mat(urls,
                     freq_bands=(0.1,1),
                     time_window=[-20,20],
                     time_point="t1"):
    """
    prepare a matix from SAC file
    Input files should be in the same or you should trans them to matrix by yourself
    npts,sampling rate
    
    detrend and taper are set to SAC's default,
    Parameters
    ----------
    urls :
        things can convey to obspy.read
    freq_bands:
        any_thing contains 2 float at[0]and [1]
    time_point : str, optional
    time_window:
        relative to time_point
    Returns
    -------
    temps:
        matrix containing signals
    norm_base:
        array containg normalization factors
    """
    try:
        signals=obspy.read(urls)
    except:
        raise(IOError,"cannot read files")
    signals.detrend();signals.taper(max_percentage=0.05,type="hann")
    if(len(freq_bands)>3 or freq_bands[0]>=freq_bands[1]):
        raise(ValueError,"freq setting fault")
    signals.filter("bandpass",
                   freqmin=freq_bands[0],
                   freqmax=freq_bands[1])
    lines=len(signals)
    sampling_rate=signals[0].stats.sampling_rate
    space=time_window[1]-time_window[0]
    
    temps=np.zeros(lines,space)
    norm_base=np.empty(len(signals))
    for i in range(0,lines):
        begin=(signals[i].stats.sac[time_point]\
            -signals[i].stats.sac['b'])*sampling_rate
        end=begin+space*sampling_rate
        if begin<0 or end>signals[i].npts:
            raise(ValueError,"Error when cutting time_window")
        temps[i]=signals[i].data[begin,end+1]
        norm_base[i]=np.sqrt(np.sum(np.square(temps[i])))
    return sampling_rate,temps,norm_base

def prepare_MMCC(mat_sig,method="fft"):
    """
    if correlate_method chose cross_correlation,
    Parameters
    ----------
    mat_sig : 
        signals stored in matrix
    method='fft' or 'iter'
    """
    line,size=np.shape(mat_sig)
    phase=np.unwrap(np.angle(hilbert(mat_sig)))
    if method=='fft':
        cig=sy_fft.fft(np.multiply(np.cos(phase),mat_sig),size*2)
        sig=sy_fft.fft(np.multiply(np.sin(phase),mat_sig),size*2)
        hig=sy_fft.fft(mat_sig,size*2)
    else:
        cig=np.multiply(np.cos(phase),mat_sig)
        sig=np.multiply(np.sin(phase),mat_sig)
        hig=mat_sig
    return cig,sig,hig


def MCC_fft(size,f1,f2,pair1,pair2):
    """
    do MCC for fft pairs

    Parameters
    ----------

    Returns
    -------
    float arrays
    """
    return sy_fft.ifft(
           f1[pair1]*f2[pair2].conj()
        ,size*2).real

def MMCC_fft(size,
        cig1,sig1,hig1,
        cig2,sig2,hig2,
        pair1,pair2):
    """
    do MMCC from fft format
    all cig should be in same size
    size is length of a signal line
    """
    return sy_fft.ifft(cig1[pair1]*cig2[pair2].conj()+\
        sig1[pair1]*sig2[pair2].conj()+\
        hig1[pair1]*hig2[pair2].conj(),\
        size*2).real/2

def MMCC_iter(
        cig1,sig1,hig1,
        cig2,sig2,hig2,
        pair1,pair2
):
    """
    do MMCC from iter format
    all cig should be in same size
    """
    return correlate(cig1[pair1],cig2[pair2],'full','direct')\
            +correlate(sig1[pair1],sig2[pair2],'full','direct')\
            +correlate(hig1[pair1],hig2[pair2],'full','direct')

def MMCC_cal_full(cig,sig,hig,method='fft'):
    """

    Parameters
    ----------
    cig,sig,hig:
        3 parts of signals
    method:
        method for correlate
    Returns
        rel_time:
            array containing reltime in points not time(s)
    -------
    """
    line,npts=np.shape(cig)
    npts=npts*2 if method=='fft' else npts
    tt=np.zeros((line,line))
    for i in range(line-1):
        for j in range(i+1,line):
            if method=='fft':
                tt[i, j] = np.argmax(
                    MMCC_fft(npts, cig, sig, hig,
                             cig, sig, hig,
                             i, j))
            else:
                tt[i,j]=np.argmax(
                    MMCC_iter(cig,sig,hig,
                        cig,sig,hig,
                        i,j))
    if method=='fft':
        (row,col)=np.where(tt>npts/2)
        tt[i,j]-=(npts+1)
    else:
        tt-=(npts+1)

    for i in range(line):
        rel_time=(-np.sum(tt[0:i+1,i])
                  +np.sum(tt[i,i+1:line]))/line
    return rel_time

def MMCC_cal_signal(cig,sig,hig,num,method='fft'):
    """
    Parameters
    ----------
    cig
    cig,sig,hig:
        3 parts of signals
    method:
        method for correlate
    Returns
        rel_time:
            array containing reltime in points not time(s)
    Returns
    -------
    rel_time:
        array containing reltime in points not time(s)
    """
    line,npts=np.shape(cig)
    npts=npts*2 if method=='fft' else npts
    rel_time=np.zeros(line)
    for i in range(line):
        if i == num:
            continue
        if method=='fft':
            rel_time[i]=np.argmax(
                MMCC_fft(npts,cig,sig,hig,
                         cig,sig,hig,
                         i,num)
            )
        else:
            rel_time[i]=np.argmax(
                MMCC_iter(cig,sig,hig,
                          cig,sig,hig,
                          i,num)
            )
        return rel_time




def MMCC_interface(url,method='fft',rel_method='full',*args):
    """

    Parameters
    ----------
    mat_signal:
        signals containing in mat
    method:
        method for correlate, fft or iter
    rel_method:
        'full' for correlate between all line
        else for correlate to one user chosed
    args:
        use when rel_method is not full
    Returns
    -------
    """
    sampling_rate,sig_mat,norm_base=prepare_sac2mat(url)
    cig,sig,hig=prepare_MMCC(sig_mat,method)
    if rel_method=='full':
        rel_time=MMCC_cal_full(cig,sig,hig,method)
    else:
        rel_time=MMCC_cal_signal(cig,sig,hig,args[0],method)
    return rel_time/sampling_rate




if __name__=='__main__':
    MMCC_interface()

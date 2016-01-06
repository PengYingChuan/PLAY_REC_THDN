from __future__ import division
import sys
#from scikits.audiolab import Sndfile
import pysoundfile_YingChuan as sf
import sounddevice as sd
from scipy.signal import blackmanharris, lfilter, filtfilt, freqz, butter, buttord, cheby2, cheb2ord
from numpy.fft import rfft, irfft
from numpy import argmax, sqrt, mean, absolute, arange, log10
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import FileDialog
sys.path.append('your_lib_path')

# method 1: absSum
def _calVolume(waveData, frameSize, overLap):
    wlen = len(waveData)
    step = frameSize - overLap
    frameNum = int(math.ceil(wlen*1.0/step))
    volume = np.zeros((frameNum,1))
    for i in range(frameNum):
        curFrame = waveData[np.arange(i*step,min(i*step+frameSize,wlen))]
        curFrame = curFrame - np.median(curFrame) # zero-justified
        volume[i] = np.sum(np.abs(curFrame))

    return volume

def _RMS_flat(a):
    """Return the root mean square of all the elements of *a*, flattened out.

    """
    return sqrt(mean(absolute(a)**2))

def _Find_range(f, x):
    """Find range between nearest local minima from peak at index x

    """
    for i in arange(x+1, len(f)):
        if f[i+1] >= f[i]:
            uppermin = i
            break
    for i in arange(x-1, 0, -1):
        if f[i] <= f[i-1]:
            lowermin = i + 1
            break
    print lowermin-20, uppermin+20
    return (lowermin, uppermin)

def _BLACKMANHARRIS_Window(signal):
    """ blackman-harris window process
        Use Windows to smooth transient data in the begin/end
    """
    windowed = signal * blackmanharris(len(signal))  # TODO Kaiser?

    return windowed

def _Removed_DC(signal):
    """Get rid of DC and window the signal

    """
    signal -= mean(signal) # TODO: Do this in the frequency domain, and take any skirts with it
    return signal

def _HP_22Hz(signal):

    B, A = butter(4, 20, btype='high')

    filtered = lfilter(B, A, signal)
    #print filtered

    return filtered, B, A

def _Truncate(signal):
    pass

def _HP_Butterworth(interval, sampling_rate, cutoff):

    nyq = sampling_rate * 0.5
    stopfreq = float(cutoff)
    cornerfreq = 0.5 * stopfreq
    Ws = cornerfreq/nyq
    Wp = stopfreq/nyq
    #N, Wn = buttord(Wp, Ws, 3, 20)   # (?)
    #N, Wn = buttord(Wp, Ws, 6, 20)   # (?)
    N, Wn = buttord(Wp, Ws, 10, 27)   # (?)
    print "The oder of HPF is: %f" %N
    b, a = butter(N, Wn, btype='high')   # should 'high' be here for bandpass?
    #b, a = butter(9, float(20.0/nyq) , btype='high')   # should 'high' be here for bandpass?
    sf = lfilter(b, a, interval)
    return sf, b, a

def _LP_Butterworth(interval, sampling_rate, cutoff, order=5):

    nyq = sampling_rate * 0.5

    stopfreq = float(cutoff)
    cornerfreq = 0.5 * stopfreq
    Ws = stopfreq/nyq
    Wp = cornerfreq/nyq
    N, Wn = buttord(Wp, Ws, 3, 20)   # (?)
    print "The oder of LPF is: %f" %N

    """
    Wp = 2 * np.pi * 100
    Ws = 2 * np.pi * 20
    Rp = 1.5
    Rs = 20
    N, Wn = buttord(Wp, Ws, Rp, Rs)   # (?)
    """

    # for hardcoded order:
    # N = order

    b, a = butter(N, Wn, btype='low')   # should 'high' be here for bandpass?
    #b, a = butter(9, float(20.0/nyq) , btype='high')   # should 'high' be here for bandpass?
    sf = lfilter(b, a, interval)
    return sf, b, a

def _NOTCH_butter_bandstop(interval, sampling_rate, cutoff, order=5):

    nyq = sampling_rate * 0.5

    stopfreq = float(cutoff) + 300.0
    cornerfreq = float(cutoff) - 300.0
    Ws = stopfreq/nyq
    Wp = cornerfreq/nyq
    N, Wn = buttord(Wp, Ws, 8, 20)   # (?)
    print N, Wn, Wp, Ws, stopfreq, cornerfreq

    b, a = butter(N, [Wp, Ws], btype='bandstop')   # should 'high' be here for bandpass?
    sf = lfilter(b, a, interval)
    return sf, b, a

def _NOTCH_Cheby2_bandstop(interval, sampling_rate, cutoff, order=5):

    nyq = sampling_rate * 0.5

    stopfreq = float(cutoff) + 300.0
    cornerfreq = float(cutoff) - 300.0
    Ws = stopfreq/nyq
    Wp = cornerfreq/nyq
    #N, Wn = cheb2ord(Wp, Ws, 3, 60)   # (?)

    N, Wn = cheb2ord([0.3, 0.5], [0.45, 0.48], 3, 160)   # (?)
    #print N, Wn, Wp, Ws, stopfreq, cornerfreq

    b, a = cheby2(N, 60, Wn, btype='stop')   # should 'high' be here for bandpass?
    sf = lfilter(b, a, interval)
    return sf, b, a

def _NOTCH_Butterworth(signal, sample_rate, f_cut):

    order = 3
    #bp_stop_Hz = np.array([0.91*f_cut, 1.1*f_cut] )
    #bp_stop_Hz = np.array([0.5*f_cut, 2*f_cut] )

    #f_cut_LOG1
    # 0=  np.log10(f_cut)
    #bp_stop_Hz = np.array([10**(f_cut_LOG10-0.2), 10**(f_cut_LOG10+0.2)] )
    bp_stop_Hz = np.array([1050, 12000] )

    B, A = butter(order, bp_stop_Hz/(sample_rate/2.0), 'bandstop')

    #print B, A, signal[0:2]
    filtered = lfilter(B, A, signal)
    #print filtered
    print "The oder of NOTCH is: %f" %order
    return filtered, B, A

def _NOTCH_NORMOLIZED(signal, sample_rate, f_center, order, Q):

    #Q = 0.5
    WDc = float(f_center) * (2*np.pi)
    c = 1 / float(np.tan(WDc / (sample_rate * 2)))
    b = 1/ Q

    if order == 2:
        # Second Order Filter
        B = [(c**2 + 1), (-2*(c**2 -1)), (c**2 + 1)]
        A = [(b * c + c**2 +1), (-2*(c**2 -1)), (-b*c + c**2 + 1)]

    elif order ==4:
        # Forth Order Filter
        n0 = (c**2 + 1)**2
        n1 = -4*(c**4 - 1)
        n2 = 2*(3*c**4 - 2*c**2 + 3)
        n3 = -4*(c**4 -1)
        n4 = (c**2 + 1)**2

        d0 = b**2 * c**2 + 2**0.5 * b * c * (c**2 + 1) + (c**2 + 1)**2
        d1 = -2 * 2**0.5 * (c**2 - 1) * (b * c + 2**0.5 * c**2 + 2**0.5)
        d2 = -2 * (b**2 * c**2 - 3 * c**4 + 2 * c**2 - 3)
        d3 = 2 * 2**0.5 * (c**2 - 1) * (b * c - 2**0.5 * (c**2 + 1))
        d4 = b**2 * c**2 - 2**0.5 * b * c * (c**2 + 1) + (c**2 + 1)**2

        B = [n0, n1, n2, n3, n4]
        A = [d0, d1, d2, d3, d4]

    elif order == 6:
        # Sixth Order Filter
        n0 = (c**2 + 1)**3
        n1 = -6 * (c**2 - 1) * (c**2 + 1)**2
        n2 = 3 * (c**2 + 1) * (5 * c**4 - 6 * c**2 + 5)
        n3 = -4 * (c**2 -1) * (5 * c**4 + 2 * c**2 + 5)
        n4 = 3 * (c**2 + 1) * (5 * c**4 - 6 * c**2 + 5)
        n5 = -6 * (c**2 - 1) * (c**2 + 1)**2
        n6 = (c**2 + 1)**3

        d0 = b**3 * c**3 + 2 * b**2 * c**2 * (c**2 + 1) + 2 * b * c *(c**2 + 1)**2 + (c**2 + 1)**3
        d1 = -2 * (c**2 - 1) * (2 * b**2 * c**2 + 4 * b * c * (c**2 + 1) + 3 * (c**2 + 1)**2)
        d2 = -(3 * b**3 * c**3 + 2 * b**2 * c**2 * (c**2 + 1) - 2 * b * c * (5 * c**4 - 6 * c**2 + 5) - 3 * (c**2 + 1)*(5 * c**4 - 6 * c**2 + 5))
        d3 = 4 * (c**2 - 1 ) * (-5 * c**4 + 2 * b**2 * c**2 - 2 * c**2 - 5)
        d4 = (3 * b**3 * c**3 - 2 * b**2 * c**2 * (c**2 + 1) - 2 *b *c * (5 * c**4 - 6 * c**2 + 5 ) + 3 * (c**2 + 1) * (5 * c**4 - 6 * c**2 + 5))
        d5 = -2 * (c**2 - 1) * (2 * b**2 * c**2 - 4 * b * c * (c**2 + 1) + 3 * (c**2 + 1)**2)
        d6 = -b**3 * c**3 + 2 * b**2 * c**2 * (c**2 + 1) - 2 * b * c * (c**2 + 1)**2 + (c**2 + 1)**3

        B = [n0, n1, n2, n3, n4, n5, n6]
        A = [d0, d1, d2, d3, d4, d5, d6]

    #print B, A, signal[0:2]
    filtered = lfilter(B, A, signal)
    #print filtered
    #print "The oder of NOTCH is: %f" %order
    return filtered, B, A

def _BANDPASS_NORMOLIZED(signal, sample_rate, f_center, order, Q):

    #Q = 0.5
    WDc = float(f_center) * (2*np.pi)
    c = 1 / float(np.tan(WDc / (sample_rate * 2)))
    b = 1/ Q

    if order == 2:
        # Second Order Filter
        n0 = b * c
        n1 = 0
        n2 = -b * c

        d0 = b * c + c**2 +1
        d1 = -2 * (c**2 - 1)
        d2 = -b * c + c**2 + 1

        B = [n0, n1, n2]
        A = [d0, d1, d2]

    elif order ==4:
        # Forth Order Filter
        n0 = b**2 * c**2
        n1 = 0
        n2 = -2 * b**2 * c**2
        n3 = 0
        n4 = b**2 * c**2

        d0 = b**2 * c**2 + 2**0.5 * b * c * (c**2 + 1) + (c**2 + 1)**2
        d1 = -2 * 2**0.5 * (c**2 - 1) * (b * c + 2**0.5 * c**2 + 2**0.5)
        d2 = -2 * (b**2 * c**2 - 3 * c**4 + 2 * c**2 - 3)
        d3 = 2 * 2**0.5 * (c**2 - 1) * (b * c - 2**0.5 * (c**2 + 1))
        d4 = b**2 * c**2 - 2**0.5 * b * c * (c**2 + 1) + (c**2 + 1)**2

        B = [n0, n1, n2, n3, n4]
        A = [d0, d1, d2, d3, d4]

    elif order == 6:
        # Sixth Order Filter
        n0 = b**3 * c**3
        n1 = 0
        n2 = -3 * b**3 * c**3
        n3 = 0
        n4 = 3 * b**3 * c**3
        n5 = 0
        n6 = -b**3 * c**3

        d0 = b**3 * c**3 + 2 * b**2 * c**2 * (c**2 + 1) + 2 * b * c * (c**2 + 1)**2 + (c**2 + 1)**3
        d1 = -2 * (c**2 - 1) * (2 * b**2 * c**2 + 4 * b * c * (c**2 + 1) + 3 * (c**2 + 1)**2)
        d2 = -(3 * b**3 * c**3 + 2 * b**2 * c**2 * (c**2 + 1) - 2 * b * c * (5 * c**4 - 6 * c**2 + 5) - 3 * (c**2 + 1) * (5 * c**4 - 6 * c**2 + 5))
        d3 = 4 * (c**2 - 1) * (-5 * c**4 + 2 * b**2 * c**2 - 2 * c**2 - 5)
        d4 = (3 * b**3 * c**3 - 2 * b**2 * c**2 * (c**2 + 1) - 2 * b * c * (5 * c**4 - 6 * c**2 + 5) + 3 * (c**2 + 1) * (5 * c**4 - 6 * c**2 + 5))
        d5 = -2 * (c**2 -1) * (2 * b**2 * c**2 - 4 *b * c * (c**2 + 1) + 3 * (c**2 + 1)**2)
        d6 = - b**3 * c**3 + 2 * b**2 * c**2 * (c**2 + 1) - 2 * b * c * (c**2 + 1)**2 + (c**2 + 1)**3

        B = [n0, n1, n2, n3, n4, n5, n6]
        A = [d0, d1, d2, d3, d4, d5, d6]

    #print B, A, signal[0:2]
    filtered = lfilter(B, A, signal)
    #print filtered
    #print "The oder of NOTCH is: %f" %order
    return filtered, B, A

def _NOTCH_Second_Equalizer(signal, sample_rate, f_cut, Q, Gain):
    """

        """
    A = 10.0**(Gain/40.0)
    w0 = 2.0*math.pi*f_cut/sample_rate
    alpha  = math.sin(w0) /2.0 /float(Q)
    b0 = 1.0 + alpha * A
    b1 = -2.0 * math.cos(w0)
    b2 = 1.0 - alpha *A
    a0 = 1.0 + alpha / A
    a1 = -2.0*math.cos(w0)
    a2 = 1 - alpha / A

    AA = [1.0, a1/a0, a2/a0]
    BB = [b0/a0, b1/a0, b2/a0]

    filtered = lfilter(BB, AA, signal)

    return filtered, BB, AA

def _Notch_Filter(signal, sample_rate, f_cut):
    #print f_cut

    wn = float(f_cut)/sample_rate

    r = 0.99
    B, A = np.zeros(3), np.zeros(3)
    A[0],A[1],A[2] = 1.0, -2.0*r*np.cos(2*np.pi*wn), r*r
    B[0],B[1],B[2] = 1.0, -2.0*np.cos(2*np.pi*wn), 1.0

    filtered = lfilter(B, A, signal, -1)

    return filtered, B, A

def _FFT(windowed, sample_rate):
    _FFT = rfft(windowed)
    _FFT_ABS = np.absolute(_FFT)
    _FFT_MAX_idex = argmax(_FFT_ABS)
    _FFT_MAX_MAG = np.nanmax(_FFT_ABS)
    FFT_FREQ = sample_rate/len(windowed) * np.arange(len(_FFT))
    FFT_MAG_LOG = 20*log10(_FFT_ABS/_FFT_MAX_MAG)
    Fc = sample_rate * (_FFT_MAX_idex / len(windowed))

    return FFT_FREQ, FFT_MAG_LOG, Fc, _FFT_MAX_idex, _FFT_ABS, _FFT_MAX_MAG

def THDN(signal, sample_rate, show_plot):
    """
                Calculate THDN
        """
    """ To removed DC offset"""
    signal = _Removed_DC(signal) # TODO: Do this in the frequency domain, and take any skirts with it

    """ To reduce first and final signal """
    windowed = _BLACKMANHARRIS_Window(signal)

    """ HPF 20Hz"""
    HP_Filtered, B_HPF, A_HPF = _HP_Butterworth(windowed, sample_rate, cutoff = 20)

    """ LPF 22KHz """
    LP_Filtered, B_LPF, A_LPF = _LP_Butterworth(HP_Filtered, sample_rate, cutoff = 22000)
    FFT_FREQ, FFT_MAG_LOG, Fc, _FFT_MAX_idex, _FFT_ABS, _FFT_MAX_MAG = _FFT(LP_Filtered, sample_rate)

    """ Band Pass Filter to get max energy of signal """
    BP_Filtered, B_BPF, A_BPF = _BANDPASS_NORMOLIZED(LP_Filtered, sample_rate, f_center = Fc, order = 6, Q = 200 )

    """ Notch Filter """
    NOTCH_Filtered, B_NOTCH, A_NOTCH = _NOTCH_NORMOLIZED(LP_Filtered, sample_rate, f_center = Fc, order = 4, Q = 0.8)

    """ Calculate FFT after notch filter """
    [NOTCH_FFT_FREQ, NOTCH_FFT_MAG_LOG, NOTCH_Fc, _NOTCH_FFT_MAX_idex,
     NOTCH_FFT_ABS, NOTCH_FFT_MAX_MAG] = _FFT(NOTCH_Filtered, sample_rate)

    total_rms = _RMS_flat(signal)
    print "Average RMS: %.2f dBFS" %(20*np.log10(total_rms / 0.70721))
    print 'Frequency: %f Hz' % Fc
    THDN_Value = _RMS_flat(NOTCH_Filtered) / _RMS_flat(BP_Filtered)
    print "THD+N:     %.7f%% or %.1f dB" % (THDN_Value * 100, 20 * log10(THDN_Value))

    if show_plot == True:

        """ Calculate FFT after windows """
        [unFiltered_FFT_FREQ, unFiltered_FFT_MAG_LOG, unFiltered_Fc,
        unFiltered_FFT_MAX_idex, unFiltered_FFT_ABS, unFiltered_FFT_MAX_MAG] = _FFT(windowed, sample_rate)

        """
         Filters Frequency Response
        """
        N = len(HP_Filtered)
        # High Pass Filter response
        w_HPF, HPF = freqz(B_HPF, A_HPF, N)
        HPF_ABS = np.absolute(HPF)
        # Low Pass Filter response
        w_LPF, LPF = freqz(B_LPF, A_LPF, N)
        LPF_ABS = np.absolute(LPF)
        # Band Pass Filter response
        w_BPF, BPF = freqz(B_BPF, A_BPF, N)
        BPF_ABS = np.absolute(BPF)
        # Notch Filter response
        w_NOTCH, NOTCH = freqz(B_NOTCH, A_NOTCH, N)
        NOTCH_ABS = np.absolute(NOTCH)



        # Plot FFT
        fig = plt.figure()
        ax1 = fig.add_subplot(241)
        ax1.set_title("Audio RAW Data", size =10)
        ax1.set_ylabel('Mag')
        ax1.plot(signal)
        plt.ylim(-1, 1)
        plt.grid(True)

        ax2 = fig.add_subplot(242)
        ax2.set_title("After Blackman-Harris Window",  size =10)
        ax2.plot(windowed)
        plt.ylim(-1, 1)
        plt.grid(True)

        ax3 = fig.add_subplot(243)
        ax3.set_title("After BandPass 20~22KHz",  size =10)
        ax3.plot(LP_Filtered)
        plt.ylim(-1, 1)
        plt.grid(True)

        ax4 = fig.add_subplot(244)
        ax4.set_title("After Notch Filter",  size =10)
        ax4.plot(NOTCH_Filtered)
        plt.ylim(-1, 1)
        plt.grid(True)

        ax5 = fig.add_subplot(245)
        ax5.set_title("FFT of Audio RAW Data", size =10)
        ax5.set_xlabel('Frequency (Hz)')
        ax5.set_ylabel('Gain (dB)')
        ax5.semilogx(unFiltered_FFT_FREQ, unFiltered_FFT_MAG_LOG)
        plt.xlim(1, sample_rate/2)
        plt.ylim(-240, 10)
        plt.grid(True)

        ax6 = fig.add_subplot(246)
        ax6.set_title("FFT of After Blackman-Harris Window",  size =10)
        ax6.set_xlabel('Frequency (Hz)')
        ax6.semilogx(unFiltered_FFT_FREQ, unFiltered_FFT_MAG_LOG)
        ax6.semilogx(w_HPF*sample_rate/(2*np.pi), 20*np.log10(HPF_ABS), 'r', label='High Pass Filter')
        ax6.semilogx(w_LPF*sample_rate/(2*np.pi), 20*np.log10(LPF_ABS), 'y', label='High Pass Filter')
        plt.xlim(1, sample_rate/2)
        plt.ylim(-240, 10)
        plt.grid(True)

        ax7 = fig.add_subplot(247)
        ax7.set_title("FFT of After BandPass 20~22KHz",  size =10)
        ax7.set_xlabel('Frequency (Hz)')
        ax7.semilogx(FFT_FREQ, FFT_MAG_LOG)
        ax7.semilogx(w_BPF*sample_rate/(2*np.pi), 20*np.log10(BPF_ABS), 'y', label='Band Pass Filter')
        ax7.semilogx(w_NOTCH*sample_rate/(2*np.pi), 20*np.log10(NOTCH_ABS), 'r', label='High Pass Filter')
        plt.xlim(1, sample_rate/2)
        plt.ylim(-240, 10)
        plt.grid(True)

        ax8= fig.add_subplot(248)
        ax8.set_title("FFT of After Notch Filter",  size =10)
        ax8.set_xlabel('Frequency (Hz)')
        ax8.semilogx(NOTCH_FFT_FREQ,  20*np.log10(NOTCH_FFT_ABS/_FFT_MAX_MAG))
        plt.xlim(1, sample_rate/2)
        plt.ylim(-240, 10)
        plt.grid(True)
        plt.show()


    return THDN_Value

def LOAD_WAVE(filename):
    """
    wave_file = Sndfile(filename, 'r')
    signal = wave_file.read_frames(wave_file.nframes)
    channels = wave_file.channels
    sample_rate = wave_file.samplerate
    """
    signal, sample_rate, channels, type = sf.read(filename,'float32') # normalizes PCM to the range[-1, 1) with 64bits floating point
    return signal, sample_rate, channels

def DATA_GEN(sample_rate, frequency, duration, volume_db, channel):

    mag = 10**(volume_db / 20)
    r = frequency / sample_rate

    def Cal_Samples(column, row):
        return (mag * (np.sin(2 * np.pi * column * r)).astype(np.float32))

    if channel == 1:
        samples = mag * (np.sin(2 * np.pi * np.arange(sample_rate * duration) * r)).astype(np.float32)
    elif channel == 2:
        samples = np.fromfunction(Cal_Samples, (sample_rate * duration, 2))

    return samples

def WAV_GEN(sample_rate, frequency, duration, volume_db, channels, filename):

    mag = 10**(volume_db / 20)
    r = frequency / sample_rate

    def Cal_Samples(column, row):
        return (mag * (np.sin(2 * np.pi * column * r)).astype(np.float32))

    samples = np.fromfunction(Cal_Samples, (sample_rate * duration, channels))
    """
    if channels == 1:
        samples = mag * (np.sin(2 * np.pi * np.arange(sample_rate * duration) * r)).astype(np.float32)
    elif channels == 2:
        samples = np.fromfunction(Cal_Samples, (sample_rate * duration, 2))
    """

    sf.write(samples, filename, samplerate=sample_rate)

def PLAYREC(data, fs, channel):
    mysound = sd.playrec(data, fs, channels = 2, blocking=True) # "blocking=True" is very important.
    return mysound

#WAV_GEN(sample_rate = 44100, frequency = 1000, duration = 2, volume_db = -1, channel = 2, filename = 'test_file5.wav')

#filename = sys.argv[1]
#filename = '1k_n1db_MONO_16bits.wav'
#signal, sample_rate, channels = LOAD_WAVE(filename)
#print 'Analyzing "' + filename + '"...'
#print signal[:,0], sample_rate, channels, len(signal)


sample_rate = 44100
frequency = 1000
duration = 2
volume_db = -16
channels = 1

signal = DATA_GEN(sample_rate, frequency, duration, volume_db, channels)

print signal

if channels == 1:
    # Monaural
    THDN_Value = THDN(signal, sample_rate, show_plot = True)
elif channels == 2:
    # Stereo
    if np.array_equal(signal[:,0],signal[:,1]):
        print '--Left and Right channels are identical:--'
        THDN_Value = THDN(signal[:,0], sample_rate)
    else:
        print '--Left channel:--'
        THDN_Value = THDN(signal[:,0], sample_rate, show_plot = True)
        print '--Right channel:--'
        THDN_Value = THDN(signal[:,1], sample_rate, show_plot = True)
else:
    # Multi-channel
    for ch_no, channel in enumerate(signal.transpose()):
        print '--Channel %d:--' % (ch_no + 1)
        THDN_Value = THDN(channel, sample_rate, show_plot = True)

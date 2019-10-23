import numpy as np

# for full rupture length   
def mw2m0(mw):
    return 10**(3. * mw / 2.  + 9.05)

def mag2len_L14(mw, ftype='scrrs'): # in MW - not sure if correct
    # assume mw is a list
    logM0 = np.log10(mw2m0(mw))
    #logM0 = np.array([logM0])

    if ftype == 'scrrs':
        b = 3.0
        a = 6.382
        rl = 10**((logM0 - a) / b) / 1000.

        
        idx = rl > 2.5
        b = 2.5
        a = 8.08
        rl = 10**((logM0 - a) / b) / 1000.

    return rl # in km


r1 = mag2len_L14(7.25)
print(r1)
import numpy as np
from apsint import apsint



def read_fgong(filename):
    """
    Read in an FGONG file.
    This function can read in FONG versions 250 and 300.
    The output gives the center first, surface last.

    Parameters
    ----------
    filename : string
        Absolute path an name of the FGONG file to read

    Returns
    -------
    starg : dict
        Global parameters of the model

    starl : dict
        Local parameters of the model
    """
    # Standard definitions and units for FGONG
    glob_pars = [('mass', 'g'),
                 ('Rphot', 'cm'),
                 ('Lphot', 'erg/s'),
                 ('Zini', None),
                 ('Xini', None),
                 ('alpha', None),
                 ('phi', None),
                 ('xi', None),
                 ('beta', None),
                 ('dp', 's'),
                 ('ddP_drr_c', None),
                 ('ddrho_drr_c', None),
                 ('Age', 'yr'),
                 ('teff', 'K'),
                 ('Gconst', 'cm3/gs2')]
    loc_pars = [('radius', 'cm'),
                ('ln(m/M)', None),
                ('Temp', 'K'),
                ('P', 'kg/m/s2'),
                ('Rho', 'g/cm3'),
                ('X', None),
                ('Lumi', 'erg/s'),
                ('opacity', 'cm2/g'),
                ('eps_nuc', None),
                ('gamma1', None),
                ('grada', None),
                ('delta', None),
                ('cp', None),
                ('free_e', None),
                ('brunt_A', None),
                ('rx', None),
                ('Z', None),
                ('R-r', 'cm'),
                ('eps_logg', None),
                ('Lg', 'erg/s'),
                ('xhe3', None),
                ('xc12', None),
                ('xc13', None),
                ('xn14', None),
                ('xo16', None),
                ('dG1_drho', None),
                ('dG1_dp', None),
                ('dG1_dY', None),
                ('xh2', None),
                ('xhe4', None),
                ('xli7', None),
                ('xbe7', None),
                ('xn15', None),
                ('xo17', None),
                ('xo18', None),
                ('xne20', None),
                ('xh1', None),
                ('na38', None),
                ('na39', None),
                ('na40', None)]

    # Start reading the file
    ff = open(filename, 'r')
    lines = ff.readlines()

    # Read file definitions from the fifth line (first four is comments)
    NN, ICONST, IVAR, IVERS = [int(i) for i in lines[4].strip().split()]
    if not ICONST == 15:
        raise ValueError('cannot interpret FGONG file: wrong ICONST')

    # Data storage
    data = []
    starg = {}

    # Read the file from the fifth line onwards
    # Change in the format for storing the numbers (February 2017):
    #  - If IVERS <= 1000, 1p5e16.9
    #  - If IVERS  > 1000, 1p,5(x,e26.18e3)
    if IVERS <= 1000:
        for line in lines[5:]:
            data.append([line[0 * 16:1 * 16], line[1 * 16:2 * 16],
                         line[2 * 16:3 * 16], line[3 * 16:4 * 16],
                         line[4 * 16:5 * 16]])
    else:
        for line in lines[5:]:
            data.append([line[0 * 27:1 * 27], line[1 * 27:2 * 27],
                         line[2 * 27:3 * 27], line[3 * 27:4 * 27],
                         line[4 * 27:5 * 27]])

    # Put the data into arrays
    data = np.ravel(np.array(data, float))
    for i in range(ICONST):
        starg[glob_pars[i][0]] = data[i]
    data = data[15:].reshape((NN, IVAR)).T

    # Reverse the profile to get center ---> surface
    data = data[:, ::-1]

    # Make it into a record array and return the data
    starl = np.rec.fromarrays(data, names=[lp[0] for lp in loc_pars])

    # Exclude the center r = 0. mesh (MESA includes it)
    if starl['radius'][0] < 1.e-14:
        starl = starl[1:]

    return starg, starl



# Main script to compute period spacing
#--------------------------------------
Gconst = 6.672320e-8
path = '/Users/au572692/gitProjects/period-spacing/'
filename = 'model-S.fgong'


# Read FGONG file
starg, starl = read_fgong(path+filename)


# Compute brunt-Vaisala frequency using A1 and A4
# Nl^2 = (G M / R^3) * A1 * A4
Nl = np.exp(starl["ln(m/M)"]) / np.power(starl["radius"] / starg["Rphot"], 3)
Nl *= starl["brunt_A"]
Nl *= Gconst * starg["mass"] / np.power(starg["Rphot"], 3)
index = (Nl < 0.0) | (starl["radius"] / starg["Rphot"] > 1.0)
Nl[index] = 1e-99
Nl[0] = 0.0
Nl = np.sqrt(Nl)


# Compute the integral using adaptive mesh
psint, repsint, nerr = apsint(starl["radius"], Nl)
ps = np.sqrt(2e0) * np.power(np.pi, 2) / psint


# Compute the integral using trapezoidal rule
try:
    psTrap = (np.sqrt(2e0) * np.power(np.pi, 2) /
              np.trapz(Nl / starl["radius"], x=starl["radius"]))
except ValueError:
    psTrap = 0.
    pass

print ('A.P.S. using adaptive mesh %16.8f, trapezoidal rule %16.8f, error %4d' 
       %(ps, psTrap, nerr))

"""
"OrbitalPhasePackage.py" is a self-contained version of "phaseFile.py"

Inputs:
    1) A list of Julian dates (an astronomer's way of expressing date and
    time with a single real number.)
    2) Details about a particular binary star.
        2a) Name: The star's name is used for a lookup table with position (Right ascension and
        declination, the sky's equivalent to latitude and longitude)
        2b) Ephemeris type and parameters.  The following types are currently implemented
        (listed along with their parameters:
        Linear:
            T1: Julian date of reference eclipse
            P0: orbital period (days)
        Quadratic:
            T1: Julian date of reference eclipse
            P0: orbital period (days) on T1
            Pdot: period derivative (dimensionless)
        Exponential:
            T1: Julian date of reference eclipse
            T0: Julian date of merger
            P0: initial orbital period (days)
            tau: exponential timescale (days)
            P1: reference period (days)

Outputs:
    For each Julian date, the following are computed:
    1) Barycentric Julian date (Julian date corrected to the center of the solar system)
    2) cycle number (number of full orbits since T1)
    3) orbital phase (fractional part of cycle number)

Usage options:
    1) For testing, run the "main" section of this script.  This will read a sample list of Julian
    dates from the file OrbitalPhasePackageJDs.csv, call the top level routine "OrbitalPhasePackage"
    to compute results and write them to a sample output file OrbitalPhasePackageOutput.csv.
    2) For production, we want to have the dates (and other parameters) specified by a user
    through a web page with results written back to that page.  The computations are still done
    by a call to the top level routine "OrbitalPhasePackage".

Public python packages used:
    astropy
    math
    numpy
    jplephem
    de423 (a specific solar system ephemeris used by jplephem)
Local code combined in bundle:
    "script name": "subroutine name": description
    jd_bjd: jd_bjd, starCat, starCoords, obsCat: convert Julian Date to barycentric Julian date
        contains a list of coordinate locations for stars
    OrbitalPhase: ExpEphemeris: creates ephemeris object
    readHeader: readHeader: reads header of input file (only called by ReadInputFile)
    wc: readData: reads data from input file (only called by ReadInputFile)

Written 12/12/2016 by L. Molnar

"""
from numpy import dot
import astropy.time as atime
import astropy.coordinates as coords
import astropy.units as u
import astropy.constants as const
import jplephem
import de423
from scipy.special import expi
from math import exp

def OrbitalPhasePackage(JD, Ephemeris, EphParams, Star):
    """Top level routine that computes BJD and orbital phase."""
    # Chreate ephemerise objects
    ephem = ExpEphemeris(Ephemeris, EphParams)
    # Convert to BJD
    bjdobj = jd_bjd(JD, Star, 'G98')
    BJD = bjdobj.value
    # Compute orbital phases
    cycles = list()
    phase = list()
    for i in range(len(JD)):
        (phasex, tk, cyclesx, P, OMC, Pdot) = ephem.GetCycle(BJD[i], 0)
        cycles.append(cyclesx)
        phase.append(phasex)
    return BJD, cycles, phase

def ReadInputFile():
    """For testing purposes, read a list of Julian dates from an input file."""
    # Open file, skip over header, read in data from first column, convert to float.
    ifile = open('OrbitalPhasePackageJDs.csv')
    (header, HasEOH) = readHeader(ifile)
    data = readData(ifile,',')
    JD = list()
    for row in data:
        JD.append(float(row[0]))
    return JD

def WriteOutputFile(Ephemeris, EphParams, Star, BJD, cycles, phase):
    """For testing purposes, write the computed values to an output file."""
    # Open file
    ofile = open('OrbitalPhasePackageOutput.csv','w')
    # Write header with information about run
    ofile.write('OrbitalPhasePackage\n'+Star+'\n')
    if Ephemeris == 'Exponential':
        ofile.write('Exponential\nT1,T0,Period (d),tau (d),P_ref (d)\n')
        ofile.write('%.6f,%.1f,%.9f,%.9f,%.9f\n' % (T1,T0,P0,tau,P1))
    elif Ephemeris == 'Quadratic':
        ofile.write('Quadratic\nT1,Period (h),Pdot\n')
        ofile.write('%.6f,%.9f,%.3e\n' % (EphParams[0],EphParams[1],EphParams[2]))
    elif Ephemeris == 'Linear':
        ofile.write('Quadratic\nT1,Period (h)\n')
        ofile.write('%.6f,%.9f,%.3e\n' % (EphParams[0],EphParams[1]))
    # Write out results
    ofile.write('BJD,cycles,Phase\nEOH\n')
    for i in range(len(BJD)):
        ofile.write('{0:13.5f},{1:10.5f},{2:.4f}\n'.format(BJD[i],cycles[i],phase[i]))
    # Close file
    ofile.close()
    return

def jd_bjd(jd, star, observatory, ra=None, dec=None):
     """Convert JD/UTC to BJD/TDB using astropy:
    Inputs:
    jd: Julian date in UTC
    star: Name of a star (coordinates will be looked up in a table)
    obs: Name of an observatory (coordinates will be looked up in a table)
    Output:
    bjdtdb: an astropy time object corrected to bjd with tdb as scale"""

     # Initialise ephemeris from jplephem
     eph = jplephem.Ephemeris(de423)

     # Source unit-vector
     if ra == None:
          src_vec = starCat(star)
     else:
          src_vec = starCoords(star,ra,dec) # directly use coords in deg.

     # Observatory location
     obs = obsCat(observatory)

     # Convert epochs to astropy.time.Time
     t = atime.Time(jd, scale='utc', format='jd', location=obs )

     # Get Earth-Moon barycenter position
     ## NB: jplephem uses Barycentric Dynamical Time, e.g. JD(TDB)
     ## and gives positions relative to solar system barycenter
     barycenter_earthmoon = eph.position('earthmoon', t.tdb.jd)

     # Get Moon position vectors
     moonvector = eph.position('moon', t.tdb.jd)

     # Compute Earth position vectors
     pos_earth = (barycenter_earthmoon - moonvector * eph.earth_share)*u.km

     # Compute and apply BJD correction
     corr = dot(pos_earth.T.value, src_vec.cartesian.xyz)/const.c
     dt = atime.TimeDelta(corr, scale='tdb', format='jd')
     new_jd = t + dt

     # Set default timescale to TDB
     bjdtdb = new_jd.tdb

     return bjdtdb

def starCat(star):
    """ Catalog of stellar coordinates
    Input: star: string name of star
    Output: SkyCoord object with icrs coords """
    # RA-hours, RA-min, RA-sec, dec-sign, dec-deg, dec-amin, dec-asec
    sCat = { 'Algol': [3., 8., 10.132, 1., 40., 57., 20.33],\
             'V0738+2950': [7., 38., 37.363, +1., 29., 50., 31.42],\
             'V859 Cyg': [19., 27., 13., 1., 28., 56., 50.],\
             'KIC 02305372': [19.,27.,57.7,1.,37.,40.,22.1],\
             'KIC 03104113': [19.,12.,9.5,1.,38.,17.,38.8],\
             'KIC 03765708': [19.,43.,42.7,1.,38.,50.,55.7],\
             'KIC 04074532': [19.,43.,4.9,1.,39.,6.,36.0],\
             'KIC 04851217': [19.,43.,20.2,1.,39.,57.,8.3],\
             'KIC 05020034': [19.,37.,12.9,1.,40.,11.,33.4],\
             'KIC 05770860': [18.,56.,45.6,1.,41.,1.,30.0],\
             'KIC 06044064': [19.,29.,24.8,1.,41.,19.,27.1],\
             'KIC 06213131': [19.,37.,12.9,1.,41.,33.,42.1],\
             'KIC 06464285': [19.,50.,41.0,1.,41.,52.,33.6],\
             'KIC 06677225': [19.,8.,33.0,1.,42.,8.,40.2],\
             'KIC 07696778': [19.,45.,9.1,1.,43.,23.,55.3],\
             'KIC 07938468': [18.,47.,46.4,1.,43.,43.,43.0],\
             'KIC 09087918': [19.,25.,56.1,1.,45.,25.,1.9],\
             'KIC 09934052': [18.,49.,43.3,1.,46.,48.,25.6],\
             'KIC 10030943': [19.,54.,12.4,1.,46.,56.,12.5],\
             'KIC 10736223': [19.,35.,23.1,1.,48.,3.,1.1],\
             'KIC 11097678': [19.,50.,28.7,1.,48.,41.,36.6],\
             'KIC 11144556': [19.,40.,28.6,1.,48.,43.,58.4],\
             'KIC 11924311': [19.,46.,39.8,1.,50.,14.,17.9],\
             'KIC 04853067': [19.,44.,53.1,1.,39.,55.,22.8],\
             'KIC 05471619': [19.,50.,48.9,1.,40.,38.,33.0],\
             'KIC 05792093': [19.,28.,48.0,1.,41.,0.,49.3],\
             'KIC 06044543': [19.,29.,56.4,1.,41.,18.,56.2],\
             'KIC 06066379': [19.,52.,31.8,1.,41.,20.,3.5],\
             'KIC 06314173': [19.,56.,57.7,1.,41.,37.,18.1],\
             'KIC 07938870': [18.,48.,44.0,1.,43.,42.,23.0],\
             'KIC 09402652': [19.,25.,6.9,1.,45.,56.,3.1],\
             'KIC 09840412': [19.,42.,10.5,1.,46.,37.,10.6],\
             'KIC 10292413': [19.,52.,29.6,1.,47.,19.,30.4],\
             'KIC 08758161': [19.,34.,28.9,1.,44.,58.,2.3],\
             'KIC 09832227': [19.,29.,16.0,1.,46.,37.,19.9],\
             'KIC 9832227': [19., 29., 15.95, 1., 46., 37., 19.88],\
             'BL Lac': [22., 02., 43.2, +1., 42., 16., 40.],\
             'Y Sex': [10., 02., 48.0, +1., 01., 05., 40.],\
             'E11': [05., 04., 48.8, +1., 53., 00., 18.],\
             'A22': [05., 14., 09.4, +1., 58., 25., 53.],\
             'A34': [05., 16., 02.5, +1., 58., 52., 50.],\
             'Test': [15., 0., 0., -1., 71., 0., 0.] }
    x = sCat[star]
    ra = (x[0] + x[1]/60. + x[2]/3600.)*(360./24.)
    dec = x[3]*(x[4] + x[5]/60. + x[6]/3600.)
    src_vec = coords.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs', distance=coords.Distance(1, u.km))
    return src_vec

def starCoords(star,ra,dec):
    """ Catalog of stellar coordinates
    Input: star: string name of star
         ra: ra (deg)
         dec: dec (deg)
    Output: SkyCoord object with icrs coords """
    src_vec = coords.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs', distance=coords.Distance(1, u.km))
    return src_vec

def obsCat(observatory):
    """ Catalog of observatory coordinates
    Input: star: string name of observatory (MPC code)
    Output: astropy location object """
    # oCat: observatory east longitude (deg), north latitude (deg), elevation (m), name
    oCat = { '645': [-105.82058, +32.779639, 2788., 'Apache Point/SDSS' ],\
             '705': [-105.82058, +32.780312, 2788., 'Apache Point 3.5 m' ],\
             'G98': [-108.65646, +35.52543, 2024., 'Calvin-Rehoboth' ],\
             'NM': [-108.65646, +35.52543, 2024., 'Calvin-Rehoboth' ],\
             'H62': [-85.5883, +42.931026, 200., 'Calvin-Grand Rapids'],\
             'GR': [-85.5883, +42.931026, 200., 'Calvin-Grand Rapids'],\
             'WIRO': [-105.977, 41.097, 2943., 'WIRO'] }
    x = oCat[observatory]
    obs = coords.EarthLocation.from_geodetic(x[0], x[1], height=x[2], ellipsoid='WGS84') # a location tuple
    return obs

class ExpEphemeris:
    """Class for dealing with light curves """
    def __init__(self, Ephemeris, EphParams):
        """Create a new instance of ExpEphemeris.
        Ephemeris: ephemeris type ('Linear', 'Quadratic', or 'Exponential')
        EphParams: list of ephemeris parameters
            None: an empty list
            Linear: BJD0, P0 (days, hours)
            Quadratic: BJD0, P0 (days, hours), Pdot
            Exponential: T1, T0, P0, tau, P1 (days*5)
                T1: time of phase zero; T0: time of exponential infinity; P0: initial period
                tau: exponential time scale; P1: period for reference linear ephemeris
            """
        self.eph = Ephemeris
        self.ephpar = EphParams
        if self.eph == 'Exponential':
            x = self.ephpar[3]/(self.ephpar[0] - self.ephpar[1])
            self.temp = -self.ephpar[3]/self.ephpar[2]*(exp(-x)/x + expi(-x))
        elif ( self.eph != 'Linear' and self.eph != 'Quadratic' and self.eph != 'None' ):
            print 'Bad news, unrecognized ephemeris type: ', self.eph
        return

    def GetCycle(self, t, EphOff):
        """INPUT:
        t: time in user units;
        EphOff: zero point offset to get to ephemeris units
        OUTPUT:
        phase: phase (cycles on [0,1])
        tk: time in ephemeris units
        cycles: real cycle count at t
        P: period(t)
        O-C: linear reference ephemeris phase - exponential ephemeris phase at t
        Pdot: period derivative(t)"""
        tk = t - EphOff # e.g., conversion from BJD to BKD
        Phase = 0.
        cycles = 0.
        P = 0.
        OMC = 0.
        Pdot = 0.
        if self.eph == 'None':
            cycles = 0.
        elif self.eph == 'Linear':
            cycles = (tk - self.ephpar[0])/(self.ephpar[1]/24.)
        elif self.eph == 'Quadratic':
            dt = tk - self.ephpar[0]
            coef = 0.5*self.ephpar[2]/(self.ephpar[1]/24.)**2
            cycles = dt/(self.ephpar[1]/24.) - dt**2*coef
        elif self.eph == 'Exponential':
            if tk < self.ephpar[1]:
                x = self.ephpar[3]/(tk - self.ephpar[1])
                cycles = self.ephpar[3]/self.ephpar[2]*(exp(-x)/x + expi(-x)) + self.temp
                P = self.ephpar[2]*exp(self.ephpar[3]/(tk - self.ephpar[1]))
                OMC = (tk - self.ephpar[0])/self.ephpar[4] - cycles
                Pdot = -self.ephpar[2]*self.ephpar[3]/(tk - self.ephpar[1])**2\
                              *exp(self.ephpar[3]/(tk - self.ephpar[1]))
        Phase = (cycles%1.)
        return Phase, tk, cycles, P, OMC, Pdot

    def GetDate(self, cycles, EphOff):
        """INPUT:
        cycles: real cycle count
        EphOff: zero point offset to get to ephemeris units
        OUTPUT:
        t: time in user units
        tk: time in ephemeris units
        P: period(t)
        O-C: linear reference ephemeris phase - exponential ephemeris phase at t
        Pdot: period derivative(t)"""
        P = 0.
        OMCA = 0.
        PdotA = 0.
        if self.eph == 'None':
            tk = 0.
        elif self.eph == 'Linear':
            tk = cycles*self.ephpar[1]/24. + self.ephpar[0]
        elif self.eph == 'Quadratic':
            temp = cycles*self.ephpar[1]/24. + self.ephpar[0]
            (phase, tkA, cyclesA, PA, OMCA, PdotA) = self.GetCycle(temp, 0.)
            tk = (cycles + OMCA)*self.ephpar[1]/24. + self.ephpar[0]
        elif self.eph == 'Exponential':
            temp = cycles*self.ephpar[4] + self.ephpar[0]
            (phase, tkA, cyclesA, PA, OMCA, PdotA) = self.GetCycle(temp, 0.)
            tk = (cycles + OMCA)*self.ephpar[4] + self.ephpar[0]
            P = self.ephpar[2]*exp(self.ephpar[3]/(tk - self.ephpar[1]))
        t = tk + EphOff
        return t, tk, P, OMCA, PdotA

def readHeader(file):
    """Read and return the header of a text file

    If no header is present, the entire file is returned

    Input: the name of a newly opened text file

    Output:
    1) a list of strings, one per header line
    2) a boolean that is True if an EOH was found

    Intended uses:
    1) Read and print the header of a file.
    2) Read past the header in preparation for reading the data.
    3) Determine the absence of a header.

    Note: Expects newline character (octal 012 = \n).
    A supplemental carriage return is OK (octal 014 = \m or \r).
    Carriage return alone is not (Mac Excel does this when asked to produce
    CSV or MSDOSCSV, but not WindowsCSV.)

    written 6/2009 by LAM
    revised 2/26/2012 by LAM to add second EOH string and to output boolean
    """
    EOH = [ '----------', 'EOH' ] # acceptable "end of header" lines
    header = list()
    while True:
        line = file.readline()
        if line[0:len(EOH[0])] == EOH[0] or \
           line[0:len(EOH[1])] == EOH[1]:   # EOH reached, header returned
            return(header, True)
        else:
            if line == "":                  # EOF reached, entire file returned
                return(header, False)
            else:                           # Another header line read
                header.append(line)

def readData(iFile,delimiter=None):
    """Read and return the data from a file, then close the file.

    Input:
    1) the file name
    2) delimiter for parsing the line.  The default delimiter is white space.

    Output:
    1) a list of lists, each of which has the data from one line.

    Intended uses:
    1) Read in the data from a text file.  Typically, this will be csv data
    and the header will already have been passed over.  However, any delimiter
    may be used, and the status of the header is up to the user.
    2) wc.py uses this to report on size statistics.
    Note: a blank final line is not counted in the tally (so that the lack of
    a trailing newline character is  unimportant).

    written 7/2009 by LAM
    revised 7/6/2009 to use readHeader function
    revised 4/21/2010 to strip the \n off the last item
    revised 2/26/2012 to use revised readHeader (which reports on EOH) and to
        give statistics on full file in absence of a header
    """
    data = list()
    line = iFile.readline()
    while line:
        if(line != ''):
            if delimiter == None:
                data.append(line.strip('\n').split())
            else:
                data.append(line.strip('\n').split(delimiter))
            line = iFile.readline()
    iFile.close()
    return(data)

if __name__ == '__main__':
    # User parameters

    Star = 'KIC 9832227' # Choose one of three preset value

    if Star == 'KIC 9832227':
        T1 = 2455688.49913 # Julian date of zero cycle
        T0 = 2459663.      # Julian date of merger
        P0 = 0.45796151    # days (initial period)
        tau = 0.113541     # days (exponential timescale)
        P1 = 0.4579515     # days (reference period)
        Ephemeris = 'Exponential'
        EphParams = [ T1, T0, P0, tau, P1 ]
    if Star == 'V1309 Sco':
        T1 = 2453118.9582  # JD (time of zero cycle), 2004 detrended model
        T0 = 2455233.5     # JD (time of merger)
        P0 = 1.4456        # days (initial period)
        tau = 15.29        # days (exponential timescale)
        P1 = 1.4349        # days (reference period)
        Ephemeris = 'Exponential'
        EphParams = [ T1, T0, P0, tau, P1 ]
    if Star == 'V859 Cyg':
        T0 = 2456805.89965 # BJD (time of zero cycle)
        P0 = 9.7200735    # hours (initial period)
        Ephemeris = 'Linear'
        EphParams = [ T0, P0 ]

    # End User parameters

    JD = ReadInputFile()
    BJD, cycles, phase = OrbitalPhasePackage(JD, Ephemeris, EphParams, Star)
    WriteOutputFile(Ephemeris, EphParams, Star, BJD, cycles, phase)

    print '\nOrbitalPhasePackage finished'

## Code to allow calls to Refprop from within Python that look like calls
## to the Refprop implementation for Matlab
## Brandon Woodland
## Original completion: 28 March 2012
## Update: 20 August 2014
## This update is designed to allow the code to work with mixtures where the
## concentrations are specified on a mass basis.  It should retain the original
## form of its API so that older code with pure fluids still works.
## 
## This code is based on the file "Refprop2.py" found on the NIST Refprop FAQ
## and written by Bruce Wernick with his information and update history of
## Refprop2.py as follows:

## REFPROP8 library
## Bruce Wernick
## info@coolit.co.za
## Last updated: 6 August 2010

#-------------------------------------------------------------------------------
#  temperature                         K
#  pressure, fugacity                  kPa
#  density                             mol/L
#  composition                         mole fraction
#  quality                             mole basis
#  enthalpy, internal energy           J/mol
#  Gibbs, Helmholtz free energy        J/mol
#  entropy, heat capacity              J/(mol-K)
#  speed of sound                      m/s
#  Joule-Thompson coefficient          K/kPa
#  d(p)/d(rho)                         kPa-L/mol
#  d2(p)/d(rho)2                       kPa-(L/mol)2
#  viscosity                           uPa-s
#  thermal conductivity                W/(m-K)
#  dipole moment                       debye
#  surface tension                     N/m
#-------------------------------------------------------------------------------


#------------------------------------------------------------------------------
"""
%       prop_req    character string showing the requested properties
%                   Each property is represented by one character:
%                           0   Refprop DLL version number
%                           A   Speed of sound [m/s]
%                           B   Volumetric expansivity (beta) [1/K]
%                           C   Cp [J/(kg K)]
%                           D   Density [kg/m^3]
%                           F   Fugacity [kPa] (returned as an array)
%                           G   Gross heating value [J/kg]
%                           H   Enthalpy [J/kg]
%                           I   Surface tension [N/m]
%                           J   Isenthalpic Joule-Thompson coeff [K/kPa]
%                           K   Ratio of specific heats (Cp/Cv) [-]
%                           L   Thermal conductivity [W/(m K)]
%                           M   Molar mass [g/mol]
%                           N   Net heating value [J/kg]
%                           O   Cv [J/(kg K)]
%                           P   Pressure [kPa]
%                           Q   Quality (vapor fraction) (kg/kg)
%                           S   Entropy [J/(kg K)]
%                           T   Temperature [K]
%                           U   Internal energy [J/kg]
%                           V   Dynamic viscosity [Pa*s]
%                           X   Liquid phase & gas phase comp.(mass frac.)
%                           Y   Heat of Vaporization [J/kg]
%                           Z   Compressibility factor
%                           $   Kinematic viscosity [cm^2/s]
%                           %   Thermal diffusivity [cm^2/s]
%                           ^   Prandtl number [-]
%                           )   Adiabatic bulk modulus [kPa]
%                           |   Isothermal bulk modulus [kPa]
%                           =   Isothermal compressibility [1/kPa]
%                           ~   Cstar [-]
%                           `   Throat mass flux [kg/(m^2 s)]
%                           +   Liquid density of equilibrium phase
%                           -   Vapor density of equilibrium phase
%
%                           E   dP/dT (along the saturation line) [kPa/K]
%                           #   dP/dT     (constant rho) [kPa/K]
%                           R   d(rho)/dP (constant T)   [kg/m^3/kPa]
%                           W   d(rho)/dT (constant p)   [kg/(m^3 K)]
%                           !   dH/d(rho) (constant T)   [(J/kg)/(kg/m^3)]
%                           &   dH/d(rho) (constant P)   [(J/kg)/(kg/m^3)]
%                           (   dH/dT     (constant P)   [J/(kg K)]
%                           @   dH/dT     (constant rho) [J/(kg K)]
%                           *   dH/dP     (constant T)   [J/(kg kPa)]
%
%       spec1           first input character:  T, P, H, D, C, R, or M
%                         T, P, H, D:  see above
%                         C:  properties at the critical point
%                         R:  properties at the triple point
%                         M:  properties at Tmax and Pmax
%                            (Note: if a fluid's lower limit is higher
%                             than the triple point, the lower limit will
%                             be returned)
%
%       value1          first input value
%
%       spec2           second input character:  P, D, H, S, U or Q
%
%       value2          second input value
%
%       substance1      file name of the pure fluid (or the first
%                       component of the mixture)
%
%       mixture1        file name of the predefined fluid mixture
%                       with the extension ".mix" included
%
%       substance2,substance3,...substanceN
%                       name of the other substances in the
%                       mixture. Up to 20 substances can be handled.
%                       Valid substance names are equal to the file names
%                       in the C:\Program Files\REFPROP\fluids\' directory.
%
%       x               vector with mass fractions of the substances
%                       in the mixture.
"""
#-------------------------------------------------------------------------------

from ctypes import *
#rp=windll.LoadLibrary("c:/program files/refprop/refprop.dll")
rp = windll.LoadLibrary('/Users/ammarbahman/Google Drive/Purdue/PhD/ACHP/REFPROP91/REFPRP64.dll')
k0=273.15

MaxComps=20
hpth=create_string_buffer('c:/program files (x86)/refprop/', 255)
hfld=create_string_buffer('', 10000)
hfm=create_string_buffer('hmx.bnc', 255)
hrf=create_string_buffer('DEF', 3)
htype=create_string_buffer('NBS', 3)
hmix=create_string_buffer('NBS', 3)
hcomp=create_string_buffer('NBS', 3)
ierr=c_long(0)
herr=create_string_buffer('', 255)
hname=create_string_buffer('', 12)
hn80=create_string_buffer('', 80)
hcas=create_string_buffer('', 12)
nc=c_long(1)
fldindex = c_long()
wm=c_double()
x=(c_double*MaxComps)()  #these arrays are declared by multiplying the desired type by a number.  This generates a new class where the optional initialization is omitted by leaving the () empty.
xmass=(c_double*MaxComps)()
xl=(c_double*MaxComps)()
xv=(c_double*MaxComps)()
xlkg=(c_double*MaxComps)()
xvkg=(c_double*MaxComps)()


#SetFluids eliminates memory leaks due to multiple calls to SETUP to change fluids
#Currently, only one mixture may be specified, and it must be specified as the first "fluid" in fluidlist. More mixtures could be supported by using a different array of mole fractions and setting uninvolved fluids' concentrations to zero.
def SetFluids(fluidlist, FluidRef='DEF'):
    
    global ierr,herr #this global statement ensures that any changes to these values are made to the ones declared outside the functions (with global scope).
    global nc,hfld,hfm,hrf,wm,x,xmass
    
    rp.SETPATHdll(byref(hpth),c_long(255)) #set the path so we don't need to include it in all the fluid names
    SETUP.fluiddict = {} #create a fluid dictionary as an attribute of the SETUP function so it will persist between calls to that function
    FluidName = ''
    nc_list = []
    for fldstr in fluidlist:
        fldargs = fldstr.split(',',1) #split the string into at most two arguments. The first contains the fluid file name or names. The second contains a list of mass concentrations of each fluid.
        FluidName += fldargs[0]+'|'
        nc_list.append(fldstr.count('|') + 1)
        if fldstr.count('|') >= 1: #then the first fluid listed as a mixture.
            SETUP.fluiddict[fldargs[0]] = 0 #the index in the fluid dictionary for the mixture is 0
    FluidName = FluidName[0:-1] #drop the '|' at the end of the last fluid name
    PureFluidList = FluidName.split('|')
    
    i = 1
    for fluid in PureFluidList:
        SETUP.fluiddict[fluid] = i
        i += 1

    hfld.value=FluidName
    hrf.value=FluidRef
    nc.value = len(PureFluidList)

    rp.SETUPdll(byref(nc),byref(hfld),byref(hfm),byref(hrf),byref(ierr),byref(herr),c_long(10000),c_long(255),c_long(3),c_long(255))

    if (ierr.value > 0) and not (ierr.value == 117):
        raise ValueError(herr.value)
    
    nc.value = nc_list[0]
    
    rp.SETNCdll(byref(nc))

def SETUP(fldstr): #fldstr is a string that contains fluids and concentrations. kwargs was used as a placeholder for a number of settings. Currently, only the fluid reference setting is supported in the call to SetFluids.
    #SETUP has to be called each time the fluid is changed.  We will call it every time so that the fluid
    #may be specified in each function call.
    """
    Old SETUP function substituted by SetFluids and SETUP @Brandon
    """
    ''''I could have used *fldargs and given fluid names and components as a tuple. But then every intance of Props()
    would have to be changed to have *fluid or multiple arguments at the end. This would also break the
    convenience functions T_hp(), T_sat(), and IsFluidType() that use fluid as the first input'''
    '''define models and initialize arrays'''

    global ierr,herr #this global statement ensures that any changes to these values are made to the ones declared outside the functions (with global scope).
    global nc,fldindex,hfld,hfm,hrf,wm,x,xmass

    if not hasattr(SETUP,'fluiddict'):
        raise AttributeError('A call to SetFluids must first be made before the PropsSI routine can be run.')
    
    fldargs = fldstr.split(',',1) #split the string into at most two arguments. The first contains the fluid file name or names. The second contains a list of mass concentrations of each fluid.
    FluidName = fldargs[0]
    
    if not FluidName in SETUP.fluiddict:
        raise KeyError('A fluid name was specified that was not in the list provided to SetFluids. Dynamic adding of fluids is not currently supported.')

    fldindex.value = SETUP.fluiddict[FluidName]
    rp.PUREFLDdll(byref(fldindex))

    if len(fldargs) > 1: #an array of mass fractions was specified
        fldargs[1] = eval(fldargs[1]) #turn the array of mass fractions into a list by evaluating the string
        if not len(fldargs[1])==nc.value:
            raise ValueError('The number of mass concentrations does not match the number of fluid components!')
        xmass.__init__(*fldargs[1]) #give all the component mass fractions in fldargs[1] to the c_double array, xmass
# for i in xmass: print i
        rp.XMOLEdll(xmass,x,byref(wm)) #convert the array of mass fractions in xmass to mole fractions and store them in x. Also store the mixture molecular weight in wm.
    else: #no array of mass fractions was specified so we have a pure fluid
        rp.WMOLdll(x,byref(wm)) #the array of mole fractions is ignored if the number of components is equal to 1.
    return





# -- GENERAL FLASH SUBROUTINES --
      
def TPFLSH(t,p):
  '''flash calculation given temperature and pressure'''
  global ierr,herr
  global xl,xv
  t=c_double(t)
  p=c_double(p)
  D,Dl,Dv=c_double(),c_double(),c_double()
  q,e,h,s=c_double(),c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
#  for i in x:  #use this in debugging to check the values in the x array because they cannot be observed in the watch window since it is a ctypes array
#      print i
#      if i <=0.01:
#          break
      
  rp.TPFLSHdll(byref(t),byref(p),x,byref(D),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(h),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return D.value,Dl.value,Dv.value,q.value,e.value*1000,h.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def TDFLSH(t,D):
  '''flash calculation given temperature and density'''
  global ierr,herr
  global xl,xv
  t=c_double(t)
  D=c_double(D)
  p=c_double()
  Dl,Dv=c_double(),c_double()
  q,e,h,s=c_double(),c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.TDFLSHdll(byref(t),byref(D),x,byref(p),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(h),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  
  #'P':0, 'Q':3, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9
  return p.value*1000,Dl.value,Dv.value,q.value,e.value*1000,h.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def PDFLSH(p,D):
  '''flash calculation given pressure and density'''
  global ierr,herr
  global xl,xv
  p=c_double(p)
  D=c_double(D)
  t=c_double()
  Dl,Dv=c_double(),c_double()
  q,e,h,s=c_double(),c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.PDFLSHdll(byref(p),byref(D),x,byref(t),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(h),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,Dl.value,Dv.value,q.value,e.value*1000,h.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def PHFLSH(p,h):
  '''flash calculation given pressure and enthalpy'''
  global ierr,herr
  global xl,xv
  p=c_double(p)
  h=c_double(h)
  t=c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  q,e,s=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.PHFLSHdll(byref(p),byref(h),x,byref(t),byref(D),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,D.value,Dl.value,Dv.value,q.value,e.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def PSFLSH(p,s):
  '''flash calculation given pressure and entropy'''
  global ierr,herr
  p=c_double(p)
  s=c_double(s)
  t=c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  q,e,h=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.PSFLSHdll(byref(p),byref(s),x,byref(t),byref(D),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(h),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,D.value,Dl.value,Dv.value,q.value,e.value*1000,h.value*1000,cv.value*1000,cp.value*1000,w.value

def PEFLSH(p,e):
  '''flash calculation given pressure and energy'''
  global ierr,herr
  p=c_double(p)
  e=c_double(e)
  t=c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  q,s,h=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.PEFLSHdll(byref(p),byref(e),x,byref(t),byref(D),byref(Dl),byref(Dv),xl,xv,byref(q),byref(h),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,D.value,Dl.value,Dv.value,q.value,h.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def THFLSH(t,h,kr=1):
  '''flash calculation given temperature and enthalpy'''
  global ierr,herr
  kr=c_long(kr)
  t=c_double(t)
  h=c_double(h)
  p=c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  q,e,s=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.THFLSHdll(byref(t),byref(h),x,byref(kr),byref(p),byref(D),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return p.value*1000,D.value,Dl.value,Dv.value,q.value,e.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def TSFLSH(t,s,kr=1):
  '''flash calculation given temperature and entropy'''
  global ierr,herr
  global xl,xv
  t=c_double(t)
  s=c_double(s)
  kr=c_long(kr)
  p=c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  q,e,h=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.TSFLSHdll(byref(t),byref(s),x,byref(kr),byref(p),byref(D),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(h),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return p.value*1000,D.value,Dl.value,Dv.value,q.value,e.value*1000,h.value*1000,cv.value*1000,cp.value*1000,w.value

def TEFLSH(t,e,kr=1):
  '''flash calculation given temperature and energy'''
  global ierr,herr
  t=c_double(t)
  e=c_double(e)
  kr=c_long(kr)
  p=c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  q,h,s=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.TEFLSHdll(byref(t),byref(e),x,byref(kr),byref(p),byref(D),byref(Dl),byref(Dv),xl,xv,byref(q),byref(h),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return p.value*1000,D.value,Dl.value,Dv.value,q.value,h.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def DHFLSH(D,h):
  '''flash calculation given density and enthalpy'''
  global ierr,herr
  D=c_double(D)
  h=c_double(h)
  t,p=c_double(),c_double()
  Dl,Dv=c_double(),c_double()
  q,e,s=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.DHFLSHdll(byref(D),byref(h),x,byref(t),byref(p),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,p.value*1000,Dl.value,Dv.value,q.value,e.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def DSFLSH(D,s):
  '''flash calculation given density and entropy'''
  global ierr,herr
  D=c_double(D)
  s=c_double(s)
  t,p=c_double(),c_double()
  Dl,Dv=c_double(),c_double(),c_double()
  q,e,h=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.DSFLSHdll(byref(D),byref(s),x,byref(t),byref(p),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(h),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,p.value*1000,Dl.value,Dv.value,q.value,e.value*1000,h.value*1000,cv.value*1000,cp.value*1000,w.value

def DEFLSH(D,e):
  '''flash calculation given density and energy'''
  global ierr,herr
  D=c_double(D)
  e=c_double(e)
  t,p=c_double(),c_double()
  Dl,Dv=c_double(),c_double()
  q,h,s=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.DEFLSHdll(byref(D),byref(e),x,byref(t),byref(p),byref(Dl),byref(Dv),xl,xv,byref(q),byref(h),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,p.value*1000,Dl.value,Dv.value,q.value,h.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def HSFLSH(h,s):
  '''flash calculation given enthalpy and entropy'''
  global ierr,herr
  h=c_double(h)
  s=c_double(s)
  t,p=c_double(),c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  q,e=c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.HSFLSHdll(byref(h),byref(s),x,byref(t),byref(p),byref(D),byref(Dl),byref(Dv),xl,xv,byref(q),byref(e),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,p.value*1000,D.value,Dl.value,Dv.value,q.value,e.value*1000,cv.value*1000,cp.value*1000,w.value

def ESFLSH(e,s):
  '''flash calculation given energy and entropy'''
  global ierr,herr
  e=c_double(e)
  s=c_double(s)
  t,p=c_double(),c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  q,h=c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.ESFLSHdll(byref(e),byref(s),x,byref(t),byref(p),byref(D),byref(Dl),byref(Dv),xl,xv,byref(q),byref(h),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,p.value*1000,D.value,Dl.value,Dv.value,q.value,h.value*1000,cv.value*1000,cp.value*1000,w.value

def TQFLSH(t,q,kq=2):
  '''flash calculation given temperature and quality'''
  global ierr,herr
  t=c_double(t)
  q=c_double(q)
  kq=c_long(kq)
  p=c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  e,h,s=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.TQFLSHdll(byref(t),byref(q),x,byref(kq),byref(p),byref(D),byref(Dl),byref(Dv),xl,xv,byref(e),byref(h),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return p.value*1000,D.value,Dl.value,Dv.value,e.value*1000,h.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value

def PQFLSH(p,q,kq=2):
  '''flash calculation given pressure and quality'''
  global ierr,herr
  p=c_double(p)
  q=c_double(q)
  kq=c_long(kq)
  t=c_double()
  D,Dl,Dv=c_double(),c_double(),c_double()
  e,h,s=c_double(),c_double(),c_double()
  cv,cp=c_double(),c_double()
  w=c_double()
  rp.PQFLSHdll(byref(p),byref(q),x,byref(kq),byref(t),byref(D),byref(Dl),byref(Dv),xl,xv,byref(e),byref(h),byref(s),byref(cv),byref(cp),byref(w),byref(ierr),byref(herr),c_long(255))
  return t.value,D.value,Dl.value,Dv.value,e.value*1000,h.value*1000,s.value*1000,cv.value*1000,cp.value*1000,w.value      

def CRITP():
  '''critical parameters'''
  global ierr,herr
  tcrit=c_double()
  pcrit=c_double()
  Dcrit=c_double()
  rp.CRITPdll(x,byref(tcrit),byref(pcrit),byref(Dcrit),byref(ierr),byref(herr),c_long(255))
  return tcrit.value,pcrit.value,Dcrit.value

def TRNPRP(t,D):
  '''transport properties of thermal conductivity and
     viscosity as functions of temperature and density
     eta--viscosity (uPa.s)
     tcx--thermal conductivity (W/m.K)'''
  global ierr,herr
  t=c_double(t)
  D=c_double(D)
  eta,tcx=c_double(),c_double()
  rp.TRNPRPdll(byref(t),byref(D),x,byref(eta),byref(tcx),byref(ierr),byref(herr),c_long(255))
  #TODO: conductivity value needs a factor 1000 to get W/mK. Maybe some unit conversion issues
  return eta.value,tcx.value*1000

def QMOLES_to_MASS(q,xv,xl):
          MW_v = c_double()
          MW_l = c_double()
          rp.WMOLdll(xv,byref(MW_v))
          rp.WMOLdll(xl,byref(MW_l))
          return MW_v.value/(q*MW_v.value + (1-q)*MW_l.value)  #conversion factor to change molar-based quality into mass-based quality
      # use QMOLES_to_MASS(output['Q'],xv,xl) as the value for the 'Q' entry
      #in the Molar_to_Mass dictionary in Props() if using a mixture
      
def PropsSI(Prop_des, Prop1, Prop1_val, Prop2, Prop2_val, fldstr):
    
    SETUP(fldstr)
    MW = wm.value
    
    if Prop_des in ['T', 'P', 'D', 'H', 'S', 'U', 'Q', 'O', 'C', 'A', 'V', 'L']:
        if Prop1=='T':
            if Prop2=='P':
                Prop2_val = Prop2_val/1000 #Pa to kPa
                props = TPFLSH(Prop1_val, Prop2_val)
                output = {'D':0, 'Q':3, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='D':
                props = TDFLSH(Prop1_val, Prop2_val/MW)
                output = {'P':0, 'Q':3, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='H':
                Prop2_val = Prop2_val/1000 #J/kg to kJ/kg
                props = THFLSH(Prop1_val, Prop2_val*MW, kr=1) #returns lower density root (if two roots exist) unless you specify kr=2
                output = {'P':0, 'D':1, 'Q':4, 'U':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='S':
                Prop2_val = Prop2_val/1000 #J/kg-K to kJ/kg-K
                props = TSFLSH(Prop1_val, Prop2_val*MW, kr=1) #returns lower density root (if two roots exist) unless you specify kr=2
                output = {'P':0, 'D':1, 'Q':4, 'U':5, 'H':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='U':
                Prop2_val = Prop2_val/1000 #J/kg to kJ/kg
                props = TEFLSH(Prop1_val, Prop2_val*MW, kr=1) #returns lower density root (if two roots exist) unless you specify kr=2
                output = {'P':0, 'D':1, 'Q':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='Q':
                props = TQFLSH(Prop1_val, Prop2_val, kq=2) #input quality is given on a MASS basis (kq=1 for MOLAR basis) - doesn't matter for pure fluids
                output = {'P':0, 'D':1, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            else: 
                raise ValueError('The second input property must be P, D, H, S, U, or Q')
        elif Prop1=='P':
            Prop1_val = Prop1_val/1000 #Pa to kPa
            if Prop2=='D':
                props = PDFLSH(Prop1_val, Prop2_val/MW)
                output = {'T':0, 'Q':3, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='H':
                Prop2_val = Prop2_val/1000 #J/kg to kJ/kg
                props = PHFLSH(Prop1_val, Prop2_val*MW)
                output = {'T':0, 'D':1, 'Q':4, 'U':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='S':
                Prop2_val = Prop2_val/1000 #J/kg-K to kJ/kg-K
                props = PSFLSH(Prop1_val, Prop2_val*MW)
                output = {'T':0, 'D':1, 'Q':4, 'U':5, 'H':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='U':
                Prop2_val = Prop2_val/1000 #J/kg to kJ/kg
                props = PEFLSH(Prop1_val, Prop2_val*MW)
                output = {'T':0, 'D':1, 'Q':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='Q':
                props = PQFLSH(Prop1_val, Prop2_val, kq=2) #input quality is given on a MASS basis (kq=1 for MOLAR basis) - doesn't matter for pure fluids
                output = {'T':0, 'D':1, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='T':
                props = TPFLSH(Prop2_val, Prop1_val) #This is called in case 'P' is the first property and 'T' is the second property entered
                output = {'D':0, 'Q':3, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            else:
                raise ValueError('The input property must be T, D, H, S, U, or Q')
        elif Prop1=='D':
            if Prop2=='H':
                Prop2_val = Prop2_val/1000 #J/kg to kJ/kg
                props = DHFLSH(Prop1_val/MW, Prop2_val*MW)
                output = {'T':0, 'P':1, 'Q':4, 'U':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='S':
                Prop2_val = Prop2_val/1000 #J/kg-K to kJ/kg-K
                props = DSFLSH(Prop1_val/MW, Prop2_val*MW)
                output = {'T':0, 'P':1, 'Q':4, 'U':5, 'H':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='U':
                Prop2_val = Prop2_val/1000 #J/kg to kJ/kg
                props = DEFLSH(Prop1_val/MW, Prop2_val*MW)
                output = {'T':0, 'P':1, 'Q':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='T':
                props = TDFLSH(Prop2_val, Prop1_val/MW)
                output = {'P':0, 'Q':3, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='P':
                Prop2_val = Prop2_val/1000 #Pa to kPa
                props = PDFLSH(Prop2_val, Prop1_val/MW)
                output = {'T':0, 'Q':3, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            else:
                raise ValueError('The second input property must be T, P, H, S, or U')
        elif Prop1=='H':
            Prop1_val = Prop1_val/1000 #J/kg to kJ/kg
            if Prop2=='S':
                Prop2_val = Prop2_val/1000 #J/kg-K to kJ/kg-K
                props = HSFLSH(Prop1_val*MW, Prop2_val*MW)
                output = {'T':0, 'P':1, 'D':2, 'Q':5, 'U':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='T':
                props = THFLSH(Prop2_val, Prop1_val*MW, kr=1) #returns lower density root (if two roots exist) unless you specify kr=2
                output = {'P':0, 'D':1, 'Q':4, 'U':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='P':
                Prop2_val = Prop2_val/1000 #Pa to kPa
                props = PHFLSH(Prop2_val, Prop1_val*MW)
                output = {'T':0, 'D':1, 'Q':4, 'U':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='D':
                props = DHFLSH(Prop2_val/MW, Prop1_val*MW)
                output = {'T':0, 'P':1, 'Q':4, 'U':5, 'S':6, 'O':7, 'C':8, 'A':9}
            else:
                raise ValueError('The second input property must be T, P, D, or S')
        elif Prop1=='U':
            Prop1_val = Prop1_val/1000 #J/kg to kJ/kg
            if Prop2=='S':
                Prop2_val = Prop2_val/1000 #J/kg-K to kJ/kg-K
                props = ESFLSH(Prop1_val*MW, Prop2_val*MW)
                output = {'T':0, 'P':1, 'D':2, 'Q':5, 'H':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='T':
                props = TEFLSH(Prop2_val, Prop1_val*MW, kr=1) #returns lower density root (if two roots exist) unless you specify kr=2
                output = {'P':0, 'D':1, 'Q':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='P':
                Prop2_val = Prop2_val/1000 #Pa to kPa
                props = PEFLSH(Prop2_val, Prop1_val*MW)
                output = {'T':0, 'D':1, 'Q':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='D':
                props = DEFLSH(Prop2_val/MW, Prop1_val*MW)
                output = {'T':0, 'P':1, 'Q':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            else:
                raise ValueError('The second input property must be T, P, D, or S')
        elif Prop1=='S':
            Prop1_val = Prop1_val/1000 #J/kg to kJ/kg-K
            if Prop2=='T':
                props = TSFLSH(Prop2_val, Prop1_val*MW, kr=1) #returns lower density root (if two roots exist) unless you specify kr=2
                output = {'P':0, 'D':1, 'Q':4, 'U':5, 'H':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='P':
                Prop2_val = Prop2_val/1000 #Pa to kPa
                props = PSFLSH(Prop2_val, Prop1_val*MW)
                output = {'T':0, 'D':1, 'Q':4, 'U':5, 'H':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='D':
                props = DSFLSH(Prop1_val/MW, Prop2_val*MW)
                output = {'T':0, 'P':1, 'Q':4, 'U':5, 'H':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='H':
                Prop2_val = Prop2_val/1000 #J/kg to kJ/kg
                props = HSFLSH(Prop2_val*MW, Prop1_val*MW)
                output = {'T':0, 'P':1, 'D':2, 'Q':5, 'U':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='U':
                Prop2_val = Prop2_val/1000 #J/kg to kJ/kg
                props = ESFLSH(Prop2_val*MW, Prop1_val*MW)
                output = {'T':0, 'P':1, 'D':2, 'Q':5, 'H':6, 'O':7, 'C':8, 'A':9}
            else:
                raise ValueError('The second input property must be T, P, D, H, or U')
        elif Prop1=='Q':
            if Prop2=='T':
                props = TQFLSH(Prop2_val, Prop1_val, kq=2) #input quality is given on a MASS basis (kq=1 for MOLAR basis) - doesn't matter for pure fluids
                output = {'P':0, 'D':1, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            elif Prop2=='P':
                Prop2_val = Prop2_val/1000 #Pa to kPa
                props = PQFLSH(Prop2_val, Prop1_val, kq=2) #input quality is given on a MASS basis (kq=1 for MOLAR basis) - doesn't matter for pure fluids
                output = {'T':0, 'D':1, 'U':4, 'H':5, 'S':6, 'O':7, 'C':8, 'A':9}
            else:
                raise ValueError('The second input property must be T or P')
        else:
            raise ValueError('The first input property must be T, P, D, H, S, U, or Q')     
    elif Prop_des=='B' or Prop_des=='E':
        props = CRITP()
        output = {'B':0, 'E':1}
    elif Prop_des=='M':
        props = MW,
        output = {'M':0}
    else:
        raise ValueError('The desired output property must be T, P, D, H, S, U, Q, O, C, A, V, L, B, E, or M')
    
    if Prop_des=='V' or Prop_des=='L':
        if Prop1=='T':
            T = Prop1_val
            if Prop2=='D':
                D = Prop2_val/MW
            else:
                D = props[output['D']]
        elif Prop1=='D':
            D = Prop1_val
            if Prop2=='T':
                T = Prop2_val
            else:
                T = props[output['T']]
        else:
            T = props[output['T']]
            D = props[output['D']]
        if Prop1=='Q':
            Q = Prop1_val
        elif Prop2=='Q':
            Q = Prop2_val
        else:
            Q = props[output['Q']]
                            
        if Q<=0 or Q>=1.0:
            props = TRNPRP(T,D)
            output = {'V':0, 'L':1}
        else: #two phase mixture -> find average property value
            if 'D' in output:
                rho_f_index = output['D'] + 1
                rho_g_index = output['D'] + 2
            elif 'P' in output:
                rho_f_index = output['P'] + 1
                rho_g_index = output['P'] + 2
            else:       #'T' is precedes saturated density values in the output list
                rho_f_index = output['T'] + 1
                rho_g_index = output['T'] + 2
                
            rho_f = props[rho_f_index]
            rho_g = props[rho_g_index]
            
            if Prop_des=='V': 
                out = 0
                output = {'V':0}
            else: 
                out = 1
                output = {'L':0}
            props = TRNPRP(T,rho_f)[out]*(1-Q) + TRNPRP(T,rho_g)[out]*Q
            props = props,  #convert props to a 1-tuple
    if Prop_des=='Q':
        q_mol_to_mass = QMOLES_to_MASS(props[output['Q']],xv,xl)
    else:
        q_mol_to_mass = 1  #a dummy value that is not used
            
    Molar_to_Mass = {'T':1, 'P':1, 'D':MW, 'Q':q_mol_to_mass, 'U':1/MW, 'H':1/MW, 'S':1/MW, 'O':1/MW, 'C':1/MW, 'A':1, 'V':1e-6, 'L':1e-3, 'B':1, 'E':1, 'M':1} #need to change the value of 'Q' if we use mixtures
    return props[output[Prop_des]]*Molar_to_Mass[Prop_des]

def pcrit(FluidName):
    SETUP(FluidName)
    return CRITP()[1]

def Tcrit(FluidName):
    SETUP(FluidName)
    return CRITP()[0]

def T_hp(FluidName, h, p, T_guess):
    return PropsSI('T','P',p,'H',h,FluidName)

def Tsat(FluidName, p, Q, T_guess):
    if Q<0 or Q>1:
        raise ValueError('Quality must be a value between 0.0 and 1.0')
    return Props('T','P',p,'Q',Q,FluidName)

def IsFluidType(FluidName, Type):
    if Type=='Brine':   #Only CoolProp handles brine so it cannot be possible for the fluid type to be 'Brine' here
        return 0
    else:
        return 1        #I don't have any routines which set the fluid type so I can't really check for any other fluid types


#TODO: Modified for Incompressible as hot fluid
def IsFluidIncompType(FluidName, Type):
    if Type=='INCOMP::T66':
        return 1
    else:
        pass
        #raise Exception('Incompressible fluid not implemented yet')
        


def UseSaturationLUT(yes_no):
    return None         #We don't have the capability to generate a lookup table of saturation properties yet

if __name__ == '__main__':
    mixture = 'r134a.fld|r245fa.fld, [0.5, 0.5]' #if concentration values do not sum to 1, they will be normalized to sum to 1 within refprop.
    fluid1 = 'r134a.fld'
    fluid2 = 'r245fa.fld'
    fluid3 = 'water.fld'
    wf = fluid1
    # fluid = 'co2.fld|acetone.fld, [0.5, 0.5]' #if a concentration is given for a single component, it is ignored

    import time
    t0 = time.time()
    
    SetFluids([mixture, fluid3, 'r245fa.fld']) #'air.ppf'])
    
    T = 60+273.15
    D = 161.37
    
    for i in range(5):
        PropsSI('H','T',T,'D',D,'r245fa.fld')
        PropsSI('H','T',T,'D',D,'r134a.fld')
        PropsSI('H','T',T,'D',D,'water.fld')
    
    #print 'h for air =', PropsSI('H','T',T,'P',100*1000,'air.ppf')
    T = 80+273.15
    P = 1000
    
    print 'mu for r134a =', PropsSI('V','T',T,'Q',0,'r134a.fld')
    print 'mu for r245fa =', PropsSI('V','T',T,'Q',0,'r245fa.fld')
    print 'mu for r134a-r245fa =', PropsSI('V','T',T,'Q',0,'r134a.fld|r245fa.fld, [0.625,0.375]')
    
    # PropsSI('H','T',T,'D',D,'r134a.fld|r245fa.fld, [0.625,0.375]')
    # print 'h for 245fa =', PropsSI('H','T',T,'D',D,'r245fa.fld')
    # print 'h for r245fa =', PropsSI('H','T',T,'D',D,'r245fa.fld')
    # print 'mu for 245fa =', PropsSI('V','T',T,'D',D,'r245fa.fld')
    # print 'mu for water =', PropsSI('V','T',T,'D',D,'water.fld')
    #
    # print 'h for r134a =', PropsSI('H','T',293.147630045,'P',571.66502251*1000,'r134a.fld')
    # print 'mu for r245fa =', PropsSI('V','T',T,'D',D,'r245fa.fld') #This line will fail strangely. Transport properties only work for the first version of the same fluid specified in the fluidlist. However, thermodynamic properties seem to work fine for both.
    # print PropsSI('V','T',60+273.15,'D',161.37,'r245fa.fld')
    # print PropsSI('Q', 'T', 400, 'P', 1000*1000, fluid)
    # print PropsSI('Q', 'T', 300, 'P', 1000*1000, fluid)
    # fluid = 'co2.fld|acetone.fld, [0.5, 0.5]'
    # Q = range(-1,12,1)
    # for val in Q:
    #   h = PropsSI('H', 'Q', float(val)/10, 'P', 1500*1000, mixture)
    #   q = PropsSI('Q','P', 1500*1000, 'H', h*1000, mixture)
    #   t = PropsSI('T','P', 1500*1000, 'H', h*1000, mixture) - 273.15
    # print 'Q = ', float(val)/10, ', h = ', h, ', q = ', q, ', T = ', t
    # print PropsSI('Q', 'T', 350.5376, 'D', 1000, fluid1)
    # print PropsSI('P', 'T', 400, 'H', 511.5*1000, fluid2)
    # print PropsSI('P', 'T', 373, 'Q', 1.0, fluid3)
    # print PropsSI('P', 'T', 400, 'U', 481.206*1000, fluid)
    # print PropsSI('P', 'T', 312.5376, 'Q', 0.5, fluid)
    # print 'testing an error here'
    # print PropsSI('H', 'P', 1000*1000, 'Q', 1.1, fluid)
    print 'Pcrit [kPa]:',pcrit(fluid1)
    # print T_hp(fluid, 511.5, 1000, 300)
    # print Tsat(fluid, 1000, 0.1, 0)
    # print PropsSI('V', 'P', 2000*1000, 'Q', 1.000, fluid)
    # print PropsSI('L', 'P', 2000*1000, 'Q', 1.000, fluid)
    # print PropsSI('E', 'T', 0, 'P', 0, fluid)
    # print PropsSI('M','T',0,'P',0, fluid)
    t1 = time.time()
    print 'time to call functions in Windows =', t1 - t0,'secs'
      
      
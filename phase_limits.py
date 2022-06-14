from parameters import Tc, pc, Tt
import math as m

def get_SVP_vap_liq(T):  #  in K
    #SVP=saturation vapor pressure
    # between the vapor and the liquid phase
    # Wagner and Pruss 2002, valid up to the critical point, with an incertitude <0.1%
    # Eq. 2.5
    th = 1. - T / Tc
    if not (T >= Tt and T <= Tc):
        raise "Outside of the valid temperature range for this saturation vapor pressure"
    else:
        Psat = Tc / T * (
                    -7.85951783 * th +
                    1.84408259 * th ** 1.5 -
                    11.7866497 * th ** 3 +
                    22.6807411 * th ** 3.5 -
                    15.9618719 * th ** 4 +
                    1.80122502 * th ** 7.5)

    return m.exp(Psat) * pc # Pascals


def get_SVP_vap_ice(T):
    # SVP=saturation vapor pressure
    # between the vapor and ice phase
    # valid from 273.16K to 190K, 2.21 in Wagner and Pruss 2002

    if (T > Tt):
        raise "Outside of the valid temperature range for this saturation vapor pressure"
    else:
        th = T / Tt
        Psat = -13.928169 * (1. - th ** (-1.5)) + 34.7078238 * (1. - th ** (-1.25))


    return m.exp(Psat) * 611.657 # Pascals, the pressure multiplying the exp(Psat) is knowingly different from Pt
    # see Wagner and Pruss 2002 and Wagner et al. 1994 for more details


def HPice(T):  # give Pmelt in PASCALS

    pmelt = 0.

    if (T >= 273.31 and T < 353.5):  # ice VI
        pmelt = (1. - 1.07476 * (1. - (T / 273.31) ** 4.6)) * 632.4  # MPa, Wagner 1994

    elif (T >= 353.5 and T < 800.):  # ice VII
        theta = T / 355.
        pmelt = 0.85 * ((theta) ** 3.47 - 1.) + 2.17  # GPa, Lin+ 2004
        pmelt = pmelt * 1000.  # MPa

    elif (T >= 800 and T <= 1273. + domega):
        pmelt = 3.3488 + 1.8696e-5 * T ** 2  # fitted from Hernandez+ 2018
        pmelt = pmelt * 1000.

    #    elif (T>1419. and T<=5000.):
    #        pmelt= 26.528+1.2261e-3*T+6.2997e-6*T**2 #fitted from Hernandez+ 2018
    #        pmelt=pmelt*1000.
    #
    #    elif (T>5000.):
    #        pmelt=190e3 #Millot+ 2018

    if pmelt > 0.:
        pmelt = pmelt * 1e6  # Pa

    return pmelt
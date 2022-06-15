from parameters import Tc, pc, Tt, domega
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
        Psat = (-13.928169 * (1. - th ** (-1.5)) +
                34.7078238 * (1. - th ** (-1.25)) )


    return m.exp(Psat) * 611.657 # Pascals, the pressure multiplying the exp(Psat) is knowingly different from Pt
    # see Wagner and Pruss 2002 and Wagner et al. 1994 for more details

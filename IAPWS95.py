from parameters import *
import math as m

#################################################################
#################################################################
#################################################################
#################################################################

def search_p(rhoi, T, ps):
    delt = rhoi / rhoc
    tau = Tc / T
    phir, dphird, dphirt, dphirtt, dphirdt, dphirdd = compute_phir(delt, tau)

    return ((1. + delt * dphird) * rhoi * Rig * T - ps) / ps  # everything in Pa


#################################################################
#################################################################
#################################################################
#################################################################

def pressure(rhoi, T):
    delt = rhoi / rhoc
    tau = Tc / T

    phir, dphird, dphirt, dphirtt, dphirdt, dphirdd = compute_phir(delt, tau)

    return (1. + delt * dphird) * rhoi * Rig * T


#################################################################
#################################################################
#################################################################
#################################################################


def compute_phi0(delt, tau):

    sum = 0.
    dsum = 0.
    dstt = 0.
    for i in range(3, 8):
        sum = sum + n0[i] * m.log(1. - m.exp(-1. * gam0[i] * tau))
        dsum = dsum + n0[i] * gam0[i] * ((1. - m.exp(-1. * gam0[i] * tau)) ** (-1.) - 1.)
        dstt = dstt + n0[i] * gam0[i] ** 2. * m.exp(-1. * gam0[i] * tau) * (1. - m.exp(-1. * gam0[i] * tau)) ** (-2.)

    phi0 = m.log(delt) + n0[0] + n0[1] * tau + n0[2] * m.log(tau) + sum
    dphi0t = n0[1] + n0[2] / tau + dsum
    dphi0tt = -1. * n0[2] / tau ** 2. - dstt

    return phi0, dphi0t, dphi0tt



#################################################################
#################################################################
#################################################################
#################################################################


def compute_phir(delt, tau):


    sum1 = 0.
    ds1d = 0.
    ds1t = 0.
    ds1tt = 0.
    ds1dd = 0.
    ds1dt = 0.

    for i in range(0, 7):
        sum1 = sum1 + n[i] * delt ** d[i] * tau ** ti[i]
        ds1d = ds1d + n[i] * d[i] * delt ** (d[i] - 1.) * tau ** ti[i]
        ds1t = ds1t + n[i] * ti[i] * delt ** d[i] * tau ** (ti[i] - 1.)
        ds1tt = ds1tt + n[i] * ti[i] * (ti[i] - 1.) * delt ** d[i] * tau ** (ti[i] - 2.)
        ds1dd = ds1dd + n[i] * d[i] * (d[i] - 1.) * delt ** (d[i] - 2.) * tau ** ti[i]
        ds1dt = ds1dt + n[i] * d[i] * ti[i] * tau ** (ti[i] - 1.) * delt ** (d[i] - 1.)



    sum2 = 0.
    ds2d = 0.
    ds2t = 0.
    ds2tt = 0.
    ds2dd = 0.
    ds2dt = 0.

    for i in range(7, 51):
        sum2 = sum2 + n[i] * delt ** (d[i]) * tau ** ti[i] * m.exp(-1. * delt ** c[i])
        ds2d = ds2d + n[i] * m.exp(-1. * delt ** c[i]) * (
                    delt ** (d[i] - 1.) * tau ** ti[i] * (d[i] - c[i] * delt ** c[i]))
        ds2t = ds2t + n[i] * ti[i] * delt ** d[i] * tau ** (ti[i] - 1.) * m.exp(-1. * delt ** c[i])
        ds2tt = ds2tt + n[i] * ti[i] * (ti[i] - 1.) * delt ** d[i] * tau ** (ti[i] - 2.) * m.exp(-1. * delt ** c[i])
        ds2dd = ds2dd + n[i] * m.exp(-1. * delt ** c[i]) * (delt ** (d[i] - 2.) * tau ** ti[i] * (
                    (d[i] - c[i] * delt ** c[i]) * (d[i] - 1. - c[i] * delt ** c[i]) - c[i] ** 2. * delt ** c[i]))
        ds2dt = ds2dt + n[i] * ti[i] * tau ** (ti[i] - 1.) * delt ** (d[i] - 1.) * (d[i] - c[i] * delt ** c[i]) * m.exp(
            -1. * delt ** c[i])



    sum3 = 0.
    ds3d = 0.
    ds3t = 0.
    ds3tt = 0.
    ds3dd = 0.
    ds3dt = 0.
    for i in range(51, 54):
        sum3 = sum3 + n[i] * delt ** d[i] * tau ** ti[i] * m.exp(
            -1. * al[i] * (delt - eps[i]) ** 2. - be[i] * (tau - gam[i]) ** 2.)
        ds3d = ds3d + n[i] * delt ** d[i] * tau ** ti[i] * m.exp(
            -1. * al[i] * (delt - eps[i]) ** 2. - be[i] * (tau - gam[i]) ** 2.) * (
                           d[i] / delt - 2. * al[i] * (delt - eps[i]))
        ds3t = ds3t + n[i] * delt ** d[i] * tau ** ti[i] * m.exp(
            -1. * al[i] * (delt - eps[i]) ** 2. - be[i] * (tau - gam[i]) ** 2.) * (
                           ti[i] / tau - 2. * be[i] * (tau - gam[i]))
        ds3tt = ds3tt + n[i] * delt ** d[i] * tau ** ti[i] * m.exp(
            -1. * al[i] * (delt - eps[i]) ** 2. - be[i] * (tau - gam[i]) ** 2.) * (
                            (ti[i] / tau - 2. * be[i] * (tau - gam[i])) ** 2. - ti[i] / tau ** 2. - 2. * be[i])

        ds3dd = ds3dd + n[i] * tau ** ti[i] * m.exp(
            -1. * al[i] * (delt - eps[i]) ** 2. - be[i] * (tau - gam[i]) ** 2.) * (
                            -2. * al[i] * delt ** d[i] + 4. * al[i] ** 2 * delt ** d[i] * (delt - eps[i]) ** 2. - 4. *
                            d[i] * al[i] * delt ** (d[i] - 1.) * (delt - eps[i]) + d[i] * (d[i] - 1.) * delt ** (
                                        d[i] - 2.))

        ds3dt = ds3dt + n[i] * delt ** d[i] * tau ** ti[i] * m.exp(
            -1. * al[i] * (delt - eps[i]) ** 2. - be[i] * (tau - gam[i]) ** 2.) * (
                            ti[i] / tau - 2. * be[i] * (tau - gam[i])) * (d[i] / delt - 2. * al[i] * (delt - eps[i]))



    sum4 = 0.
    ds4d = 0.
    ds4t = 0.
    ds4tt = 0.
    ds4dd = 0.
    ds4dt = 0.
    for i in range(54, 56):
        psi = m.exp(-1. * C[i] * (delt - 1.) ** 2. - D[i] * (tau - 1.) ** 2.)  # ok
        theta = (1. - tau) + A[i] * ((delt - 1.) ** 2.) ** (1. / (2. * bet[i]))  # ok
        Delta = theta ** 2 + B[i] * ((delt - 1.) ** 2.) ** (a[i])  # ok

        dDeltad = (delt - 1.) * (
                    A[i] * theta * 2. / bet[i] * ((delt - 1.) ** 2.) ** (1. / (2. * bet[i]) - 1.) + 2. * B[i] * a[i] * (
                        (delt - 1.) ** 2.) ** (a[i] - 1.))
        dDeltbid = b[i] * Delta ** (b[i] - 1.) * dDeltad
        dDeltbit = -2. * theta * b[i] * Delta ** (b[i] - 1.)
        dpsid = -2. * C[i] * (delt - 1.) * psi
        dpsit = -2. * D[i] * (tau - 1.) * psi

        dpsitt = (2. * D[i] * (tau - 1.) ** 2. - 1.) * 2. * D[i] * psi
        dpsidd = (2. * C[i] * (delt - 1) ** 2. - 1.) * 2. * C[i] * psi
        dpsidt = 4. * C[i] * D[i] * (delt - 1.) * (tau - 1.) * psi

        dDeldd = 1. / (delt - 1.) * dDeltad + (delt - 1.) ** 2 * (
                    4. * B[i] * a[i] * (a[i] - 1.) * ((delt - 1.) ** 2.) ** (a[i] - 2.) + 2. * A[i] ** 2. * (
                        1. / bet[i]) ** 2. * (((delt - 1.) ** 2.) ** (1. / (2. * bet[i]) - 1.)) ** 2 + A[
                        i] * theta * 4. / bet[i] * (1. / (2. * bet[i]) - 1.) * ((delt - 1.) ** 2.) ** (
                                1. / (2. * bet[i]) - 2.))
        dDelbitt = 2. * b[i] * Delta ** (b[i] - 1.) + 4. * theta ** 2 * b[i] * (b[i] - 1.) * Delta ** (b[i] - 2.)
        dDelbidd = b[i] * (Delta ** (b[i] - 1.) * dDeldd + (b[i] - 1.) * Delta ** (b[i] - 2.) * dDeltad ** 2.)
        dDelbidt = -1. * A[i] * b[i] * 2. / bet[i] * Delta ** (b[i] - 1.) * (delt - 1.) * ((delt - 1.) ** 2.) ** (
                    1. / (2. * bet[i]) - 1.) - 2. * theta * b[i] * (b[i] - 1.) * Delta ** (b[i] - 2.) * dDeltad

        sum4 = sum4 + n[i] * Delta ** (b[i]) * delt * psi
        ds4t = ds4t + n[i] * delt * (dDeltbit * psi + Delta ** b[i] * dpsit)
        ds4d = ds4d + n[i] * (Delta ** b[i] * (psi + delt * dpsid) + dDeltbid * delt * psi)
        ds4tt = ds4tt + n[i] * delt * (dDelbitt * psi + 2. * dDeltbit * dpsit + Delta ** b[i] * dpsitt)

        ds4dd = ds4dd + n[i] * (Delta ** b[i] * (2. * dpsid + delt * dpsidd) + 2. * dDeltbid * (
                    psi + delt * dpsid) + dDelbidd * delt * psi)
        ds4dt = ds4dt + n[i] * (Delta ** b[i] * (dpsit + delt * dpsidt) + delt * dDeltbid * dpsit + dDeltbit * (
                    psi + delt * dpsid) + dDelbidt * delt * psi)

    phir = sum1 + sum2 + sum3 + sum4
    dphird = ds1d + ds2d + ds3d + ds4d
    dphirt = ds1t + ds2t + ds3t + ds4t
    dphirtt = ds1tt + ds2tt + ds3tt + ds4tt
    dphirdd = ds1dd + ds2dd + ds3dd + ds4dd
    dphirdt = ds1dt + ds2dt + ds3dt + ds4dt

    return phir, dphird, dphirt, dphirtt, dphirdt, dphirdd
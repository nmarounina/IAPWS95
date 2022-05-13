


#################################################################
#################################################################
#################################################################
#################################################################

def search_p(rhoi, T, ps):
    delt = rhoi / rhoc
    tau = Tc / T
    phir, dphird, dphirt, dphirtt, dphirdt, dphirdd = CALC_phir1(delt, tau)

    return ((1. + delt * dphird) * rhoi * Rig * T - ps) / ps  # everything in Pa


#################################################################
#################################################################
#################################################################
#################################################################

def pressure(rhoi, T):
    delt = rhoi / rhoc
    tau = Tc / T

    phir, dphird, dphirt, dphirtt, dphirdt, dphirdd = CALC_phir1(delt, tau)

    return (1. + delt * dphird) * rhoi * Rig * T


#################################################################
#################################################################
#################################################################
#################################################################


def CALC_phi01(delt, tau):
    n0 = nu.zeros(8)
    gam = nu.zeros(8)

    n0[0] = -8.32044648201;
    gam[0] = 0.
    n0[1] = 6.6832105268;
    gam[1] = 0.
    n0[2] = 3.00632;
    gam[2] = 0.
    n0[3] = 0.012436;
    gam[3] = 1.28728967
    n0[4] = 0.97315;
    gam[4] = 3.53734222
    n0[5] = 1.27950;
    gam[5] = 7.74073708
    n0[6] = 0.96956;
    gam[6] = 9.24437796
    n0[7] = 0.24873;
    gam[7] = 27.5075105

    sum = 0.
    dsum = 0.
    dstt = 0.
    for i in range(3, 8):
        sum = sum + n0[i] * m.log(1. - m.exp(-1. * gam[i] * tau))
        dsum = dsum + n0[i] * gam[i] * ((1. - m.exp(-1. * gam[i] * tau)) ** (-1.) - 1.)
        dstt = dstt + n0[i] * gam[i] ** 2. * m.exp(-1. * gam[i] * tau) * (1. - m.exp(-1. * gam[i] * tau)) ** (-2.)

    phi0 = m.log(delt) + n0[0] + n0[1] * tau + n0[2] * m.log(tau) + sum
    dphi0t = n0[1] + n0[2] / tau + dsum
    dphi0tt = -1. * n0[2] / tau ** 2. - dstt

    return phi0, dphi0t, dphi0tt


#################################################################
#################################################################
#################################################################
#################################################################


def CALC_phir1(delt, tau):
    # ---------------------INITIALIZATIONS

    c = nu.zeros(56)
    d = nu.zeros(56)
    ti = nu.zeros(56)
    n = nu.zeros(56)
    al = nu.zeros(56)
    be = nu.zeros(56)
    gam = nu.zeros(56)
    eps = nu.zeros(56)
    a = nu.zeros(56)
    b = nu.zeros(56)
    B = nu.zeros(56)
    C = nu.zeros(56)
    D = nu.zeros(56)
    A = nu.zeros(56)
    bet = nu.zeros(56)
    # ------------------------------------

    c[0] = 0.;
    d[0] = 1.;
    ti[0] = -0.5;
    n[0] = 0.12533547935523e-1;
    c[1] = 0.;
    d[1] = 1.;
    ti[1] = 0.875;
    n[1] = 0.78957634722828e1;
    c[2] = 0.;
    d[2] = 1.;
    ti[2] = 1.;
    n[2] = -0.87803203303561e1;
    c[3] = 0.;
    d[3] = 2.;
    ti[3] = 0.5;
    n[3] = 0.31802509345418e0;
    c[4] = 0.;
    d[4] = 2.;
    ti[4] = 0.75;
    n[4] = -0.26145533859358e0;
    c[5] = 0.;
    d[5] = 3.;
    ti[5] = 0.375;
    n[5] = -0.78199751687981e-2;
    c[6] = 0.;
    d[6] = 4.;
    ti[6] = 1.;
    n[6] = 0.88089493102134e-2;

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

    c[7] = 1.;
    d[7] = 1.;
    ti[7] = 4.;
    n[7] = -0.66856572307965
    c[8] = 1.;
    d[8] = 1.;
    ti[8] = 6.;
    n[8] = 0.20433810950965
    c[9] = 1.;
    d[9] = 1.;
    ti[9] = 12.;
    n[9] = -0.66212605039687e-4  # >>>> ti
    c[10] = 1.;
    d[10] = 2.;
    ti[10] = 1.;
    n[10] = -0.19232721156002
    c[11] = 1.;
    d[11] = 2.;
    ti[11] = 5.;
    n[11] = -0.25709043003438
    c[12] = 1.;
    d[12] = 3.;
    ti[12] = 4.;
    n[12] = 0.16074868486251
    c[13] = 1.;
    d[13] = 4.;
    ti[13] = 2.;
    n[13] = -0.40092828925807e-1
    c[14] = 1.;
    d[14] = 4.;
    ti[14] = 13.;
    n[14] = 0.39343422603254e-6
    c[15] = 1.;
    d[15] = 5.;
    ti[15] = 9.;
    n[15] = -0.75941377088144e-5
    c[16] = 1.;
    d[16] = 7.;
    ti[16] = 3.;
    n[16] = 0.56250979351888e-3
    c[17] = 1.;
    d[17] = 9.;
    ti[17] = 4.;
    n[17] = -0.15608652257135e-4
    c[18] = 1.;
    d[18] = 10.;
    ti[18] = 11.;
    n[18] = 0.11537996422951e-8
    c[19] = 1.;
    d[19] = 11.;
    ti[19] = 4.;
    n[19] = 0.36582165144204e-6
    c[20] = 1.;
    d[20] = 13.;
    ti[20] = 13.;
    n[20] = -0.13251180074668e-11
    c[21] = 1.;
    d[21] = 15.;
    ti[21] = 1.;
    n[21] = -0.62639586912454e-9
    c[22] = 2.;
    d[22] = 1.;
    ti[22] = 7.;
    n[22] = -0.10793600908932
    c[23] = 2.;
    d[23] = 2.;
    ti[23] = 1.;
    n[23] = 0.17611491008752e-1
    c[24] = 2.;
    d[24] = 2.;
    ti[24] = 9.;
    n[24] = 0.22132295167546
    c[25] = 2.;
    d[25] = 2.;
    ti[25] = 10.;
    n[25] = -0.40247669763528
    c[26] = 2.;
    d[26] = 3.;
    ti[26] = 10.;
    n[26] = 0.58083399985759
    c[27] = 2.;
    d[27] = 4.;
    ti[27] = 3.;
    n[27] = 0.49969146990806e-2
    c[28] = 2.;
    d[28] = 4.;
    ti[28] = 7.;
    n[28] = -0.31358700712549e-1
    c[29] = 2.;
    d[29] = 4.;
    ti[29] = 10.;
    n[29] = -0.74315929710341
    c[30] = 2.;
    d[30] = 5.;
    ti[30] = 10.;
    n[30] = 0.47807329915480
    c[31] = 2.;
    d[31] = 6.;
    ti[31] = 6.;
    n[31] = 0.20527940895948e-1
    c[32] = 2.;
    d[32] = 6.;
    ti[32] = 10.;
    n[32] = -0.13636435110343
    c[33] = 2.;
    d[33] = 7.;
    ti[33] = 10.;
    n[33] = 0.14180634400617e-1
    c[34] = 2.;
    d[34] = 9.;
    ti[34] = 1.;
    n[34] = 0.83326504880713e-2  # >>>>?
    c[35] = 2.;
    d[35] = 9.;
    ti[35] = 2.;
    n[35] = -0.29052336009585e-1
    c[36] = 2.;
    d[36] = 9.;
    ti[36] = 3.;
    n[36] = 0.38615085574206e-1
    c[37] = 2.;
    d[37] = 9.;
    ti[37] = 4.;
    n[37] = -0.20393486513704e-1
    c[38] = 2.;
    d[38] = 9.;
    ti[38] = 8.;
    n[38] = -0.16554050063734e-2
    c[39] = 2.;
    d[39] = 10.;
    ti[39] = 6.;
    n[39] = 0.19955571979541e-2
    c[40] = 2.;
    d[40] = 10.;
    ti[40] = 9.;
    n[40] = 0.15870308324157e-3
    c[41] = 2.;
    d[41] = 12.;
    ti[41] = 8.;
    n[41] = -0.16388568342530e-4
    c[42] = 3.;
    d[42] = 3.;
    ti[42] = 16.;
    n[42] = 0.43613615723811e-1
    c[43] = 3.;
    d[43] = 4.;
    ti[43] = 22.;
    n[43] = 0.34994005463765e-1
    c[44] = 3.;
    d[44] = 4.;
    ti[44] = 23.;
    n[44] = -0.76788197844621e-1
    c[45] = 3.;
    d[45] = 5.;
    ti[45] = 23.;
    n[45] = 0.22446277332006e-1
    c[46] = 4.;
    d[46] = 14.;
    ti[46] = 10.;
    n[46] = -0.62689710414685e-4
    c[47] = 6.;
    d[47] = 3.;
    ti[47] = 50.;
    n[47] = -0.55711118565645e-9
    c[48] = 6.;
    d[48] = 6.;
    ti[48] = 44.;
    n[48] = -0.19905718354408
    c[49] = 6.;
    d[49] = 6.;
    ti[49] = 46.;
    n[49] = 0.31777497330738
    c[50] = 6.;
    d[50] = 6.;
    ti[50] = 50.;
    n[50] = -0.11841182425981

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

    c[51] = 0.;
    d[51] = 3.;
    ti[51] = 0.;
    n[51] = -0.31306260323435e2;
    al[51] = 20.;
    be[51] = 150;
    gam[51] = 1.21;
    eps[51] = 1.;
    c[52] = 0.;
    d[52] = 3.;
    ti[52] = 1.;
    n[52] = 0.31546140237781e2;
    al[52] = 20.;
    be[52] = 150;
    gam[52] = 1.21;
    eps[52] = 1.;
    c[53] = 0.;
    d[53] = 3.;
    ti[53] = 4.;
    n[53] = -0.25213154341695e4;
    al[53] = 20.;
    be[53] = 250;
    gam[53] = 1.25;
    eps[53] = 1.;

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

    a[54] = 3.5;
    b[54] = 0.85;
    B[54] = 0.2;
    n[54] = -0.14874640856724;
    C[54] = 28.;
    D[54] = 700;
    A[54] = 0.32;
    bet[54] = 0.3;
    a[55] = 3.5;
    b[55] = 0.95;
    B[55] = 0.2;
    n[55] = 0.31806110878444;
    C[55] = 32.;
    D[55] = 800;
    A[55] = 0.32;
    bet[55] = 0.3;

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
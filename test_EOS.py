from parameters import M_h2o
import IAPWS95


point=IAPWS95.DataPoint( 700, 22.5e6 )
print( " %s %f %s" %("Density =",  point.rho*M_h2o, "kg.m-3" ))
print(" %s %f %s" %("Entropy =" , point.s, "J.mol-1.K-1"))
print(" %s %f %s" %("Heat capacity =" , point.Cp, "J.mol-1.K-1"))
print(" %s %f %s" %("Internal enerjy =" , point.u, "J.mol-1"))
from parameters import M_h2o, Rig
import IAPWS95
from scipy import optimize

point=IAPWS95.DataPoint( 700, 22.5e6 )
print( " %s %f %s" %("Density =",  point.rho*M_h2o, "kg.m-3" ))
print(" %s %f %s" %("Entropy =" , point.s, "J.mol-1.K-1"))
print(" %s %f %s" %("Heat capacity =" , point.Cp, "J.mol-1.K-1"))
print(" %s %f %s" %("Internal energy =" , point.u, "J.mol-1"))



#Benchmark, see table 6.6 of Wagner and Pruss (2002)
#T=500 K and rho=838.025 kg.m-3 :
def func_search_for_p(p):

    point=IAPWS95.DataPoint( 500., p )
    #print( "%e %f" %(p, point.rho*M_h2o) )

    return (point.rho*M_h2o - 838.025)/838.025

pressure=optimize.bisect(func_search_for_p, 1e7, 2.63e8)



point_benchmark=IAPWS95.DataPoint( 500., pressure)

print("\n")
print( "Benchmark:")
print( "T=", point_benchmark.T, "rho=", point_benchmark.rho*M_h2o)
print("\n")
print("Ideal-gas part and derivatives:")
print( point_benchmark.phi0 )
print( point_benchmark.dphi0t )
print( point_benchmark.dphi0tt )

print("\n")
print("Residual part and derivatives:")
print( point_benchmark.phir )
print( point_benchmark.dphird )
print( point_benchmark.dphirt )
print( point_benchmark.dphirdd )
print( point_benchmark.dphirtt )
print( point_benchmark.dphirdt )

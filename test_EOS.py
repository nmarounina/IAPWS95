from parameters import M_h2o
import IAPWS95


point=IAPWS95.DataPoint( 700, 22.5e6 )
print("Density =" , point.rho*M_h2o, "kg.m-3" )
# IAPWS95
IAPWS95 is the pure water equation of state (EoS), detailed in Wagner and Pruß (2002). 
This EoS is valid for fluid (vapor, liquid, and supercritical) phases of water, from the melting line and up to 1273 K and 1 GPa.

This program takes as input the temperature (in K) and pressure (in Pa) of pure water and provides several thermodynamic properties of pure water fluid phases.

## Usage:
To use this EoS, import it in a python script and then initiate a DataPoint object:

```python
import IAPWS95

my_point = IAPWS95.DataPoint(my_temperature, my_pressure)
```
This call will generate an object containing the following list of properties:

```python
my_point.T # my_temperature, provided as input
my_point.p # my_pressure, provided as input
my_point.rho # density, mol.m^{-3}
my_point.Cp  # heat capacity at constant pressure J.mol^{-1}.K^{-1}
my_point.s  # entropy J.mol^{-1}.K^{-1}
my_point.u  # internal energy J.mol^{-1}
my_point.dsdp  # the derivative of entropy with pressure at constant volume J.mol^{-1}.K^{-1}.Pa^{-1}
```
On the saturation vapor pressure line, this program returns the properties of the vapor phase.

The call of the IAPWS95 DataPoint object is further illustrated in the test_EOS.py routine, that also contains a benchmark used to test the validity of the EoS

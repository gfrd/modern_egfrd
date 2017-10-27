#!/usr/bin/env python

# --------------------------------------------------------------------------------------------------------------------------------
#
# python -O calcspectrum <power_rec.dat>
# 
# where power_rec.dat is the output file of the eGFRD reaction-recorder
# that should contain 2 reactions (1=association, 2=dissociation)
#
# See: Kaizu K, de Ronde W, Paijmans J, Takahashi K, Tostevin F, ten Wolde PR (2014) The Berg-Purcell Limit Revisited Biophys J, 106:976-985.
#      doi ( https://dx.doi.org/10.1016/j.bpj.2013.12.030 )
#
# Since Python is a very slow languague:
# I recommend to use the calcspectrum.cpp instead (it does the same thing in 1/15 the time)
#
# --------------------------------------------------------------------------------------------------------------------------------

import sys
import time
import math
import cmath
import numpy

# --------------------------------------------------------------------------------------------------------------------------------

def calc_powerspectrum(filename) :

    fourier_size = 1 << 15
    output_size = 1 << 8
    average_size = fourier_size / output_size
    print("Fourier Size: {}, Output Size: {}".format(fourier_size, output_size))

    # Take time
    start = time.clock()

    # pre-process file
    count = 0
    offset = 0
    with open(filename) as fp:
        for line in fp:
            if "Time" in line :
                offset = count + 1; count = 0 ; continue
            count = count + 1
    final_time = float(line.split()[0])
    print("Events: {}, Last: {:06.3f}".format(count, final_time))

    omega_max = 10.0
    omega_min = math.log10(2.0 * math.pi * omega_max / final_time)
    weight = (omega_max - omega_min) / fourier_size

    # Calculate omega table
    print("Omega Min: {:06.3f}, Max: {:06.3f}".format(omega_min, omega_max))
    omega = numpy.empty(fourier_size)
    for i in range(0, fourier_size) :
        omega[i] = math.pow(10, (i + 1.0) * weight + omega_min)

    # Calculate fourier
    fourier = numpy.zeros(fourier_size, dtype=complex)
    with open(filename) as fp:
        
        for i in range(0,offset) :
            fp.readline()

        for i in range(0, count) :
            
            line = fp.readline()
            split = line.split()
            etime = float(split[0])
            flip_state = 1 if int(split[1]) == 1 else -1
            for j in range(0, fourier_size) :
                fourier[j] += cmath.rect(flip_state / omega[j], omega[j] * etime)
    
    # Calculate power-spectrum
    power = numpy.empty(fourier_size)
    for i in range(0, fourier_size) :
        power[i] = (fourier[i].real**2 + fourier[i].imag**2) / final_time

    # Reduce number of points by averaging
    omega_avg = numpy.empty(output_size)
    power_avg = numpy.empty(output_size)
    for i in range(0, output_size) :

        sum = 0;
        for j in range(0, average_size) :
            sum += power[average_size * i + j]
        power_avg[i] = sum / average_size
        omega_avg[i] = omega[average_size * i + average_size / 2]

    # Calculate intrinsic and effective
    D = 1e-12
    r = 5e-9
    ka = 9.16639E-19
    kd = 220.8
    kD = 4 * math.pi * r * 2 * D
    kon = 1 / (1 / ka + 1 / kD)
    koff = 1 / (1 / kd + 1 / kD)
    c = 0.4E-3 * 6.022140857E23
    mu_int = ka * c + kd
    mu_eff = kon * c + koff
    intrinsic = numpy.empty(output_size)
    effective = numpy.empty(output_size)
    for i in range (0, output_size) :
        intrinsic[i] = 0.5 * mu_int / (math.pow(mu_int, 2) + math.pow(omega_avg[i], 2))
        effective[i] = 0.5 * mu_eff / (math.pow(mu_eff, 2) + math.pow(omega_avg[i], 2))

    # Print elapsed calculation time
    elapsed = time.clock()
    print("Calculation Time: {:10.3f} sec".format(elapsed - start))

    # Print output table
    print("\n{:20}\t{:20}\t{:20}\t{:20}".format("omega", "power", "intrinsic", "effective"))
    for i in range(0, output_size) :
        print("{:20.12e}\t{:20.12e}\t{:20.12e}\t{:20.12e}".format(omega_avg[i], power_avg[i], intrinsic[i], effective[i]))

# --------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__' :
    calc_powerspectrum(sys.argv[1])

# --------------------------------------------------------------------------------------------------------------------------------

#!/usr/bin/python
import sys
import numpy as np
from numpy import linalg as LA

def print_header_Xis():
    return '''
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.2.1                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      Xis;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
'''

def print_header_weights():
    return '''
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.2.1                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      weights;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
'''

def print_help():
    print '''
Usage exampels:

1. Using half-range Gasuu-Hermit quaderature points:
    setDV GH 408.1 28
where 408 is the most probable molecular speed, and 28 is 
the number of discrete velocities in each velocity space directions.

OR 

2. Using compound Newton-Cotes rule:
  setDV NC 1600.0 41
where 1600 is max discrete velocity in each velocity space directions,
and 41 is number of directions velocity in each velocity space directions.
NOTE: The number should be like 4*Z+1 where (Z=1,2,3.....). 
'''

def save_file(Xis, weights):
    print("Writting constant/Xis and constant/weights...\n")
    with open('constant/Xis','w') as f_Xi:
        f_Xi.write(print_header_Xis())
        f_Xi.write(str(len(Xis)))
        f_Xi.write("\n(\n")
        for i in Xis:
            f_Xi.write("% 18.15e\n" % i)
        f_Xi.write(");")

    with open('constant/weights','w') as f_weight:
        f_weight.write(print_header_weights())
        f_weight.write(str(len(weights)))
        f_weight.write("\n(\n")
        for i in weights:
            f_weight.write("% 18.15e\n" % i)
        f_weight.write(");")
    print("Writting Done!\n")

def dvGH(C,N2):
    N = N2/2

    a = np.zeros(N)
    b = np.zeros(N)
    a[0] = 1.0/np.sqrt(np.pi)
    a[1] = 2.0/np.sqrt(np.pi)/(np.pi-2.0)
    b[1] = a[0]/( a[0] + a[1])/2.0

    for i in range(2,N):
        b[i] = (i-1)+1.0/2.0-b[i-1]-a[i-1]**2
        a[i] = (i**2/4.0/b[i]-b[i-1]-1.0/2)/a[i-1]-a[i-1]

    J = np.diag(a) + np.diag(np.sqrt(b[1:N]),1) \
       + np.diag(np.sqrt(b[1:N]),-1)

    v,V = LA.eig(J)

    w = V[0,:]*V[0,:]*np.sqrt(np.pi)/2.0

    vw = np.transpose(np.vstack((v,w)))
    vw = vw[vw[:,0].argsort()]
    v = vw[:,0]
    w = vw[:,1]

    Xis = np.hstack((-np.flipud(v),v))
    weights = np.hstack((np.flipud(w),w))
    weights = weights*np.exp(Xis**2)*C
    Xis = Xis*C
    return (Xis, weights)
    
def dvNC(xiMax,N):
    nXi = N
    xiMin = -xiMax
    dv = (xiMax - xiMin)/(nXi-1)
    nBy4 = (nXi-1)/4

    Xis = np.zeros(nXi)
    weights = np.zeros(nXi)

    for i in range(nXi):
        Xis[i] = xiMin + dv*i
         
    for i in range(nBy4):
        weights[4*i+0] = 14.0/90*4
        weights[4*i+1] = 32.0/90*4
        weights[4*i+2] = 12.0/90*4
        weights[4*i+3] = 32.0/90*4

    weights[0] = 7.0/90*4
    weights[nXi-1] = 7.0/90*4

    for i in range(nXi):
        weights[i] = dv*weights[i]

    return (Xis, weights)

if __name__ == '__main__':
    if len(sys.argv) != 4 :
        print_help()
        exit()
    elif sys.argv[1] == 'GH':
        if int(sys.argv[3])%2 != 0 or int(sys.argv[3]) < 7:
            print("ERROR!  Number of discrete velocity should be even, and at least 6.")
            print_help()
            exit()

        Xis, weights = dvGH(float(sys.argv[2]), int(sys.argv[3]))
        save_file(Xis, weights)
    elif sys.argv[1] == 'NC':
        if int(sys.argv[3])%4 != 1:
            print("The last number should be 4*Z+1, where Z = 1,2,3,...\n")
            print_help()
            exit()
        Xis, weights = dvNC(float(sys.argv[2]), int(sys.argv[3]))
        save_file(Xis, weights)
    else:
        print_help()
        exit()

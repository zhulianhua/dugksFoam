#include "pwlInterp2.C"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;
int main(int argc, char* argv[])
{
    int nx, ny;
    nx = atoi(argv[1]);
    ny = nx;
    vector<double> xd(nx);
    vector<double> yd(ny);
    vector<double> zd(nx*ny);

    vector<double> xi(4);
    vector<double> yi(4);

    for(int i=0; i<xd.size(); ++i)
        xd[i] = 2.0*M_PI/((double)xd.size()-1)*i;
    for(int i=0; i<yd.size(); ++i)
        yd[i] = 2.0*M_PI/((double)yd.size()-1)*i;

    for(int j=0; j<yd.size(); ++j)
        for(int i=0; i<xd.size(); ++i)
            zd[j*nx+i] = sin(xd[i])*cos(yd[j]);
    
    xi[0] = 4;
    yi[0] = 2;
    xi[1] = 2.40;
    yi[1] = 0.48;
    xi[2] = 2.8;
    yi[2] = 1.4;
    xi[3] = 4.0;
    yi[3] = 0.98;

    pwlInterp2 itp(xd, yd, zd);
    vector<double> zi = itp.interp(xi,yi);

    //cout << "zi[0] = " << zi[0] << endl;
    //cout << "za[0] = " << sin(xi[0])*cos(yi[0]) <<  endl;
    //
    cout.precision(6);
    cout<<scientific<<endl;

    for(int i =0; i<4; i++)
        cout << xi[i] << "\t" << yi[i] << "\t" <<  zi[i] << "\t" << sin(xi[i])*cos(yi[i]) << endl;;
    return 0;
}

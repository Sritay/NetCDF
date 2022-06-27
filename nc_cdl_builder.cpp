/**********************************************************************
*    This program is to build cdlpart files from a plaintext LAMMPS
*    dump file. The cdlpart files can then be transformed to a single 
*    NetCDF file using standard NetCDF utilities, such as ncgen.
*    To compile: g++ nc_cdl_builder.cpp
*    To run: ./a.out
*
**********************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <ctime>
using namespace std;

int main()
{
    string line, fName;
    int id, n, t, tTime, nAtoms, yStep, tSkip, currentTimeStep, typ;
    double x, y, z, deltaT, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r, vAVG = 0;
    deltaT = 2.0;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    int nTimeSteps = 2; // CHANGE, use command line: grep -o 'TIMESTEP' dump_file | wc -l

    // open text dump file
    ifstream data("dump_file", ios::in);

    ofstream file("vel.dat", ios::out);
    ofstream velo("velocities.cdlpart");
    ofstream ids("id.cdlpart");
    ofstream types("type.cdlpart");
    ofstream positions("positions.cdlpart");
    ofstream stresses("tau.cdlpart");
    ofstream times("time.cdlpart");

    ids << "id = " << endl;
    types << "type = " << endl;
    velo << "velocities = ";
    positions << "coordinates = ";
    stresses << "c_dstress = ";
    times << "time = " << endl;

    clock_t tstart = clock();
    clock_t tend;

    for (t = 0; t < nTimeSteps; t++)
    {
        for (n = 1; n < 10; n++)
        {
            if (n == 4)
            {
                data >> nAtoms;
            }

            if (n == 2)
            {
                data >> currentTimeStep;

                cout << "currentTimeStep = " << currentTimeStep
                     << "; t = " << t << " [ "
                     << 100 * float(t + 1) / float(nTimeSteps)
                     << "% ]" << endl;
                times << "," << currentTimeStep;
            }

            getline(data, line);
        }

        vAVG = 0;

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3 >> tau4 >> tau5 >> tau6;
            ids << "," << id;
            types << "," << typ;
            positions << ",\n"<< x << "," << y << "," << z;
            velo << ",\n"<< vx << "," << vy << "," << vz;
            stresses << ",\n"<< tau1 << "," << tau2 << "," << tau3 << "," << tau4 << "," << tau5 << "," << tau6;
        }

        getline(data, line);

    }

    ids << "; " << endl;
    types << "; " << endl;
    velo << "; " << endl;
    positions << "; " << endl;
    stresses << ";\n } " << endl;
    times << "; " << endl;

    tend = clock();

    cout << "elapsed time: " << double(tend - tstart) / CLOCKS_PER_SEC << "s\n"
         << endl;

    return 0;
}

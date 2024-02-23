#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "PSToolboxBaseEdge.h"
#include "EpanetReader.h"
#include "Connector.h"
#include "SCP.h"

using namespace std;


double getPN10Delta(double D)
{
    int N = 9;
    double Dn[N] = {0.09, 0.11, 0.14, 0.18, 0.20, 0.225, 0.28, 0.315, 0.4};
    double Delta[N] = {5.4e-3, 6.6e-3, 8.3e-3, 10.7e-3, 11.9e-3, 13.4e-3, 16.6e-3, 18.7e-3, 23.7e-3};

    int idx = -1;
    for(int i = 0; i < N; i++)
    {
        if(D < Dn[i])
        {
            idx = i;
            break;
        }
    }
    if(idx == 0) return Delta[0]; //below first
    if(idx == -1) return Delta[N-1]; //above last

    return (Delta[idx-1]*(Dn[idx]-D) + Delta[idx]*(D - Dn[idx-1])) / (Dn[idx+1] - Dn[idx]);
}

void calculatePropagationVelocity(EpanetReader & reader)
{
    //go through all the pipes
    for (int i = 0; i < reader.pipes.size(); i++)
    {
        double r = reader.pipes[i].Roughness;
        double D = reader.pipes[i].Diameter;
        double Ef = 2.18e9; //El. mod. of water
        double Ec, Delta = -1;

        if(r == 100) //KPE PN10 pipes
        {
            Ec = 300.0e6; //El. mod. of pipe
            if(D == 150) Delta = 9.5e-3;
            if(D == 173) Delta = 10.7e-3;
        }
        if(r == 101) //ACPN10 pipes
        {
            Ec = 19613.3e6;
            //hogy jön ki a Delta?
        }
        if(r == 150 && (reader.pipes[i].ID == "23" || reader.pipes[i].ID == "24")) //KPE PN16 pipes
        {
            Ec = 300.0e6; //El. mod. of pipe
            Delta = 50.8e-3;
        }
        if(r == 150) //KPE PN16 pipes
        {
            Ec = 1400.0e6; //El. mod. of pipe
            Delta = getPN10Delta(D);
        }

        if(Delta == -1)
        {
            cout << "Pipe data not found: " <<  reader.pipes[i].ID <<endl;
        }

        double Er = 1.0/(1.0/Ef + D/(Delta * Ec));

        reader.pipes[i].SpeedOfSound = sqrt(Er/1000.0);
        reader.pipes[i].Delta = Delta;

        cout << "Pipe " << i << " ID: " << reader.pipes[i].ID << endl;
        cout << "\tRoughness= " << reader.pipes[i].Roughness << endl;
        cout << "\t       D = " << reader.pipes[i].Diameter << endl;
        cout << "\t   Delta = " << reader.pipes[i].Delta << endl;
        cout << "\t       a = " << reader.pipes[i].SpeedOfSound << endl << endl;
    }
    
}


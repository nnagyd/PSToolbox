#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "PSToolboxBaseEdge.h"
#include "Connector.h"
#include "SCP.h"

using namespace std;

struct JunctionReader
{
    string ID;
    int Elev;
    int Head;
    double Demand;
    int Pattern;
    vector<int> idxPipe;
    vector<bool> end;
    int type; // 0:junction, 1:reservoir
};

struct ReservoirReader
{
    string ID;
    int Head;
    int Pattern;
};

struct TankReader
{
    string ID;
    int Elevation;
    int InitLevel;
    int MinLevel;
    double MaxLevel;
    double Diameter;
    double MinVol;
    string VolCurve;
    string Overflow;
};

struct PipeReader
{
    string ID;
    string Node1;
    string Node2;
    double Length;
    double Diameter;
    double Roughness;
    double Delta;
    double SpeedOfSound;
    double MinorLoss;
};

struct PumpReader
{
    string ID;
    string Node1;
    string Node2;
    string Parameters;
};

struct ValveReader
{
    string ID;
    string Node1;
    string Node2;
    double Diameter;
    string Type;
    int Setting;
    double MinorLoss;
};

class EpanetReader
{
public:
    vector<JunctionReader> junctions;
    vector<ReservoirReader> reservoirs;
    vector<TankReader> tanks;
    vector<PipeReader> pipes;
    vector<PumpReader> pumps;
    vector<ValveReader> valves;

    vector<PSToolboxBaseEdge *> edges;
    vector<Connector *> cons;
    vector<int> con_at_edge_start;
    vector<int> con_at_edge_end;

    void readFromFile(const std::string &filename);
    void convertToRunner();
    void convertToRunner2();
    int nextPipeAtNode(int idx);
    string getOtherNodeOfPipe(int idxPipe, string Node);
    void printEdgesAndCons();
    void unifyJunctions();
    int findNodeByID(const std::string ID);
    int findReservoirByID(const std::string ID);
    bool checkIfIncludedInVector(vector<int> list, int elem);
    vector<int> findConnectingPipes(const std::string node);
    JunctionReader parseJunction(const std::string &line);
    ReservoirReader parseReservoir(const std::string &line);
    TankReader parseTank(const std::string &line);
    PipeReader parsePipe(const std::string &line);
    PumpReader parsePump(const std::string &line);
    ValveReader parseValve(const std::string &line);
};
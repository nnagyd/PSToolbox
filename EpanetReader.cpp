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

void EpanetReader::readFromFile(const std::string &filename)
{
    ifstream file(filename);
    string line;
    string section;

    while (std::getline(file, line))
    {
        // cout << "len=" << line.length() << "\t";
        if (line.find("[JUNCTIONS]") != string::npos)
        {
            section = "JUNCTIONS";
        }
        else if (line.find("[RESERVOIRS]") != string::npos)
        {
            section = "RESERVOIRS";
        }
        else if (line.find("[TANKS]") != string::npos)
        {
            section = "TANKS";
        }
        else if (line.find("[PIPES]") != string::npos)
        {
            section = "PIPES";
        }
        else if (line.find("[PUMPS]") != string::npos)
        {
            section = "PUMPS";
        }
        else if (line.find("[VALVES]") != string::npos)
        {
            section = "VALVES";
            cout << "Valve found!\n";
        }
        else if (line.length() <= 1 || line[0] == ';')
        {
            cout << "Line skipped"
                    << "\n";
            continue; // Skip empty lines or comments
        }
        else
        {
            // cout << "Section = " << section<< "\n";
            if (section == "JUNCTIONS")
            {
                JunctionReader junction = parseJunction(line);
                cout << "Junction ID = " << junction.ID << "\telev = " << junction.Elev << "\tDemand = " << junction.Demand << "\t pattern = " << junction.Pattern << "\n";

                junctions.push_back(junction);
            }
            else if (section == "RESERVOIRS")
            {
                ReservoirReader reservoir = parseReservoir(line);
                cout << "Reservoir ID = " << reservoir.ID << "\thead = " << reservoir.Head << "\tPattern = " << reservoir.Pattern << "\n";

                reservoirs.push_back(reservoir);
            }
            else if (section == "TANKS")
            {
                TankReader tank = parseTank(line);
                cout << "Tank ID = " << tank.ID << "\tCurve = " << tank.VolCurve << "\t Overflow: " << tank.Overflow << "\n";

                tanks.push_back(tank);
            }
            else if (section == "PIPES")
            {
                PipeReader pipe = parsePipe(line);
                if (pipe.ID == "ClosedPipe - exclude")
                {
                    cout << "Closed PipeReader\n";
                }
                else
                {
                    cout << "Pipe ID = " << pipe.ID << "\tNode1 = " << pipe.Node1 << "\t Node2: " << pipe.Node2 << "\n";

                    pipes.push_back(pipe);
                }
            }
            else if (section == "PUMPS")
            {
                PumpReader pump = parsePump(line);
                cout << "Pump ID = " << pump.ID << "\tNode1 = " << pump.Node1 << "\t Node2: " << pump.Node2 << "\tPars: " << pump.Parameters << "\n";

                pumps.push_back(pump);
            }

            else if (section == "VALVES")
            {
                ValveReader valve = parseValve(line);
                cout << "Valve ID = " << valve.ID << "\tNode1 = " << valve.Node1 << "\t Node2: " << valve.Node2 << "\n";

                valves.push_back(valve);
            }
        }
    }
}

void EpanetReader::convertToRunner()
{
    for (int i = 0; i < pipes.size(); i++) // go through all the pipes
    {
        edges.push_back(new SCP(pipes[i].ID, pipes[i].Node1, pipes[i].Node2, 1000, 1300, pipes[i].Length, 0.001 * pipes[i].Diameter, 0.02, 0, 0, false));
    }

    vector<int> openEnds;
    vector<int> connectedEnds;
    vector<int> connectedFronts;
    openEnds.push_back(0);

    // string startingNode = junctions[0].ID;
    // vector<int> connectingEdges = findConnectingEdges(startingNode);
    // cout << "Node connects " << connectingEdges.size() << " edges\n";

    do
    {
        // start of step print
        cout << "\n---------- STEP ------------\nOpen ends: ";
        for (int i = 0; i < openEnds.size(); i++)
        {
            cout << junctions[openEnds[i]].ID << ",";
        }
        cout << "Connected ends: ";
        for (int i = 0; i < connectedEnds.size(); i++)
        {
            cout << junctions[connectedEnds[i]].ID << ",";
        }
        cout << "Connected fronts: ";
        for (int i = 0; i < connectedFronts.size(); i++)
        {
            cout << junctions[connectedFronts[i]].ID << ",";
        }
        cout << "\n\n";

        int actualNodeIdx = openEnds.back();
        string actualNode = junctions[actualNodeIdx].ID;
        vector<int> connectingPipes = findConnectingPipes(actualNode);

        cout << "Connecting pipes to node " << actualNode << " are ";
        for (int i = 0; i < connectingPipes.size(); i++)
        {
            int idx = connectingPipes[i];
            cout << pipes[idx].ID << " (id=" << idx << ") and ";
        }
        cout << "thats it\n";

        if (connectingPipes.size() == 1)
        {
            // create the connection -> boundary - pipe
            int idxPipe = connectingPipes[0];
            bool isFront = !checkIfIncludedInVector(connectedEnds, findNodeByID(getOtherNodeOfPipe(idxPipe, actualNode)));
            double demand = junctions[actualNodeIdx].Demand;
            cons.push_back(new Connector(edges[idxPipe], isFront, "Pressure", 1.e5, demand, true));
            cout << "!!! Connect pipe-BC: " << pipes[idxPipe].ID << "," << isFront << "\n\n";

            // the BC must be applied properly e.g. Tank

            // refresh open and connected end lists
            openEnds.pop_back();
            connectedEnds.push_back(actualNodeIdx);

            // new open ends through the connected pipe
            string otherNode = getOtherNodeOfPipe(idxPipe, actualNode);
            int nextNodeIdx = findNodeByID(otherNode);
            cout << "Next Node: " << otherNode << "\n";

            // check if already included
            if (!checkIfIncludedInVector(connectedEnds, nextNodeIdx) && !checkIfIncludedInVector(openEnds, nextNodeIdx))
            {
                openEnds.push_back(nextNodeIdx);
            }
        }

        if (connectingPipes.size() == 2)
        {
            // create the connection -> boundary - pipe
            cout << "Actual Node Idx = " << actualNodeIdx << "\n";
            cout << "Actual Node = " << actualNode << "\n";

            vector<bool> isFronts;
            int nextNodeIdx;

            for (int i = 0; i < 2; i++)
            {
                int idx = connectingPipes[i];
                bool isFront = !checkIfIncludedInVector(connectedEnds, findNodeByID(getOtherNodeOfPipe(idx, actualNode)));

                isFronts.push_back(isFront);
                cout << "Pipe " << pipes[idx].ID << " (id=" << idx << ") is front = " << isFront << "\n";

                if (isFront)
                {
                    string otherNode = getOtherNodeOfPipe(idx, actualNode);
                    nextNodeIdx = findNodeByID(otherNode);
                    cout << "Next Node: " << otherNode << "\n";
                }
            }

            // check validity
            bool valid = isFronts[0] != isFronts[1];
            cout << "Valid? = " << valid << "\n";
            // what if not valid?
            /*if(!valid)
            {
                openEnds.pop_back();
                continue;
            }*/
            if (!valid)
            {
                isFronts[0] = !isFronts[0];
            }

            // create connection -> pipe - pipe
            double demand = junctions[actualNodeIdx].Demand;
            cons.push_back(new Connector(edges[connectingPipes[0]], isFronts[0], edges[connectingPipes[1]], isFronts[1], demand, true));
            cout << "!!! Connect pipe-pipe: " << pipes[connectingPipes[0]].ID << "," << isFronts[0] << "," << pipes[connectingPipes[1]].ID << "," << isFronts[1] << "\n\n";

            // refresh open and connected end lists
            openEnds.pop_back();
            connectedEnds.push_back(actualNodeIdx);
            if (!checkIfIncludedInVector(connectedEnds, nextNodeIdx) && !checkIfIncludedInVector(openEnds, nextNodeIdx))
            {
                openEnds.push_back(nextNodeIdx);
            }
        }

        if (connectingPipes.size() == 3)
        {
            // create the connection -> pipe - pipe
            cout << "Actual Node Idx = " << actualNodeIdx << "\n";
            cout << "Actual Node = " << actualNode << "\n";

            vector<bool> isFronts;
            vector<int> nextNodeIdx;

            for (int i = 0; i < 3; i++)
            {
                int idx = connectingPipes[i];
                bool isFront = !checkIfIncludedInVector(connectedEnds, findNodeByID(getOtherNodeOfPipe(idx, actualNode)));

                isFronts.push_back(isFront);
                cout << "Pipe " << pipes[idx].ID << " (id=" << idx << ") is front = " << isFront << "\n";

                if (isFront)
                {
                    string otherNode = getOtherNodeOfPipe(idx, actualNode);
                    int nextNodeIdxTmp = findNodeByID(otherNode);
                    cout << "Next Node: " << otherNode << "\n";

                    nextNodeIdx.push_back(nextNodeIdxTmp);
                }
            }

            // create connection -> pipe - pipe
            double demand = junctions[actualNodeIdx].Demand;
            // cons.push_back(new Connector(edges[connectingPipes[0]],isFronts[0],edges[connectingPipes[1]],isFronts[1],demand,true));
            cout << "!!! Connect pipe-pipe-pipe: " << pipes[connectingPipes[0]].ID << "," << isFronts[0] << "," << pipes[connectingPipes[1]].ID << "," << isFronts[1] << "," << pipes[connectingPipes[2]].ID << "," << isFronts[2] << "\n\n";

            // refresh open and connected end lists
            openEnds.pop_back();
            connectedEnds.push_back(actualNodeIdx);
            for (int i = 0; i < nextNodeIdx.size(); i++)
            {
                if (!checkIfIncludedInVector(connectedEnds, nextNodeIdx[i]) && !checkIfIncludedInVector(openEnds, nextNodeIdx[i]))
                {
                    openEnds.push_back(nextNodeIdx[i]);
                }
            }
        }
    } while (openEnds.size() > 0);
}

void EpanetReader::convertToRunner2()
{
    // allow the uniform treatment of junctionsans reservoirs
    unifyJunctions();
    cout << "Junctions unified!\n";

    for (int i = 0; i < pipes.size(); i++) // go through all the pipes
    {
        double a = 1300;
        edges.push_back(new SCP(pipes[i].ID, pipes[i].Node1, pipes[i].Node2, 1000, a, pipes[i].Length, 0.001 * pipes[i].Diameter, 0.02, 0, 0, false));
    }
    cout << "SCP pipes ready!\n";

    /*for(int i = 0; i < junctions.size(); i++) //go through all the nodes to reset
    {
        for(int j = 0; j < 3; j++)
        {
            junctions[i].idxPipe[j] = -1;
            junctions[i].end[j] = false;
        }
    }*/

    cout << "Njunctions = " << junctions.size() << "\n";

    // go through all the pipes and find their respective nodes
    for (int i = 0; i < pipes.size(); i++)
    {
        // find the node index of pipe ends
        int idx1 = findNodeByID(pipes[i].Node1); // idx is start (0)
        int idx2 = findNodeByID(pipes[i].Node2); // idx is end (1)

        junctions[idx1].idxPipe.push_back(i);
        junctions[idx1].end.push_back(true);
        junctions[idx2].idxPipe.push_back(i);
        junctions[idx2].end.push_back(false);

        con_at_edge_start.push_back(idx1);
        con_at_edge_end.push_back(idx2);
    }
    cout << "Pipe end nodes found!\n";

    // create the connectivities
    for (int i = 0; i < junctions.size(); i++)
    {
        int size = junctions[i].idxPipe.size();
        if (size == 1)
        {
            // Connect pipe-BC
            int idx1 = junctions[i].idxPipe[0];
            bool end1 = junctions[i].end[0];

            // check type
            if (junctions[i].type == 0) //"free" end
            {
                double demand = junctions[i].Demand;
                cons.push_back(new Connector(edges[idx1], end1, "Velocity", 0, demand, false));
                cout << "Node " << junctions[i].ID << "\t call Connector(" << pipes[idx1].ID << "," << end1 << ",Velocity,0," << demand << ")\n";
            }
            if (junctions[i].type == 1) // reservoir
            {
                double head = junctions[i].Head;
                double demand = 0;
                cons.push_back(new Connector(edges[idx1], end1, "Pressure", head, demand, false));
                cout << "Node " << junctions[i].ID << "\t call Connector(" << pipes[idx1].ID << "," << end1 << ",Pressure,"<< head << "," << demand << ")\n";
            }
        }

        if (size == 2)
        {
            // Connect pipe-pipe
            int idx1 = junctions[i].idxPipe[0];
            int idx2 = junctions[i].idxPipe[1];
            bool end1 = junctions[i].end[0];
            bool end2 = junctions[i].end[1];
            double demand = junctions[i].Demand;
            cons.push_back(new Connector(edges[idx1], end1, edges[idx2], end2, demand, false));
            cout << "Node " << junctions[i].ID << "\t call Connector(" << pipes[idx1].ID << "," << end1 << "," << pipes[idx2].ID << "," << end2 << "," << demand << ")\n";
        }
        if (size == 3)
        {
            // Connect pipe-pipe-pipe
            int idx1 = junctions[i].idxPipe[0];
            int idx2 = junctions[i].idxPipe[1];
            int idx3 = junctions[i].idxPipe[2];
            bool end1 = junctions[i].end[0];
            bool end2 = junctions[i].end[1];
            bool end3 = junctions[i].end[2];
            double demand = junctions[i].Demand;
            // cons.push_back(new Connector(edges[idx1],end1,edges[idx2],end2,demand,false));
            cout << "Node " << junctions[i].ID << "\t call Connector(" << pipes[idx1].ID << "," << end1 << "," << pipes[idx2].ID << "," << end2 << "," << pipes[idx3].ID << "," << end3  << "," << demand << ")\n";
        }
    }
}

int EpanetReader::nextPipeAtNode(int idx)
{
    if (junctions[idx].idxPipe[0] == -1)
        return 0;
    if (junctions[idx].idxPipe[1] == -1)
        return 1;
    if (junctions[idx].idxPipe[2] == -1)
        return 2;
    return -1;
}

string EpanetReader::getOtherNodeOfPipe(int idxPipe, string Node)
{
    if (pipes[idxPipe].Node1 == Node)
        return pipes[idxPipe].Node2;
    else
        return pipes[idxPipe].Node1;
}

void EpanetReader::printEdgesAndCons()
{
    for (int i = 0; i < edges.size(); i++)
    {
        cout << edges[i]->Info();
    }
}

void EpanetReader::unifyJunctions()
{
    for (int i = 0; i < junctions.size(); i++)
    {
        junctions[i].type = 0;
    }

    for (int i = 0; i < reservoirs.size(); i++)
    {
        cout << "Reservoir " << reservoirs[i].ID << "\n";
        JunctionReader J;
        J.type = 1;
        J.Head = reservoirs[i].Head;
        J.ID = reservoirs[i].ID;
        junctions.push_back(J);
    }

    cout << "Njunctions = " << junctions.size() << "\n";
}

int EpanetReader::findNodeByID(const std::string ID)
{
    for (int i = 0; i < junctions.size(); i++)
    {
        if (junctions[i].ID == ID)
        {
            return i;
        }
    }
    return -1;
}

int EpanetReader::findReservoirByID(const std::string ID)
{
    for (int i = 0; i < reservoirs.size(); i++)
    {
        if (reservoirs[i].ID == ID)
        {
            return i;
        }
    }
    return -1;
}

bool EpanetReader::checkIfIncludedInVector(vector<int> list, int elem)
{
    for (int i = 0; i < list.size(); i++)
    {
        if (list[i] == elem)
            return true;
    }
    return false;
}

vector<int> EpanetReader::findConnectingPipes(const std::string node)
{
    vector<int> connectingPipes;
    for (int i = 0; i < pipes.size(); i++)
    {
        if (node == pipes[i].Node1 || node == pipes[i].Node2)
        {
            // cout << "Found at " << i << "\n";
            connectingPipes.push_back(i);
        }
    }
    return connectingPipes;
}

JunctionReader EpanetReader::parseJunction(const std::string &line)
{
    istringstream iss(line);
    JunctionReader junction;
    iss >> junction.ID >> junction.Elev >> junction.Demand >> junction.Pattern;
    return junction;
}

ReservoirReader EpanetReader::parseReservoir(const std::string &line)
{
    istringstream iss(line);
    ReservoirReader reservoir;
    iss >> reservoir.ID >> reservoir.Head >> reservoir.Pattern;
    return reservoir;
}

TankReader EpanetReader::parseTank(const std::string &line)
{
    istringstream iss(line);
    TankReader tank;
    iss >> tank.ID >> tank.Elevation >> tank.InitLevel >> tank.MinLevel >> tank.MaxLevel >> tank.Diameter >> tank.MinVol >> tank.VolCurve >> tank.Overflow;
    return tank;
}

PipeReader EpanetReader::parsePipe(const std::string &line)
{
    istringstream iss(line);
    PipeReader pipe;
    string Status;
    iss >> pipe.ID >> pipe.Node1 >> pipe.Node2 >> pipe.Length >> pipe.Diameter >> pipe.Roughness >> pipe.MinorLoss >> Status;

    if (Status == "Open")
    {
        return pipe;
    }
    else
    {
        pipe.ID = "ClosedPipe - exclude";
        return pipe;
    }
}

PumpReader EpanetReader::parsePump(const std::string &line)
{
    istringstream iss(line);
    PumpReader pump;
    string tmp;
    iss >> pump.ID >> pump.Node1 >> pump.Node2 >> tmp >> pump.Parameters;
    return pump;
}

ValveReader EpanetReader::parseValve(const std::string &line)
{
    istringstream iss(line);
    ValveReader valve;
    string tmp;
    iss >> valve.ID >> valve.Node1 >> valve.Node2 >> valve.Diameter >> valve.Type >> valve.Setting >> valve.MinorLoss;
    return valve;
}
#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include "EpanetReader.h"
#include "Nyiregyhaza.h"
#include "SCP.h"
#include "Connector.h"
#include "Reservoir.h"
#include "Valve.h"
#include "Valve_with_Absorber.h"
#include "PSToolboxBaseEdge.h"
#include "PSToolboxRunner.h"
#include "my_tools.h"


using namespace std;

int main(int argc, char **argv) {
    EpanetReader reader;
    string location = "test11.inp";
    reader.readFromFile(location);
    cout << "------- File Read Complete ---------\n";
    calculatePropagationVelocity(reader);
    cout << "------- Velocities calculated ---------\n";
    //reader.convertToRunner2();
    cout << "------- Runner ready ---------\n";
    //reader.printEdgesAndCons();
}

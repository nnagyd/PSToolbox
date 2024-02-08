#include <vector>
#include "PSToolboxBaseEdge.h"
#include "Connector.h"

class PSToolboxRunner
{
  public:
    vector<PSToolboxBaseEdge *> edges;
    vector<Connector *> connectors;
    PSToolboxRunner(vector<PSToolboxBaseEdge *>&, vector<Connector *>&);
};

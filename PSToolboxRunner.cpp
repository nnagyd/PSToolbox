#include "PSToolboxRunner.h"
#include "PSToolboxBaseEdge.h"
#include "Connector.h"

PSToolboxRunner::PSToolboxRunner(vector<PSToolboxBaseEdge *>& _e, vector<Connector *>& _c){
  edges = _e;
  connectors = _c;
};

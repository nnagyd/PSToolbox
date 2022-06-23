#pragma once

#include "Reservoir.h"
#include "LWP.h"
#include "SCP.h"
#include "Valve.h"

// @file
//! Class to contain Connectors
/*! blabla
*/
class Connector
{
public:
  Connector();
  ~Connector();
  void Connector_LWP_Reservoir_and_Pipe_Front(double t_target, 
    Reservoir *r, LWP *p, double& pf, double& Tf);  
  
  void Connector_LWP_Pipe_Back_and_Valve(double t_target, 
    LWP *p, Valve *v, double p_downstream, double& pb, double& Tb);

  void Connector_SCP_Reservoir_and_Pipe_Front(double t_target, 
    Reservoir *r, SCP *p, double rho, double a, double& pf);
  bool Connector_SCP_Pipe_Back_and_Valve(double t_target, 
    SCP *p, Valve *v, double rho, double a, double p_downstream, double& pb);

  void Connector_SCP_Pipes(double t_target, 
    SCP* p1, bool is_front1,
    SCP* p2, bool is_front2,
    double mpout,
    double& p, double& v1, double& v2);

  void Connector_SCP_Pipes(double t_target, 
    SCP *p1, bool is_front1,
    SCP *p2, bool is_front2,
    SCP *p3, bool is_front3,
    double mpout,
    double& p, double& v1, double& v2, double& v3);

  void Connector_LWP_Pipes(double t_target, 
    LWP *p1, LWP *p2, 
    double& pL, double& pR, 
    double& TL, double& TR,
    double& rhoL, double& rhoR,
    double& vL, double& vR);

private:
  double Connector_LWP_Reservoir_and_Pipe_Front_fun(double rho, double beta, Reservoir* r1, LWP* p1);
  double signed_sqrt(double x);
};

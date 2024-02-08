#ifndef PSToolboxBaseEdge_H
#define PSToolboxBaseEdge_H
#include <iostream>
#include <vector>
using namespace std;

class PSToolboxBaseEdge
{
  protected:
    string name; 
   unsigned int type;

  public:
    PSToolboxBaseEdge(const string name);
    virtual ~PSToolboxBaseEdge();

    string Get_name(){return name;};
    unsigned int Get_type(){return type;};
    virtual double Get_dt() = 0;
    virtual double Get_tnext()=0;
};
#endif

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
    ~PSToolboxBaseEdge(){};

    string Get_name(){return name;};
    unsigned int Get_type(){return type;};
    virtual double Get_t()=0;
    virtual double Get_dt()=0;
    virtual double Get_tnext()=0;
    virtual void Ini(int)=0;
    virtual void Ini()=0;

    virtual string Info()=0;
    
    virtual void Set_BC_Left(string type, double val)=0;
    virtual void Set_BC_Right(string type, double val)=0;
    virtual void Step()=0;
};
#endif

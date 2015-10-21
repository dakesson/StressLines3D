// Gmsh - Copyright (C) 1997-2012 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#ifndef _ELASTICITY_SOLVER_H_
#define _ELASTICITY_SOLVER_H_

#include <map>
#include <string>
#include "SVector3.h"
#include "dofManager.h"
#include "functionSpace.h"

template <class scalar> class simpleFunction;
class GModel;
class PView;
class groupOfElements;

struct LagrangeMultiplierField {
  int _tag;
  groupOfElements *g;
  double _tau;
  SVector3 _d;
  simpleFunction<double> _f;
  LagrangeMultiplierField() : _tag(0), g(0){}
};

struct elasticField {
  int _tag; // tag for the dofManager
  groupOfElements *g; // support for this field
  double _E, _nu; // specific elastic datas (should be somewhere else)
  elasticField () : _tag(0), g(0){}
};

struct BoundaryCondition
{
  int _tag; // tag for the dofManager
  enum location{UNDEF, ON_VERTEX, ON_EDGE, ON_FACE, ON_VOLUME};
  location onWhat; // on vertices or elements
  groupOfElements *g; // support for this BC
  BoundaryCondition() : _tag(0), onWhat(UNDEF), g(0) {}
};

struct dirichletBC : public BoundaryCondition
{
  int _comp; // component
  simpleFunction<double> *_f;
  dirichletBC ():BoundaryCondition(), _comp(0), _f(0){}
};

struct neumannBC  : public BoundaryCondition
{
  simpleFunction<SVector3> *_f;
  neumannBC () : BoundaryCondition(), _f(NULL){}
};
// an elastic solver ...
class elasticitySolver
{
 public:
  GModel *pModel;
  int _dim, _tag;
  dofManager<double> *pAssembler;
  FunctionSpace<SVector3> *LagSpace;
  FunctionSpace<double> *LagrangeMultiplierSpace;

  // young modulus and poisson coefficient per physical
  std::vector<elasticField> elasticFields; 

  std::vector<LagrangeMultiplierField> LagrangeMultiplierFields;
  // neumann BC
  std::vector<neumannBC> allNeumann;
  // dirichlet BC
  std::vector<dirichletBC> allDirichlet;

 public:
  elasticitySolver(int tag) : _tag(tag), pAssembler(0), LagSpace(0), LagrangeMultiplierSpace(0) {}

  elasticitySolver(GModel *model, int tag);

  void addDirichletBC (int dim, int entityId, int component, double value);
  void addNeumannBC (int dim, int entityId, const std::vector<double> value);
  void addElasticDomain (int tag, double e, double nu);

  virtual ~elasticitySolver()
  {
    if (LagSpace) delete LagSpace;
    if (LagrangeMultiplierSpace) delete LagrangeMultiplierSpace;
    if (pAssembler) delete pAssembler;
  }
  void assemble (linearSystem<double> *lsys);
  void readInputFile(const std::string &meshFileName);
  void read(const std::string s) {readInputFile(s.c_str());}
  virtual void setMesh(const std::string &meshFileName);
  void solve();
  void postSolve();
  void exportKb();
  void getSolutionOnElement(MElement *el, fullMatrix<double> &sol);
  virtual PView *buildDisplacementView(const std::string postFileName);
  virtual PView *buildStressesView(const std::string postFileName);
  virtual PView *buildLagrangeMultiplierView(const std::string posFileName);
  virtual PView *buildElasticEnergyView(const std::string postFileName);
  virtual PView *buildVonMisesView(const std::string postFileName);
  // std::pair<PView *, PView*> buildErrorEstimateView
  //   (const std::string &errorFileName, double, int);
  // std::pair<PView *, PView*> buildErrorEstimateView
  //   (const std::string &errorFileName, const elasticityData &ref, double, int);
};

#endif

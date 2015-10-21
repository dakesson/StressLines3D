// Gmsh - Copyright (C) 1997-2012 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#ifndef _MESH_GFACE_DELAUNAY_INSERTIONFACE_H_
#define _MESH_GFACE_DELAUNAY_INSERTIONFACE_H_

#include "MTriangle.h"
#include "MQuadrangle.h"
#include "STensor3.h"
#include <list>
#include <set>
#include <map>

class GModel;
class GFace;
class BDS_Mesh;
class BDS_Point;

void buildMetric(GFace *gf, double *uv, double *metric);
int inCircumCircleAniso(GFace *gf, double *p1, double *p2, double *p3, 
                        double *p4, double *metric);
int inCircumCircleAniso(GFace *gf, MTriangle *base, const double *uv, 
                        const double *metric, const std::vector<double> &Us,
                        const std::vector<double> &Vs);
void circumCenterMetric(double *pa, double *pb, double *pc, const double *metric,
                        double *x, double &Radius2);
void circumCenterMetric(MTriangle *base, const double *metric,
                        const std::vector<double> &Us, 
                        const std::vector<double> &Vs,
                        double *x, double &Radius2);
bool circumCenterMetricInTriangle(MTriangle *base, const double *metric,
                                  const std::vector<double> &Us,
                                  const std::vector<double> &Vs);
bool invMapUV(MTriangle *t, double *p,
              const std::vector<double> &Us, const std::vector<double> &Vs,
              double *uv, double tol);

class MTri3
{
 protected :
  bool deleted;
  double circum_radius;
  MTriangle *base;
  MTri3 *neigh[3];

 public :
  static int radiusNorm; // 2 is euclidian norm, -1 is infinite norm  
  bool isDeleted() const { return deleted; }
  void forceRadius(double r) { circum_radius = r; }
  double getRadius() const { return circum_radius; }

  MTri3(MTriangle *t, double lc, SMetric3 *m = 0, const std::vector<double> *Us = 0, const std::vector<double> *Vs = 0, GFace *gf = 0);
  inline MTriangle *tri() const { return base; }
  inline void  setNeigh(int iN , MTri3 *n) { neigh[iN] = n; }
  inline MTri3 *getNeigh(int iN ) const { return neigh[iN]; }
  int inCircumCircle(const double *p) const;
  inline int inCircumCircle(double x, double y) const
  {
    const double p[2] = {x, y};
    return inCircumCircle(p);
  }
  inline int inCircumCircle(const MVertex * v) const
  {
    return inCircumCircle(v->x(), v->y());
  }
  inline void setDeleted(bool d){ deleted = d; }
  inline bool assertNeigh() const
  {
    if(deleted) return true;
    for(int i = 0; i < 3; i++)
      if(neigh[i] && (neigh[i]->isNeigh(this) == false)) return false;
    return true;
  }
  inline bool isNeigh(const MTri3 *t) const
  {
    for(int i = 0; i < 3; i++)
      if(neigh[i] == t) return true;
    return false;
  }
};

/*        2
  3 +------------+ 2
    |            |
    |            |
  3 |            | 1
    |            |
    |            |
    +------------+
  0       0        1 
	  
   We require that quads are oriented
   We'd like to walk into the quad mesh and 
   create sheets
   
   We start from a given quad and one edge
   We give a path as an array of int's
   If path[i] == 1 --> 
      
 */

class MQua4
{
 protected :
  MQuadrangle *base;
  MQua4 *neigh[4];
 public :
 MQua4(MQuadrangle *q) : base(q) {
    neigh[0] = neigh[1] = neigh[2] = neigh[3] = NULL;}
  inline MQuadrangle *qua() const { return base; }
  inline void  setNeigh(int iN , MQua4 *n) { neigh[iN] = n; }
  inline MQua4 *getNeigh(int iN ) const { return neigh[iN]; }
  inline int getEdge(MVertex *v1, MVertex *v2){
    for (int i=0;i<4;i++){
      MEdge e = base->getEdge(i);
      if (e.getVertex(0) == v1 && e.getVertex(1) == v2)return i;
      if (e.getVertex(0) == v2 && e.getVertex(1) == v1)return i;      
    }
    return -1;
  }
};


class compareTri3Ptr
{
 public:
  inline bool operator () (const MTri3 *a, const MTri3 *b)  const
  {
    if(a->getRadius() > b->getRadius()) return true;
    if(a->getRadius() < b->getRadius()) return false;
    return a<b;
  }
};

void connectQuads(std::vector<MQua4*> &);
void connectTriangles(std::list<MTri3*> &);
void connectTriangles(std::vector<MTri3*> &);
void connectTriangles(std::set<MTri3*,compareTri3Ptr> &AllTris);
void bowyerWatson(GFace *gf, int MAXPNT= 1000000000);
void bowyerWatsonFrontal(GFace *gf);
void bowyerWatsonFrontalLayers(GFace *gf, bool quad);
void buildBackGroundMesh (GFace *gf);

struct edgeXface
{
  MVertex *v[2];
  MTri3 * t1;
  int i1;
  edgeXface(MTri3 *_t, int iFac) : t1(_t), i1(iFac)
  {
    v[0] = t1->tri()->getVertex(iFac == 0 ? 2 : iFac-1);
    v[1] = t1->tri()->getVertex(iFac);
    std::sort(v, v + 2);
  }
  inline bool operator < ( const edgeXface &other) const
  {
    if(v[0] < other.v[0]) return true;
    if(v[0] > other.v[0]) return false;
    if(v[1] < other.v[1]) return true;
    return false;
  }
};

struct edgeXquad
{
  MVertex *v[2];
  MQua4 * t1;
  int i1;
  edgeXquad(MQua4 *_t, int iFac) : t1(_t), i1(iFac)
  {
    v[0] = t1->qua()->getVertex(iFac);
    v[1] = t1->qua()->getVertex((iFac+1)%4);
    std::sort(v, v + 2);
  }
  inline bool operator < ( const edgeXquad &other) const
  {
    if(v[0] < other.v[0]) return true;
    if(v[0] > other.v[0]) return false;
    if(v[1] < other.v[1]) return true;
    return false;
  }
};


#endif

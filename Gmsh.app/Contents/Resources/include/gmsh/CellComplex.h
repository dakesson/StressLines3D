// Gmsh - Copyright (C) 1997-2012 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.
//
// Contributed by Matti Pellikka <matti.pellikka@tut.fi>.

#ifndef _CELLCOMPLEX_H_
#define _CELLCOMPLEX_H_

#include <map>
#include <string.h>
#include <set>
#include <algorithm>
#include <queue>
#include <string>
#include "Cell.h"
#include "MElement.h"
#include "GModel.h"

class Cell;
class BdInfo;

class CellComplex
{
 private:

  GModel* _model;

  // sorted containers of unique cells in this cell complex
  // one for each dimension
  std::set<Cell*, Less_Cell> _cells[4];

  // original cells of this cell complex
  std::set<Cell*, Less_Cell> _ocells[4];

  // new cells created during reductions
  std::vector<Cell*> _newcells;
  std::vector<Cell*> _removedcells;

  int _dim;
  bool _simplicial;
  bool _saveorig;

  int _deleteCount;

  bool _reduced;

  // for constructor
  bool _insertCells(std::vector<MElement*>& elements, int domain);
  bool _removeCells(std::vector<MElement*>& elements, int domain);

  bool _immunizeCells(std::vector<MElement*>& elements);

  // enqueue cells in queue if they are not there already
  void enqueueCells(std::map<Cell*, short int, Less_Cell>& cells,
		    std::queue<Cell*>& Q, std::set<Cell*, Less_Cell>& Qset);

  // insert/remove a cell from this cell complex
  void removeCell(Cell* cell, bool other=true);
  void insertCell(Cell* cell);

  // queued coreduction
  int coreduction(Cell* startCell, bool omit,
		  std::vector<Cell*>& omittedCells);

 public:
  CellComplex(GModel* model,
	      std::vector<MElement*>& domainElements,
	      std::vector<MElement*>& subdomainElements,
              std::vector<MElement*>& nondomainElements,
              std::vector<MElement*>& nonsubdomainElements,
              std::vector<MElement*>& immuneElements,
              bool saveOriginalComplex=true);
  ~CellComplex();


  GModel* getModel() const { return _model; }
  int getDim() { return _dim; }
  bool simplicial() { return _simplicial; }

  // get the number of certain dimensional cells
  int getSize(int dim, bool orig=false){
    if(!orig) return _cells[dim].size();
    else return _ocells[dim].size(); }

  // get dim-dimensional cells
  // domain = 0: cells in domain relative to subdomain
  // domain = 1: cells in domain
  // domain = 2: cells in subdomain
  void getCells(std::set<Cell*, Less_Cell>& cells, int dim, int domain=0);
  //std::set<Cell*, Less_Cell> getOrigCells(int dim){ return _ocells[dim]; }

  // iterator for the cells of same dimension
  typedef std::set<Cell*, Less_Cell>::iterator citer;

  // iterators to the first and last cells of certain dimension
  citer firstCell(int dim, bool orig=false) {
    return orig ? _ocells[dim].begin() : _cells[dim].begin(); }
  citer lastCell(int dim, bool orig=false) {
    return orig ? _ocells[dim].end() : _cells[dim].end(); }

  // true if cell complex has given cell
  bool hasCell(Cell* cell, bool orig=false);

  // check whether two cells both belong to subdomain or if neither one does
  bool inSameDomain(Cell* c1, Cell* c2) const {
    return (c1->getDomain() ==  c2->getDomain()); }

  // remove cells in subdomain from this cell complex
  void removeSubdomain();

  // (co)reduction of this cell complex
  // removes (co)reduction pairs of cell of dimension dim and dim-1
  int reduction(int dim, bool omit, std::vector<Cell*>& omittedCells);
  int coreduction(int dim, bool omit, std::vector<Cell*>& omittedCells);

  // Cell combining for reduction and coreduction
  int combine(int dim);
  int cocombine(int dim);

  // check whether all boundary cells of a cell has this cell
  // as coboundary cell and vice versa
  // also check whether all (co)boundary cells of a cell
  // belong to this cell complex
  bool coherent();

  // full (co)reduction of this cell complex (all dimensions, with combining)
  // (with highest dimensional cell omitting?)
  int reduceComplex(bool docombine=true, bool omit=true);
  int coreduceComplex(bool docombine=true, bool omit=true);

  bool isReduced() const { return _reduced; }

  int eulerCharacteristic() {
    return getSize(0) - getSize(1) + getSize(2) - getSize(3); }
  void printEuler() {
    printf("Euler characteristic: %d. \n", eulerCharacteristic()); }

  // restore the cell complex to its original state before (co)reduction
  bool restoreComplex();

  // print the vertices of cells of certain dimension
  void printComplex(int dim);


  // experimental
  int saveComplex(std::string filename);
  int loadComplex(std::string filename);
};

#endif

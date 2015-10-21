// Gmsh - Copyright (C) 1997-2012 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#ifndef _GMSH_H_
#define _GMSH_H_

#include <string>
#include "GmshMessage.h"

int GmshInitialize(int argc=0, char **argv=0);
int GmshSetMessageHandler(GmshMessage *callback);
int GmshSetBoundingBox(double xmin, double xmax,
                       double ymin, double ymax, 
                       double zmin, double zmax);
int GmshSetOption(std::string category, std::string name, std::string value, int index=0);
int GmshSetOption(std::string category, std::string name, double value, int index=0);
int GmshSetOption(std::string category, std::string name, unsigned int value, int index=0);
int GmshGetOption(std::string category, std::string name, std::string &value, int index=0);
int GmshGetOption(std::string category, std::string name, double &value, int index=0);
int GmshGetOption(std::string category, std::string name, unsigned int &value, int index=0);
int GmshMergeFile(std::string fileName);
int GmshWriteFile(std::string fileName);
int GmshFinalize();
int GmshBatch();
int GmshFLTK(int argc=0, char **argv=0);

#endif

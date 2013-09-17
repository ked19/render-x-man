#ifndef _RENDER_WIRE_H
#define _RENDER_WIRE_H

#include "TRIModel.h"
#include "define.h"
#include "myMath.h"

#include <windows.h>
#include <vector>
//#include <GL/gl.h>
#include <GL/glew.h>
//#include "glee.h"

using namespace std;

typedef struct 
{
	DATA aP0[3];
	DATA aP1[3];
} Line;

void RenderWire(const TRIModel &canModel, DATA clipR);

#endif
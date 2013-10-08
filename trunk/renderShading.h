#ifndef _RENDER_SHADING_H
#define _RENDER_SHADING_H

#include "TRIModel.h"
#include "shader.h"
#include "myMath.h"

#include <windows.h>
//#include <GL/gl.h>
#include <GL/glew.h>
//#include "glee.h"

void RenderShading(const TRIModel &worldModel);

typedef enum {WIRE, FLAT_SH, GOURAUD_SH, PHONG_SH, VOLUME} RENDER_TYPE;
extern RENDER_TYPE eRenderMethod;

extern Shader *pFlatSh;
extern Shader *pGouraudSh;
extern Shader *pPhongSh;

#endif
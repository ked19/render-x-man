#include "renderShading.h"

RENDER_TYPE eRenderMethod = FLAT_SH;

extern DATA viewScl;

Shader *pFlatSh = 0;
Shader *pGouraudSh = 0;
Shader *pPhongSh = 0;

void RenderShading(const TRIModel &worldModel)
{
	float fVs = (float)viewScl;
	float aaLightPos[][4] = {{ 1.f/fVs,  1.f/fVs,  -1.f/fVs,  1.f},		
							 {-1.f/fVs, -1.f/fVs,  -1.f/fVs,  1.f},		
							 {-1.f/fVs, -1.f/fVs,  1.f/fVs,  1.f}};
	float aaLightAmb[][4] = {{0.01f, 0.01f, 0.01f, 1.f},		
							 {0.01f, 0.01f, 0.01f, 1.f},		
							 {0.01f, 0.01f, 0.01f, 1.f}};
	float aaLightDif[][4] = {{0.6f, 0.6f, 0.6f, 1.f},		
							 {0.6f, 0.6f, 0.6f, 1.f},		
							 {0.6f, 0.6f, 0.6f, 1.f}};
	float aaLightSpc[][4] = {{0.2f, 0.2f, 0.2f, 1.f},		
							 {0.2f, 0.2f, 0.2f, 1.f},		
							 {0.2f, 0.2f, 0.2f, 1.f}};

	Shader *pShNow = 0;
	if(eRenderMethod == FLAT_SH)
	{
		pShNow = pFlatSh;
	}
	else if(eRenderMethod == GOURAUD_SH)
	{
		pShNow = pGouraudSh;
	}
	else if(eRenderMethod == PHONG_SH)
	{
		pShNow = pPhongSh;
	}
	else {}

	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);

	glShadeModel(GL_SMOOTH);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);
	for(unsigned i=0; i<3; i++)
	{
		glLightfv(GL_LIGHT0+i, GL_POSITION, aaLightPos[i]);
		glLightfv(GL_LIGHT0+i, GL_AMBIENT,  aaLightAmb[i]);
		glLightfv(GL_LIGHT0+i, GL_DIFFUSE,  aaLightDif[i]);
		glLightfv(GL_LIGHT0+i, GL_SPECULAR, aaLightSpc[i]);
		glLightf(GL_LIGHT0+i, GL_SPOT_EXPONENT, 20.f);
		glEnable(GL_LIGHT0+i);
	}
	glPopMatrix();

	glUniform3f(pShNow->GetUniLoc("v3Eye"), 0, 0, 0);
	glUniform1i(pShNow->GetUniLoc("lightNum"), 3);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	const vector<Triangle> &vTri = worldModel.triangleList;
	for(unsigned i=0; i<vTri.size(); i++)
	{
		
		const int *pFC = vTri[i].foreColor;
		//const int *pBC = vTri[i].backColor;		
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		glColor3ub(pFC[0], pFC[1], pFC[2]);
		
		double aNAvg[3] = {0};
		for(unsigned v=0; v<3; v++)
		{
			myMath.Add3V(aNAvg, vTri[i].normal[v]);
		}
		myMath.Normal3V(aNAvg);

		double aVAvg[3] = {0};
		for(unsigned v=0; v<3; v++)
		{
			myMath.Add3V(aVAvg, vTri[i].vertex[v]);
		}
		myMath.Mul3V(aVAvg, 1.F/3.F);

		if(eRenderMethod == FLAT_SH)
		{
			glUniform3f(pShNow->GetUniLoc("v3Pos"), (float)aVAvg[0], (float)aVAvg[1], (float)aVAvg[2]);
			glUniform3f(pShNow->GetUniLoc("v3Norm"), (float)aNAvg[0], (float)aNAvg[1], (float)aNAvg[2]);
		}
		else {}
		
		glBegin(GL_TRIANGLES);
			glNormal3dv(vTri[i].normal[0]);		
			glVertex3dv(vTri[i].vertex[0]);
			glNormal3dv(vTri[i].normal[1]);		
			glVertex3dv(vTri[i].vertex[1]);
			glNormal3dv(vTri[i].normal[2]);		
			glVertex3dv(vTri[i].vertex[2]);
		glEnd();
	}

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);
}
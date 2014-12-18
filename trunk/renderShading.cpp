#include "renderShading.h"

RENDER_TYPE eRenderMethod = FLAT_SH;

extern DATA viewScl;

Shader *pFlatSh = 0;
Shader *pGouraudSh = 0;
Shader *pPhongSh = 0;

void RenderShading(const TRIModel &worldModel, unsigned tfNo)
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

	const vector<Triangle> &vTri = worldModel.triangleList;
	//*******************************************
	// texture mapping
	//*******************************************
	static bool bTx = false;
	static Tex2D *pTxF = 0;
	static float aTf[8];
	static float aTf0[8];
	static float aTf1[8];
	if (!bTx) {
		bTx = true;

		MyImg *pImg = imgIO.Read("amanda_peet.jpg");
		Vect3D<unsigned> dim = pImg->GetDim();
		static float *pF = new float[dim.m_x * dim.m_y * dim.m_z];
		pImg->CopyTo_zFirst(pF);
		for (unsigned i = 0; i < dim.m_x * dim.m_y * dim.m_z; i++) {
			pF[i] /= 255.F;
		}
		pTxF = new Tex2D(pF, dim.m_x, dim.m_y, dim.m_z);

		//ifstream ifs("amanda_peet_tf.txt");
		ifstream ifs("amanda_peet_tf0.txt");
		for (unsigned i = 0; i < 8; i++) {
			ifs >> aTf[i];
			cout << aTf[i] << " ";
			if (i <= 3) {
				//aTf[i] /= (dim.m_x - 1);
			} else {
				//aTf[i] /= (dim.m_y - 1);
			}
		}
		//getchar();
		ifs.close();

		ifstream ifs0("amanda_peet_tf0.txt");
		for (unsigned i = 0; i < 8; i++) {
			ifs0 >> aTf0[i];
			cout << aTf0[i] << " ";
			if (i <= 3) {
				//aTf[i] /= (dim.m_x - 1);
			} else {
				//aTf[i] /= (dim.m_y - 1);
			}
		}
		//getchar();
		ifs0.close();

		//ifstream ifs1("amanda_peet_tf1.txt");
		ifstream ifs1("amanda_peet_tf0.txt");
		for (unsigned i = 0; i < 8; i++) {
			ifs1 >> aTf1[i];
			cout << aTf1[i] << " ";
			if (i <= 3) {
				//aTf[i] /= (dim.m_x - 1);
			} else {
				//aTf[i] /= (dim.m_y - 1);
			}
		}
		//getchar();
		ifs1.close();

		/*
		for (unsigned i = 0; i < vTri.size(); i++) {
			for (unsigned v = 0; v < 3; v++) {
				DATA xx = vTri[i].vertex[v][0] * aTf[0] + vTri[i].vertex[v][1] * aTf[1]  + vTri[i].vertex[v][2] * aTf[2] + aTf[3]; 
				DATA yy = vTri[i].vertex[v][0] * aTf[4] + vTri[i].vertex[v][1] * aTf[5]  + vTri[i].vertex[v][2] * aTf[6] + aTf[7]; 
				cout << xx << " " << yy << " " << vTri[i].vertex[v][0] << " " << vTri[i].vertex[v][1] << endl;
			}
		}
		*/
	} else {}
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, pTxF->GetTexID());
	glUniform1i(pShNow->GetUniLoc("smpFace"), 0);

	glUniform3f(pShNow->GetUniLoc("v3Cnt"), -1.0119, -1.0119, -1.0119);
	glUniform1i(pShNow->GetUniLoc("tfNo"), tfNo);
	glUniform4f(pShNow->GetUniLoc("v4TfX"), aTf[0], aTf[1], aTf[2], aTf[3]);
	glUniform4f(pShNow->GetUniLoc("v4TfY"), aTf[4], aTf[5], aTf[6], aTf[7]);
	glUniform4f(pShNow->GetUniLoc("v4Tf0X"), aTf0[0], aTf0[1], aTf0[2], aTf0[3]);
	glUniform4f(pShNow->GetUniLoc("v4Tf0Y"), aTf0[4], aTf0[5], aTf0[6], aTf0[7]);
	glUniform4f(pShNow->GetUniLoc("v4Tf1X"), aTf1[0], aTf1[1], aTf1[2], aTf1[3]);
	glUniform4f(pShNow->GetUniLoc("v4Tf1Y"), aTf1[4], aTf1[5], aTf1[6], aTf1[7]);

	glUniform3f(pShNow->GetUniLoc("v3Eye"), 0, 0, 0);
	glUniform1i(pShNow->GetUniLoc("lightNum"), 3);

	float sym = (bSymetry) ? 1.f : 0;
	glUniform1f(pShNow->GetUniLoc("symetry"), sym), 

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//const vector<Triangle> &vTri = worldModel.triangleList;
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
#include <windows.h>
#include <GL/glew.h>

#include "shader.h"
#include "myWinApp.h"

TRIModel *pTriModel = 0;
TRIModel canonicalModel;

MyVolume *pVolModel = 0;

DATA clipRatio = 0.8F;

DATA aObjCnt[3];
DATA aaObjVec[][4] = {{1.F, 0, 0, 0}, {0, 1.F, 0, 0}, {0, 0, 1.F, 0}};
DATA aViewCnt[] = {0, 0, -1.5F};
DATA viewScl;
DATA viewWide = 0.5F;
DATA steDis = 0.1F; //0.15F; //0.06F; //0.02F;

DATA trsStep = 0.001F;
DATA rotStep = 0.5F * D2R;
DATA zoomStep = 0.01F;
DATA vWideStep = 0.001F;

Vect2D<int> trsAnchor(0, 0);
Vect2D<int> rotAnchor(0, 0);
Vect2D<int> sclAnchor(0, 0);
DATA aEyePos[] = {0, 0, 0};
DATA aaEyeVec[][4] = {{1.F, 0, 0, 0}, {0, 1.F, 0, 0}, {0, 0, 1.F, 0}};

bool bEyes = false;
bool bParal = false;
bool bCoord = false;

DATA disAtt = 1.F;
DATA briScl = 1.F;
DATA smpNum = 500;
DATA opaScl = 0.2F;
DATA mskThrd = 0; //0.5F;

bool bMIP = false;
DATA occScl = 0.2;

bool bXClip = false;	bool bXFront = false;	DATA xPlane = 0;
bool bYClip = false;	bool bYFront = false;	DATA yPlane = 0;
bool bZClip = false;	bool bZFront = false;	DATA zPlane = 0;

bool bStereo = false;
bool bFusion = false;

void Unstereo();
Vect2D<int> mainWinPos(0, 0);	Vect2D<int> mainWinLen(0, 0);
Vect2D<int> glWinPos(0, 0);		Vect2D<int> glWinLen(0, 0);

bool bContour = false;
float *pFrame = 0;
Tex2D *pTexFrame = 0;	

extern bool bMap;

//*************************************************************************************************

void GetContour(Tex2D* &pTexContour)
{
	static Layer *pLyrFrame = 0;
	static Mtx *pMtxGray = 0;
	static Mtx *apMtxGauss[] = {0, 0};

	Vect2D<unsigned> dimFrame = pTexContour->GetDim(); 

	if(!pLyrFrame ||
		pLyrFrame->GetDim().m_x<dimFrame.m_x || 
		pLyrFrame->GetDim().m_y<dimFrame.m_y)
	{
		delete pLyrFrame;	pLyrFrame = new Layer(dimFrame.m_x, dimFrame.m_y, 3);
		delete pMtxGray;	pMtxGray = new Mtx(dimFrame.m_x, dimFrame.m_y);	
	}
	else {}
	if(!apMtxGauss[0])	apMtxGauss[0] = new Mtx(20, 20);	else {}
	if(!apMtxGauss[1])	apMtxGauss[1] = new Mtx(20, 20);	else {}

	pTexContour->GetCell(*pLyrFrame);

	lyrOp.gray.Gen(*pMtxGray, *pLyrFrame);
	
	DATA coeSharp = 1.F;
	mtxOp.DoG.GenNPR(*pMtxGray, apMtxGauss, 2.F, coeSharp);
	mtxOp.sub.Gen(1.F, *pMtxGray);

	for(unsigned y=0; y<dimFrame.m_y; y++)
	{
		for(unsigned x=0; x<dimFrame.m_x; x++)
		{
			float val = (float)pMtxGray->CellVal(x, y);
			pTexContour->SetCell(val, x, y, 0);
			pTexContour->SetCell(val, x, y, 1);
			pTexContour->SetCell(val, x, y, 2);
		}
	}
	glEnable(GL_TEXTURE_2D);
	pTexContour->Update();
}

//*************************************************************************************************

MyGlWindow::MyGlWindow(int x, int y, int w, int h, const char *l)
:Fl_Gl_Window(x, y, w, h, l)
{}

int MyGlWindow::handle(int evnt)
{	
	if(evnt == FL_PUSH)
	{
		int buttn = Fl::event_button();
		if(buttn == 1)
		{
			trsAnchor.m_x = Fl::event_x();
			trsAnchor.m_y = h()-1 - Fl::event_y();
			redraw();
		}
		else if(buttn == 3)
		{
			int xIn = Fl::event_x();
			int yIn = h()-1 - Fl::event_y();
			rotAnchor.m_x = xIn;
			rotAnchor.m_y = yIn;
			redraw();
		}
		else {}
		return 1;
	}
	else if(evnt == FL_DRAG)
	{
		int buttn = Fl::event_button();
		if(buttn == 1)
		{
			int xIn = Fl::event_x();
			int yIn = h()-1 - Fl::event_y();
			
			DATA aX[3], aY[3];
			for(unsigned i=0; i<3; i++)
			{
				aX[i] = -aaEyeVec[0][i] * (xIn - trsAnchor.m_x) * trsStep;
				aY[i] = -aaEyeVec[1][i] * (yIn - trsAnchor.m_y) * trsStep;
			}
			if(bEyes)
			{
				myMath.Add3V(aEyePos, aX);
				myMath.Add3V(aEyePos, aY);
			}
			else
			{
				myMath.Sub3V(aViewCnt, aX);
				myMath.Sub3V(aViewCnt, aY);
			}

			trsAnchor.m_x = xIn;
			trsAnchor.m_y = yIn;
			redraw();
		}
		else if(buttn == 3)
		{
			int xIn = Fl::event_x();
			int yIn = h()-1 - Fl::event_y();

			DATA aRV[3];
			int xDiff = xIn - rotAnchor.m_x;
			int yDiff = yIn - rotAnchor.m_y;
			for(unsigned i=0; i<3; i++)
			{
				aRV[i] = xDiff*aaEyeVec[0][i] + yDiff*aaEyeVec[1][i]; 
			}
			DATA aRAxis[3];
			myMath.Cross3V(aRAxis, aRV, aaEyeVec[2]);
			myMath.Normal3V(aRAxis);

			DATA aDiff[] = {xDiff, yDiff};
			DATA ang = myMath.Len2V(aDiff) * rotStep;

			if(bEyes)
			{
				for(unsigned i=0; i<3; i++)
				{
					DATA aVOld[4];
					myMath.Copy4V(aVOld, aaEyeVec[i]);
					myMath.RotateVec(aaEyeVec[i], aVOld, ang, aRAxis);
				}
			}
			else
			{
				for(unsigned i=0; i<3; i++)
				{
					DATA aVOld[4];
					myMath.Copy4V(aVOld, aaObjVec[i]);
					myMath.RotateVec(aaObjVec[i], aVOld, ang, aRAxis);
				}
			}

			rotAnchor.m_x = xIn;
			rotAnchor.m_y = yIn;
			redraw();
		}
		return 1;
	}
	else if(evnt == FL_RELEASE)
	{}
	else if(evnt == FL_MOUSEWHEEL)
	{
		DATA scrollY = Fl::event_dy() * zoomStep;
		DATA sY2 = Fl::event_dy() * vWideStep;
		DATA aDiff[3];
		for(unsigned i=0; i<3; i++)
		{
			aDiff[i] =  scrollY * aaEyeVec[2][i];
		}
		if(bEyes)
		{
			myMath.Add3V(aEyePos, aDiff);
			viewWide -= sY2;
		}
		else
		{
			myMath.Sub3V(aViewCnt, aDiff);
			viewWide += sY2;
		}
		redraw();
	}
	else if(evnt == FL_FOCUS)
	{
		return 1;
	}
	else if(evnt == FL_KEYDOWN)
	{
		char key = Fl::event_key();
		if(key == 'n')
		{
			Unstereo();
			redraw();
			return 1;
		}
		else if(key == 'v')
		{
			DATA aVec[12];
			for(unsigned j=0; j<3; j++)
			{
				for(unsigned i=0; i<4; i++)
				{
					aVec[j*4+i] = aaObjVec[j][i];
				}
			}
			ofstream out("out.bin", ios::binary);
			out.write((char*)aVec, sizeof(DATA)*12);
			out.close();

			redraw();
			return 1;
		}
		else if(key == 'r')
		{
			DATA aVec[12];
			ifstream inF("out.bin", ios::binary);
			inF.read((char*)aVec, sizeof(DATA)*12);
			inF.close();

			for(unsigned j=0; j<3; j++)
			{
				for(unsigned i=0; i<4; i++)
				{
					aaObjVec[j][i] = aVec[j*4+i];
				}
			}

			redraw();
			return 1;
		}
		else {}
		return 0;
	}
	else
	{
		return Fl_Widget::handle(evnt);
	}
	return 0;
}

//*************************************************************************************************

void RotateObj()
{
	DATA aBasis[] = {aaObjVec[0][0], aaObjVec[0][1], aaObjVec[0][2], aaObjVec[0][3],
					 aaObjVec[1][0], aaObjVec[1][1], aaObjVec[1][2], aaObjVec[1][3],
					 aaObjVec[2][0], aaObjVec[2][1], aaObjVec[2][2], aaObjVec[2][3],
					 0, 0, 0, 1.F};
	glMultMatrixd(aBasis);
}

void CreatePerspectivePMtx(DATA aaPM[][4], DATA left, DATA right, DATA bottom, DATA top, DATA nearZ, DATA farZ)
{
	DATA xDiff = right - left;
	DATA yDiff = top - bottom;
	DATA zDiff = farZ - nearZ;
	DATA divXD = 1.F / xDiff;
	DATA divYD = 1.F / yDiff;
	DATA divZD = 1.F / zDiff;

	aaPM[0][0] = 2.F*nearZ*divXD;	aaPM[0][1] = 0;					aaPM[0][2] = (right+left)*divXD;		aaPM[0][3] = 0;
	aaPM[1][0] = 0;					aaPM[1][1] = 2.F*nearZ*divYD;	aaPM[1][2] = (top+bottom)*divYD;		aaPM[1][3] = 0;
	aaPM[2][0] = 0;					aaPM[2][1] = 0;					aaPM[2][2] = -(farZ+nearZ)*divZD;		aaPM[2][3] = -2.F*farZ*nearZ*divZD;
	aaPM[3][0] = 0;					aaPM[3][1] = 0;					aaPM[3][2] = -1.F;						aaPM[3][3] = 0;
}
void CreateParallelPMtx(DATA aaPM[][4], DATA left, DATA right, DATA bottom, DATA top, DATA nearZ, DATA farZ)
{
	DATA xDiff = right - left;
	DATA yDiff = top - bottom;
	DATA zDiff = farZ - nearZ;
	DATA divXD = 1.F / xDiff;
	DATA divYD = 1.F / yDiff;
	DATA divZD = 1.F / zDiff;

	aaPM[0][0] = 2.F*divXD;		aaPM[0][1] = 0;				aaPM[0][2] = 0;				aaPM[0][3] = -(right+left)*divXD;
	aaPM[1][0] = 0;				aaPM[1][1] = 2.F*divYD;		aaPM[1][2] = 0;				aaPM[1][3] = -(top+bottom)*divYD;
	aaPM[2][0] = 0;				aaPM[2][1] = 0;				aaPM[2][2] = 2.F*divZD;		aaPM[2][3] = -(farZ+nearZ)*divZD;
	aaPM[3][0] = 0;				aaPM[3][1] = 0;				aaPM[3][2] = 0;				aaPM[3][3] = 1.F;
}

void TransformCanonical(TRIModel &canModel, const TRIModel &inModel)
{
	DATA aVM[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, aVM);
	DATA aaViewMtx[][4] = {{aVM[0], aVM[4], aVM[8],  aVM[12]},
						   {aVM[1], aVM[5], aVM[9],  aVM[13]},
						   {aVM[2], aVM[6], aVM[10], aVM[14]},
						   {aVM[3], aVM[7], aVM[11], aVM[15]}};

	DATA aPM[16];
	glGetDoublev(GL_PROJECTION_MATRIX, aPM);
	DATA aaProjMtx[][4] = {{aPM[0], aPM[4], aPM[8],  aPM[12]},
						   {aPM[1], aPM[5], aPM[9],  aPM[13]},
						   {aPM[2], aPM[6], aPM[10], aPM[14]},
						   {aPM[3], aPM[7], aPM[11], aPM[15]}};
	
	DATA aaTransfMtx[4][4];
	myMath.Mul4MM(aaTransfMtx, aaProjMtx, aaViewMtx);
	
	canModel.triangleList.clear();
	for(unsigned t=0; t<inModel.triangleList.size(); t++)
	{
		Triangle triIn = inModel.triangleList[t];
		canModel.triangleList.push_back(triIn);
	}

	for(unsigned t=0; t<canModel.triangleList.size(); t++)
	{
		Triangle &triCan = canModel.triangleList[t];
		for(unsigned v=0; v<3; v++)
		{
			DATA vOld[4] = {0, 0, 0, 1.F};
			myMath.Copy3V(vOld, triCan.vertex[v]);

			DATA vNew[4];
			myMath.Mul4MV(vNew, aaTransfMtx, vOld);
			DATA div = 1.F / vNew[3];
			for(unsigned d=0; d<3; d++)
			{
				triCan.vertex[v][d] = vNew[d] * div;
			}
		}

		const Triangle &triIn = inModel.triangleList[t];
		for(unsigned v=0; v<3; v++)
		{
			DATA nOld[4] = {0, 0, 0, 1.F};
			myMath.Copy3V(nOld, triIn.vertex[v]);
			myMath.Add3V(nOld, triIn.normal[v]);

			DATA nNew[4];
			myMath.Mul4MV(nNew, aaTransfMtx, nOld);
			DATA div = 1.F / nNew[3];
			for(unsigned d=0; d<3; d++)
			{
				triCan.normal[v][d] = nNew[d] * div;
				triCan.normal[v][d] -= triCan.vertex[v][d];
			}
			myMath.Normal3V(triCan.normal[v]);
		}
	}
}

void GetNearCPlane(DATA aaCPlane[][4], DATA l, DATA r, DATA b, DATA t, DATA n)
{
	DATA aX[] = {l, r, r, l};
	DATA aY[] = {b, b, t, t};
	for(unsigned i=0; i<4; i++)
	{
		for(unsigned j=0; j<3; j++)
		{
			aaCPlane[i][j] = aEyePos[j] + aX[i]*aaEyeVec[0][j] + aY[i]*aaEyeVec[1][j] + n*1.001F*aaEyeVec[2][j]; 
		}
		aaCPlane[i][3] = 1.F;
	}
}

void MyGlWindow::draw()
{
	static int bShInit = false;
	if(!bShInit)
	{
		bShInit = true;
		GLenum err = glewInit();

		pFlatSh = new Shader("flatSh.vert", "flatSh.frag");
		pFlatSh->Load();

		pGouraudSh = new Shader("GouraudSh.vert", "GouraudSh.frag");
		pGouraudSh->Load();

		pPhongSh = new Shader("PhongSh.vert", "PhongSh.frag");
		pPhongSh->Load();

		pRCastingSh = new Shader("rayCasting.vert", "rayCasting.frag");
		pRCastingSh->Load();

		pRc2Sh = new Shader("rayCasting.vert", "rayCasting.frag");
		pRc2Sh->Load();

		pLengthImgSh = new Shader("lengthImg.vert", "lengthImg.frag");
		pLengthImgSh->Load();

		pStereoSh = new Shader("stereo.vert", "stereo.frag");
		pStereoSh->Load();

		pFusionSh = new Shader("fusion.vert", "fusion.frag");
		pFusionSh->Load();
	}
	else {}

	static bool bStereoInit = false;
	if(bStereo && !bStereoInit)
	{
		bStereoInit = true;
		for(unsigned i=0; i<2; i++)
		{
			//delete apTexRcast[i];	apTexRcast[i] = new Tex2D(w(), h(), 4);
		}
	}
	else {}

	glUseProgram(0);
	glClearColor(0, 0.2F, 0, 1.F);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(!pTriModel && !pVolModel)
	{
		return;
	}
	else {}

	unsigned runNo = 0;
	unsigned fusNo = 0;
	while((!bStereo && runNo<1) || 
		  ( bStereo && runNo<2))
	{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//**********************************************
	// my transformation
	//**********************************************
	// my projection matrix
	DATA aaNearCPlane[4][4];
	DATA aaPM[4][4];
	DATA r = (DATA)h() / w();
	if(bParal)
	{
		CreateParallelPMtx(aaPM, -viewWide, viewWide, -viewWide*r, viewWide*r, -0.1F, -3.F);
		GetNearCPlane(aaNearCPlane, -viewWide, viewWide, -viewWide*r, viewWide*r, -0.1F);
	}
	else
	{
		CreatePerspectivePMtx(aaPM, -0.2F, 0.2F, -0.2F*r, 0.2F*r, 0.3F, 3.F);
		GetNearCPlane(aaNearCPlane, -0.2F, 0.2F, -0.2F*r, 0.2F*r, -0.3F);
	}
	DATA aProjMtx[] = {aaPM[0][0], aaPM[1][0], aaPM[2][0], aaPM[3][0],
					   aaPM[0][1], aaPM[1][1], aaPM[2][1], aaPM[3][1],
					   aaPM[0][2], aaPM[1][2], aaPM[2][2], aaPM[3][2],
					   aaPM[0][3], aaPM[1][3], aaPM[2][3], aaPM[3][3]};
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(aProjMtx);
	// end of my projection matrix

	// my view matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	DATA aEyeNow[3];
	if(!bStereo)
	{
		aEyeNow[0] = aEyePos[0];
		aEyeNow[1] = aEyePos[1];
		aEyeNow[2] = aEyePos[2];
		//glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	}
	else if(runNo == 0)
	{
		aEyeNow[0] = aEyePos[0] - steDis*aaEyeVec[0][0];
		aEyeNow[1] = aEyePos[1] - steDis*aaEyeVec[0][1];
		aEyeNow[2] = aEyePos[2] - steDis*aaEyeVec[0][2];
		//glColorMask(GL_FALSE, GL_TRUE, GL_TRUE, GL_TRUE);
	}
	else if(runNo == 1)
	{
		aEyeNow[0] = aEyePos[0] + steDis*aaEyeVec[0][0];
		aEyeNow[1] = aEyePos[1] + steDis*aaEyeVec[0][1];
		aEyeNow[2] = aEyePos[2] + steDis*aaEyeVec[0][2];
		//glColorMask(GL_TRUE, GL_FALSE, GL_FALSE, GL_TRUE);
		//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}
	else
	{
		assert(0);
	}
	gluLookAt(aEyeNow[0], aEyeNow[1], aEyeNow[2],
			  aEyeNow[0]-aaEyeVec[2][0], aEyeNow[1]-aaEyeVec[2][1], aEyeNow[2]-aaEyeVec[2][2],
			  //aViewCnt[0], aViewCnt[1], aViewCnt[2],
			  aaEyeVec[1][0], aaEyeVec[1][1], aaEyeVec[1][2]);
	//gluLookAt(aEyePos[0], aEyePos[1], aEyePos[2], 
	//		  aEyePos[0]-aaEyeVec[2][0], aEyePos[1]-aaEyeVec[2][1], aEyePos[2]-aaEyeVec[2][2], 
	//		  aaEyeVec[1][0], aaEyeVec[1][1], aaEyeVec[1][2]);

	glTranslated(aViewCnt[0], aViewCnt[1], aViewCnt[2]);
	RotateObj();
	glScaled(viewScl, viewScl, viewScl);
	glTranslated(-aObjCnt[0], -aObjCnt[1], -aObjCnt[2]);
	// end of my view matrix

	//**********************************************
	// end of my transformation
	//**********************************************

	glViewport(0, 0, w(), h());
	if(eRenderMethod == WIRE)
	{
		TransformCanonical(canonicalModel, *pTriModel);
		
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glUseProgram(0);
		RenderWire(canonicalModel, clipRatio);

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
	}
	else if(eRenderMethod == FLAT_SH)
	{
		glUseProgram(pFlatSh->GetProg());
		RenderShading(*pTriModel);
	}
	else if(eRenderMethod == GOURAUD_SH)
	{
		glUseProgram(pGouraudSh->GetProg());
		RenderShading(*pTriModel);
	}
	else if(eRenderMethod == PHONG_SH)
	{
		glUseProgram(pPhongSh->GetProg());
		RenderShading(*pTriModel);
	}
	else if(eRenderMethod == VOLUME)
	{
		pVolModel->Draw(aaNearCPlane, w(), h(), runNo, fusNo, mskThrd);
	}
	else
	{}

	glUseProgram(0);
	//glPopMatrix();
	
	if(bContour)
	{
		if(!pTexFrame)	
		{
			pFrame = new float[w()*h()*3];
			pTexFrame = new Tex2D(pFrame, w(), h(), 3);	
		}
		else {}
		
		int aWinLow[] = {0, 0};
		pTexFrame->ReadFrame(aWinLow);
		glFlush();

		GetContour(pTexFrame);
		for(unsigned i=0; i<(unsigned)(w()*h()*3); i++)
		{
			pFrame[i] *= 0.5f;
		}

		glMatrixMode(GL_PROJECTION);	glPushMatrix();		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);		glPushMatrix();		glLoadIdentity();

		glRasterPos2f(-1.f, -1.f);
		glDrawPixels(w(), h(), GL_RGB, GL_FLOAT, pFrame);

		glMatrixMode(GL_MODELVIEW);		glPopMatrix();	
		glMatrixMode(GL_PROJECTION);	glPopMatrix();	

		/*
		glMatrixMode(GL_PROJECTION);	glPushMatrix();		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);		glPushMatrix();		glLoadIdentity();
	
		glDisable(GL_DEPTH_TEST);
		float aCoorLow[] = {-1.f, -1.f};
		pTexFrame->DrawRec2D(aCoorLow, 2.f, 2.f, false);
		
		glMatrixMode(GL_MODELVIEW);		glPopMatrix();	
		glMatrixMode(GL_PROJECTION);	glPopMatrix();	
		*/
		if(!bStereo || runNo==1)
		{
			bContour = false;
		}
		else {}
	}
	else {}
	

	/*
	if(bStereo)
	{
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, apTexRcast[runNo]->GetTexID());
		glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, w(), h(), 0);
		glDisable(GL_TEXTURE_2D);
	}
	else {}
	*/
	/*
	if(bStereo && runNo==1)
	{
		glUseProgram(pStereoSh->GetProg());
		//glUseProgram(pFusionSh->GetProg());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, apTexRcast[0]->GetTexID());
		glUniform1i(pStereoSh->GetUniLoc("smpImgR"), 0);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, apTexRcast[1]->GetTexID());
		glUniform1i(pStereoSh->GetUniLoc("smpImgB"), 1);

		glUniform2f(pStereoSh->GetUniLoc("v2WinDim"), (float)w(), (float)h());
		
		glMatrixMode(GL_PROJECTION);	glPushMatrix();		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);		glPushMatrix();		glLoadIdentity();
	
		glDisable(GL_DEPTH_TEST);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glBegin(GL_QUADS);
			glVertex2f(-1.f, -1.f);
			glVertex2f(1.f, -1.f);
			glVertex2f(1.f, 1.f);
			glVertex2f(-1.f, 1.f);
		glEnd();
		
		glMatrixMode(GL_MODELVIEW);		glPopMatrix();	
		glMatrixMode(GL_PROJECTION);	glPopMatrix();	
	}
	else {}
	*/

	if(bCoord)
	{
		DATA boLen = 0.5F / viewScl;
		glBegin(GL_LINES);
			glColor4f(1.f, 0, 0, 1.f);
			glVertex3d(aObjCnt[0],		 aObjCnt[1],	     aObjCnt[2]);
			glVertex3d(aObjCnt[0]+boLen, aObjCnt[1],	     aObjCnt[2]);
			glColor4f(0, 1.f, 0, 1.f);
			glVertex3d(aObjCnt[0],		  aObjCnt[1],		aObjCnt[2]);
			glVertex3d(aObjCnt[0],		  aObjCnt[1]+boLen, aObjCnt[2]);
			glColor4f(0, 0, 1.f, 1.f);
			glVertex3d(aObjCnt[0],		  aObjCnt[1],	     aObjCnt[2]);
			glVertex3d(aObjCnt[0],		  aObjCnt[1],	     aObjCnt[2]+boLen);
		glEnd();
	}

	glFlush();

	if(bFusion)
	{
		fusNo++;
		if(fusNo == 2)
		{
			fusNo = 0;
			runNo++;
		}
		else {}
	}
	else
	{
		runNo++;
	}
	}
}

//*************************************************************************************************

void LoadModel(Fl_Widget *w, void *v)
{
	Fl_File_Chooser *pFChooser = new Fl_File_Chooser("../visualIR/cluster_25d/", "(*.{tri,tf})", Fl_File_Chooser::SINGLE, 0);
	pFChooser->show();
	while(pFChooser->shown())
	{
		Fl::wait();
	}

	if(!pFChooser->value())
	{
		return;
	}
	else {}
	string type = pFChooser->value();
	type = type.substr( type.find_last_of(".")+1 );
	if(type.compare("tf") == 0)
	{
		delete pVolModel;
		pVolModel = new MyVolume(pFChooser->value());

		aObjCnt[0] = 0.5F;
		aObjCnt[1] = 0.5F;
		aObjCnt[2] = 0.5F;
		viewScl = 1.F;

		Fl_Group *pGroup = (Fl_Group*)v;
		int childNum = pGroup->children();
		for(unsigned i=0; (int)i<childNum; i++)
		{
			Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
			pChild->clear();
			pChild->deactivate();
		}

		eRenderMethod = VOLUME;
		Fl_Round_Button *pVol = (Fl_Round_Button*)pGroup->child(4);
		pVol->activate();
		pVol->set();

		pRender_glWin->redraw();
		return;
	}
	else {}

	//**********************************************

	if(eRenderMethod == VOLUME)
	{
		Fl_Group *pGroup = (Fl_Group*)v;
		int childNum = pGroup->children();
		for(unsigned i=0; (int)i<childNum; i++)
		{
			Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
			pChild->activate();
			pChild->clear();
		}
		pGroup->child(4)->deactivate();

		eRenderMethod = FLAT_SH;
		Fl_Round_Button *pFlat = (Fl_Round_Button*)pGroup->child(1);
		pFlat->set();
	}
	else {}

	delete pTriModel;
	pTriModel = new TRIModel;
	pTriModel->loadFromFile(pFChooser->value());
	delete pFChooser;

	DATA aaBoundBox[3][2];
	for(unsigned i=0; i<3; i++)
	{
		aaBoundBox[i][0] = 1e10;	
		aaBoundBox[i][1] = -1e10;
	}
	for(unsigned t=0; t<pTriModel->triangleList.size(); t++)
	{
		Triangle &tri = pTriModel->triangleList[t];
		for(unsigned v=0; v<3; v++)
		{
			for(unsigned d=0; d<3; d++)
			{
				if(tri.vertex[v][d] < aaBoundBox[d][0])
				{
					aaBoundBox[d][0] = tri.vertex[v][d];
				}
				else {}

				if(tri.vertex[v][d] > aaBoundBox[d][1])
				{
					aaBoundBox[d][1] = tri.vertex[v][d];
				}
				else {}
			}
		}
	}
	for(unsigned i=0; i<3; i++)
	{
		aObjCnt[i] = (aaBoundBox[i][0] + aaBoundBox[i][1]) * 0.5F;
	}

	DATA maxLen = 0;
	for(unsigned i=0; i<3; i++)
	{
		DATA len = aaBoundBox[i][1] - aaBoundBox[i][0];
		if(len > maxLen)
		{
			maxLen = len;
		}
		else {}
	}
	viewScl = 1.F / maxLen;

	pRender_glWin->redraw();
}

//*************************************************************************************************

void ReloadTF(Fl_Widget *w, void *v)
{
	if(!pVolModel)	return;
	else {}

	pVolModel->ReloadTF(mskThrd);
	pRender_glWin->redraw();
}

//*************************************************************************************************
//
//*************************************************************************************************

void SetWire(Fl_Widget *w, void *v)
{
	Fl_Group *pGroup = (Fl_Group*)v;
	int childNum = pGroup->children();
	for(unsigned i=0; (int)i<childNum; i++)
	{
		Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
		pChild->clear();
	}

	eRenderMethod = WIRE;
	Fl_Round_Button *pWire = (Fl_Round_Button*)w;
	pWire->set();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetFlatShading(Fl_Widget *w, void *v)
{
	Fl_Group *pGroup = (Fl_Group*)v;
	int childNum = pGroup->children();
	for(unsigned i=0; (int)i<childNum; i++)
	{
		Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
		pChild->clear();
	}

	eRenderMethod = FLAT_SH;
	Fl_Round_Button *pShading = (Fl_Round_Button*)w;
	pShading->set();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetGouraudShading(Fl_Widget *w, void *v)
{
	Fl_Group *pGroup = (Fl_Group*)v;
	int childNum = pGroup->children();
	for(unsigned i=0; (int)i<childNum; i++)
	{
		Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
		pChild->clear();
	}

	eRenderMethod = GOURAUD_SH;
	Fl_Round_Button *pShading = (Fl_Round_Button*)w;
	pShading->set();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetPhongShading(Fl_Widget *w, void *v)
{
	Fl_Group *pGroup = (Fl_Group*)v;
	int childNum = pGroup->children();
	for(unsigned i=0; (int)i<childNum; i++)
	{
		Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
		pChild->clear();
	}

	eRenderMethod = PHONG_SH;
	Fl_Round_Button *pShading = (Fl_Round_Button*)w;
	pShading->set();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void SetEyes(Fl_Widget *w, void *v)
{
	Fl_Group *pGroup = (Fl_Group*)v;
	int childNum = pGroup->children();
	for(unsigned i=0; (int)i<childNum; i++)
	{
		Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
		pChild->clear();
	}

	bEyes = true;
	Fl_Round_Button *pEyes = (Fl_Round_Button*)w;
	pEyes->set();
}
void SetObj(Fl_Widget *w, void *v)
{
	Fl_Group *pGroup = (Fl_Group*)v;
	int childNum = pGroup->children();
	for(unsigned i=0; (int)i<childNum; i++)
	{
		Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
		pChild->clear();
	}

	bEyes = false;
	Fl_Round_Button *pObj = (Fl_Round_Button*)w;
	pObj->set();
}

void SetParal(Fl_Widget *w, void *v)
{
	Fl_Group *pGroup = (Fl_Group*)v;
	int childNum = pGroup->children();
	for(unsigned i=0; (int)i<childNum; i++)
	{
		Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
		pChild->clear();
	}

	bParal = true;
	Fl_Round_Button *pParal = (Fl_Round_Button*)w;
	pParal->set();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetPersp(Fl_Widget *w, void *v)
{
	Fl_Group *pGroup = (Fl_Group*)v;
	int childNum = pGroup->children();
	for(unsigned i=0; (int)i<childNum; i++)
	{
		Fl_Round_Button *pChild = (Fl_Round_Button*)pGroup->child(i);
		pChild->clear();
	}

	bParal = false;
	Fl_Round_Button *pPersp = (Fl_Round_Button*)w;
	pPersp->set();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetCoord(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pCoord = (Fl_Check_Button*)w;
	bCoord = (pCoord->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

//*************************************************************************************************
//
//*************************************************************************************************

void ChangeDisAtt(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pDisAtt = (Fl_Value_Slider*)w;
	double val = pDisAtt->value();
	disAtt = val;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void ChangeOpaScl(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pOpaScl = (Fl_Value_Slider*)w;
	double val = pOpaScl->value();
	opaScl = val;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void ChangeBriScl(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pBriScl = (Fl_Value_Slider*)w;
	double val = pBriScl->value();
	briScl = val;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void ChangeSmpNum(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pSmpNum = (Fl_Value_Slider*)w;
	double val = pSmpNum->value();
	smpNum = val;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void ChangeMskThrd(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pMskThrd = (Fl_Value_Slider*)w;
	double val = pMskThrd->value();
	mskThrd = val;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void Relight(Fl_Widget *w, void *v)
{
	pVolModel->MakeShadow(mskThrd);

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void SetMIP(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pM = (Fl_Check_Button*)w;
	bMIP = (pM->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void ChangeOccScl(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pOccScl = (Fl_Value_Slider*)w;
	double val = pOccScl->value();
	occScl = val;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void SetXClipping(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBX = (Fl_Check_Button*)w;
	bXClip = (pBX->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetXFront(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBF = (Fl_Check_Button*)w;
	bXFront = (pBF->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void ChangeXPlane(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pXPlane = (Fl_Value_Slider*)w;
	double val = pXPlane->value();
	xPlane = val;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetYClipping(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBY = (Fl_Check_Button*)w;
	bYClip = (pBY->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetYFront(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBF = (Fl_Check_Button*)w;
	bYFront = (pBF->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void ChangeYPlane(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pYPlane = (Fl_Value_Slider*)w;
	double val = pYPlane->value();
	yPlane = val;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetZClipping(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBZ = (Fl_Check_Button*)w;
	bZClip = (pBZ->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void SetZFront(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBF = (Fl_Check_Button*)w;
	bZFront = (pBF->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
void ChangeZPlane(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pZPlane = (Fl_Value_Slider*)w;
	double val = pZPlane->value();
	zPlane = val;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

//*************************************************************************************************
//
//*************************************************************************************************

void RunStereo(Fl_Widget *w, void *v)
{
	mainWinPos.m_x = pMain_win->x();	mainWinPos.m_y = pMain_win->y();
	mainWinLen.m_x = pMain_win->w();	mainWinLen.m_y = pMain_win->h();
	glWinPos.m_x = pRender_glWin->x();	glWinPos.m_y = pRender_glWin->y();
	glWinLen.m_x = pRender_glWin->w();	glWinLen.m_y = pRender_glWin->h();

	//pMain_win->fullscreen();
	//pRender_glWin->fullscreen();
	bStereo = true;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void p(void*)
{
	pRender_glWin->redraw();
	Fl::repeat_timeout(0, p);
}
void RunStFusion(Fl_Widget *w, void *v)
{
	mainWinPos.m_x = pMain_win->x();	mainWinPos.m_y = pMain_win->y();
	mainWinLen.m_x = pMain_win->w();	mainWinLen.m_y = pMain_win->h();
	glWinPos.m_x = pRender_glWin->x();	glWinPos.m_y = pRender_glWin->y();
	glWinLen.m_x = pRender_glWin->w();	glWinLen.m_y = pRender_glWin->h();

	//pMain_win->fullscreen();
	//pRender_glWin->fullscreen();
	bStereo = true;
	bFusion = true;

	Fl::focus(pRender_glWin);
	//pRender_glWin->redraw();

	wglSwapIntervalEXT(1);
	Fl::add_timeout(0.1F, p);
}
void Unstereo()
{
	pMain_win->fullscreen_off(mainWinPos.m_x, mainWinPos.m_y, mainWinLen.m_x, mainWinLen.m_y);
	pRender_glWin->fullscreen_off(glWinPos.m_x, glWinPos.m_y, glWinLen.m_x, glWinLen.m_y);
	bStereo = false;
	bFusion = false;
	wglSwapIntervalEXT(0);
	Fl::remove_timeout(p);
	bMap = false;
}

//*************************************************************************************************
//
//*************************************************************************************************

void fps(void*)
{
	static unsigned numF = 0;
	static time_t strT;
	static time_t endT;

	if(numF == 0)
	{
		strT = time(0);
	}
	else {}
	pRender_glWin->redraw();
	numF++;
	if(numF == 100)
	{
		endT = time(0);
		cout << "fps = " << (DATA)numF/(endT-strT) << endl;
		numF = 0;
	}
	else {}

	Fl::repeat_timeout(0, fps);
}
void CountFps(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pFps = (Fl_Check_Button*)w;
	wglSwapIntervalEXT(0);
	if(pFps->value())
	{
		cout << "counting fps.." << endl;
		Fl::add_timeout(0.1F, fps);
	}
	else
	{
		cout << "disable fps" << endl;
		Fl::remove_timeout(fps);
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

void SetContour(Fl_Widget *w, void *v)
{
	bContour = true;	

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void SetDivShow(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBDiv = (Fl_Check_Button*)w;
	bDivShow = (pBDiv->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void SetRandDis(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pRnd = (Fl_Check_Button*)w;
	bRndDis = (pRnd->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

//***********************************************

void ChangeDivThrd(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pDiv = (Fl_Value_Slider*)w;
	divThrd = pDiv->value();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

//***********************************************

void SetBrightness(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBBrightness = (Fl_Check_Button*)w;
	bBrightness_div = (pBBrightness->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void ChangeBrightness(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pBrightness = (Fl_Value_Slider*)w;
	brightnessVal = pBrightness->value();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

//***********************************************

void SetContrast(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBContrast = (Fl_Check_Button*)w;
	bContrast_div = (pBContrast->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void ChangeContrast(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pContrast = (Fl_Value_Slider*)w;
	contrastVal = pContrast->value();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();	
}

//***********************************************

void SetDim(Fl_Widget *w, void *v)
{
	Fl_Check_Button *pBDim = (Fl_Check_Button*)w;
	bDim_div = (pBDim->value())? true: false;

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}

void ChangeDim(Fl_Widget *w, void *v)
{
	Fl_Value_Slider *pDim = (Fl_Value_Slider*)w;
	dimVal = pDim->value();

	Fl::focus(pRender_glWin);
	pRender_glWin->redraw();
}
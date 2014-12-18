#ifndef _MYWIN_APP_H
#define _MYWIN_APP_H

#include "TRIModel.h"
#include "define.h"
#include "myMath.h"
#include "renderWire.h"
#include "renderShading.h"
#include "myVolume.h"
#include "texture.h"

#include <windows.h>
#include <GL/glew.h>
#include <FL/Fl_Gl_Window.h>
#include <FL/Fl_File_Chooser.h>
#include <FL/Fl_Round_Button.h>
#include <FL/Fl_Check_Button.h>
#include <FL/Fl_Value_Slider.h>
#include <GL/wglew.h>
#include <cassert>
#include <ctime>

class MyGlWindow: public Fl_Gl_Window
{
public:
	MyGlWindow(int x, int y, int w, int h, const char *l=0);
	void draw();
	int handle(int evnt);

private:
};

extern Fl_Double_Window *pMain_win;
extern MyGlWindow *pRender_glWin;
extern DATA disAtt;
extern DATA opaScl;
extern DATA briScl;
extern DATA smpNum;
extern DATA localL;

extern bool bMIP;
extern DATA occScl;

extern bool bXClip;	extern bool bXFront;	extern DATA xPlane;
extern bool bYClip;	extern bool bYFront;	extern DATA yPlane;
extern bool bZClip;	extern bool bZFront;	extern DATA zPlane;

extern bool bStereo;
extern bool bFusion;
//extern float *apRcast[2];	
//extern Tex2D *apTexRcast[2];
//extern float *apRfuse[2];	
//extern Tex2D *apTexRfuse[2];
//extern float *apImg[2];		
//extern Tex2D *apTexImg[2];

extern bool bDivShow;
extern DATA divThrd;

extern Vect3D<DATA> dimRatio;
extern Vect3D<DATA> pG;
extern Vect3D<DATA> pR;

extern bool bSymetry;

void LoadModel	(Fl_Widget *w, void *v);
void ReloadTF	(Fl_Widget *w, void *v);

void SetEyes			(Fl_Widget *w, void *v);
void SetObj				(Fl_Widget *w, void *v);
void SetParal			(Fl_Widget *w, void *v);
void SetPersp			(Fl_Widget *w, void *v);
void SetCoord			(Fl_Widget *w, void *v);
void SetWire			(Fl_Widget *w, void *v);
void SetFlatShading		(Fl_Widget *w, void *v);
void SetGouraudShading	(Fl_Widget *w, void *v);
void SetPhongShading	(Fl_Widget *w, void *v);

void ChangeDisAtt		(Fl_Widget *w, void *v);
void ChangeOpaScl		(Fl_Widget *w, void *v);
void ChangeBriScl		(Fl_Widget *w, void *v);
void ChangeSmpNum		(Fl_Widget *w, void *v);

void ChangeLocalL		(Fl_Widget *w, void *v);
void Relight			(Fl_Widget *w, void *v);

void SetMIP				(Fl_Widget *w, void *v);
void ChangeOccScl		(Fl_Widget *w, void *v);

void CountFps			(Fl_Widget *w, void *v);

void SetXClipping		(Fl_Widget *w, void *v);
void SetYClipping		(Fl_Widget *w, void *v);
void SetZClipping		(Fl_Widget *w, void *v);
void SetDClipping		(Fl_Widget *w, void *v);
void SetXFront			(Fl_Widget *w, void *v);
void SetYFront			(Fl_Widget *w, void *v);
void SetZFront			(Fl_Widget *w, void *v);
void SetDFront			(Fl_Widget *w, void *v);
void ChangeXPlane		(Fl_Widget *w, void *v);
void ChangeYPlane		(Fl_Widget *w, void *v);
void ChangeZPlane		(Fl_Widget *w, void *v);
void ChangeDPlane		(Fl_Widget *w, void *v);

//*************************************************************************************************

void RunStereo			(Fl_Widget *w, void *v);

void SetDivShow			(Fl_Widget *w, void *v);
void SetRandDis			(Fl_Widget *w, void *v);
void SetAODiv			(Fl_Widget *w, void *v);
void SetDisDiv			(Fl_Widget *w, void *v);
void ChangeDivThrd		(Fl_Widget *w, void *v);

void SetKey				(Fl_Widget *w, void *v);		void ChangeKey			(Fl_Widget *w, void *v);
void SetFeature			(Fl_Widget *w, void *v);		void ChangeFeaScl		(Fl_Widget *w, void *v);
void SetDim				(Fl_Widget *w, void *v);		void ChangeDim			(Fl_Widget *w, void *v);

void RunStFusion		(Fl_Widget *w, void *v);
void SetContour			(Fl_Widget *w, void *v);

#endif
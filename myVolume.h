#ifndef _MY_VOLUME_H
#define _MY_VOLUME_H

#include "define.h"
#include "myMath.h"
#include "shader.h"
#include "volumeData.h"
#include "matrixOperation.h"
#include "myWinApp.h"
#include "layerOperation.h"
#include "texture.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <GL/glew.h>
#include <ctime>

using namespace std;

typedef struct
{
	DATA vMin;
	DATA vMax;
} BLCK;

class MyVolume
{
public:
	MyVolume(string f);
	~MyVolume();

	void GetColor();
	void ReloadTF(DATA mskThd);

	void MakeShadow(DATA mskThd);
	//void Relight();

	void Draw(DATA aaCPlane[][4], unsigned w, unsigned h, unsigned runNo, unsigned fusNo, DATA mskThd);

private:
	void BuildBlock(unsigned step);
	void DecideActBlock();

	Vect3D<unsigned> m_vDim;
	Vect3D<DATA> m_vStep;

	float *m_pData;
	float m_valMin, m_valMax;
	Tex3D *m_pTexData;

	float *m_pColor;
	float *m_pId;
	unsigned m_colDim;
	Tex1D *m_pTexId;

	float *m_pPreint;
	Vect2D<unsigned> m_vPreintDim;
	Tex2D *m_pTexPreint;

	BLCK *m_pBlock;
	float *m_pBAct;
	Tex3D *m_pTexAct;
	//GLuint m_texAct;
	Vect3D<unsigned> m_vBlckDim;
	Vect3D<float> m_vBlckRatio;
	unsigned m_blckStep;

	string m_fName;
	VolumeFile *m_pVolF;

	bool m_bShdInit;
	Layer *m_pLyrShdInit;
	float *m_pShadow;	Tex3D *m_pTexShadow;
	float *m_pShLocal;	Tex3D *m_pTexShLocal;

	float *m_pMsk;
	Tex3D *m_pTexMsk;
};

extern Shader *pLengthImgSh;
extern Shader *pRCastingSh;
extern Shader *pRc2Sh;
extern Shader *pStereoSh;
extern Shader *pFusionSh;

//typedef enum {AO, DISTANCE} DIVISION_TYPE;
//extern DIVISION_TYPE eDivMethod;
extern bool bDivShow;
extern bool bRndDis;
extern DATA divThrd;

extern bool bKey_div;
extern DATA keyVal;

extern bool bFeature_div;
extern DATA featureVal;

extern bool bDim_div;
extern DATA dimVal;

#endif
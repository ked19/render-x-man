#include "myVolume.h"

extern DATA aObjCnt[3];
extern DATA aaObjVec[4][4];
extern DATA aEyePos[3];
extern DATA aaEyeVec[4][4];
extern DATA aViewCnt[3];
extern DATA steDis;
extern DATA viewScl;

Shader *pLengthImgSh = 0;
Shader *pRCastingSh = 0;
Shader *pRc2Sh = 0;
Shader *pStereoSh = 0;
Shader *pFusionSh = 0;

//DIVISION_TYPE eDivMethod = DISTANCE;
bool bDivShow = false;
bool bRndDis = true;
DATA divThrd = 2.F;

bool bBrightness_div = true;
DATA brightnessVal = 0.15F; 

bool bContrast_div = true;
DATA contrastVal = 1.5F;

bool bDim_div = true;
DATA dimVal = 0.3F;

MyVolume::MyVolume(string f)
	:m_vDim(0, 0, 0), m_vStep(0, 0, 0), m_pData(0), m_pTexData(0) 
	,m_vPreintDim(0, 0) ,m_pPreint(0), m_pTexPreint(0)
	,m_pVolF(0)
	,m_vBlckDim(0, 0, 0), m_blckStep(0), m_pBlock(0), m_vBlckRatio(0, 0, 0)
	,m_pBAct(0), m_pTexAct(0)
	,m_colDim(0), m_pColor(0), m_pId(0), m_pTexId(0)
	,m_bShdInit(false), m_pLyrShdInit(0)
	,m_pShadow(0), m_pTexShadow(0)
	,m_pMsk(0), m_pTexMsk(0)
{
	m_fName = f;
	m_fName = m_fName.substr(0, m_fName.find_last_of('.'));
	cout << m_fName << endl;

	m_pVolF = new VolumeFile(m_fName+".vol");
	const VolumeData &rVD = m_pVolF->GetOrgVolVal();
	m_vStep = rVD.GetStep();
	m_vDim = rVD.GetDim();

	unsigned tSize = m_vDim.m_x * m_vDim.m_y * m_vDim.m_z;
	m_pData = new float[tSize];
	rVD.GetLyrVal().CopyTo_zLast(m_pData);

	myMath.Range(m_valMin, m_valMax, m_pData, tSize);
	cout << "valMax: " << m_valMax << endl;
	cout << "valMin: " << m_valMin << endl;

	float invDis = 1.f / (m_valMax - m_valMin);
	for(unsigned i=0; i<tSize; i++)
	{
		m_pData[i] = (m_pData[i] - m_valMin) * invDis;
	}

	//*******************************************

	GetColor();

	BuildBlock(8);
	DecideActBlock();

	//*******************************************

	m_pTexData = new Tex3D(m_pData, m_vDim.m_x, m_vDim.m_y, m_vDim.m_z, 1, false, true); 
	
	m_pTexAct = new Tex3D(m_pBAct, m_vBlckDim.m_x, m_vBlckDim.m_y, m_vBlckDim.m_z, 1, false, false);
}

MyVolume::~MyVolume()
{
	delete []m_pData;		delete m_pTexData;
	delete []m_pColor;		
	delete []m_pId;			delete m_pTexId;
	delete []m_pPreint;		delete m_pTexPreint;
	delete []m_pShadow;		delete m_pTexShadow;
	delete []m_pMsk;		delete m_pTexMsk;
	delete []m_pBAct;		delete m_pTexAct;
	delete m_pVolF;
	delete m_pLyrShdInit;
}

//*************************************************************************************************

void MyVolume::ReloadTF(DATA mskThd)
{
	GetColor();

	//BuildBlock(8);
	DecideActBlock();

	MakeShadow(mskThd);

	m_pTexAct->Update(0, 0, 0, m_vBlckDim.m_x, m_vBlckDim.m_y, m_vBlckDim.m_z);
	//m_pTexAct = new Tex3D(m_pBAct, m_vBlckDim.m_x, m_vBlckDim.m_y, m_vBlckDim.m_z, 1, false, false);
}

//*************************************************************************************************

void MyVolume::BuildBlock(unsigned step)
{
	m_blckStep = step;

	m_vBlckDim.m_x = m_vDim.m_x / step;
	if(m_vDim.m_x%step != 0)
	{
		m_vBlckDim.m_x++;
	}
	else {}
	//m_vBlckRatio.m_x = (float)m_vDim.m_x / (m_vBlckDim.m_x * step);
	m_vBlckRatio.m_x = (m_vDim.m_x-1.f) / (m_vBlckDim.m_x*step-1.f);

	m_vBlckDim.m_y = m_vDim.m_y / step;
	if(m_vDim.m_y%step != 0)
	{
		m_vBlckDim.m_y++;
	}
	else {}
	//m_vBlckRatio.m_y = (float)m_vDim.m_y / (m_vBlckDim.m_y * step);
	m_vBlckRatio.m_y = (m_vDim.m_y-1.f) / (m_vBlckDim.m_y*step-1.f);

	m_vBlckDim.m_z = m_vDim.m_z / step;
	if(m_vDim.m_z%step != 0)
	{
		m_vBlckDim.m_z++;
	}
	else {}
	//m_vBlckRatio.m_z = (float)m_vDim.m_z / (m_vBlckDim.m_z * step);
	m_vBlckRatio.m_z = (m_vDim.m_z-1.f) / (m_vBlckDim.m_z*step-1.f);

	unsigned size = m_vBlckDim.m_x * m_vBlckDim.m_y * m_vBlckDim.m_z;
	delete []m_pBlock;
	m_pBlock = new BLCK[size];

	for(unsigned z=0; z<m_vBlckDim.m_z; z++)
	{
		for(unsigned y=0; y<m_vBlckDim.m_y; y++)
		{
			for(unsigned x=0; x<m_vBlckDim.m_x; x++)
			{
				DATA vMin = 1e10;
				DATA vMax = -1e10;
				for(int w=-1; w<=(int)step; w++)
				{
					int zLoc = z*step + w;
					zLoc = (int)(zLoc*m_vBlckRatio.m_z + 0.5f);
					if(zLoc>=(int)m_vDim.m_z || zLoc<0)
					{
						continue;
					}
					else {}

					for(int v=-1; v<=(int)step; v++)
					{
						int yLoc = y*step + v;
						yLoc = (int)(yLoc*m_vBlckRatio.m_y + 0.5f);
						if(yLoc>=(int)m_vDim.m_y || yLoc<0)
						{
							continue;
						}
						else {}

						for(int u=-1; u<=(int)step; u++)
						{
							int xLoc = x*step + u;
							xLoc = (int)(xLoc*m_vBlckRatio.m_x + 0.5f);
							if(xLoc>=(int)m_vDim.m_x || xLoc<0)
							{
								continue;
							}
							else {}

							unsigned vPos = zLoc*m_vDim.m_y*m_vDim.m_x + yLoc*m_vDim.m_x + xLoc;
							DATA val = m_pData[vPos];
							if(val < vMin)
							{
								vMin = val;
							}
							else {}
							if(val > vMax)
							{
								vMax = val;
							}
							else {}
						} // u
					} // v
				} // w

				unsigned bPos = z*m_vBlckDim.m_y*m_vBlckDim.m_x + y*m_vBlckDim.m_x + x;
				m_pBlock[bPos].vMin = vMin;
				m_pBlock[bPos].vMax = vMax;
			} // x
		} // y
	} // z
}

void MyVolume::DecideActBlock()
{
	assert(m_pBlock);
	assert(m_pColor);

	DATA *pOpaAcc = new DATA[m_colDim];
	pOpaAcc[0] = m_pColor[3];
	for(unsigned i=1; i<m_colDim; i++)
	{
		pOpaAcc[i] = pOpaAcc[i-1] + m_pColor[i*4+3];	
	}

	unsigned bSize = m_vBlckDim.m_x * m_vBlckDim.m_y * m_vBlckDim.m_z;
	delete []m_pBAct;
	m_pBAct = new float[bSize];

	for(unsigned i=0; i<bSize; i++)
	{
		unsigned minIdx = (unsigned)(m_pBlock[i].vMin * (m_colDim-1.F) + 0.5F);
		unsigned maxIdx = (unsigned)(m_pBlock[i].vMax * (m_colDim-1.F) + 0.5F);

		DATA opa = pOpaAcc[maxIdx] - pOpaAcc[minIdx];
		if(opa >= 1e-8 || m_pColor[minIdx*4+3]>0)
		{
			m_pBAct[i] = 1.f;
		}
		else
		{
			m_pBAct[i] = 0;
		}
	}

	delete []pOpaAcc;
}

//*************************************************************************************************

void GenPreint(Layer &lyrPreint, float *pColor, unsigned dimColor)
{
	Vect3D<unsigned> dimPreint = lyrPreint.GetDim();
	unsigned intNum = 1000;

	cout << "pre-integration." << endl;
	for(unsigned y=0; y<dimPreint.m_y; y++)
	{
		cout << y << " ";
		DATA sY = 1.F * y / (dimPreint.m_y-1);
		for(unsigned x=0; x<dimPreint.m_x; x++)
		{
			DATA sX = 1.F * x / (dimPreint.m_x-1);
			DATA aCPreint[] = {0, 0, 0, 0};
			for(unsigned ii=0; ii<intNum; ii++)
			{
				DATA sNow = sX + (sY-sX)*ii/(intNum-1.F);
				sNow *= (dimColor - 1.F);
				if(sNow < 0)
				{
					sNow = 0;
				}

				unsigned idxL = (unsigned)sNow; 
				unsigned idxH = idxL + 1;
				float aColor[4];
				for(unsigned i=0; i<3; i++)
				{
					DATA cL = pColor[idxL*4+i] * pColor[idxL*4+3];
					DATA cH = pColor[idxH*4+i] * pColor[idxH*4+3];
					aColor[i] = (float)myMath.Interpolate_linear(cL, cH, idxH-sNow);
				}
				DATA aL = pColor[idxL*4+3];
				DATA aH = pColor[idxH*4+3];
				if(idxH-sNow>1.f || idxH-sNow<0)
					cout << idxH << " " << sNow << ",  ";   
				aColor[3] = (float)myMath.Interpolate_linear(aL, aH, idxH-sNow); 

				aCPreint[0] = aCPreint[0] + (1.F-aCPreint[3])*aColor[0]/intNum; 
				aCPreint[1] = aCPreint[1] + (1.F-aCPreint[3])*aColor[1]/intNum; 
				aCPreint[2] = aCPreint[2] + (1.F-aCPreint[3])*aColor[2]/intNum; 
				aCPreint[3] = aCPreint[3] + (1.F-aCPreint[3])*aColor[3]/intNum; 
			}

			lyrPreint.CellRef(x, y, 0) = aCPreint[0];
			lyrPreint.CellRef(x, y, 1) = aCPreint[1];
			lyrPreint.CellRef(x, y, 2) = aCPreint[2];
			lyrPreint.CellRef(x, y, 3) = aCPreint[3];
		}
	}
	cout << "ok" << endl;
}

void MyVolume::GetColor()
{
	string strTf = m_fName + ".tf";
	ifstream inF(strTf.c_str(), ios::binary);
	inF.read((char*)&m_colDim, sizeof(unsigned));
	delete []m_pColor;
	m_pColor = new float[m_colDim*4];
	inF.read((char*)m_pColor, sizeof(float)*m_colDim*4);
	inF.close();

	//*******************************************

	string strId = m_fName + ".id";
	ifstream inF2(strId.c_str(), ios::binary);
	inF2.read((char*)&m_colDim, sizeof(unsigned));
//cout << m_colDim << endl;
	delete []m_pId;
	m_pId = new float[m_colDim];
	inF2.read((char*)m_pId, sizeof(float)*m_colDim);
for(unsigned i=0; i<m_colDim; i++)
{
	cout << m_pId[i] << " ";
}
	inF2.close();

	m_pTexId = new Tex1D(m_pId, m_colDim, 1, false, false, false); 

	//*******************************************

	bool bRead = false; //false;
	string volName = m_fName.substr(m_fName.find_last_of("/")+1);
	string preFName = volName+"_out.pre"; 
	if(!bRead)
	{
		Layer lyrPreint(m_colDim, m_colDim, 4); //(512, 512, 4);
		GenPreint(lyrPreint, m_pColor, m_colDim);

		Vect3D<unsigned> dimPreint = lyrPreint.GetDim();
		m_pPreint = new float[dimPreint.m_x*dimPreint.m_y*4];
		m_vPreintDim.m_x = dimPreint.m_x;
		m_vPreintDim.m_y = dimPreint.m_y;
		for(unsigned y=0; y<dimPreint.m_y; y++)
		{
			for(unsigned x=0; x<dimPreint.m_x; x++)
			{
				unsigned loc = y*dimPreint.m_x + x;
				m_pPreint[loc*4] = (float)lyrPreint.CellVal(x, y, 0);
				m_pPreint[loc*4+1] = (float)lyrPreint.CellVal(x, y, 1);
				m_pPreint[loc*4+2] = (float)lyrPreint.CellVal(x, y, 2);
				m_pPreint[loc*4+3] = (float)lyrPreint.CellVal(x, y, 3);
			}
		}

		ofstream preF(preFName.c_str(), ios::binary);
		preF.write((char*)&m_vPreintDim.m_x, sizeof(unsigned));
		preF.write((char*)&m_vPreintDim.m_y, sizeof(unsigned));
		preF.write((char*)m_pPreint, sizeof(float)*m_vPreintDim.m_x*m_vPreintDim.m_y*4);
		preF.close();
	}
	else
	{
		ifstream preF(preFName.c_str(), ios::binary);
		preF.read((char*)&m_vPreintDim.m_x, sizeof(unsigned));
		preF.read((char*)&m_vPreintDim.m_y, sizeof(unsigned));
		m_pPreint = new float[m_vPreintDim.m_x*m_vPreintDim.m_y*4];
		preF.read((char*)m_pPreint, sizeof(float)*m_vPreintDim.m_x*m_vPreintDim.m_y*4);
		preF.close();
	}

	m_pTexPreint = new Tex2D(m_pPreint, m_vPreintDim.m_x, m_vPreintDim.m_y, 4); 
}

//*************************************************************************************************

void valTex2d(float aVal[], float *pTex, unsigned w, unsigned h, unsigned ch, unsigned x, unsigned y)
{
	unsigned loc = (y*w + x) * ch;
	for(unsigned i=0; i<ch; i++)
	{
		aVal[i] = pTex[loc + i];
	}
}
void valTex3d(float aVal[], float *pTex, unsigned w, unsigned h, unsigned d, unsigned ch, unsigned x, unsigned y, unsigned z)
{
	unsigned loc = (z*w*h + y*w + x) * ch;
	for(unsigned i=0; i<ch; i++)
	{
		aVal[i] = pTex[loc + i];
	}
}

void InterpolateTex2d(float aVal[], float *pTex, unsigned w, unsigned h, unsigned ch, float xIn, float yIn)
{
	float x = xIn * (w-1);
	float y = yIn * (h-1);

	unsigned xL = (unsigned)x;
	unsigned yB = (unsigned)y;
	float aLbVal[4];
	valTex2d(aLbVal, pTex, w, h, ch, xL, yB);
	if(xL>=w-1 || yB>=h-1)
	{
		for(unsigned i=0; i<ch; i++)
		{
			aVal[i] = aLbVal[i];
		}
		return;
	}
	else {}

	float aRbVal[4], aLtVal[4], aRtVal[4];
	valTex2d(aRbVal, pTex, w, h, ch, xL+1, yB);
	valTex2d(aLtVal, pTex, w, h, ch, xL,   yB+1);
	valTex2d(aRtVal, pTex, w, h, ch, xL+1, yB+1);

	DATA ratioL = x - xL;
	DATA ratioB = y - yB;
	for(unsigned i=0; i<ch; i++)
	{
		DATA lv = aLbVal[i];
		DATA rb = aRbVal[i];
		DATA lt = aLtVal[i];
		DATA rt = aRtVal[i];
		aVal[i] = (float)myMath.Interpolate_linear(lv, rb, lt, rt, ratioL, ratioB);
	}
}
void InterpolateTex3d(float aVal[], float *pTex, unsigned w, unsigned h, unsigned d, unsigned ch, float xIn, float yIn, float zIn)
{
	float x = xIn * (w-1);
	float y = yIn * (h-1);
	float z = zIn * (d-1);

	unsigned xL = (unsigned)x;
	unsigned yB = (unsigned)y;
	unsigned zN = (unsigned)z;
	float aLbnVal[4];
	valTex3d(aLbnVal, pTex, w, h, d, ch, xL, yB, zN);
	if(xL>=w-1 || yB>=h-1 || zN>=d-1)
	{
		for(unsigned i=0; i<ch; i++)
		{
			aVal[i] = aLbnVal[i];
		}
		return;
	}
	else {}

	float aRbnVal[4], aLtnVal[4], aRtnVal[4], aLbfVal[4], aRbfVal[4], aLtfVal[4], aRtfVal[4];
	valTex3d(aRbnVal, pTex, w, h, d, ch, xL+1, yB,   zN);
	valTex3d(aLtnVal, pTex, w, h, d, ch, xL,   yB+1, zN);
	valTex3d(aRtnVal, pTex, w, h, d, ch, xL+1, yB+1, zN);
	valTex3d(aLbfVal, pTex, w, h, d, ch, xL,   yB,   zN+1);
	valTex3d(aRbfVal, pTex, w, h, d, ch, xL+1, yB,   zN+1);
	valTex3d(aLtfVal, pTex, w, h, d, ch, xL,   yB+1, zN+1);
	valTex3d(aRtfVal, pTex, w, h, d, ch, xL+1, yB+1, zN+1);

	DATA ratioL = x - xL;
	DATA ratioB = y - yB;
	DATA ratioN = z - zN;
	for(unsigned i=0; i<ch; i++)
	{
		DATA lbn = aLbnVal[i];
		DATA rbn = aRbnVal[i];
		DATA ltn = aLtnVal[i];
		DATA rtn = aRtnVal[i];
		DATA lbf = aLbfVal[i];
		DATA rbf = aRbfVal[i];
		DATA ltf = aLtfVal[i];
		DATA rtf = aRtfVal[i];
		aVal[i] = (float)myMath.Interpolate_linear(lbn, rbn, ltn, rtn, lbf, rbf, ltf, rtf, ratioL, ratioB, ratioN);
	}
}

Layer *ReadMsk(string fName)
{
	string strMsk = fName + ".msk";
	ifstream mskF(strMsk.c_str(), ios::binary);
	cout << "read " << fName << ".msk" << endl;

	Vect3D<unsigned> dimMsk(0, 0, 0);
	mskF.read((char*)&dimMsk.m_x, sizeof(unsigned));
	mskF.read((char*)&dimMsk.m_y, sizeof(unsigned));
	mskF.read((char*)&dimMsk.m_z, sizeof(unsigned));

	Layer *pLyrMsk = new Layer(dimMsk.m_x, dimMsk.m_y, dimMsk.m_z);
	for(unsigned z=0; z<dimMsk.m_z; z++)
	{
		for(unsigned y=0; y<dimMsk.m_y; y++)
		{
			for(unsigned x=0; x<dimMsk.m_x; x++)
			{
				float v;
				mskF.read((char*)&v, sizeof(float));
				pLyrMsk->CellRef(x, y, z) = v;
			}
		}
	}
	mskF.close();
	cout << "read mask ok." << endl;
	return pLyrMsk;
}
/*
void ReadMsk(float *pMsk, string fName, Vect3D<unsigned> dimShd)
{
		ifstream mskF(fName+".msk", ios::binary);
		cout << "read " << fName << ".msk" << endl;

		Vect3D<unsigned> dimMskIn(0, 0, 0);
		mskF.read((char*)&dimMskIn.m_x, sizeof(unsigned));
		mskF.read((char*)&dimMskIn.m_y, sizeof(unsigned));
		mskF.read((char*)&dimMskIn.m_z, sizeof(unsigned));

		Layer *pLyrMskIn = new Layer(dimMskIn.m_x, dimMskIn.m_y, dimMskIn.m_z);
		for(unsigned z=0; z<dimMskIn.m_z; z++)
		{
			for(unsigned y=0; y<dimMskIn.m_y; y++)
			{
				for(unsigned x=0; x<dimMskIn.m_x; x++)
				{
					float v;
					mskF.read((char*)&v, sizeof(float));
					pLyrMskIn->CellRef(x, y, z) = (DATA)v;
				}
			}
		}
		mskF.close();

		Layer *pLyrMskScl = new Layer(dimShd.m_x, dimShd.m_y, dimShd.m_z);
		lyrOp.scaleDim.Gen(*pLyrMskScl, *pLyrMskIn);
		delete pLyrMskIn;

		for(unsigned z=0; z<dimShd.m_z; z++)
		{
			for(unsigned y=0; y<dimShd.m_y; y++)
			{
				for(unsigned x=0; x<dimShd.m_x; x++)
				{
					unsigned no = z*dimShd.m_x*dimShd.m_y + y*dimShd.m_x + x;
					pMsk[no] = (float)pLyrMskScl->CellVal(x, y, z);
				}
			}
		}
		delete pLyrMskScl;

		cout << "read mask ok." << endl;
}
*/

void GetLocZ(Vect3D<float> &vLoc, Vect3D<float> vPlane, Vect3D<unsigned> dimP)
{
	vLoc.m_x = vPlane.m_x / (dimP.m_x-1.f);
	vLoc.m_y = vPlane.m_y / (dimP.m_y-1.f);
	vLoc.m_z = vPlane.m_z / (dimP.m_z-1.f);
}
void GetLocY(Vect3D<float> &vLoc, Vect3D<float> vPlane, Vect3D<unsigned> dimP)
{
	vLoc.m_x = vPlane.m_x / (dimP.m_x-1.f);
	vLoc.m_y = vPlane.m_z / (dimP.m_z-1.f);
	vLoc.m_z = vPlane.m_y / (dimP.m_y-1.f);
}
void GetLocX(Vect3D<float> &vLoc, Vect3D<float> vPlane, Vect3D<unsigned> dimP)
{
	vLoc.m_x = vPlane.m_z / (dimP.m_z-1.f);
	vLoc.m_y = vPlane.m_x / (dimP.m_x-1.f);
	vLoc.m_z = vPlane.m_y / (dimP.m_y-1.f);
}

unsigned GetShdLocZ(Vect3D<unsigned> vPlane, Vect3D<unsigned> dim)
{
	return vPlane.m_z*dim.m_x*dim.m_y + vPlane.m_y*dim.m_x + vPlane.m_x;
}
unsigned GetShdLocY(Vect3D<unsigned> vPlane, Vect3D<unsigned>dim)
{
	return vPlane.m_y*dim.m_x*dim.m_y + vPlane.m_z*dim.m_x + vPlane.m_x;
}
unsigned GetShdLocX(Vect3D<unsigned> vPlane, Vect3D<unsigned>dim)
{
	return vPlane.m_y*dim.m_x*dim.m_y + vPlane.m_x*dim.m_x + vPlane.m_z;
}

void MyVolume::MakeShadow(DATA mskThd)
{
	unsigned shdLen = 200;
	unsigned avgLen = 1;
	DATA devG = 1.8F;

	Vect3D<DATA> dimRatio(m_vDim.m_x*m_vStep.m_x,
						  m_vDim.m_y*m_vStep.m_y,
						  m_vDim.m_z*m_vStep.m_z);	
	DATA ratioMax = FindMax(dimRatio);
	Mult(dimRatio, 1.F/ratioMax);

	Vect3D<unsigned> dimShd((unsigned)(shdLen*dimRatio.m_x + 0.5F),
							(unsigned)(shdLen*dimRatio.m_y + 0.5F),
							(unsigned)(shdLen*dimRatio.m_z + 0.5F));
	cout << "dimShd: " << dimShd.m_x << " " << dimShd.m_y << " " << dimShd.m_z << endl;

	unsigned size = dimShd.m_x * dimShd.m_y * dimShd.m_z;
	if(!m_pShadow)
	{
		delete []m_pShadow;
		m_pShadow = new float[size*3];
	}
	else {}

	bool bRead = false; //false; 
	string volName = m_fName.substr(m_fName.find_last_of("/")+1);
	string shdFName = volName+"_out.shd";
	Vect3D<unsigned> dimMsk(0, 0, 0);
	if(!bRead)
	{
		time_t stT = time(0);

		for(unsigned i=0; i<size*3; i++)
		{
			m_pShadow[i] = 1.f;
		}

		//ReadMsk(m_pMsk, m_fName, dimShd);
		Layer *pLyrMsk = ReadMsk(m_fName);

		dimMsk = pLyrMsk->GetDim();
		delete []m_pMsk;
		m_pMsk = new float[dimMsk.m_x * dimMsk.m_y * dimMsk.m_z];
		pLyrMsk->CopyTo_zLast(m_pMsk);

		Layer *pLyrSM = new Layer(dimShd.m_x, dimShd.m_y, dimShd.m_z);
		lyrOp.scaleDim.Gen(*pLyrSM, *pLyrMsk);
		delete pLyrMsk;
		float *pSM = new float[dimShd.m_x * dimShd.m_y * dimShd.m_z];
		pLyrSM->CopyTo_zLast(pSM);
		delete pLyrSM;

		unsigned maxLen = (dimShd.m_x>dimShd.m_y)? dimShd.m_x: dimShd.m_y;
		if(dimShd.m_z > maxLen)	{maxLen = dimShd.m_z;}
		else {}
		Mtx *pMtxP = new Mtx(maxLen, maxLen);
		Mtx *pMtxN = new Mtx(maxLen, maxLen);
		Mtx *pMtxF = new Mtx(maxLen, maxLen);

		Mtx *pMtxG = new Mtx(avgLen, avgLen);
		mtxOp.Gauss.Gen(*pMtxG, devG, devG, true);

		//***************************************

		unsigned aPlaneX[] = {dimShd.m_x, dimShd.m_x, dimShd.m_y};
		unsigned aPlaneY[] = {dimShd.m_y, dimShd.m_z, dimShd.m_z};
		unsigned aPlaneZ[] = {dimShd.m_z, dimShd.m_y, dimShd.m_x};
		void (*GetLoc[])(Vect3D<float>&, Vect3D<float>, Vect3D<unsigned>) = {GetLocZ, GetLocY, GetLocX};
		unsigned (*GetShdLoc[])(Vect3D<unsigned>, Vect3D<unsigned>) = {GetShdLocZ, GetShdLocY, GetShdLocX};
		for(unsigned i=0; i<3; i++)
		{
			cout << "basis " << i << ":" << endl;
			Vect3D<unsigned> dimP(aPlaneX[i], aPlaneY[i], aPlaneZ[i]);
			for(unsigned z=0; z<aPlaneZ[i]; z++)
			{
				cout << z << " ";
				float zN = (aPlaneZ[i]-1.f) / aPlaneZ[i] * z;
				float zF = (aPlaneZ[i]-1.f) / aPlaneZ[i] * (z+1.f);
				for(unsigned y=0; y<aPlaneY[i]; y++)
				{ 
					for(unsigned x=0; x<aPlaneX[i]; x++)
					{
						Vect3D<float> vNPlane(x, y, zN);
						Vect3D<float> vNLoc(0, 0, 0);
						GetLoc[i](vNLoc, vNPlane, dimP);
						//cout << "n: " << vNLoc.m_x << " " << vNLoc.m_y << " " << vNLoc.m_z << endl;
						float nnVal;
						InterpolateTex3d(&nnVal, m_pData, m_vDim.m_x, m_vDim.m_y, m_vDim.m_z, 1, 
										vNLoc.m_x, vNLoc.m_y, vNLoc.m_z);
									
						Vect3D<float> vFPlane(x, y, zF);
						Vect3D<float> vFLoc(0, 0, 0);
						GetLoc[i](vFLoc, vFPlane, dimP);
						//cout << "f: " << vFLoc.m_x << " " << vFLoc.m_y << " " << vFLoc.m_z << endl;
						float ffVal;
						InterpolateTex3d(&ffVal, m_pData, m_vDim.m_x, m_vDim.m_y, m_vDim.m_z, 1, 
										vFLoc.m_x, vFLoc.m_y, vFLoc.m_z);

						float aCol[4];
						InterpolateTex2d(aCol, m_pPreint, m_vPreintDim.m_x, m_vPreintDim.m_y, 4, 
										nnVal, ffVal);

						pMtxP->CellRef(x, y) = aCol[3];
						pMtxN->CellRef(x, y) = nnVal;
						pMtxF->CellRef(x, y) = ffVal;
					} // x
				} // y

				for(unsigned y=0; y<aPlaneY[i]; y++)
				{
					for(unsigned x=0; x<aPlaneX[i]; x++)
					{
						//int yB = (int)y - avgLen/2;	int yT = yB + avgLen - 1;
						//int xL = (int)x - avgLen/2;	int xR = xL + avgLen - 1;
						//if(yB < 0)	yB = 0;	else {}		if(yT >= (int)dimP.m_y)	yT = dimP.m_y - 1;	else {}
						//if(xL < 0)	xL = 0;	else {}		if(xR >= (int)dimP.m_x)	xR = dimP.m_x - 1;	else {}

						DATA tVal = 0;
						DATA nVal = 0;
						DATA fVal = 0;
						DATA avgSize = 0;
						//for(unsigned yy=(unsigned)yB; yy<=(unsigned)yT; yy++)
						for(unsigned yy=0; yy<avgLen; yy++)
						{
							int yLoc = (int)y - avgLen/2 + yy;
							if(yLoc<0 || yLoc>=(int)dimP.m_y) continue;
							else {}
							//for(unsigned xx=(unsigned)xL; xx<=(unsigned)xR; xx++)
							for(unsigned xx=0; xx<avgLen; xx++)
							{
								int xLoc = (int)x - avgLen/2 + xx;
								if(xLoc<0 || xLoc>=(int)dimP.m_x) continue;
								else {}

								tVal += pMtxP->CellVal(xLoc, yLoc)*pMtxG->CellVal(xx, yy);
								nVal += pMtxN->CellVal(xLoc, yLoc)*pMtxG->CellVal(xx, yy);
								fVal += pMtxF->CellVal(xLoc, yLoc)*pMtxG->CellVal(xx, yy);
								avgSize += pMtxG->CellVal(xx, yy);
							} // xx
						} // yy
						tVal /= avgSize;
						nVal /= avgSize;
						fVal /= avgSize;
//if(nVal > 0.1)
//cout << nVal << " " << fVal << ",  ";

						Vect3D<unsigned> locNow(x, y, z);
						unsigned shdLoc = GetShdLoc[i](locNow, dimShd);
						
						float preVal = 0;
						if(z > 0)
						{
							Vect3D<unsigned> locPre(x, y, z-1);
							unsigned shdPre = GetShdLoc[i](locPre, dimShd);
							preVal = m_pShadow[shdPre*3+2-i];
						}
						else {}
						
						//float nMsk = m_pId[(unsigned)(tSum*m_colDim+0.5f)];
						//float fMsk = nMsk;
						float nMsk = m_pId[(unsigned)(nVal*m_colDim+0.5f)];
						float fMsk = m_pId[(unsigned)(fVal*m_colDim+0.5f)];
						float mskVal = 1.f;
						if(pSM[shdLoc] == 0)	mskVal = 0;
						else {}
						if(pSM[shdLoc]==nMsk||pSM[shdLoc]==fMsk)	mskVal = 0;
						else {}
						m_pShadow[shdLoc*3+2-i] = preVal + (float)tVal*mskVal; //aCol[3]*mskVal;						
					} // x
				} // y
			} // z
			cout << endl;
		} // i	
		delete []pSM;
		delete pMtxP;
		delete pMtxN;
		delete pMtxF;

		float shdMax = 0;
		for(unsigned z=0; z<dimShd.m_z; z++)
		{
			for(unsigned y=0; y<dimShd.m_y; y++)
			{
				for(unsigned x=0; x<dimShd.m_x; x++)
				{
					unsigned loc = z*dimShd.m_x*dimShd.m_y + y*dimShd.m_x + x;
					for(unsigned c=0; c<3; c++)
					{
						if(m_pShadow[loc*3+c] > shdMax)
						{
							shdMax = m_pShadow[loc*3+c];
						}
						else {}
					}
				}
			}
		}
		cout << "shdMax: " << shdMax << endl;
		for(unsigned i=0; i<size; i++)
		{
			m_pShadow[i*3]  /= shdMax;
			m_pShadow[i*3+1] /= shdMax;
			m_pShadow[i*3+2] /= shdMax;
		}

		time_t endT = time(0);
		double diffT = difftime(endT, stT);
		cout << "time: " << diffT << " secs." << endl;

		ofstream shdF(shdFName.c_str(), ios::binary);
		shdF.write((char*)&dimShd.m_x, sizeof(unsigned));
		shdF.write((char*)&dimShd.m_y, sizeof(unsigned));
		shdF.write((char*)&dimShd.m_z, sizeof(unsigned));
		shdF.write((char*)m_pShadow, sizeof(float)*dimShd.m_x*dimShd.m_y*dimShd.m_z*3);
		shdF.close();
	} // bRead
	else
	{
		Layer *pLyrMsk = ReadMsk(m_fName);
		dimMsk = pLyrMsk->GetDim();
		delete []m_pMsk;
		m_pMsk = new float[dimMsk.m_x * dimMsk.m_y * dimMsk.m_z];
		pLyrMsk->CopyTo_zLast(m_pMsk);
		delete pLyrMsk;

		ifstream shdF(shdFName.c_str(), ios::binary);
		shdF.read((char*)&dimShd.m_x, sizeof(unsigned));
		shdF.read((char*)&dimShd.m_y, sizeof(unsigned));
		shdF.read((char*)&dimShd.m_z, sizeof(unsigned));
		shdF.read((char*)m_pShadow, sizeof(float)*dimShd.m_x*dimShd.m_y*dimShd.m_z*3);
		shdF.close();
	}

	m_pTexShadow = new Tex3D(m_pShadow, dimShd.m_x, dimShd.m_y, dimShd.m_z, 3);
	m_pTexMsk = new Tex3D(m_pMsk, dimMsk.m_x, dimMsk.m_y, dimMsk.m_z, 1, false, false, false);
}

//*************************************************************************************************

void CreateParallelPMtx_light(DATA aaPM[][4], DATA left, DATA right, DATA bottom, DATA top, DATA nearZ, DATA farZ)
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

void RotateObj_light()
{
	DATA aBasis[] = {aaObjVec[0][0], aaObjVec[0][1], aaObjVec[0][2], aaObjVec[0][3],
					 aaObjVec[1][0], aaObjVec[1][1], aaObjVec[1][2], aaObjVec[1][3],
					 aaObjVec[2][0], aaObjVec[2][1], aaObjVec[2][2], aaObjVec[2][3],
					 0, 0, 0, 1.F};
	glMultMatrixd(aBasis);
}

void GetCont(Mtx &mtxIn)
{
	static Mtx *apMtxGauss[] = {0, 0};

	Vect2D<unsigned> dim = mtxIn.GetDim(); 

	if(!apMtxGauss[0])	apMtxGauss[0] = new Mtx(20, 20);	else {}
	if(!apMtxGauss[1])	apMtxGauss[1] = new Mtx(20, 20);	else {}
	
	DATA coeSharp = 1.F;
	mtxOp.DoG.GenNPR(mtxIn, apMtxGauss, 2.F, coeSharp);
	mtxOp.sub.Gen(1.F, mtxIn);
}

bool bMap = false;
void MyVolume::Draw(DATA aaCPlane[][4], unsigned w, unsigned h, unsigned runNo, unsigned fusNo, DATA mskThd)
{
	static unsigned winW = 0;
	static unsigned winH = 0;
	bool bText = false;
	if(winW!=w || winH!= h)
	{
		bText = true;
	}
	else {}

	static float *pStrPos = 0;					static Tex2D *pTexStrPos = 0;
	static float *pLengthImg = 0;				static Tex2D *pTexLengthImg = 0;
	static float *apRcast[] = {0, 0};			static Tex2D *apTexRcast[] = {0, 0};
	static float *apRfuse[] = {0, 0};			static Tex2D *apTexRfuse[] = {0, 0};
	static float *apImg[] = {0, 0, 0, 0, 0};	static Tex2D *apTexImg[] = {0, 0, 0, 0, 0};
	static MyImg *pImgW = 0; 
	static Layer *pLyrFrame = 0;				static Mtx *pMtxGray = 0;	
	if(bText)
	{
		winW = w;
		winH = h;
		unsigned texSize = w * h * 4;
		
		delete []pStrPos;	pStrPos = new float[texSize];
		delete pTexStrPos;	pTexStrPos = new Tex2D(pStrPos, w, h, 4);

		delete []pLengthImg;	pLengthImg = new float[texSize];
		delete pTexLengthImg;	pTexLengthImg = new Tex2D(pLengthImg, w, h, 4);
		
		for(unsigned i=0; i<2; i++)
		{
			delete []apRcast[i];	apRcast[i] = new float[texSize];
			delete apTexRcast[i];	apTexRcast[i] = new Tex2D(apRcast[i], w, h, 4);

			delete []apRfuse[i];	apRfuse[i] = new float[texSize];
			delete apTexRfuse[i];	apTexRfuse[i] = new Tex2D(apRfuse[i], w, h, 4);
		}
		for(unsigned i=0; i<5; i++)
		{
			delete []apImg[i];		apImg[i] = new float[texSize];
			delete apTexImg[i];		apTexImg[i] = new Tex2D(apImg[i], w, h, 4);
		}

		pImgW = new  MyImg(w, h, 4);
		
		pLyrFrame = new Layer(w, h, 4);
		pMtxGray = new Mtx(w, h);
		
		srand( (unsigned)time(0) );
	}
	else {}

	string volName = m_fName.substr(m_fName.find_last_of('/')+1);
	//cout << volName << endl;
	bool bReadSte = false;
	static bool bFlip = false;
	static bool bRIni = false;
	if(bReadSte)
	{
		if(!bRIni)
		{
			for(unsigned i=0; i<2; i++)
			{
				string str;
				stringstream ss;
				ss.clear();
				ss << i;
				ss >> str;

				imgIO.Read(*pImgW, volName+"_"+str+".bmp");
				lyrOp.mul.Gen(*pImgW, INV_255);
				apTexImg[i]->Read(*pImgW);
				apTexImg[i]->Update();
			}
			imgIO.Read(*pImgW, volName+"_trad.bmp");
			lyrOp.mul.Gen(*pImgW, INV_255);
			apTexImg[2]->Read(*pImgW);
			apTexImg[2]->Update();
			for(unsigned i=3; i<=4; i++)
			{
				string str;
				stringstream ss;
				ss.clear();
				ss << i;
				ss >> str;

				imgIO.Read(*pImgW, volName+"_"+str+".bmp");
				lyrOp.mul.Gen(*pImgW, INV_255);
				apTexImg[i]->Read(*pImgW);
				apTexImg[i]->Update();
			}

			bRIni = true;
		}
		else {}

		if(bStereo)
		{
			bFlip = bFusion;
			bFusion = false;
			bStereo = false;	
		}
		else {}
		bMap = true;
	}
	else {}

	static int mapNo = 0;
	if(bMap)
	{
		glUseProgram(pFusionSh->GetProg());

		//glActiveTexture(GL_TEXTURE0);
		//glBindTexture(GL_TEXTURE_2D, apTexImg[mapNo]->GetTexID());
		//glUniform1i(pFusionSh->GetUniLoc("smpImg"), 0);
		glUniform2f(pFusionSh->GetUniLoc("v2WinDim"), (float)w, (float)h);

		if(bFlip)
		{
			glActiveTexture(GL_TEXTURE0);
			if(!bReadSte || bBrightness_div)	glBindTexture(GL_TEXTURE_2D, apTexImg[0]->GetTexID());
			else								glBindTexture(GL_TEXTURE_2D, apTexImg[3]->GetTexID());
			glUniform1i(pFusionSh->GetUniLoc("smpImg0"), 0);
			glActiveTexture(GL_TEXTURE1);
			if(!bReadSte || bBrightness_div)	glBindTexture(GL_TEXTURE_2D, apTexImg[1]->GetTexID());
			else								glBindTexture(GL_TEXTURE_2D, apTexImg[4]->GetTexID());
			glUniform1i(pFusionSh->GetUniLoc("smpImg1"), 1);
		}
		else
		{	glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, apTexImg[2]->GetTexID());
			glUniform1i(pFusionSh->GetUniLoc("smpImg0"), 0);
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, apTexImg[2]->GetTexID());
			glUniform1i(pFusionSh->GetUniLoc("smpImg1"), 1);
		}
		
		static int no = 0;
		no = (no + 1) % 2;
		static float rndRatio; 
		if(!bRndDis)
		{
			rndRatio = (no==0)? 0.8f: 0.2f;
		}
		else
		{
			rndRatio = (no==0)? 1.f: 0;
		}
		glUniform1f(pFusionSh->GetUniLoc("rndRatio"), rndRatio);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glDisable(GL_DEPTH_TEST);
		glBegin(GL_QUADS);
			glVertex3dv(aaCPlane[0]);
			glVertex3dv(aaCPlane[1]);
			glVertex3dv(aaCPlane[2]);
			glVertex3dv(aaCPlane[3]);
		glEnd();
		glPopMatrix();

		mapNo++;
		if(mapNo >= 2)
		{
			mapNo = 0;
		}
		else {}
		return;
	}
	else {}

	//*********************************************************************************************
	// photon shooting
	//*********************************************************************************************

	if(!m_pShadow)
	{
		MakeShadow(mskThd);
	}
	else {}

	Vect3D<DATA> dimRatio(m_vDim.m_x*m_vStep.m_x,
						  m_vDim.m_y*m_vStep.m_y,
						  m_vDim.m_z*m_vStep.m_z);	
	DATA ratioMax = FindMax(dimRatio);
	Mult(dimRatio, 1.F/ratioMax);

	Vect3D<float> vScl((float)dimRatio.m_x, (float)dimRatio.m_y, (float)dimRatio.m_z);
	if(bStereo)
	{
		float sScl = 1.f; //0.8f;
		vScl.m_x *= sScl;
		vScl.m_y *= sScl;
		vScl.m_z *= sScl;
	}
	else {}
	glTranslatef(0.5f, 0.5f, 0.5f);
	glScalef(vScl.m_x, vScl.m_y, vScl.m_z);
	glTranslatef(-0.5f, -0.5f, -0.5f);

	unsigned aaCubeIdx[][4] = {{0, 1, 2, 3}, {1, 5, 6, 2}, {5, 4, 7, 6}, {0, 4, 7, 3}, {3, 2, 6, 7}, {4, 5, 1, 0}};
	Vect3D<DATA> blckStLen(1.F/m_vBlckDim.m_x, 1.F/m_vBlckDim.m_y, 1.F/m_vBlckDim.m_z);

	glUseProgram(0);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glClear(GL_DEPTH_BUFFER_BIT);

	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
	
	unsigned blckPos = 0;
	DATA zN = 0;
	DATA zF = blckStLen.m_z;
	for(unsigned z=0; z<m_vBlckDim.m_z; z++)
	{
		DATA yB = 0;
		DATA yT = blckStLen.m_y;
		for(unsigned y=0; y<m_vBlckDim.m_y; y++)
		{
			DATA xL = 0;
			DATA xR = blckStLen.m_x;
			for(unsigned x=0; x<m_vBlckDim.m_x; x++)
			{	
				if(m_pBAct[blckPos] != 0)
				{
					DATA aaV[][4] = {{xL, yB, zN, 1.F}, {xR, yB, zN, 1.F}, {xR, yT, zN, 1.F}, {xL, yT, zN, 1.F},
									 {xL, yB, zF, 1.F}, {xR, yB, zF, 1.F}, {xR, yT, zF, 1.F}, {xL, yT, zF, 1.F}};

					glBegin(GL_QUADS);
					for(unsigned i=0; i<6; i++)
					{
						glColor4dv(aaV[ aaCubeIdx[i][0] ]);
						glVertex3dv(aaV[ aaCubeIdx[i][0] ]);
						glColor4dv(aaV[ aaCubeIdx[i][1] ]);
						glVertex3dv(aaV[ aaCubeIdx[i][1] ]);
						glColor4dv(aaV[ aaCubeIdx[i][2] ]);
						glVertex3dv(aaV[ aaCubeIdx[i][2] ]);
						glColor4dv(aaV[ aaCubeIdx[i][3] ]);
						glVertex3dv(aaV[ aaCubeIdx[i][3] ]);
					}
					glEnd();
				}
				//else {}
				blckPos++;
			
				xL = xR;
				xR += blckStLen.m_x;
			} // x
			yB = yT;
			yT += blckStLen.m_y;
		} // y
		zN = zF;
		zF += blckStLen.m_z;
	} // z

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, pTexStrPos->GetTexID()); //texStrPos);
	glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, w, h, 0);
	glDisable(GL_TEXTURE_2D);
	
	//*******************************************
	// varify the result
	//*******************************************
	
	/*
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	gluLookAt(aEyePos[0], aEyePos[1], aEyePos[2], 
			  aEyePos[0]-aaEyeVec[2][0], aEyePos[1]-aaEyeVec[2][1], aEyePos[2]-aaEyeVec[2][2], 
			  aaEyeVec[1][0], aaEyeVec[1][1], aaEyeVec[1][2]);

	glUseProgram(0);
	
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texStrPos);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glDisable(GL_DEPTH_TEST);
	//glColor4f(0, 0, 1.f, 1.f);
	glBegin(GL_QUADS);
		glTexCoord2f(0, 0);
		glVertex3dv(aaCPlane[0]);
		glTexCoord2f(1.f, 0);
		glVertex3dv(aaCPlane[1]);
		glTexCoord2f(1.f, 1.f);
		glVertex3dv(aaCPlane[2]);
		glTexCoord2f(0, 1.f);
		glVertex3dv(aaCPlane[3]);
	glEnd();

	glPopMatrix();
	glDisable(GL_TEXTURE_2D);
	return;
	*/

	//*********************************************************************************************
	// ray length image
	//*********************************************************************************************

	glUseProgram(pLengthImgSh->GetProg());

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, pTexStrPos->GetTexID()); //texStrPos);
	glUniform1i(pLengthImgSh->GetUniLoc("smpStrPos"), 0);
	glUniform2f(pLengthImgSh->GetUniLoc("v2WinDim"), (float)w, (float)h);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_GREATER);

	blckPos = 0;
	zN = 0;
	zF = blckStLen.m_z;
	for(unsigned z=0; z<m_vBlckDim.m_z; z++)
	{
		DATA yB = 0;
		DATA yT = blckStLen.m_y;
		for(unsigned y=0; y<m_vBlckDim.m_y; y++)
		{
			DATA xL = 0;
			DATA xR = blckStLen.m_x;
			for(unsigned x=0; x<m_vBlckDim.m_x; x++)
			{	
				if(m_pBAct[blckPos] != 0)
				{
					DATA aaV[][4] = {{xL, yB, zN, 1.F}, {xR, yB, zN, 1.F}, {xR, yT, zN, 1.F}, {xL, yT, zN, 1.F},
									 {xL, yB, zF, 1.F}, {xR, yB, zF, 1.F}, {xR, yT, zF, 1.F}, {xL, yT, zF, 1.F}};

					glBegin(GL_QUADS);
					for(unsigned i=0; i<6; i++)
					{
						glColor4dv(aaV[ aaCubeIdx[i][0] ]);
						glVertex3dv(aaV[ aaCubeIdx[i][0] ]);
						glColor4dv(aaV[ aaCubeIdx[i][1] ]);
						glVertex3dv(aaV[ aaCubeIdx[i][1] ]);
						glColor4dv(aaV[ aaCubeIdx[i][2] ]);
						glVertex3dv(aaV[ aaCubeIdx[i][2] ]);
						glColor4dv(aaV[ aaCubeIdx[i][3] ]);
						glVertex3dv(aaV[ aaCubeIdx[i][3] ]);
					}
					glEnd();
				}
				//else {}
				blckPos++;
			
				xL = xR;
				xR += blckStLen.m_x;
			} // x
			yB = yT;
			yT += blckStLen.m_y;
		} // y
		zN = zF;
		zF += blckStLen.m_z;
	} // z

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, pTexLengthImg->GetTexID());
	glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, w, h, 0);
	glDisable(GL_TEXTURE_2D);

	//*********************************************************************************************
	// ray casting
	//*********************************************************************************************

	Shader *pSh = 0;
	if(!bFusion)
	{
		pSh = pRCastingSh;
		if(bDivShow)
		{
			pSh = pRc2Sh;
		}
		else
		{
			pSh = pRCastingSh;
		}
	}
	else if(fusNo == 0) //if(fusNo == 1)
	{
		pSh = pRCastingSh;
	}
	else if(fusNo == 1) //if(fusNo == 0)
	{
		pSh = pRc2Sh;
	}
	else
	{
		assert(0);
	}

	glUseProgram(pSh->GetProg());

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, pTexStrPos->GetTexID()); //texStrPos);
	glUniform1i(pSh->GetUniLoc("smpStrPos"), 0);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, pTexLengthImg->GetTexID());
	glUniform1i(pSh->GetUniLoc("smpLengthImg"), 1);

	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_3D, m_pTexData->GetTexID());
	glUniform1i(pSh->GetUniLoc("smpData"), 3);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_3D, m_pTexAct->GetTexID());
	glUniform1i(pSh->GetUniLoc("smpAct"), 2);

	glActiveTexture(GL_TEXTURE5);
	glBindTexture(GL_TEXTURE_3D, m_pTexShadow->GetTexID());
	glUniform1i(pSh->GetUniLoc("smpPhot"), 5);

	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_2D, m_pTexPreint->GetTexID());
	glUniform1i(pSh->GetUniLoc("smpColor"), 4);

	glActiveTexture(GL_TEXTURE6);
	glBindTexture(GL_TEXTURE_1D, m_pTexId->GetTexID());
	glUniform1i(pSh->GetUniLoc("smpId"), 6);
	
	glActiveTexture(GL_TEXTURE7);
	glBindTexture(GL_TEXTURE_3D, m_pTexMsk->GetTexID());
	glUniform1i(pSh->GetUniLoc("smpMsk"), 7);

	glUniform2f(pSh->GetUniLoc("v2WinDim"), (float)w, (float)h);
	glUniform3f(pSh->GetUniLoc("v3BlckRatio"), m_vBlckRatio.m_x, m_vBlckRatio.m_y, m_vBlckRatio.m_z);

	glUniform1f(pSh->GetUniLoc("disAtt"), (float)disAtt);
	if(pSh != pRc2Sh)	glUniform1f(pSh->GetUniLoc("opaScl"), (float)opaScl);
	else				glUniform1f(pSh->GetUniLoc("opaScl"), (float)(opaScl*divThrd));
	glUniform1f(pSh->GetUniLoc("briScl"), (float)briScl);
	glUniform1f(pSh->GetUniLoc("smpNum"), (float)smpNum);
	//glUniform1f(pSh->GetUniLoc("mskThrd"), (float)mskThrd);

	float fMIP = (bMIP)? 1.f: 0;
	glUniform1f(pSh->GetUniLoc("bMIP"), fMIP);
	glUniform1f(pSh->GetUniLoc("occScl"), (float)occScl);

	float fXClip = (bXClip)? 1.f: 0;
	float fXFront= (bXFront)? 1.f: 0;
	glUniform1f(pSh->GetUniLoc("bXClip"), fXClip);
	glUniform1f(pSh->GetUniLoc("bXFront"), fXFront);
	glUniform1f(pSh->GetUniLoc("xPlane"), (float)xPlane);
	float fYClip = (bYClip)? 1.f: 0;
	float fYFront= (bYFront)? 1.f: 0;
	glUniform1f(pSh->GetUniLoc("bYClip"), fYClip);
	glUniform1f(pSh->GetUniLoc("bYFront"), fYFront);
	glUniform1f(pSh->GetUniLoc("yPlane"), (float)yPlane);
	float fZClip = (bZClip)? 1.f: 0;
	float fZFront= (bZFront)? 1.f: 0;
	glUniform1f(pSh->GetUniLoc("bZClip"), fZClip);
	glUniform1f(pSh->GetUniLoc("bZFront"), fZFront);
	glUniform1f(pSh->GetUniLoc("zPlane"), (float)zPlane);
	//glUniform3f(pSh->GetUniLoc("v3DimRatio"), (float)dimRatio.m_x, (float)dimRatio.m_y, (float)dimRatio.m_z);

	float bBG = (bStereo)? 0: 1.f; 
	glUniform1f(pSh->GetUniLoc("bBG"), bBG);

	static float aView[3];
	static float aViewX[3];
	static float aViewY[3];
	if(runNo == 0)
	{		
		aView[0] = (float)(aaEyeVec[2][0]*aaObjVec[0][0] + aaEyeVec[2][1]*aaObjVec[0][1] + aaEyeVec[2][2]*aaObjVec[0][2]);
		aView[1] = (float)(aaEyeVec[2][0]*aaObjVec[1][0] + aaEyeVec[2][1]*aaObjVec[1][1] + aaEyeVec[2][2]*aaObjVec[1][2]);
		aView[2] = (float)(aaEyeVec[2][0]*aaObjVec[2][0] + aaEyeVec[2][1]*aaObjVec[2][1] + aaEyeVec[2][2]*aaObjVec[2][2]);

		aViewX[0] = (float)(aaEyeVec[0][0]*aaObjVec[0][0] + aaEyeVec[0][1]*aaObjVec[0][1] + aaEyeVec[0][2]*aaObjVec[0][2]);
		aViewX[1] = (float)(aaEyeVec[0][0]*aaObjVec[1][0] + aaEyeVec[0][1]*aaObjVec[1][1] + aaEyeVec[0][2]*aaObjVec[1][2]);
		aViewX[2] = (float)(aaEyeVec[0][0]*aaObjVec[2][0] + aaEyeVec[0][1]*aaObjVec[2][1] + aaEyeVec[0][2]*aaObjVec[2][2]);

		aViewY[0] = (float)(aaEyeVec[1][0]*aaObjVec[0][0] + aaEyeVec[1][1]*aaObjVec[0][1] + aaEyeVec[1][2]*aaObjVec[0][2]);
		aViewY[1] = (float)(aaEyeVec[1][0]*aaObjVec[1][0] + aaEyeVec[1][1]*aaObjVec[1][1] + aaEyeVec[1][2]*aaObjVec[1][2]);
		aViewY[2] = (float)(aaEyeVec[1][0]*aaObjVec[2][0] + aaEyeVec[1][1]*aaObjVec[2][1] + aaEyeVec[1][2]*aaObjVec[2][2]);
	}
	else {}
	glUniform3f(pSh->GetUniLoc("v3View1"),  (float)-aView[0],  (float)-aView[1],  (float)-aView[2]);
	glUniform3f(pSh->GetUniLoc("v3ViewX1"), (float)-aViewX[0], (float)-aViewX[1], (float)-aViewX[2]);
	//glUniform3f(pSh->GetUniLoc("v3ViewY1"), (float)-aViewY[0], (float)-aViewY[1], (float)-aViewY[2]);

	if(pSh == pRc2Sh) {}
	else {}
	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glDisable(GL_DEPTH_TEST);
	glColor4f(0, 0, 1.f, 1.f);
	glBegin(GL_QUADS);
		glVertex3dv(aaCPlane[0]);
		glVertex3dv(aaCPlane[1]);
		glVertex3dv(aaCPlane[2]);
		glVertex3dv(aaCPlane[3]);
	glEnd();

	static DATA avgV = 0;
	if(bStereo)
	{
		if(!bFusion)
		{			
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, apTexRcast[runNo]->GetTexID());
			glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, w, h, 0);
			glDisable(GL_TEXTURE_2D);
		}
		else 
		{
			int aWinLow[] = {0, 0};
			Tex2D *pTexRef = (fusNo==0)? apTexRfuse[runNo]: apTexRcast[runNo];
			pTexRef->ReadFrame(aWinLow);
			Vect2D<unsigned> dim = pTexRef->GetDim();		

			if(fusNo == 0)
			{
				avgV = 0;	
				DATA wSum = 0;
				for(unsigned y=0; y<dim.m_y; y++)
				{
					for(unsigned x=0; x<dim.m_x; x++)
					{						
						DATA aRgba[] = {
							pTexRef->GetCell(x, y, 0),
							pTexRef->GetCell(x, y, 1),
							pTexRef->GetCell(x, y, 2),
							pTexRef->GetCell(x, y, 3)};

						DATA val = myMath.RGB2Gray(aRgba);
						DATA ww = aRgba[3];
						avgV += log( val*255.F + 1.F) * ww;
						wSum += ww;
					} // x
				} // y
				avgV /= wSum;
				cout << avgV << " avg0" << endl;
				avgV = (exp(avgV) - 1.F) / 255.F; 
				cout << avgV << " normal0" << endl;	
			} // funNo
			else {}

			//***********************************
			
			for(unsigned y=0; y<dim.m_y; y++)
			{
				for(unsigned x=0; x<dim.m_x; x++)
				{
					DATA aRgba[] = {
						pTexRef->GetCell(x, y, 0),
						pTexRef->GetCell(x, y, 1),
						pTexRef->GetCell(x, y, 2),
						pTexRef->GetCell(x, y, 3)};
					
					for(unsigned i=0; i<3; i++)
					{
						DATA val = aRgba[i];

						if(fusNo == 1)
						{
							if(bContrast_div) {val = pow(val, contrastVal);}
							else {}

							if(bDim_div) {val *= dimVal;}
							else {}
						}
						else {}

						if(bBrightness_div)
						{ 
							DATA valOld = val;
							val = brightnessVal * val / avgV;
							val = (val + valOld*valOld/4.F) / (1.F+val);
							val *= 1.5F;
						}
						else {}
						
						if(val < 0)			{val = 0;}
						else if(val > 1.F)	{val = 1.F;}
						else {}

						aRgba[i] = val;
					}

					for(unsigned c=0; c<3; c++)
					{
						pTexRef->SetCell((float)aRgba[c], x, y, c);
					}
				}
			}
			pTexRef->Update();

			string strFus;
			stringstream ss;
			ss.clear();
			ss << fusNo;
			ss >> strFus;
						
			string strRun;
			ss.clear();
			ss << runNo;
			ss >> strRun;

			string fName = volName + "_" + strFus + "_" + strRun + ".bmp";

			pTexRef->GetCell(*pImgW);
			lyrOp.mul.Gen(*pImgW, 255.F);
			imgIO.Write(fName, *pImgW);					
		} // bFusion
	}
	else {} // bStereo
	
	if(bStereo && runNo==1)
	{
		Tex2D *apTex[] = {0, 0};
		if(!bFusion)			{apTex[0] = apTexRcast[0];	apTex[1] = apTexRcast[1];}
		else if(fusNo == 0)		{apTex[0] = apTexRfuse[0];	apTex[1] = apTexRfuse[1];}
		else if(fusNo == 1)		{apTex[0] = apTexRcast[0];	apTex[1] = apTexRcast[1];}
		else {}

		glUseProgram(pStereoSh->GetProg());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, apTex[0]->GetTexID());
		glUniform1i(pStereoSh->GetUniLoc("smpImgR"), 0);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, apTex[1]->GetTexID());
		glUniform1i(pStereoSh->GetUniLoc("smpImgB"), 1);

		glUniform2f(pStereoSh->GetUniLoc("v2WinDim"), (float)w, (float)h);
		
		glDisable(GL_DEPTH_TEST);
		glColor4f(0, 0, 0, 1.f);
		glBegin(GL_QUADS);
		glVertex3dv(aaCPlane[0]);
		glVertex3dv(aaCPlane[1]);
		glVertex3dv(aaCPlane[2]);
		glVertex3dv(aaCPlane[3]);
		glEnd();
	}
	else {}

	if(bStereo && !bFusion && runNo==1)
	{
		int aWinLow[] = {0, 0};
		apTexImg[2]->ReadFrame(aWinLow);

		apTexImg[2]->GetCell(*pImgW);
		lyrOp.mul.Gen(*pImgW, 255.F);
		imgIO.Write(volName+"_trad.bmp", *pImgW);
	}
	else {}

	if(bFusion && runNo==1)
	{	
		int aWinLow[] = {0, 0};
		apTexImg[fusNo]->ReadFrame(aWinLow);
		
		if(fusNo == 1)
		{			
			if(bRndDis)
			{
				Vect2D<unsigned> dim = apTexImg[fusNo]->GetDim();
				//for(unsigned y=0; y<dim.m_y; y+=2)
				for(unsigned y=0; y<dim.m_y; y++)
				{
					//if(y%4 == 0)
					if(y%2 == 0)
					//for(unsigned yy=y; yy<=y+1; yy++)
					for(unsigned yy=y; yy<=y; yy++)
					{
						if(yy == dim.m_y)  break;
						for(unsigned x=0; x<dim.m_x; x++)
						{
							for(unsigned c=0; c<4; c++)
							{
								float tmp = apTexImg[0]->GetCell(x, yy, c);
								apTexImg[0]->SetCell(apTexImg[1]->GetCell(x, yy, c), x, yy, c);
								apTexImg[1]->SetCell(tmp, x, yy, c);
							}
						}						
					}
				}
			}
			else {}
			apTexImg[0]->Update();
			apTexImg[1]->Update();

			for(unsigned i=0; i<2; i++)
			{
				string str;
				stringstream ss;
				ss.clear();
				if(bBrightness_div)		ss << i;
				else					ss << (i+3);
				ss >> str;

				apTexImg[i]->GetCell(*pImgW);
				lyrOp.mul.Gen(*pImgW, 255.F);
				imgIO.Write(volName+"_"+str+".bmp", *pImgW);
			}
		}
		else {}

		if(fusNo == 1)
		{
			bFlip = true;
			bStereo = false;
			bFusion = false;
			bMap = true;
		}
		else {}		
	}
	else {}	

	glPopMatrix();
}
#include "renderWire.h"

vector<Line> drawLine;

int ComputeIntersect_one(DATA aPInt[], DATA aPFrom[], DATA aPTo[], DATA clipR, unsigned idxS=0)
{
	DATA aaFunc[][4] = {{1.F,   0,   0, -clipR},
						{1.F,   0,   0,  clipR},
						{  0, 1.F,   0, -clipR},
						{  0, 1.F,   0,  clipR},
						{  0,   0, 1.F, -clipR},
				        {  0,   0, 1.F,  clipR}};
	
	for(unsigned i=idxS; i<6; i++)
	{
		bool bInt = myMath.IntersectLinePlane3D(aPInt, aPFrom, aPTo, aaFunc[i]);
		if(bInt)
		{
			return i;
		}
		else {}
	}
	return -1;
}

bool ComputeIntersect_two(DATA aaPInt[][3], DATA aPFrom[], DATA aPTo[], DATA clipR)
{
	int idx = 0;
	for(unsigned i=0; i<2; i++)
	{
		if(idx < 6)
		{
			return false;
		}

		while(idx < 6)
		{
			idx = ComputeIntersect_one(aaPInt[0], aPFrom, aPTo, clipR, idx);
			if(idx == -1)
			{
				return false;
			}
			else {}
			idx++;

			if(aaPInt[i][0]>=-clipR && aaPInt[i][0]<=clipR &&
			   aaPInt[i][1]>=-clipR && aaPInt[i][1]<=clipR &&
			   aaPInt[i][2]>=-clipR && aaPInt[i][2]<=clipR)
			{
				break;
			}
		}
	}

	return true;
}

unsigned TestRecord(DATA aP[], DATA clipR)
{
	unsigned r = 0;
	if(aP[2] > clipR)	r++;
	else {}
	r = r << 1;

	if(aP[2] < -clipR)	r++;
	else {}
	r = r << 1;

	if(aP[1] > clipR)	r++;
	else {}
	r = r << 1;

	if(aP[1] < -clipR)	r++;
	else {}
	r = r << 1;

	if(aP[0] > clipR)	r++;
	else {}
	r = r << 1;

	if(aP[0] < -clipR)	r++;	
	else {}
	r = r << 1;

	return r;
}
void DeterminClipCase(vector<Line> &lineGen, DATA aPA[], DATA aPB[], DATA clipR)
{
	unsigned rA = TestRecord(aPA, clipR);
	unsigned rB = TestRecord(aPB, clipR);

	if(rA==0 && rB==0)
	{
		Line ln;
		myMath.Copy3V(ln.aP0, aPA);
		myMath.Copy3V(ln.aP1, aPB);
		lineGen.push_back(ln);
	}
	else
	{
		if(rA & rB)
		{}
		else 
		{
			if(rA==0 || rB==0)
			{
				DATA aPFrom[3], aPTo[3];
				if(rA == 0)
				{
					myMath.Copy3V(aPFrom, aPA);
					myMath.Copy3V(aPTo, aPB);
				}
				else
				{
					myMath.Copy3V(aPFrom, aPB);
					myMath.Copy3V(aPTo, aPA);
				}

				DATA aPInt[3];
				DATA minDis = 1e10;
				int idx = 0;
				while(idx < 6)
				{
					DATA aPTmp[3];
					idx = ComputeIntersect_one(aPTmp, aPFrom, aPTo, clipR, idx);
					if(idx == -1)
					{
						break;
					}
					else {}
					idx++;
					
					DATA aDis[3];
					myMath.Copy3V(aDis, aPFrom);
					myMath.Sub3V(aDis, aPTmp);
					DATA dis = aDis[0]*aDis[0] + aDis[1]*aDis[1] + aDis[2]*aDis[2];
					if(dis < minDis)
					{
						minDis = dis;
						myMath.Copy3V(aPInt, aPTmp);
					}
				}

				Line ln;
				myMath.Copy3V(ln.aP0, aPFrom);
				myMath.Copy3V(ln.aP1, aPInt);
				lineGen.push_back(ln);
			}
			else
			{
				DATA aaPInt[2][3];
				bool bGet = ComputeIntersect_two(aaPInt, aPA, aPB, clipR);
				if(bGet)
				{
					Line ln;
					myMath.Copy3V(ln.aP0, aaPInt[0]);
					myMath.Copy3V(ln.aP1, aaPInt[1]);
					lineGen.push_back(ln);
				}
				else {}
			}
		}
	}
}

void ClipLines(vector<Line> &lineGen, const TRIModel &canModel, DATA clipR)
{
	unsigned aAIdx[] = {0, 1, 2};
	unsigned aBIdx[] = {1, 2, 0};

	lineGen.clear();
	for(unsigned t=0; t<canModel.triangleList.size(); t++)
	{
		Triangle tri = canModel.triangleList[t];
		for(unsigned i=0; i<3; i++)
		{
			DeterminClipCase(lineGen, tri.vertex[aAIdx[i]], tri.vertex[aBIdx[i]], clipR);
		}
	}
}

void RenderWire(const TRIModel &canModel, DATA clipR)
{
	ClipLines(drawLine, canModel, clipR);

	glColor4f(1.f, 1.f, 1.f, 1.f);
	glBegin(GL_LINES);
	for(unsigned i=0; i<drawLine.size(); i++)
	{		
		glVertex3dv(drawLine[i].aP0);	
		glVertex3dv(drawLine[i].aP1);				
	}
	glEnd();

	glColor4f(1.f, 0, 0, 1.f);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_QUADS);
		glVertex2d(-clipR, -clipR);
		glVertex2d( clipR, -clipR);
		glVertex2d( clipR,  clipR);
		glVertex2d(-clipR,  clipR);
	glEnd();
}
#include "main.h"

#include <GLFW/glfw3.h>


void main()
{
	pMain_win = make_window();
	pMain_win->show();
	pRender_glWin->show();

	Fl::run();
}

/*
const Vect2D<unsigned> dimFbo(400, 400);

void main() 
{
	bool bDebug = true;

	int fErr = glfwInit();
	MyAssert(fErr == GL_TRUE);
	if (bDebug) {
		cout << "glfwInit ok" << endl;
	} else {}
	
	GLFWwindow* pWin = glfwCreateWindow(dimFbo.m_x, dimFbo.m_y, "", 0, 0);
	if (bDebug) {
		cout << "glfwCreateWindow ok" << endl;
	} else {}

	MyAssert(pWin != 0);
	glfwMakeContextCurrent(pWin);
	glfwHideWindow(pWin);

	//glewExperimental = GL_TRUE;
	GLenum err = glewInit();
	MyAssert(err == GLEW_OK);

	pFlatSh = new Shader("flatSh.vert", "flatSh.frag");
	pFlatSh->Load();
	
	// RGBA 2D texture, 24 bit depth texture, 
	GLuint txColor;
	glGenTextures(1, &txColor);
	glBindTexture(GL_TEXTURE_2D, txColor);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	// 0 means reserve texture memory, but texels are undefined
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, dimFbo.m_x, dimFbo.m_y, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);

	GLuint frameBf;
	glGenFramebuffers(1, &frameBf);
	glBindFramebuffer(GL_FRAMEBUFFER, frameBf);
	// attach 2D texture to this FBO
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, txColor, 0);

	GLuint depthBf;
	glGenRenderbuffers(1, &depthBf);
	glBindRenderbuffer(GL_RENDERBUFFER, depthBf);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, dimFbo.m_x, dimFbo.m_y);
	// attach depth buffer to FBO
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBf);

	// does the GPU support current FBO configuration?
	GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	MyAssert(status == GL_FRAMEBUFFER_COMPLETE);

	TRIModel *pTriModel = 0;
	delete pTriModel;
	pTriModel = new TRIModel;
	pTriModel->loadFromFile("D:\\working\\renderXman\\faith_all.tri");

	// rendering
	glClearColor(0.5, 0.5, 0.5, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// my projection matrix
	DATA r = (DATA)dimFbo.m_y / dimFbo.m_x;
	DATA dd = 8.F;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glFrustum(-0.02F, 0.02F, -0.02F * r, 0.02F * r, 0.03F, 2.F);
	glFrustum(-0.013F, 0.013F, -0.013F * r, 0.013F * r, 0.03F, 2.F);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//DATA aENow[] = {aEPos[0], aEPos[1], aEPos[2]};
	//DATA aENow[] = {-0.94F, 3.40F, -3.59F + 1.F};
	DATA aENow[] = {-1.F, 3.5F, -3.59F + 1.F};

	gluLookAt(aENow[0], aENow[1], aENow[2],
			  aENow[0], aENow[1], aENow[2] - 1.F,
			  0, 1.F, 0);

	
	//glTranslated(aVCnt[0], aVCnt[1], aVCnt[2]);
	//RotateObj();
	//glScaled(viewScl, viewScl, viewScl);
	//glTranslated(-aObjCnt[0], -aObjCnt[1], -aObjCnt[2]);
	

	glViewport(0, 0, dimFbo.m_x, dimFbo.m_y);

	glEnable(GL_DEPTH_TEST);
	glUseProgram(pFlatSh->GetProg());
	RenderShading(*pTriModel, 0);
	glDisable(GL_DEPTH_TEST);

	GLubyte *pPixel = new GLubyte [dimFbo.m_x * dimFbo.m_y * 4];
	glReadPixels(0, 0, dimFbo.m_x, dimFbo.m_y, GL_RGBA, GL_UNSIGNED_BYTE, pPixel);

	Layer lyrOut(dimFbo.m_x, dimFbo.m_y, 3);
	Vect3D<unsigned> dimOut = lyrOut.GetDim();
	for (unsigned y = 0; y < dimFbo.m_y; y++) {
		for (unsigned x = 0; x < dimFbo.m_x; x++) {
			for (unsigned c = 0; c < dimOut.m_z; c++) {
				lyrOut.CellRef(x, y, c) = pPixel[(y * dimFbo.m_x + x) * 4 + c];
			}
		}
	}
	imgIO.Write("out.jpg", MyImg(lyrOut));
}
*/
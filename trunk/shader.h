#ifndef _SHADER_H
#define _SHADER_H

#include <windows.h>
#include <GL/glew.h>
//#include "glee.h"
#include <iostream>
#include <fstream>

using namespace std;

class Shader
{
public:
	Shader(const string &vFname, const string &fFname);
	~Shader();

	bool Load();

	GLuint GetProg();
	GLint GetUniLoc(char *pName);
	GLint GetAttrLoc(char *pName);

private:
	void Dump(GLuint handle);
	bool Compile(GLuint shader, const string &fname);

	GLuint m_prog;
	GLuint m_vertShader;
	GLuint m_fragShader;
	bool m_bReady;
	bool m_bInvalid;

	string m_vertFname;
	string m_fragFname;
};

#endif
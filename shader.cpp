#include "shader.h"

Shader::Shader(const string &vFname, const string &fFname)
	:m_vertFname(vFname), m_fragFname(fFname),
	m_prog(0), m_vertShader(0), m_fragShader(0), m_bReady(false), m_bInvalid(false)
{}

Shader::~Shader()
{}

//************************************************

void Shader::Dump(GLuint handle)
{
	char* buffer = 0;
    int length = 0;

    // reset the error code
    glGetError();

    // this function accepts both shaders AND programs, so do a type check.
    if(glIsProgram(handle))
    {
		glGetProgramiv(handle, GL_INFO_LOG_LENGTH, &length);
		if(glGetError() != GL_NO_ERROR)
		{
			return;
		}
		else if(length <= 0)
		{
			return;
		}
		else {}

		buffer = new char[length];
        glGetProgramInfoLog(handle, length, 0, buffer);
        if(glGetError() != GL_NO_ERROR)
		{
			delete []buffer;
            return;
		}
		else {}
    }
    else if(glIsShader(handle))
    {
        glGetShaderiv(handle, GL_INFO_LOG_LENGTH, &length);
        if(glGetError() != GL_NO_ERROR)
		{
            return;
		}
		else if(length <= 0)
		{
			return;
		}

        buffer = new char[length];
        glGetShaderInfoLog(handle, length, 0, buffer);
        if(glGetError() != GL_NO_ERROR)
		{
			delete []buffer;
            return;
		}
		else {}
    }
	else {}

    if(buffer)
    {
        cout << buffer << endl;
        delete []buffer;
    }
	else {}
}

//***********************************************

bool Shader::Compile(GLuint shader, const string &fname)
{
	cout << "compiling " << fname.c_str() << "..." << endl;

	ifstream fp(fname.c_str(), ios::binary);
	if(fp.bad())
	{
		cout << "unable to open file" << endl;
		return false;
	}
	else {}

	fp.seekg(0, ios::end);
 	unsigned size = (unsigned)fp.tellg();
	char *pBuffer = new char[size+1];

	fp.seekg(0, ios::beg);
	fp.read(pBuffer, size);
	pBuffer[size] = '\0';
	//cout << pBuffer << endl;
	
	glShaderSource(shader, 1, (const GLchar**)&pBuffer, 0);
	delete []pBuffer;

    int success;
    glCompileShader(shader);
    Dump(shader);
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if(success)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Shader::Load()
{
	if(!glUseProgram || m_bInvalid)
	{
		return false;
	}
	else {}

	if(m_bReady)
	{
		glUseProgram(m_prog);
	}
	else {}

	//**********************************************

	m_vertShader = glCreateShader(GL_VERTEX_SHADER);
    if(!Compile(m_vertShader, m_vertFname)) 
	{
		cout << "unable to compile " << m_vertFname.c_str() << endl;
        m_bInvalid = true;
        return false;
    }
	else {}

	m_fragShader = glCreateShader(GL_FRAGMENT_SHADER);
    if(!Compile(m_fragShader, m_fragFname)) 
	{
		cout << "unable to compile " << m_fragFname.c_str() << endl;
        m_bInvalid = true;
        return false;
    }
	else {}

	m_prog = glCreateProgram();
	glAttachShader(m_prog, m_vertShader);
	glAttachShader(m_prog, m_fragShader);

	//**********************************************

    int success;
	cout << "linking " << m_vertFname.c_str() << " with " << m_fragFname.c_str() << "..." << endl;
	glLinkProgram(m_prog);
	glGetProgramiv(m_prog, GL_LINK_STATUS, &success);
    if(!success) 
	{
		cout << "failed link" << endl;
        Dump(m_prog);
        m_bInvalid = true;
        return false;
    }
	else {}

	cout << "validatin..." << endl;
    glValidateProgram(m_prog);
    glGetProgramiv(m_prog, GL_VALIDATE_STATUS, &success);
    if(!success) 
	{
		cout << "failed validation" << endl;
        Dump(m_prog);
        m_bInvalid = true;
        return false;
    }
	else {}

	//**********************************************

    Dump(m_prog);
    glUseProgram(m_prog);
    m_bReady = true;
	cout << "load shader ok" << endl;
	return true;
}

//*************************************************************************************************

GLuint Shader::GetProg()
{
	return m_prog;
}

GLint Shader::GetUniLoc(char *pName)
{
	GLint loc = glGetUniformLocation(m_prog, pName);
    if(loc == -1)
	{
        cout << "no such uniform named " << pName << endl;
	}
    return loc;
}

GLint Shader::GetAttrLoc(char *pName)
{
	GLint loc = glGetAttribLocation(m_prog, pName);
	if(loc == -1)
	{
		cout << "no such attribute named " << pName << endl;
	}
	return loc;
}
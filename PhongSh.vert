uniform vec3 v3Eye;
uniform int lightNum;

varying vec3 v3N;
varying vec3 v3StdPos;

void main(void)
{
    v3N = normalize(gl_NormalMatrix * gl_Normal);

    vec4 v4StdPos = gl_ModelViewMatrix * gl_Vertex;
    v3StdPos = v4StdPos.xyz; 

    gl_FrontColor = gl_Color;
    gl_Position = ftransform(); 
}

uniform vec4 v4TfX;
uniform vec4 v4TfY;
uniform vec4 v4Tf0X;
uniform vec4 v4Tf0Y;
uniform vec4 v4Tf1X;
uniform vec4 v4Tf1Y;
uniform vec3 v3Cnt;
uniform int tfNo;
uniform float symetry;

uniform vec3 v3Eye;
uniform int lightNum;
uniform vec3 v3Pos;
uniform vec3 v3Norm;

varying vec2 v2TxF;
varying vec2 v2TxF0;
varying vec2 v2TxF1;

void main(void)
{
    vec3 v3StdN = normalize(gl_NormalMatrix * v3Norm);

    vec4 v4Pos = vec4(v3Pos.xyz, 1.0);
    vec4 v4StdPos = gl_ModelViewMatrix * v4Pos;
    vec3 v3StdPos = v4StdPos.xyz; 

    vec3 v3StdEye = v3Eye;
    vec3 v3V_E = normalize(v3StdEye - v3StdPos);

    vec4 v4X = v4TfX;
    vec4 v4Y = v4TfY;
    if (tfNo == 1) {
        v4X = v4Tf1X;
        v4Y = v4Tf1Y;
    }

    if (symetry == 1.0 && gl_Vertex.x < v3Cnt.x) {
        v2TxF0.x = (v3Cnt.x * 2 - gl_Vertex.x) * v4X.r + 
                   gl_Vertex.y * v4X.g + 
                   gl_Vertex.z * v4X.b + 
                   v4X.a;
        v2TxF0.y = (v3Cnt.x * 2 - gl_Vertex.x) * v4Y.r + 
                   gl_Vertex.y * v4Y.g + 
                   gl_Vertex.z * v4Y.b + 
                   v4Y.a;
    } else {
        v2TxF0.x = gl_Vertex.x * v4X.r + 
                   gl_Vertex.y * v4X.g + 
                   gl_Vertex.z * v4X.b + 
                   v4X.a;
        v2TxF0.y = gl_Vertex.x * v4Y.r + 
                   gl_Vertex.y * v4Y.g + 
                   gl_Vertex.z * v4Y.b + 
                   v4Y.a;
    }

    v2TxF.x = v2TxF.x / 332.0;
    v2TxF.y = v2TxF.y / 419.0;
    v2TxF0.x = v2TxF0.x / 332.0;
    v2TxF0.y = v2TxF0.y / 419.0;
    v2TxF1.x = v2TxF1.x / 332.0;
    v2TxF1.y = v2TxF1.y / 419.0;


    gl_FrontColor = vec4(0.0, 0.0, 0.0, 1.0);
    for(int i=0; i<lightNum; i++)
    {
        vec4 v4StdLPos = gl_LightSource[i].position;

        vec4 v4Ambi = gl_LightSource[i].ambient * gl_Color;

        vec3 v3V_L = normalize(v4StdLPos.xyz - v3StdPos);
        float cosNL = max(dot(v3StdN.xyz, v3V_L), 0.0);
        vec4 v4Diff = gl_LightSource[i].diffuse * cosNL * gl_Color;

        //vec3 v3V_E = normalize(v3StdEye - v3StdPos);
        vec3 v3HV = normalize((v3V_L + v3V_E)); 
        float cosH = max(dot(v3StdN.xyz, v3HV), 0.0);
        vec4 v4Spec = gl_LightSource[i].specular * pow(cosH, gl_LightSource[i].spotExponent);
	
        gl_FrontColor.rgb += v4Ambi.rgb + v4Diff.rgb + v4Spec.rgb;
    }
    gl_Position = ftransform();
}

/*
void main(void)
{
    vec3 v3StdN = normalize(gl_NormalMatrix * v3Norm);

    vec4 v4Pos = vec4(v3Pos.xyz, 1.0);
    vec4 v4StdPos = gl_ModelViewMatrix * v4Pos;
    vec3 v3StdPos = v4StdPos.xyz; 

    vec3 v3StdEye = v3Eye;
 
    gl_FrontColor = vec4(0.0, 0.0, 0.0, 1.0);
    for(int i=0; i<lightNum; i++)
    {
        vec4 v4StdLPos = gl_LightSource[i].position;

        vec4 v4Ambi = gl_LightSource[i].ambient * gl_Color;

        vec3 v3V_L = normalize(v4StdLPos.xyz - v3StdPos);
        float cosNL = max(dot(v3StdN.xyz, v3V_L), 0.0);
        vec4 v4Diff = gl_LightSource[i].diffuse * cosNL * gl_Color;

        vec3 v3V_E = normalize(v3StdEye - v3StdPos);
        vec3 v3HV = normalize((v3V_L + v3V_E)); 
        float cosH = max(dot(v3StdN.xyz, v3HV), 0.0);
        vec4 v4Spec = gl_LightSource[i].specular * pow(cosH, gl_LightSource[i].spotExponent);
	
        gl_FrontColor.rgb += v4Ambi.rgb + v4Diff.rgb + v4Spec.rgb;
    }
    gl_Position = ftransform();
}
*/

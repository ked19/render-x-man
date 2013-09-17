uniform vec3 v3Eye;
uniform int lightNum;
uniform vec3 v3Pos;
uniform vec3 v3Norm;

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

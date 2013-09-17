uniform vec3 v3Eye;
uniform int lightNum;

varying vec3 v3N;
varying vec3 v3StdPos;

void main(void)
{	
    vec3 v3StdEye = v3Eye;
    vec3 v3StdN = normalize(v3N);

    gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
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
	
        gl_FragColor.rgb += v4Ambi.rgb + v4Diff.rgb + v4Spec.rgb;
    }
    gl_FragColor = min(gl_FragColor, 1.0);
}

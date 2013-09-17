uniform sampler2D smpStrPos; 
uniform vec2 v2WinDim;

void main(void)
{
    vec2 v2TexPos = vec2(gl_FragCoord.x/v2WinDim.x, 
           	         gl_FragCoord.y/v2WinDim.y);
    vec4 v4StrPosVal = texture2D(smpStrPos, v2TexPos);

    vec3 v3Vect = gl_Color.rgb - v4StrPosVal.rgb;
    gl_FragColor.rgb = v3Vect * v4StrPosVal.a;
    gl_FragColor.a =  length(v4StrPosVal.rgb);
    gl_FragColor.rgb = normalize(gl_FragColor.rgb);

    gl_FragColor.rgb = (gl_FragColor.rgb + 1.0) * 0.5;
}

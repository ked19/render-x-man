//uniform sampler2D smpImg; 
uniform sampler2D smpImg0; 
uniform sampler2D smpImg1; 
uniform float rndRatio;
uniform vec2 v2WinDim;

void main(void)
{
    vec2 v2TexPos = vec2(gl_FragCoord.x/v2WinDim.x, 
           	         gl_FragCoord.y/v2WinDim.y);
    //vec4 v4Tmp = texture2D(smpImg, v2TexPos);
    //gl_FragColor = v4Tmp;
    
    vec4 v4Tmp0 = texture2D(smpImg0, v2TexPos);
    vec4 v4Tmp1 = texture2D(smpImg1, v2TexPos);
    gl_FragColor = v4Tmp0*rndRatio + (1.0-rndRatio)*v4Tmp1;
}

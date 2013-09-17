uniform sampler2D smpImgR; 
uniform sampler2D smpImgB; 
uniform vec2 v2WinDim;

void main(void)
{
    vec2 v2Coord = gl_FragCoord.xy;
    vec2 v2TexPos = vec2(gl_FragCoord.x/v2WinDim.x, 
           	         gl_FragCoord.y/v2WinDim.y);
    vec4 v4ImgR = texture2D(smpImgR, v2TexPos);
    vec4 v4ImgB = texture2D(smpImgB, v2TexPos);
    v4ImgR.rgb = v4ImgR.rgb + (1.0-v4ImgR.a) * vec3(0.0, 0.2, 0.4);
    v4ImgB.rgb = v4ImgB.rgb + (1.0-v4ImgB.a) * vec3(0.0, 0.2, 0.4);

//if(mod(gl_FragCoord.y, 2.0) < 1.0)
//gl_FragColor.rgb = v4ImgR.rgb;
//else
//gl_FragColor.rgb = v4ImgB.rgb;
//gl_FragColor.a = 1.0;

    gl_FragColor.rgb = vec3(v4ImgR.g*0.7+v4ImgR.b*0.3, v4ImgB.g, v4ImgB.b);
    gl_FragColor.r = pow(gl_FragColor.r, 1.0/1.5);
    
    gl_FragColor.a = 1.0;
}

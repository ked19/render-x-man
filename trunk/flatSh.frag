uniform sampler2D smpFace;
varying vec2 v2TxF;
varying vec2 v2TxF0;
varying vec2 v2TxF1;

void main(void)
{	
    vec4 v4FColor = texture2D(smpFace, v2TxF);
    vec4 v4FColor0 = texture2D(smpFace, v2TxF0);
    vec4 v4FColor1 = texture2D(smpFace, v2TxF1);

    //gl_FragColor = gl_Color;
    gl_FragColor = vec4(v4FColor0.rgb, 1.0);
    //gl_FragColor = vec4(v4FColor.rgb * length(gl_Color.rgb), 1.0);
    //gl_FragColor = vec4(gl_Color.rgb, 1.F);

	//gl_FragColor = vec4(v4FColor0.rgb * 0.8 + gl_Color.rgb * 0.2, 1.0);
}

/*
void main(void)
{	
    gl_FragColor = gl_Color;
}
*/
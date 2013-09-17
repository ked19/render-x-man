uniform sampler2D smpStrPos;
uniform sampler2D smpLengthImg;
uniform sampler3D smpAct;
uniform sampler3D smpData;
uniform sampler2D smpColor;
uniform sampler1D smpId;
uniform sampler3D smpMsk;
uniform sampler3D smpPhot;
uniform vec2 v2WinDim;
uniform vec3 v3BlckRatio;
uniform vec3 v3View1;
uniform vec3 v3ViewX1;
//uniform vec3 v3ViewY1;
//uniform vec3 v3DimRatio;

uniform float disAtt;
uniform float opaScl;
uniform float briScl;
//uniform float mskThrd;
uniform float smpNum;

uniform float bMIP;
uniform float occScl;

uniform float bXClip;  uniform float bXFront;  uniform float xPlane;
uniform float bYClip;  uniform float bYFront;  uniform float yPlane;
uniform float bZClip;  uniform float bZFront;  uniform float zPlane;

uniform float bBG;

#ifdef fusion_2
uniform float opaEnd;
#endif

void main(void)
{
    #ifdef DEPTH_TEST
    vec3 v3G = vec3(0.5, 0.7, 0.5);
    vec3 v3R = vec3(0.7, 0.5, 0.5);
    float disLen = 0.005;
    #endif

    vec2 v2PlaneTxt = vec2(gl_FragCoord.x/v2WinDim.x, 
           	           gl_FragCoord.y/v2WinDim.y);
    vec4 v4StrPos = texture2D(smpStrPos, v2PlaneTxt);
    vec4 v4LengthImg = texture2D(smpLengthImg, v2PlaneTxt);
    v4LengthImg.rgb = v4LengthImg.rgb*2.0 - 1.0;

    vec3 v3View  = v3View1;	v3View	= normalize( v3View );
    vec3 v3ViewX = v3ViewX1; 	v3ViewX = normalize( v3ViewX );
    //vec3 v3ViewY = v3ViewY1; 	v3ViewY = normalize( v3ViewY );

    vec4 v4Rel = vec4(0.0, 0.0, 0.0, 0.0);
    vec3 v3Pos = v4StrPos.xyz; 
    vec3 v3PPre;
    vec3 v3BlckPos = v3Pos.xyz / v3BlckRatio.xyz;
    for(int i=1; i<=smpNum; i++)
    {
        if(v3Pos.x<0.0 || v3Pos.x>1.0 ||
           v3Pos.y<0.0 || v3Pos.y>1.0 ||
           v3Pos.z<0.0 || v3Pos.z>1.0) 
        {
            break;
        }
        #ifdef fusion_2
        if(v4Rel.a > opaEnd)
        {
            break;
        }
        #endif

        #ifdef DEPTH_TEST
        vec3 v3GDif = (v3Pos - v3G) * v3DimRatio;
        vec3 v3RDif = (v3Pos - v3R) * v3DimRatio;
        float lenG = length(v3GDif);
        float lenR = length(v3RDif);
        float bSetG = 0.0;
        float bSetR = 0.0;
        if(lenG < disLen)  bSetG = 1.0;
        if(lenR < disLen)  bSetR = 1.0;
        #endif

        vec4 v4Act = texture3D(smpAct, v3BlckPos); 
        #ifdef DEPTH_TEST
        if(v4Act.r>=0.5 || bSetG>0.5 || bSetR>0.5)
        #else
        if(v4Act.r >= 0.5) 
        #endif
        {
            float val = texture3D(smpData, v3Pos).r;
            float valPre = texture3D(smpData, v3PPre).r;
            vec4 v4Col = texture2D(smpColor, vec2(valPre, val));
            v4Col.rgb /= v4Col.a;

            //if(v4Col.a == 0.0)
            //{
            //    v4Col.rgb = vec3(0.0, 0.0, 0.0);
            //}
float eg = 0.0001;
//v4Col.rgb = step(eg, v4Col.a)*v4Col.rgb + (1.0-step(eg, v4Col.a))*vec3(0.0, 0.0, 0.0);
v4Col.rgb = step(eg, v4Col.a)*v4Col.rgb;

            vec3 v3Shadow = texture3D(smpPhot, v3Pos).rgb;
            float msk = texture3D(smpMsk, v3Pos).r;
            vec3 v3EndShd = vec3( texture3D(smpPhot, vec3(1.0, v3Pos.y, v3Pos.z)).r,
                                  texture3D(smpPhot, vec3(v3Pos.x, 1.0, v3Pos.z)).g,
                                  texture3D(smpPhot, vec3(v3Pos.x, v3Pos.y, 1.0)).b );
            vec3 v3ClipShd = vec3( texture3D(smpPhot, vec3(xPlane, v3Pos.y, v3Pos.z)).r,
                                   texture3D(smpPhot, vec3(v3Pos.x, yPlane, v3Pos.z)).g,
                                   texture3D(smpPhot, vec3(v3Pos.x, v3Pos.y, zPlane)).b );
            
            vec3 v3Dis = v3Pos;
            vec3 v3EndDis = vec3(1.0, 1.0, 1.0);

//if(msk < 0.01 )
//v4Col.a = 0.0;
float eg2 = 0.01;
//v4Col.a = step(eg2, msk)*v4Col.a + (1.0-step(eg2, msk))*0.0;
v4Col.a = step(eg2, msk)*v4Col.a;

float tId = texture1D(smpId, val).r;
if(msk == tId)
v4Col.a = 0.0;

v4Col.a *= opaScl;
            
            v4Col.a = bXClip * ( bXFront*(1.0-step(xPlane, v3Pos.x))*v4Col.a + (1.0-bXFront)*(1.0-step(v3Pos.x, xPlane))*v4Col.a ) +
                      (1.0-bXClip) * v4Col.a;
            v3EndShd.r = bXClip * ( bXFront*v3ClipShd.r + (1.0-bXFront)*(v3EndShd.r-v3ClipShd.r) ) +
                         (1.0-bXClip) * v3EndShd.r;
            //v3Shadow.r = bXClip * ( (1.0-bXFront)*(v3Shadow.r-v3ClipShd.r) + bXFront*v3Shadow.r ) +
            //             (1.0-bXClip) * v3Shadow.r;  
            v3Shadow.r -= bXClip * ( (1.0-bXFront)*(v3ClipShd.r));  

            v4Col.a = bYClip * ( bYFront*(1.0-step(yPlane, v3Pos.y))*v4Col.a + (1.0-bYFront)*(1.0-step(v3Pos.y, yPlane))*v4Col.a ) +
                      (1.0-bYClip) * v4Col.a;
            v3EndShd.g = bYClip * ( bYFront*v3ClipShd.g + (1.0-bYFront)*(v3EndShd.g-v3ClipShd.g) ) +
                         (1.0-bYClip) * v3EndShd.g;
            //v3Shadow.g = bYClip * ( (1.0-bYFront)*(v3Shadow.g-v3ClipShd.g) + bYFront*v3Shadow.g ) +
            //             (1.0-bYClip) * v3Shadow.g;  
            v3Shadow.g -= bYClip * ( (1.0-bYFront)*(v3ClipShd.g));  

            v4Col.a = bZClip * ( bZFront*(1.0-step(zPlane, v3Pos.z))*v4Col.a + (1.0-bZFront)*(1.0-step(v3Pos.z, zPlane))*v4Col.a ) +
                      (1.0-bZClip) * v4Col.a;
            v3EndShd.b = bZClip * ( bZFront*v3ClipShd.b + (1.0-bZFront)*(v3EndShd.b-v3ClipShd.b) ) +
                         (1.0-bZClip) * v3EndShd.b;
            //v3Shadow.b = bZClip * ( (1.0-bZFront)*(v3Shadow.b-v3ClipShd.b) + bZFront*v3Shadow.b ) +
            //             (1.0-bZClip) * v3Shadow.b;  
            v3Shadow.b -= bZClip * ( (1.0-bZFront)*(v3ClipShd.b));  

            vec3 v3ShNeg = vec3(v3EndShd.r - v3Shadow.r, 
                                v3EndShd.g - v3Shadow.g, 
                                v3EndShd.b - v3Shadow.b); 
               
            float att = disAtt * 20.0;  
            vec3 v3ExpPos = vec3(exp(-v3Shadow.r*att), exp(-v3Shadow.g*att), exp(-v3Shadow.b*att));
            vec3 v3ExpNeg = vec3(exp(-v3ShNeg.r*att), exp(-v3ShNeg.g*att), exp(-v3ShNeg.b*att));
   
            float hSq3 = sqrt(3.0) / 2.0;
            float sXKey = hSq3*v3ViewX.x + 0.5*v3View.x;	//sXKey /= v3DimRatio.x;
            float sYKey = hSq3*v3ViewX.y + 0.5*v3View.y;	//sYKey /= v3DimRatio.y;
            float sZKey = hSq3*v3ViewX.z + 0.5*v3View.z;	//sZKey /= v3DimRatio.z;
            //float sumKey = abs(sXKey) + abs(sYKey) + abs(sZKey);
            //sXKey = sXKey / sumKey;
            //sYKey = sYKey / sumKey;
            //sZKey = sZKey / sumKey;
            float shXKey = step(0.0, sXKey)*v3ExpPos.r*(sXKey) + (1.0-step(0.0, sXKey))*v3ExpNeg.r*(-sXKey);
            float shYKey = step(0.0, sYKey)*v3ExpPos.g*(sYKey) + (1.0-step(0.0, sYKey))*v3ExpNeg.g*(-sYKey);
            float shZKey = step(0.0, sZKey)*v3ExpPos.b*(sZKey) + (1.0-step(0.0, sZKey))*v3ExpNeg.b*(-sZKey);
            float shadowValKey = (shXKey + shYKey + shZKey); 
            
            float sXFil = -hSq3*v3ViewX.x + 0.5*v3View.x;	//sXFil /= v3DimRatio.x;
            float sYFil = -hSq3*v3ViewX.y + 0.5*v3View.y;	//sYFil /= v3DimRatio.y;
            float sZFil = -hSq3*v3ViewX.z + 0.5*v3View.z;	//sZFil /= v3DimRatio.z;
            //float sumFil = abs(sXFil) + abs(sYFil) + abs(sZFil);
            //sXFil = sXFil / sumFil;
            //sYFil = sYFil / sumFil;
            //sZFil = sZFil / sumFil;
            float shXFil = step(0.0, sXFil)*v3ExpPos.r*(sXFil) + (1.0-step(0.0, sXFil))*v3ExpNeg.r*(-sXFil);
            float shYFil = step(0.0, sYFil)*v3ExpPos.g*(sYFil) + (1.0-step(0.0, sYFil))*v3ExpNeg.g*(-sYFil);
            float shZFil = step(0.0, sZFil)*v3ExpPos.b*(sZFil) + (1.0-step(0.0, sZFil))*v3ExpNeg.b*(-sZFil);
            float shadowValFil = (shXFil + shYFil + shZFil); 

            float sXBak = -v3View.x;	//sXBak /= v3DimRatio.x;
            float sYBak = -v3View.y;	//sYBak /= v3DimRatio.y;
            float sZBak = -v3View.z;	//sZBak /= v3DimRatio.z;
            //float sumBak = abs(sXBak) + abs(sYBak) + abs(sZBak);
            //sXBak = sXBak / sumBak;
            //sYBak = sYBak / sumBak;
            //sZBak = sZBak / sumBak;
            float shXBak = step(0.0, sXBak)*v3ExpPos.r*(sXBak) + (1.0-step(0.0, sXBak))*v3ExpNeg.r*(-sXBak);
            float shYBak = step(0.0, sYBak)*v3ExpPos.g*(sYBak) + (1.0-step(0.0, sYBak))*v3ExpNeg.g*(-sYBak);
            float shZBak = step(0.0, sZBak)*v3ExpPos.b*(sZBak) + (1.0-step(0.0, sZBak))*v3ExpNeg.b*(-sZBak);
            float shadowValBak = (shXBak + shYBak + shZBak); 

            float sXHed = v3View.x;	//sXHed /= v3DimRatio.x;
            float sYHed = v3View.y;	//sYHed /= v3DimRatio.y;
            float sZHed = v3View.z;	//sZHed /= v3DimRatio.z;
            //float sumHed = abs(sXHed) + abs(sYHed) + abs(sZHed);
            //sXHed = sXHed / sumHed;
            //sYHed = sYHed / sumHed;
            //sZHed = sZHed / sumHed;
            float shXHed = step(0.0, sXHed)*v3ExpPos.r*(sXHed) + (1.0-step(0.0, sXHed))*v3ExpNeg.r*(-sXHed);
            float shYHed = step(0.0, sYHed)*v3ExpPos.g*(sYHed) + (1.0-step(0.0, sYHed))*v3ExpNeg.g*(-sYHed);
            float shZHed = step(0.0, sZHed)*v3ExpPos.b*(sZHed) + (1.0-step(0.0, sZHed))*v3ExpNeg.b*(-sZHed);
            float shadowValHed = (shXHed + shYHed + shZHed);
            
            float shadowVal = (shadowValKey/3.0 + shadowValFil/3.0 + shadowValBak/2.0 + shadowValHed); // / (1.0/3.0 + 1.0/3.0 + 1.0/2.0 + 1.0);
            shadowVal = shadowVal / (1/3.0 + 1/3.0 + 1/2.0 + 1.0);
            //shadowVal = min(shadowVal, 1.0);
            //shadowVal = max(shadowVal, 0.01);
            v4Col.rgb *= shadowVal; 

            #ifdef DEPTH_TEST
            if(bSetG > 0.5)  v4Col = vec4(0.0, 1.0, 0.0, 1.0);
            if(bSetR > 0.5)  v4Col = vec4(1.0, 0.0, 0.0, 1.0);
            #endif

            float shadowValViw = 1.0;
            v4Rel.rgb = (1.0-bMIP)*(v4Rel.rgb + (1.0-v4Rel.a)*v4Col.rgb*v4Col.a) + bMIP*(v4Rel.rgb + v4Col.rgb*shadowValViw*v4Col.a);
            v4Rel.a = (1.0-bMIP)*(v4Rel.a + (1.0-v4Rel.a)*v4Col.a) + bMIP*(v4Rel.a + shadowValViw*v4Col.a);

            if(bMIP>0.5 && v4Rel.a>occScl*5.0)
            {
                break;
            }
            
            v4Rel.a = (1.0-bMIP)*min(1.0, v4Rel.a) + bMIP*v4Rel.a;
        }

        v3PPre = v3Pos;
        v3Pos += v4LengthImg.xyz / smpNum;
        v3BlckPos = v3Pos.xyz / v3BlckRatio.xyz;
    }

    if(v4Rel.a > 1.0)
    {
        v4Rel.rgb /= v4Rel.a;
        v4Rel.a = 1.0;
    }

v4Rel.rgb = v4Rel.rgb * briScl * 255.0 / 200.0;
v4Rel.rgb = v4Rel.rgb / (v4Rel.rgb+1.0);

    //if(bBG > 0.5)
    //{
    //    v4Rel.rgb = v4Rel.rgb + (1.0 - v4Rel.a) * vec3(0, 0.2, 0.4) ;
    //}
    v4Rel.rgb += bBG*(1.0-v4Rel.a)*vec3(0.0, 0.2, 0.4);
    gl_FragColor = v4Rel;
}

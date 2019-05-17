//==============================================================================
// Base P5 Fragment Shader
//==============================================================================
#ifdef GL_ES
precision mediump float; // set precision level if available
#endif
//==============================================================================
// Global Constants
const float twoPi = 6.283185307179586;
//==============================================================================
// input variables
uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;
//==============================================================================
// varying: these have all come from the vertex shader
varying vec3 var_vertPos;
varying vec4 var_vertCol;
varying vec3 var_vertNormal;
varying vec2 var_vertTexCoord;
varying vec4 v_color;
//==============================================================================
// functions
vec3 getPixelColour(vec2 pos)
{
    float dis = distance(pos, vec2(0.0)); // distance from centre
    return vec3(sin(twoPi * dis * time / 2.0), 0.2, sin(dis));
}
//==============================================================================
void main( void )
{
    vec2 pos = ((gl_FragCoord.xy / resolution.xy) * 2.0) - 1.0; // scale screen co-ordinates to -1 < (x, y) < 1
    pos.x *= resolution.x / resolution.y;
    gl_FragColor = vec4(var_vertPos.z/100.0,0.0,0.0,1.0);
}

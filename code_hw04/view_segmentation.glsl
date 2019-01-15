#version 150
// Copyright 2015, Christopher J. Foster and the other displaz contributors.
// Use of this code is governed by the BSD-style license found in LICENSE.txt

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform mat4 modelViewProjectionMatrix;

//------------------------------------------------------------------------------
#if defined(VERTEX_SHADER)

uniform float pointRadius = 0.1;   //# uiname=Point Radius; min=0.001; max=10
uniform float trimRadius = 1000000;//# uiname=Trim Radius; min=1; max=1000000
uniform int showSegIDzero = 0;     //# uiname=Show unsegmented; enum=No|Yes
uniform int selectionMode = 0;     //# uiname=Selection; enum=All segments|One segment
uniform int selectedSegment = 1;   //# uiname=Segment ID; min=0; max=1000
uniform float minPointSize = 0;
uniform float maxPointSize = 400.0;
// Point size multiplier to get from a width in projected coordinates to the
// number of pixels across as required for gl_PointSize
uniform float pointPixelScale = 0;
uniform vec3 cursorPos = vec3(0);
uniform int fileNumber = 0;
in float intensity;
in vec3 position;
in vec3 color;
in int returnNumber;
in int segment_id;
in int numberOfReturns;
in int pointSourceId;
in int classification;
in float plane_distance;

flat out float modifiedPointRadius;
flat out float pointScreenSize;
flat out vec3 pointColor;
flat out int markerShape;

void main()
{
    vec4 p = modelViewProjectionMatrix * vec4(position,1.0);
    float r = length(position.xy - cursorPos.xy);
    float trimFalloffLen = min(5, trimRadius/2);
    float trimScale = min(1, (trimRadius - r)/trimFalloffLen);
    modifiedPointRadius = pointRadius * trimScale;
    pointScreenSize = clamp(2*pointPixelScale*modifiedPointRadius / p.w, minPointSize, maxPointSize);
    markerShape = 1;

    // Set color based on segment_id
    int s = segment_id % 63+1;
    if (segment_id == 0) pointColor = vec3(0.000, 0.000, 0.000);
    else if (s == 1) pointColor = vec3(0.004, 0.000, 0.404);
    else if (s == 2) pointColor = vec3(0.835, 1.000, 0.000);
    else if (s == 3) pointColor = vec3(1.000, 0.000, 0.337);
    else if (s == 4) pointColor = vec3(0.620, 0.000, 0.557);
    else if (s == 5) pointColor = vec3(0.055, 0.298, 0.631);
    else if (s == 6) pointColor = vec3(1.000, 0.898, 0.008);
    else if (s == 7) pointColor = vec3(0.000, 0.373, 0.224);
    else if (s == 8) pointColor = vec3(0.000, 1.000, 0.000);
    else if (s == 9) pointColor = vec3(0.584, 0.000, 0.227);
    else if (s == 10) pointColor = vec3(1.000, 0.576, 0.494);
    else if (s == 11) pointColor = vec3(0.643, 0.141, 0.000);
    else if (s == 12) pointColor = vec3(0.000, 0.082, 0.267);
    else if (s == 13) pointColor = vec3(0.569, 0.816, 0.796);
    else if (s == 14) pointColor = vec3(0.384, 0.055, 0.000);
    else if (s == 15) pointColor = vec3(0.420, 0.408, 0.510);
    else if (s == 16) pointColor = vec3(0.000, 0.000, 1.000);
    else if (s == 17) pointColor = vec3(0.000, 0.490, 0.710);
    else if (s == 18) pointColor = vec3(0.416, 0.510, 0.424);
    else if (s == 19) pointColor = vec3(0.000, 0.682, 0.494);
    else if (s == 20) pointColor = vec3(0.761, 0.549, 0.624);
    else if (s == 21) pointColor = vec3(0.745, 0.600, 0.439);
    else if (s == 22) pointColor = vec3(0.000, 0.561, 0.612);
    else if (s == 23) pointColor = vec3(0.373, 0.678, 0.306);
    else if (s == 24) pointColor = vec3(1.000, 0.000, 0.000);
    else if (s == 25) pointColor = vec3(1.000, 0.000, 0.965);
    else if (s == 26) pointColor = vec3(1.000, 0.008, 0.616);
    else if (s == 27) pointColor = vec3(0.408, 0.239, 0.231);
    else if (s == 28) pointColor = vec3(1.000, 0.455, 0.639);
    else if (s == 29) pointColor = vec3(0.588, 0.541, 0.910);
    else if (s == 30) pointColor = vec3(0.596, 1.000, 0.322);
    else if (s == 31) pointColor = vec3(0.655, 0.341, 0.251);
    else if (s == 32) pointColor = vec3(0.004, 1.000, 0.996);
    else if (s == 33) pointColor = vec3(1.000, 0.933, 0.910);
    else if (s == 34) pointColor = vec3(0.996, 0.537, 0.000);
    else if (s == 35) pointColor = vec3(0.741, 0.776, 1.000);
    else if (s == 36) pointColor = vec3(0.004, 0.816, 1.000);
    else if (s == 37) pointColor = vec3(0.733, 0.533, 0.000);
    else if (s == 38) pointColor = vec3(0.459, 0.267, 0.694);
    else if (s == 39) pointColor = vec3(0.647, 1.000, 0.824);
    else if (s == 40) pointColor = vec3(1.000, 0.651, 0.996);
    else if (s == 41) pointColor = vec3(0.467, 0.302, 0.000);
    else if (s == 42) pointColor = vec3(0.478, 0.278, 0.510);
    else if (s == 43) pointColor = vec3(0.149, 0.204, 0.000);
    else if (s == 44) pointColor = vec3(0.000, 0.278, 0.329);
    else if (s == 45) pointColor = vec3(0.263, 0.000, 0.173);
    else if (s == 46) pointColor = vec3(0.710, 0.000, 1.000);
    else if (s == 47) pointColor = vec3(1.000, 0.694, 0.404);
    else if (s == 48) pointColor = vec3(1.000, 0.859, 0.400);
    else if (s == 49) pointColor = vec3(0.565, 0.984, 0.573);
    else if (s == 50) pointColor = vec3(0.494, 0.176, 0.824);
    else if (s == 51) pointColor = vec3(0.741, 0.827, 0.576);
    else if (s == 52) pointColor = vec3(0.898, 0.435, 0.996);
    else if (s == 53) pointColor = vec3(0.871, 1.000, 0.455);
    else if (s == 54) pointColor = vec3(0.000, 1.000, 0.471);
    else if (s == 55) pointColor = vec3(0.000, 0.608, 1.000);
    else if (s == 56) pointColor = vec3(0.000, 0.392, 0.004);
    else if (s == 57) pointColor = vec3(0.000, 0.463, 1.000);
    else if (s == 58) pointColor = vec3(0.522, 0.663, 0.000);
    else if (s == 59) pointColor = vec3(0.000, 0.725, 0.090);
    else if (s == 60) pointColor = vec3(0.471, 0.510, 0.192);
    else if (s == 61) pointColor = vec3(0.000, 1.000, 0.776);
    else if (s == 62) pointColor = vec3(1.000, 0.431, 0.255);
    else if (s == 63) pointColor = vec3(0.910, 0.369, 0.745);

    // show only one segment mode
    if (selectionMode == 1)
    {
        if (selectedSegment != segment_id)
            markerShape = -1;
    }
    // show unsegmented points or not
    if (segment_id==0) {
        if(showSegIDzero==1)
            markerShape = 1;
        else
            markerShape = -1;
    }

    // Ensure zero size points are discarded.  The actual minimum point size is
    // hardware and driver dependent, so set the markerShape to discarded for
    // good measure.
    if (pointScreenSize <= 0)
    {
        pointScreenSize = 0;
        markerShape = -1;
    }
    else if (pointScreenSize < 1)
    {
        // Clamp to minimum size of 1 to avoid aliasing with some drivers
        pointScreenSize = 1;
    }
    gl_PointSize = pointScreenSize;
    gl_Position = p;
}


//------------------------------------------------------------------------------
#elif defined(FRAGMENT_SHADER)

uniform float markerWidth = 0.3;

flat in float modifiedPointRadius;
flat in float pointScreenSize;
flat in vec3 pointColor;
flat in int markerShape;

out vec4 fragColor;

// Limit at which the point is rendered as a small square for antialiasing
// rather than using a specific marker shape
const float pointScreenSizeLimit = 2;
const float sqrt2 = 1.414213562;

void main()
{
    if (markerShape < 0) // markerShape == -1: discarded.
        discard;
    // (markerShape == 0: Square shape)
#   ifndef BROKEN_GL_FRAG_COORD
    gl_FragDepth = gl_FragCoord.z;
#   endif
    if (markerShape > 0 && pointScreenSize > pointScreenSizeLimit)
    {
        float w = markerWidth;
        if (pointScreenSize < 2*pointScreenSizeLimit)
        {
            // smoothly turn on the markers as we get close enough to see them
            w = mix(1, w, pointScreenSize/pointScreenSizeLimit - 1);
        }
        vec2 p = 2*(gl_PointCoord - 0.5);
        if (markerShape == 1) // shape: .
        {
            float r = length(p);
            if (r > 1)
                discard;
#           ifndef BROKEN_GL_FRAG_COORD
            gl_FragDepth += projectionMatrix[3][2] * gl_FragCoord.w*gl_FragCoord.w
                            // TODO: Why is the factor of 0.5 required here?
                            * 0.5*modifiedPointRadius*sqrt(1-r*r);
#           endif
        }
        else if (markerShape == 2) // shape: o
        {
            float r = length(p);
            if (r > 1 || r < 1 - w)
                discard;
        }
        else if (markerShape == 3) // shape: x
        {
            w *= 0.5*sqrt2;
            if (abs(p.x + p.y) > w && abs(p.x - p.y) > w)
                discard;
        }
        else if (markerShape == 4) // shape: +
        {
            w *= 0.5;
            if (abs(p.x) > w && abs(p.y) > w)
                discard;
        }
    }
    fragColor = vec4(pointColor, 1);
}

#endif

